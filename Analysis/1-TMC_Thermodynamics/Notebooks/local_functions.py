import csv
import json
import math
import os
import re
from collections import Counter, deque
from html import escape
from pathlib import Path

import numpy as np
from IPython.display import HTML, IFrame, SVG, display
from matplotlib import colormaps
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize, to_hex
from matplotlib.ticker import MaxNLocator
from rdkit import Chem, RDLogger
from rdkit.Chem import Draw, rdMolDescriptors
from scope.operations.graphs import compute_topological_distances
from scope.thermodynamics import get_Hvib, get_Svib

RDLogger.DisableLog('rdApp.warning')

# Function to load systems correctly, from paths included in the .csv output files
def resolve_system_path(path_from_csv, data_path):
    """Resolve a stored system path against the configured local Data directory."""
    path = os.path.abspath(os.path.expanduser(str(path_from_csv)))
    if os.path.isfile(path):
        return path

    normalized = str(path_from_csv).replace("\\", "/")
    if "/Data/" in normalized:
        normalized = normalized.split("/Data/", 1)[1]
    elif normalized.startswith("Data/"):
        normalized = normalized[len("Data/"):]
    return os.path.join(data_path, normalized.lstrip("/"))

def get_CH_volume(coord: list):
    import scipy
    from scipy.spatial import ConvexHull, QhullError
    coords = np.asarray(coord, dtype=float)
    if coords.shape[0] < 4: return np.nan
    try:
        ## Scipy returns a "hull" class. Here, we only report the volume
        return np.round(float(ConvexHull(coords).volume),2)
    except QhullError:
        return np.nan

#####################
#### CSV Related ####
#####################
def read_csv(path):
    with open(path, newline='') as f:
        return list(csv.DictReader(f))

def count_csv_entries(csv_path, column, value):
    count = 0
    with open(csv_path, newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            if row[column] == value:
                count += 1
    return count

def get_csv_entries(csv_path, column, value):
    matches = []
    with open(csv_path, newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            if row[column] == value:
                matches.append(row)
    return matches 

def filter_rows(csv_path, condition, first: bool=False):
    matches = []
    with open(csv_path, newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            if condition(row):
                matches.append(row)
                if first: return row
    return matches

def sum_csv_entries(csv_path, column):
    total = 0
    with open(csv_path, newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            value = row[column].strip()
            if value:
                total += int(value)
    return total

def get_values(csv_path, column):
    values = []
    with open(csv_path, newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            values.append(float(row[column]))
    return values

def write_csv_rows(path, rows, fieldnames=None):
    if not rows and fieldnames is None:
        raise ValueError('fieldnames are required when rows is empty')
    fieldnames = fieldnames or list(rows[0])
    with open(path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, extrasaction='ignore')
        writer.writeheader()
        writer.writerows(rows)

def sanitize_filename(text):
    text = str(text).strip().replace('1/2', '12').replace('/', '_').replace('\\', '_')
    text = re.sub(r'[^A-Za-z0-9._-]+', '_', text)
    return re.sub(r'_+', '_', text).strip('._-')

#################################
#### Optimized State Selection ##
#################################
def find_optimized_minimum(source):
    """Return the valid optimized minimum, preferring b3lyp_opt over corr_opt."""
    found, state = source.find_state('b3lyp_opt')
    if found and state.check_minimum():
        return True, state

    found, state = source.find_state('corr_opt')
    if found and state.check_minimum():
        return True, state

    return False, None

def get_hs_ls_minimum_states(system):
    found_hs, hs_mol = system.find_source('ref_hs_mol')
    found_ls, ls_mol = system.find_source('ref_ls_mol')
    if not (found_hs and found_ls):
        raise ValueError('Missing HS or LS reference molecule')

    found_hs_state, hs_state = find_optimized_minimum(hs_mol)
    found_ls_state, ls_state = find_optimized_minimum(ls_mol)
    if not (found_hs_state and found_ls_state):
        raise ValueError('Missing valid HS or LS optimized minimum')
    return hs_state, ls_state

#################################
#### Vibrational Analysis #######
#################################
def assess_vnm_mapping(mappings: list) -> dict:
    """Summarize the quality of a complete VNM mapping."""
    if len(mappings) == 0:
        raise ValueError('ASSESS_VNM_MAPPING: Mapping cannot be empty')

    overlaps              = np.asarray([mapping['overlap'] for mapping in mappings], dtype=float)
    frequency_differences = np.asarray([mapping['frequency_difference'] for mapping in mappings], dtype=float)
    costs                 = np.asarray([mapping['cost'] for mapping in mappings], dtype=float)

    return {
        'mapping_quality': float(np.exp(-np.mean(costs))),
        'mean_cost': float(np.mean(costs)),
        'mean_overlap': float(np.mean(overlaps)),
        'mean_squared_overlap': float(np.mean(overlaps**2)),
        'median_overlap': float(np.median(overlaps)),
        'minimum_overlap': float(np.min(overlaps)),
        'overlap_10th_percentile': float(np.percentile(overlaps, 10)),
        'fraction_overlap_above_0.8': float(np.mean(overlaps >= 0.8)),
        'fraction_overlap_above_0.5': float(np.mean(overlaps >= 0.5)),
        'frequency_mae': float(np.mean(frequency_differences)),
        'frequency_rmse': float(np.sqrt(np.mean(frequency_differences**2))),
        'maximum_frequency_difference': float(np.max(frequency_differences)),
    }

def get_fcs_internal_basis(coords, masses):
    """Return an orthonormal basis for internal motion of an atomic group."""
    coords = np.asarray(coords, dtype=float)
    masses = np.asarray(masses, dtype=float)
    if coords.ndim != 2 or coords.shape[1] != 3:
        raise ValueError('GET_FCS_INTERNAL_BASIS: Coordinates must have shape (natoms, 3)')
    if masses.shape != (len(coords),):
        raise ValueError('GET_FCS_INTERNAL_BASIS: Masses must match the coordinates')

    sqrt_masses = np.sqrt(masses)
    center_of_mass = np.average(coords, axis=0, weights=masses)
    centered_coords = coords - center_of_mass

    rigid_vectors = []
    for axis in np.eye(3):
        rigid_vectors.append((sqrt_masses[:, None] * axis).reshape(-1))
    for axis in np.eye(3):
        rotation = np.cross(axis, centered_coords)
        rigid_vectors.append((sqrt_masses[:, None] * rotation).reshape(-1))

    rigid_vectors = np.asarray(rigid_vectors).T
    left_vectors, singular_values, _ = np.linalg.svd(rigid_vectors, full_matrices=True)
    tolerance = max(rigid_vectors.shape) * np.finfo(float).eps * singular_values[0]
    rigid_rank = np.count_nonzero(singular_values > tolerance)
    if rigid_rank != 6:
        raise ValueError(f'GET_FCS_INTERNAL_BASIS: Expected six rigid directions, found {rigid_rank}')
    return left_vectors[:, rigid_rank:]

def project_vnms_on_fcs(vnms, coords, fcs_indices):
    """Project normalized molecular VNMs onto all and internal FCS motion."""
    if not vnms:
        raise ValueError('PROJECT_VNMS_ON_FCS: No VNMs were provided')

    fcs_indices = np.asarray(fcs_indices, dtype=int)
    masses = np.asarray(vnms[0].masses, dtype=float)[fcs_indices]
    internal_basis = get_fcs_internal_basis(np.asarray(coords)[fcs_indices], masses)
    fcs_weights = []
    internal_weights = []

    for vnm in vnms:
        mode = vnm.mass_weight_mode(permanent=False)
        norm = np.linalg.norm(mode)
        if not np.isfinite(norm) or norm <= 0.0:
            raise ValueError('PROJECT_VNMS_ON_FCS: VNM has an invalid norm')
        fcs_mode = (mode / norm)[fcs_indices].reshape(-1)
        coefficients = internal_basis.T @ fcs_mode
        fcs_weights.append(np.sum(fcs_mode**2))
        internal_weights.append(np.sum(coefficients**2))

    return np.asarray(fcs_weights), np.asarray(internal_weights), internal_basis

def get_mode_entropies(vnms, temperature):
    """Return one harmonic vibrational-entropy contribution per VNM."""
    entropies = []
    for vnm in vnms:
        entropy = get_Svib([vnm.freq_cm], temperature, freq_units='cm', outunits='kj')
        entropies.append(entropy.value * 1000.0)
    return np.asarray(entropies)

def get_mode_enthalpies(vnms, temperature):
    """Return one harmonic vibrational-enthalpy contribution per VNM."""
    enthalpies = []
    for vnm in vnms:
        enthalpy = get_Hvib([vnm.freq_cm], temperature, freq_units='cm', outunits='kj')
        enthalpies.append(enthalpy.value)
    return np.asarray(enthalpies)

def get_ligand_metal_weights(vnms, ligand_indices, metal_index):
    """Partition each VNM's atomic participation among ligands and the metal."""
    ligand_weights = []
    metal_weights = []

    for vnm in vnms:
        atomic_participation = vnm.get_atomic_participation()
        this_ligand_weights = []
        for indices in ligand_indices:
            this_ligand_weights.append(np.sum(atomic_participation[indices]))
        ligand_weights.append(this_ligand_weights)
        metal_weights.append(atomic_participation[metal_index])

    return np.asarray(ligand_weights), np.asarray(metal_weights)

def bin_mode_contributions(frequencies, contributions, frequency_bins):
    """Sum mode-resolved values inside frequency windows."""
    return np.histogram(frequencies, bins=frequency_bins, weights=contributions)[0]

def get_weighted_frequency_metrics(hs_vnms, ls_vnms, hs_weights, ls_weights):
    """Return pair-weighted state means and their mapped HS-minus-LS shift."""
    hs_frequencies = np.asarray([vnm.freq_cm for vnm in hs_vnms], dtype=float)
    ls_frequencies = np.asarray([vnm.freq_cm for vnm in ls_vnms], dtype=float)
    hs_weights = np.asarray(hs_weights, dtype=float)
    ls_weights = np.asarray(ls_weights, dtype=float)

    if not (hs_frequencies.shape == ls_frequencies.shape == hs_weights.shape == ls_weights.shape):
        raise ValueError('GET_WEIGHTED_FREQUENCY_METRICS: Frequencies and weights must have equal shapes')

    valid = (
        np.isfinite(hs_frequencies)
        & np.isfinite(ls_frequencies)
        & np.isfinite(hs_weights)
        & np.isfinite(ls_weights)
        & (hs_frequencies > 0.0)
        & (ls_frequencies > 0.0)
    )
    if not np.any(valid):
        raise ValueError('GET_WEIGHTED_FREQUENCY_METRICS: No positive mapped mode pairs were found')

    hs_frequencies = hs_frequencies[valid]
    ls_frequencies = ls_frequencies[valid]
    hs_weights = hs_weights[valid]
    ls_weights = ls_weights[valid]
    pair_weights = 0.5 * (hs_weights + ls_weights)

    if np.sum(hs_weights) <= 0.0 or np.sum(ls_weights) <= 0.0 or np.sum(pair_weights) <= 0.0:
        raise ValueError('GET_WEIGHTED_FREQUENCY_METRICS: Participation weights must have a positive sum')

    hs_mean_frequency = float(np.average(hs_frequencies, weights=pair_weights))
    ls_mean_frequency = float(np.average(ls_frequencies, weights=pair_weights))
    return {
        'HS_mean_frequency': hs_mean_frequency,
        'LS_mean_frequency': ls_mean_frequency,
        'dnu_mapped': hs_mean_frequency - ls_mean_frequency,
        'num_mode_pairs': int(np.count_nonzero(valid)),
    }

def get_atomic_frequency_shifts(hs_vnms, ls_vnms):
    """Return the participation-weighted mapped HS-minus-LS shift for each atom."""
    if len(hs_vnms) != len(ls_vnms) or len(hs_vnms) == 0:
        raise ValueError('GET_ATOMIC_FREQUENCY_SHIFTS: Mapped VNM sets must have equal nonzero lengths')

    hs_frequencies = np.asarray([vnm.freq_cm for vnm in hs_vnms], dtype=float)
    ls_frequencies = np.asarray([vnm.freq_cm for vnm in ls_vnms], dtype=float)
    hs_participation = np.asarray([vnm.get_atomic_participation() for vnm in hs_vnms])
    ls_participation = np.asarray([vnm.get_atomic_participation() for vnm in ls_vnms])

    if hs_participation.shape != ls_participation.shape:
        raise ValueError('GET_ATOMIC_FREQUENCY_SHIFTS: Atomic participation arrays must have equal shapes')
    if hs_participation.ndim != 2 or hs_participation.shape[0] != len(hs_frequencies):
        raise ValueError('GET_ATOMIC_FREQUENCY_SHIFTS: Expected one atom-participation row per VNM')

    valid = (
        np.isfinite(hs_frequencies)
        & np.isfinite(ls_frequencies)
        & (hs_frequencies > 0.0)
        & (ls_frequencies > 0.0)
    )
    if not np.any(valid):
        raise ValueError('GET_ATOMIC_FREQUENCY_SHIFTS: No positive mapped mode pairs were found')

    pair_participation = 0.5 * (hs_participation[valid] + ls_participation[valid])
    frequency_shifts = hs_frequencies[valid] - ls_frequencies[valid]
    total_participation = np.sum(pair_participation, axis=0)
    weighted_shifts = np.sum(pair_participation * frequency_shifts[:, None], axis=0)

    return np.divide(
        weighted_shifts,
        total_participation,
        out=np.full(total_participation.shape, np.nan),
        where=total_participation > 0.0,
    )

def plot_representative_frequency_shifts(rows, output_path=None):
    """Compare ligand-, FCS-, and Fe-weighted shifts for representative systems."""
    components = [
        ('All ligands', 'dnu_all_ligands[cm-1]'),
        ('FCS', 'dnu_FCS[cm-1]'),
        ('Fe atom', 'dnu_metal[cm-1]'),
    ]
    colors = ['#d95f02', '#1b9e77']
    y_positions = np.arange(len(components))
    bar_height = 0.34
    fig, ax = plt.subplots(figsize=(6.4, 3.0), tight_layout=True)

    for index, (row, color) in enumerate(zip(rows, colors)):
        shifts = [float(row[column]) for _, column in components]
        offset = (index - (len(rows) - 1) / 2) * bar_height
        ax.barh(
            y_positions + offset,
            shifts,
            height=0.9 * bar_height,
            color=color,
            edgecolor='black',
            linewidth=0.6,
            label=f"{row['refcode']} ({row['group']})",
        )

    ax.set_yticks(y_positions, [label for label, _ in components])
    ax.set_xlabel(r'$\Delta\bar{\nu}_G$ [cm$^{-1}$]')
    ax.axvline(0.0, color='black', linewidth=0.9)
    ax.grid(axis='x', alpha=0.25)
    ax.legend(
        frameon=False,
        loc='lower center',
        bbox_to_anchor=(0.5, 1.02),
        ncol=len(rows),
    )
    ax.invert_yaxis()
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_color('black')
        spine.set_linewidth(1.0)

    if output_path:
        fig.savefig(output_path, dpi=400, bbox_inches='tight')
    plt.show()
    return fig, ax

def plot_vnm_mapping_summary(rows, output_path=None):
    """Plot compact distributions of the principal VNM-mapping metrics."""
    metrics = [
        ('mapping_quality', 'Mapping quality'),
        ('mean_overlap', 'Mean mode overlap'),
        ('frequency_mae', r'Mean $|\Delta\nu|$ [cm$^{-1}$]'),
    ]
    fig, axs = plt.subplots(1, len(metrics), figsize=(10, 3.0), tight_layout=True)

    for ax, (column, label) in zip(axs, metrics):
        values = np.asarray([
            float(row[column]) for row in rows
            if row.get(column, '') not in ('', None) and np.isfinite(float(row[column]))
        ])
        if not len(values):
            raise ValueError(f'PLOT_VNM_MAPPING_SUMMARY: No values found for {column}')
        num_bins = min(25, max(8, int(np.sqrt(len(values)))))
        median = np.median(values)
        ax.hist(values, bins=num_bins, color='#547f86', edgecolor='black', linewidth=0.5)
        ax.axvline(median, color='#b43c2f', linestyle='--', linewidth=1.3, label=f'Median: {median:.2f}')
        ax.set_xlabel(label)
        ax.set_ylabel('Systems')
        ax.legend(frameon=False, fontsize=8)

    if output_path:
        fig.savefig(output_path, dpi=400, bbox_inches='tight')
    plt.show()
    return fig, axs

def plot_binned_vnm_group_contributions(
    values_by_metric,
    group_labels,
    bin_labels,
    contribution_label,
    output_path=None,
    show_titles=True,
):
    """Plot mean frequency-binned VNM contributions for chemical groups."""
    metrics = list(values_by_metric)
    num_columns = 2
    num_rows = int(np.ceil(len(metrics) / num_columns))
    fig, axs = plt.subplots(num_rows, num_columns, figsize=(10, 3.1 * num_rows), constrained_layout=True)
    axs = np.asarray(axs).reshape(-1)
    colors = colormaps.get_cmap('turbo')(np.linspace(0.05, 0.95, len(bin_labels)))
    y_positions = np.arange(len(group_labels))

    for ax, metric in zip(axs, metrics):
        values = np.nan_to_num(np.asarray(values_by_metric[metric], dtype=float), nan=0.0)
        left_positive = np.zeros(len(group_labels))
        left_negative = np.zeros(len(group_labels))

        for bin_index, (bin_label, color) in enumerate(zip(bin_labels, colors)):
            positive = np.where(values[:, bin_index] > 0, values[:, bin_index], 0.0)
            negative = np.where(values[:, bin_index] < 0, values[:, bin_index], 0.0)
            ax.barh(y_positions, positive, left=left_positive, color=color, edgecolor='white', linewidth=0.5, label=bin_label)
            ax.barh(y_positions, negative, left=left_negative, color=color, edgecolor='white', linewidth=0.5)
            left_positive += positive
            left_negative += negative

        ax.set_yticks(y_positions, group_labels)
        ax.set_xlabel(contribution_label)
        if show_titles:
            ax.set_title(metric)
        ax.axvline(0, color='black', linewidth=0.8)
        ax.grid(axis='x', alpha=0.2)

    for ax in axs[len(metrics):]:
        ax.set_visible(False)

    handles, labels = axs[0].get_legend_handles_labels()
    fig.legend(
        handles,
        labels,
        title='Frequency window [cm$^{-1}$]',
        frameon=False,
        loc='center left',
        bbox_to_anchor=(1.01, 0.5),
    )
    if output_path:
        fig.savefig(output_path, dpi=400, bbox_inches='tight')
    plt.show()
    return fig, axs

def binned_vibrational_values(freqs_cm, temperature, freq_bins_cm):
    freqs_cm = np.asarray(freqs_cm, dtype=float)
    abs_freqs_cm = np.abs(freqs_cm)
    entropy = []
    enthalpy = []
    mode_counts = []
    for low, high in zip(freq_bins_cm[:-1], freq_bins_cm[1:]):
        # Bin by magnitude, but retain each sign so SCOPE applies its default
        # imaginary-frequency treatment.
        bin_freqs = freqs_cm[(abs_freqs_cm >= low) & (abs_freqs_cm < high)]
        entropy.append(1000 * get_Svib(bin_freqs, temperature, freq_units='cm', outunits='kj').value)
        enthalpy.append(get_Hvib(bin_freqs, temperature, freq_units='cm', outunits='kj').value)
        mode_counts.append(np.count_nonzero(bin_freqs > 0))

    entropy = np.asarray(entropy)
    enthalpy = np.asarray(enthalpy)
    mode_counts = np.asarray(mode_counts, dtype=int)
    total_modes = mode_counts.sum()
    mode_fractions = mode_counts / total_modes if total_modes else np.zeros_like(mode_counts, dtype=float)
    entropy_per_mode = np.divide(
        entropy,
        mode_counts,
        out=np.full_like(entropy, np.nan, dtype=float),
        where=mode_counts > 0,
    )
    return entropy, enthalpy, mode_counts, mode_fractions, entropy_per_mode

def build_vibrational_bin_rows(
    systems,
    thermo_by_refcode,
    temperature,
    freq_bins_cm,
    entropy_method,
    enthalpy_method,
):
    rows = []
    failures = []

    for system in systems:
        refcode = system.refcode_wo_digits
        try:
            thermo_row = thermo_by_refcode[refcode]
            hs_state, ls_state = get_hs_ls_minimum_states(system)
            hs_entropy, hs_enthalpy, hs_mode_counts, hs_mode_fractions, hs_entropy_per_mode = binned_vibrational_values(
                hs_state.freqs_cm, temperature, freq_bins_cm
            )
            ls_entropy, ls_enthalpy, ls_mode_counts, ls_mode_fractions, ls_entropy_per_mode = binned_vibrational_values(
                ls_state.freqs_cm, temperature, freq_bins_cm
            )

            for i, (low, high) in enumerate(zip(freq_bins_cm[:-1], freq_bins_cm[1:])):
                rows.append({
                    'refcode': refcode,
                    'denticity': thermo_row['denticity'],
                    'homoleptic': thermo_row['homoleptic'] == 'True',
                    'HS_minimum_state': hs_state.name,
                    'LS_minimum_state': ls_state.name,
                    'temperature[K]': temperature,
                    'entropy_method': entropy_method,
                    'enthalpy_method': enthalpy_method,
                    'bin_low[cm-1]': low,
                    'bin_high[cm-1]': high,
                    'bin_label[cm-1]': f'{int(low)}-{int(high)}',
                    'HS_Svib[J·mol-1·K-1]': hs_entropy[i],
                    'LS_Svib[J·mol-1·K-1]': ls_entropy[i],
                    'dSvib[J·mol-1·K-1]': hs_entropy[i] - ls_entropy[i],
                    'HS_Hvib[kJ·mol-1]': hs_enthalpy[i],
                    'LS_Hvib[kJ·mol-1]': ls_enthalpy[i],
                    'dHvib[kJ·mol-1]': hs_enthalpy[i] - ls_enthalpy[i],
                    'HS_n_modes': hs_mode_counts[i],
                    'LS_n_modes': ls_mode_counts[i],
                    'd_n_modes': hs_mode_counts[i] - ls_mode_counts[i],
                    'HS_mode_fraction': hs_mode_fractions[i],
                    'LS_mode_fraction': ls_mode_fractions[i],
                    'd_mode_fraction': hs_mode_fractions[i] - ls_mode_fractions[i],
                    'HS_mean_S_per_mode[J·mol-1·K-1]': (
                        hs_entropy_per_mode[i] if np.isfinite(hs_entropy_per_mode[i]) else None
                    ),
                    'LS_mean_S_per_mode[J·mol-1·K-1]': (
                        ls_entropy_per_mode[i] if np.isfinite(ls_entropy_per_mode[i]) else None
                    ),
                    'd_mean_S_per_mode[J·mol-1·K-1]': (
                        hs_entropy_per_mode[i] - ls_entropy_per_mode[i]
                        if np.isfinite(hs_entropy_per_mode[i]) and np.isfinite(ls_entropy_per_mode[i])
                        else None
                    ),
                })
        except Exception as exc:
            failures.append((refcode, repr(exc)))

    return rows, failures

def mean_vibrational_bin_values(rows, value_column, denticities, bin_labels):
    available_denticities = {str(row['denticity']) for row in rows}
    ordered = [str(denticity) for denticity in denticities if str(denticity) in available_denticities]
    means = np.full((len(ordered), len(bin_labels)), np.nan, dtype=float)

    for i, denticity in enumerate(ordered):
        for j, bin_label in enumerate(bin_labels):
            values = []
            for row in rows:
                if str(row['denticity']) != denticity or row['bin_label[cm-1]'] != bin_label:
                    continue
                if row.get(value_column, '') in ('', None):
                    continue
                value = float(row[value_column])
                if np.isfinite(value):
                    values.append(value)
            if values:
                means[i, j] = np.mean(values)

    return means, ordered

def plot_frequency_contributions(values, denticities, bin_labels, xlabel, output_path=None):
    values = np.nan_to_num(np.asarray(values, dtype=float), nan=0.0)
    fig, ax = plt.subplots(figsize=(7.2, 3.4), tight_layout=True)
    color_positions = np.linspace(0.05, 0.95, len(bin_labels)) if len(bin_labels) > 1 else [0.5]
    colors = colormaps.get_cmap('turbo')(color_positions)
    y_pos = np.arange(len(denticities))
    left_pos = np.zeros(len(denticities))
    left_neg = np.zeros(len(denticities))

    for j, (label, color) in enumerate(zip(bin_labels, colors)):
        positive = np.where(values[:, j] > 0, values[:, j], 0)
        negative = np.where(values[:, j] < 0, values[:, j], 0)
        ax.barh(y_pos, positive, left=left_pos, color=color, edgecolor='white', linewidth=0.6, label=f'{label} cm$^{{-1}}$')
        ax.barh(y_pos, negative, left=left_neg, color=color, edgecolor='white', linewidth=0.6)
        left_pos += positive
        left_neg += negative

    ax.set_yticks(y_pos, [f'd={denticity}' for denticity in denticities])
    ax.set_xlabel(xlabel)
    ax.axvline(0, color='black', linewidth=0.9)
    ax.legend(title='Frequency window', frameon=False, fontsize=8, title_fontsize=9, bbox_to_anchor=(1.02, 1), loc='upper left')
    ax.grid(axis='x', alpha=0.25)
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_color('black')
        spine.set_linewidth(1.0)

    if output_path:
        fig.savefig(output_path, dpi=400)
    plt.show()
    return fig, ax

def plot_hs_ls_frequency_profiles(
    hs_values,
    ls_values,
    denticities,
    bin_labels,
    ylabel,
    output_path=None,
):
    hs_values = np.asarray(hs_values, dtype=float)
    ls_values = np.asarray(ls_values, dtype=float)
    x = np.arange(len(bin_labels))
    colors = colormaps.get_cmap('tab10')(np.linspace(0, 0.6, len(denticities)))
    bar_width = 0.8 / len(denticities)
    fig, axs = plt.subplots(1, 2, figsize=(10, 3.6), sharey=True, constrained_layout=True)

    for ax, values, title in zip(axs, [hs_values, ls_values], ['HS', 'LS']):
        for i, (denticity, color) in enumerate(zip(denticities, colors)):
            offset = (i - (len(denticities) - 1) / 2) * bar_width
            ax.bar(
                x + offset, values[i], width=0.9 * bar_width,
                color=color, edgecolor='black', linewidth=0.4,
                label=f'd={denticity}',
            )
        ax.set_title(title)
        ax.set_xticks(x, bin_labels, rotation=45, ha='right')
        ax.set_xlabel(r'Frequency window [cm$^{-1}$]')
        ax.grid(axis='y', alpha=0.25)
        for spine in ax.spines.values():
            spine.set_visible(True)
            spine.set_color('black')
            spine.set_linewidth(1.0)

    axs[0].set_ylabel(ylabel)
    handles, labels = axs[0].get_legend_handles_labels()
    fig.legend(
        handles, labels, frameon=False, ncol=len(denticities),
        loc='upper center', bbox_to_anchor=(0.5, 1.08),
    )

    if output_path:
        fig.savefig(output_path, dpi=400, bbox_inches='tight')
    plt.show()
    return fig, axs

#################################
#### Structural Correlations ####
#################################
def first_float(value):
    if value in (None, '', 'None', 'nan', 'NaN', 'error'):
        return None
    if isinstance(value, (int, float)):
        return float(value)
    text = str(value).strip().strip('[]')
    if not text:
        return None
    try:
        return float(text.replace(',', ' ').split()[0])
    except (IndexError, ValueError):
        return None

def safe_diff(a, b):
    return a - b if a is not None and b is not None else None

def safe_abs(value):
    return abs(value) if value is not None else None

def safe_ratio(a, b):
    return a / b if a is not None and b not in (None, 0) else None

def base_refcode(refcode):
    text = str(refcode)
    return text[:6] if len(text) >= 6 and text[:6].isalpha() else text

def sorted_denticities(rows):
    return sorted(d for d in {row['denticity'] for row in rows} if d != 'unknown')

def finite_pairs(subset, feature, target):
    return [
        (row[feature], row[target])
        for row in subset
        if row[feature] is not None and row[target] is not None
        and math.isfinite(row[feature]) and math.isfinite(row[target])
    ]

def pearson_correlation(subset, feature, target='dHelec'):
    pairs = finite_pairs(subset, feature, target)
    if len(pairs) < 3:
        return None, len(pairs)

    x = np.asarray([pair[0] for pair in pairs], dtype=float)
    y = np.asarray([pair[1] for pair in pairs], dtype=float)
    if np.std(x) == 0 or np.std(y) == 0:
        return None, len(pairs)
    return float(np.corrcoef(x, y)[0, 1]), len(pairs)

def correlation_table(subset, label, target='dHelec', feature_list=()):
    table = []
    for feature in feature_list:
        correlation, num_entries = pearson_correlation(subset, feature, target=target)
        if correlation is not None:
            table.append({
                'group': label,
                'target': target,
                'feature': feature,
                'r': correlation,
                'abs_r': abs(correlation),
                'n': num_entries,
            })
    return sorted(table, key=lambda row: row['abs_r'], reverse=True)

def print_correlation_table(table, limit=None):
    print(f"{'target':<8} {'group':<14} {'feature':<14} {'r':>8} {'n':>5}")
    print('-' * 57)
    selected_rows = table if limit is None else table[:limit]
    for row in selected_rows:
        print(f"{row['target']:<8} {row['group']:<14} {row['feature']:<14} {row['r']:>8.3f} {row['n']:>5}")

def fit_line(x, y):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    slope, intercept = np.polyfit(x, y, 1)
    x_fit = np.linspace(x.min(), x.max(), 200)
    y_fit = slope * x_fit + intercept
    y_pred = slope * x + intercept
    r2 = 1 - np.sum((y - y_pred) ** 2) / np.sum((y - np.mean(y)) ** 2)
    return x_fit, y_fit, r2

def scatter_by_group(rows, feature, target='dHelec', feature_labels=None, target_labels=None):
    feature_labels = feature_labels or {}
    target_labels = target_labels or {}
    fig, ax = plt.subplots(figsize=(5, 4), tight_layout=True)
    groups = [
        ('homoleptic', [row for row in rows if row['homoleptic']], 'tab:blue'),
        ('heteroleptic', [row for row in rows if not row['homoleptic']], 'tab:orange'),
    ]

    for label, subset, color in groups:
        pairs = finite_pairs(subset, feature, target)
        if not pairs:
            continue
        correlation, num_entries = pearson_correlation(subset, feature, target=target)
        label_text = (
            f'{label}: r={correlation:.2f}, n={num_entries}'
            if correlation is not None else f'{label}: n={num_entries}'
        )
        ax.scatter(
            [pair[0] for pair in pairs], [pair[1] for pair in pairs],
            s=22, alpha=0.75, edgecolor='black', linewidth=0.3,
            color=color, label=label_text,
        )

    ax.set_xlabel(feature_labels.get(feature, feature))
    ax.set_ylabel(target_labels.get(target, target))
    ax.legend(frameon=False, fontsize=9)
    ax.grid(alpha=0.25)
    return fig, ax

def scatter_by_denticity(rows, feature, target='dHelec', feature_labels=None, target_labels=None):
    feature_labels = feature_labels or {}
    target_labels = target_labels or {}
    colors = {'1': 'purple', '2': 'green', '3': 'blue', '6': 'orange', 'unknown': 'pink'}
    fig, ax = plt.subplots(figsize=(5, 4), tight_layout=True)

    for denticity in sorted_denticities(rows):
        subset = [row for row in rows if row['denticity'] == denticity]
        pairs = finite_pairs(subset, feature, target)
        if not pairs:
            continue
        correlation, num_entries = pearson_correlation(subset, feature, target=target)
        label_text = (
            f'dent={denticity}: r={correlation:.2f}, n={num_entries}'
            if correlation is not None else f'dent={denticity}: n={num_entries}'
        )
        ax.scatter(
            [pair[0] for pair in pairs], [pair[1] for pair in pairs],
            s=22, alpha=0.75, edgecolor='black', linewidth=0.3,
            color=colors.get(denticity, 'gray'), label=label_text,
        )

    ax.set_xlabel(feature_labels.get(feature, feature))
    ax.set_ylabel(target_labels.get(target, target))
    ax.legend(frameon=False, fontsize=8)
    ax.grid(alpha=0.25)
    return fig, ax

def scatter_feature_by_denticity(
    rows,
    feature,
    target='dHelec',
    min_points_for_fit=6,
    feature_labels=None,
    target_labels=None,
):
    feature_labels = feature_labels or {}
    target_labels = target_labels or {}
    colors = {'1': 'purple', '2': 'green', '3': 'blue', '6': 'orange', 'unknown': 'pink'}
    fig, ax = plt.subplots(figsize=(5.5, 4), tight_layout=True)

    for denticity in sorted_denticities(rows):
        subset = [row for row in rows if row['denticity'] == denticity]
        pairs = finite_pairs(subset, feature, target)
        if not pairs:
            continue
        x = [pair[0] for pair in pairs]
        y = [pair[1] for pair in pairs]
        correlation, num_entries = pearson_correlation(subset, feature, target=target)
        label_text = (
            f'dent={denticity}: r={correlation:.2f}, n={num_entries}'
            if correlation is not None else f'dent={denticity}: n={num_entries}'
        )
        color = colors.get(denticity, 'gray')
        ax.scatter(
            x, y, s=24, alpha=0.75, edgecolor='black',
            linewidth=0.3, color=color, label=label_text,
        )
        if correlation is not None and num_entries >= min_points_for_fit and abs(correlation) > 0.5:
            x_fit, y_fit, _ = fit_line(x, y)
            ax.plot(x_fit, y_fit, color=color, linewidth=1.2, alpha=0.8)

    ax.set_xlabel(feature_labels.get(feature, feature))
    ax.set_ylabel(target_labels.get(target, target))
    ax.legend(frameon=False, fontsize=8)
    ax.grid(alpha=0.25)
    return fig, ax

def draw_correlation_panel(ax, spec, rows, plot_groups):
    feature, target, x_scale, xlabel, ylabel = spec

    for group_label, subset, color in plot_groups:
        pairs = finite_pairs(subset, feature, target)
        if not pairs:
            continue
        x = x_scale * np.asarray([pair[0] for pair in pairs], dtype=float)
        y = np.asarray([pair[1] for pair in pairs], dtype=float)
        ax.scatter(
            x, y, s=18, alpha=0.55, color=color, edgecolor='black',
            linewidth=0.25, label=group_label,
        )

    all_pairs = finite_pairs(rows, feature, target)
    x_all = x_scale * np.asarray([pair[0] for pair in all_pairs], dtype=float)
    y_all = np.asarray([pair[1] for pair in all_pairs], dtype=float)
    x_fit, y_fit, _ = fit_line(x_all, y_all)
    ax.plot(x_fit, y_fit, color='black', linewidth=1.4, zorder=5)

    correlation, num_entries = pearson_correlation(rows, feature, target=target)
    ax.text(
        0.97, 0.96, f'$r={correlation:.2f}$\n$N={num_entries}$', transform=ax.transAxes,
        ha='right', va='top', fontsize=9,
        bbox=dict(boxstyle='round', facecolor='white', edgecolor='black', alpha=0.85),
    )
    ax.set_xlabel(xlabel, fontsize=11)
    ax.set_ylabel(ylabel, fontsize=11)
    ax.tick_params(axis='both', direction='out', labelsize=9)
    ax.grid(alpha=0.20)
    for spine in ax.spines.values():
        spine.set_color('black')
        spine.set_linewidth(1.0)
    return ax.get_legend_handles_labels()

def build_family_structural_data(rows, membership_rows, statistics_rows):
    structural_by_id = {base_refcode(row['id']): row for row in rows}
    statistics_by_family = {
        (str(row['Denticity']), row['Family']): row for row in statistics_rows
    }
    family_structural_rows = []

    for membership in membership_rows:
        if str(membership['Main Distribution']).strip().lower() != 'true':
            continue
        refcode = base_refcode(membership['Refcode'])
        structural = structural_by_id.get(refcode)
        if structural is None:
            continue
        denticity = str(membership['Denticity'])
        family = membership['Family']
        family_structural_rows.append({
            'id': refcode,
            'Denticity': int(float(denticity)),
            'Family': family,
            'Family Label': f'd={int(float(denticity))}, {family}',
            'HS_FeN': structural['HS_FeN'],
            'HS_Epsylon': structural['HS_Epsylon'],
            'LS_FeN': structural['LS_FeN'],
            'LS_Epsylon': structural['LS_Epsylon'],
        })

    grouped_rows = {}
    for row in family_structural_rows:
        grouped_rows.setdefault((str(row['Denticity']), row['Family']), []).append(row)

    family_structural_statistics = []
    for (denticity, family), family_rows in grouped_rows.items():
        electronic = statistics_by_family[(denticity, family)]
        record = {
            'Denticity': int(denticity),
            'Family': family,
            'Num Main-Distribution Entries': len(family_rows),
            'dHelec Base [kJ/mol]': float(electronic['dHelec Base [kJ/mol]']),
            'dHelec Width [kJ/mol]': float(electronic['dHelec Width [kJ/mol]']),
        }
        for state in ['HS', 'LS']:
            for parameter, label, unit in [('FeN', 'dFeN', 'Angstrom'), ('Epsylon', 'Sigma', 'deg')]:
                values = np.asarray([row[f'{state}_{parameter}'] for row in family_rows], dtype=float)
                record[f'{state} {label} Median [{unit}]'] = float(np.median(values))
                record[f'{state} {label} IQR [{unit}]'] = float(np.quantile(values, 0.75) - np.quantile(values, 0.25))
        family_structural_statistics.append(record)

    return family_structural_rows, family_structural_statistics

######################
#### Ligand Split ####
######################
def molecule_graph_with_current_metal_connectivity(molecule, debug: int=0):
    molecule_graph = molecule.get_graph(debug=debug).copy()
    metal_indices = {metal.get_parent_index('molecule') for metal in molecule.metals}

    # Bond objects can retain metal contacts inconsistent with the current madjmat.
    molecule_graph.remove_edges_from([
        (first, second)
        for first, second in molecule_graph.edges()
        if first in metal_indices or second in metal_indices
    ])
    for metal_index in metal_indices:
        for atom_index, connected in enumerate(molecule.madjmat[metal_index]):
            if connected and atom_index != metal_index:
                distance = np.linalg.norm(
                    np.asarray(molecule.atoms[metal_index].coord, dtype=float)
                    - np.asarray(molecule.atoms[atom_index].coord, dtype=float)
                )
                molecule_graph.add_edge(metal_index, atom_index, distance=float(distance))
    return molecule_graph

def get_metal_bound_ligand_atoms(ligand):
    ligand_indices = ligand.get_parent_indices('molecule')
    molecule       = ligand.get_parent('molecule')
    metal_index    = molecule.metals[0].get_parent_index('molecule')
    return {
        ligand_index
        for ligand_index, molecule_index in enumerate(ligand_indices)
        if molecule.madjmat[metal_index, molecule_index] > 0
    }

def metal_connected_ring_core_atoms(mol, metal_bound_atoms):
    if not metal_bound_atoms:
        return set()

    rings = [set(ring) for ring in mol.GetRingInfo().AtomRings()]
    selected_ring_atoms = set()
    seen = set()

    for start in range(len(rings)):
        if start in seen:
            continue

        stack = [start]
        ring_system_atoms = set()

        while stack:
            ring_id = stack.pop()
            if ring_id in seen:
                continue

            seen.add(ring_id)
            ring_system_atoms.update(rings[ring_id])

            for other_id, other_ring in enumerate(rings):
                if other_id not in seen and rings[ring_id] & other_ring:
                    stack.append(other_id)

        if ring_system_atoms & metal_bound_atoms:
            selected_ring_atoms.update(ring_system_atoms)

    return selected_ring_atoms

def donor_backbone_core_atoms(mol, metal_bound_atoms, monodentate_radius=2):
    donors = sorted(metal_bound_atoms)
    if not donors:
        return set()

    if len(donors) == 1:
        backbone = {donors[0]}
        frontier = {donors[0]}
        for _ in range(monodentate_radius):
            next_frontier = {
                nbr.GetIdx()
                for atom_idx in frontier
                for nbr in mol.GetAtomWithIdx(atom_idx).GetNeighbors()
                if nbr.GetAtomicNum() != 1 and nbr.GetIdx() not in backbone
            }
            backbone.update(next_frontier)
            frontier = next_frontier
        return backbone

    backbone = set(donors)
    for i, first in enumerate(donors):
        for second in donors[i + 1:]:
            backbone.update(Chem.GetShortestPath(mol, first, second))
    return {idx for idx in backbone if mol.GetAtomWithIdx(idx).GetAtomicNum() != 1}

def mol_fragment_from_atom_indices(mol, atom_indices):
    atom_indices = sorted(set(int(idx) for idx in atom_indices))
    if not atom_indices:
        return Chem.MolFromSmiles('')

    index_map = {}
    fragment = Chem.RWMol()
    for old_idx in atom_indices:
        atom = Chem.Atom(mol.GetAtomWithIdx(old_idx))
        # Prevent RDKit from saturating cut bonds with new implicit hydrogens.
        atom.SetNoImplicit(True)
        index_map[old_idx] = fragment.AddAtom(atom)

    atom_set = set(atom_indices)
    for bond in mol.GetBonds():
        begin = bond.GetBeginAtomIdx()
        end = bond.GetEndAtomIdx()
        if begin in atom_set and end in atom_set:
            fragment.AddBond(index_map[begin], index_map[end], bond.GetBondType())
            new_bond = fragment.GetBondBetweenAtoms(index_map[begin], index_map[end])
            new_bond.SetIsAromatic(bond.GetIsAromatic())

    fragment = fragment.GetMol()
    fragment.UpdatePropertyCache(strict=False)
    Chem.GetSymmSSSR(fragment)
    return fragment

def fragment_composition_label(fragment):
    counts = {}
    for atom in fragment.GetAtoms():
        if atom.GetAtomicNum() == 1:
            continue
        symbol = atom.GetSymbol()
        counts[symbol] = counts.get(symbol, 0) + 1

    return ''.join(f'{symbol}{counts[symbol] if counts[symbol] > 1 else ""}' for symbol in sorted(counts)) or 'H'

def topology_graph_label(mol):
    heavy_atoms = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomicNum() != 1]
    index_map = {}
    topology = Chem.RWMol()

    for old_idx in heavy_atoms:
        atom = Chem.Atom(mol.GetAtomWithIdx(old_idx).GetAtomicNum())
        atom.SetFormalCharge(0)
        atom.SetIsAromatic(False)
        atom.SetNoImplicit(True)
        index_map[old_idx] = topology.AddAtom(atom)

    heavy_atom_set = set(heavy_atoms)
    for bond in mol.GetBonds():
        begin = bond.GetBeginAtomIdx()
        end = bond.GetEndAtomIdx()
        if begin in heavy_atom_set and end in heavy_atom_set:
            topology.AddBond(index_map[begin], index_map[end], Chem.BondType.SINGLE)

    topology = topology.GetMol()
    topology.UpdatePropertyCache(strict=False)
    return Chem.MolToSmiles(topology, canonical=True)

def fragment_family_label(mol):
    fragments = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=False)
    labels = sorted(fragment_composition_label(fragment) for fragment in fragments)
    compact_labels = []
    for label in labels:
        if compact_labels and compact_labels[-1][0] == label:
            compact_labels[-1][1] += 1
        else:
            compact_labels.append([label, 1])
    return '.'.join(
        f'{count}x[{label}]' if count > 1 else f'1x[{label}]'
        for label, count in compact_labels
    )

def donor_degree_label(core_mol, non_ring_only=False):
    ring_atoms = {idx for ring in core_mol.GetRingInfo().AtomRings() for idx in ring}
    degrees = []
    for atom in core_mol.GetAtoms():
        if not atom.HasProp('_MetalDonor') or not atom.GetBoolProp('_MetalDonor'):
            continue
        if non_ring_only and atom.GetIdx() in ring_atoms:
            continue
        degrees.append(atom.GetIntProp('_LigandEffectiveDegree'))
    return '.'.join(f'N{degree}' for degree in sorted(degrees))

def hybrid_ring_site_label(core_mol):
    def restricted_distance(allowed_atoms, starts, target):
        queue = deque((start, 0) for start in starts)
        seen = set(starts)
        while queue:
            atom_idx, distance = queue.popleft()
            if atom_idx == target:
                return distance
            for neighbour in core_mol.GetAtomWithIdx(atom_idx).GetNeighbors():
                neighbour_idx = neighbour.GetIdx()
                if neighbour_idx in allowed_atoms and neighbour_idx not in seen:
                    seen.add(neighbour_idx)
                    queue.append((neighbour_idx, distance + 1))
        return -1

    rings = [set(ring) for ring in core_mol.GetRingInfo().AtomRings()]
    ring_systems = []
    seen = set()
    for start in range(len(rings)):
        if start in seen:
            continue
        stack = [start]
        atoms = set()
        while stack:
            ring_id = stack.pop()
            if ring_id in seen:
                continue
            seen.add(ring_id)
            atoms.update(rings[ring_id])
            for other_id, other_ring in enumerate(rings):
                if other_id not in seen and rings[ring_id] & other_ring:
                    stack.append(other_id)
        ring_systems.append(atoms)

    labels = []
    for atoms in ring_systems:
        donors = {
            idx for idx in atoms
            if core_mol.GetAtomWithIdx(idx).HasProp('_MetalDonor')
            and core_mol.GetAtomWithIdx(idx).GetBoolProp('_MetalDonor')
        }
        if not donors:
            continue
        attachments = []
        for atom_idx in atoms:
            atom = core_mol.GetAtomWithIdx(atom_idx)
            if any(
                neighbour.GetIdx() not in atoms and neighbour.GetAtomicNum() != 1
                for neighbour in atom.GetNeighbors()
            ):
                attachments.append(restricted_distance(atoms, donors, atom_idx))
        composition = fragment_composition_label(mol_fragment_from_atom_indices(core_mol, atoms))
        labels.append(
            f"{composition}:donors={len(donors)}:distance={','.join(map(str, sorted(attachments)))}"
        )
    return '|'.join(sorted(labels))

def core_family_label(core_mol):
    if core_mol is None or core_mol.GetNumAtoms() == 0:
        return 'non-ring coordination core'

    ring_atoms = {idx for ring in core_mol.GetRingInfo().AtomRings() for idx in ring}
    heavy_atoms = {atom.GetIdx() for atom in core_mol.GetAtoms() if atom.GetAtomicNum() != 1}

    if not ring_atoms:
        return f'acyclic[{topology_graph_label(core_mol)};degrees={donor_degree_label(core_mol)}]'

    if heavy_atoms - ring_atoms:
        ring_mol = mol_fragment_from_atom_indices(core_mol, ring_atoms)
        return (
            f'hybrid[scaffold={fragment_family_label(ring_mol)};'
            f'sites={hybrid_ring_site_label(core_mol)};'
            f'eff. connectivity={donor_degree_label(core_mol, non_ring_only=True)}]'
        )

    return fragment_family_label(core_mol)

def split_core_substituent(ligand, mols: bool=False, debug: int=0):
    lig_mol, smiles = get_rdkit_info(ligand)
    ligand_indices  = ligand.get_parent_indices('molecule')
    molecule        = ligand.get_parent('molecule')
    metal           = molecule.metals[0]
    metal_index     = metal.get_parent_index('molecule')
    molecule_graph  = molecule_graph_with_current_metal_connectivity(molecule, debug=debug)
    top_distances   = compute_topological_distances(molecule_graph, metal_index)
    lig_distances   = np.array([top_distances[idx] for idx in ligand_indices])
    donor_atoms     = get_metal_bound_ligand_atoms(ligand)

    ring_atoms = set(metal_connected_ring_core_atoms(lig_mol, donor_atoms))
    if ring_atoms and donor_atoms <= ring_atoms:
        core_heavy_atoms = ring_atoms
    else:
        backbone_atoms = set(donor_backbone_core_atoms(lig_mol, donor_atoms))
        if ring_atoms:
            backbone_atoms = {
                atom_idx for atom_idx in backbone_atoms
                if not lig_mol.GetAtomWithIdx(atom_idx).IsInRing()
            }
        core_heavy_atoms = ring_atoms | backbone_atoms
    core_h_atoms    = {atom.GetIdx() for atom in lig_mol.GetAtoms() if atom.GetAtomicNum() == 1 and any(nbr.GetIdx() in core_heavy_atoms for nbr in atom.GetNeighbors())}
    core_atoms      = np.array(sorted(core_heavy_atoms | core_h_atoms), dtype=int)
    subs_atoms      = np.array(sorted(set(range(lig_mol.GetNumAtoms())) - set(core_atoms)), dtype=int)

    if mols:
        core_mol = mol_fragment_from_atom_indices(lig_mol, core_atoms)
        core_index_map = {old_idx: new_idx for new_idx, old_idx in enumerate(core_atoms)}
        for atom_idx in donor_atoms & set(core_atoms):
            core_atom = core_mol.GetAtomWithIdx(core_index_map[atom_idx])
            ligand_atom = lig_mol.GetAtomWithIdx(atom_idx)
            core_atom.SetBoolProp('_MetalDonor', True)
            core_atom.SetIntProp(
                '_LigandEffectiveDegree',
                int(ligand_atom.GetTotalDegree() + max(0, -ligand_atom.GetFormalCharge())),
            )
        subs_mol = mol_fragment_from_atom_indices(lig_mol, subs_atoms)
        return core_mol, subs_mol
    else:
        return core_atoms, subs_atoms

########################################
#### Distance Weighted Fingerprints ####
########################################
def get_weights(ligand, alpha=0.25, debug: int=0):
    def distance_weight(distance, alpha=0.25):
        return np.round(float(np.exp(-alpha * (float(distance) - 1.0))), 4)

    ligand_indices  = ligand.get_parent_indices('molecule')
    molecule        = ligand.get_parent('molecule')
    metal           = molecule.metals[0]
    metal_index     = metal.get_parent_index('molecule')
    molecule_graph  = molecule_graph_with_current_metal_connectivity(molecule, debug=debug)

    top_distances   = compute_topological_distances(molecule_graph, metal_index)
    distances       = np.array([top_distances[idx] for idx in ligand_indices])
    weights         = np.array([distance_weight(distance, alpha=alpha) for distance in distances])
    return weights, distances

def atom_centered_weighted_morgan_fp(mol, atom_weights, centers=None, radius=2, n_bits=1024, normalize=True):
    if centers is None:
        centers = range(mol.GetNumAtoms())

    fp = np.zeros(n_bits, dtype=float)
    contribution_rows = []

    for center_idx in centers:
        center_idx = int(center_idx)
        center_weight = float(atom_weights[center_idx])
        bit_info = {}
        rdMolDescriptors.GetMorganFingerprintAsBitVect(
            mol,
            radius,
            nBits=n_bits,
            fromAtoms=[center_idx],
            bitInfo=bit_info,
        )

        for bit_id, environments in bit_info.items():
            contribution = center_weight * len(environments)
            fp[bit_id] += contribution
            contribution_rows.append({
                'center_atom_idx': center_idx,
                'center_symbol': mol.GetAtomWithIdx(center_idx).GetSymbol(),
                'center_weight': round(center_weight, 4),
                'bit_id': int(bit_id),
                'num_environments_for_bit': len(environments),
                'contribution': round(float(contribution), 4),
                'environment_radii': sorted({int(env_radius) for _, env_radius in environments}),
            })

    norm = np.linalg.norm(fp)
    if normalize and norm > 0:
        fp = fp / norm

    return fp, contribution_rows

def get_ligand_descriptor_type2(ligand, ligand_mol, alpha=0.25, radius=2, n_bits=1024, debug: int=0):
    weights, distances = get_weights(ligand, alpha=alpha, debug=debug)
    if len(weights) != ligand_mol.GetNumAtoms():
        raise ValueError(f'Ligand/RDKit atom mismatch: {len(weights)} weights for {ligand_mol.GetNumAtoms()} atoms')

    core_mol, _                               = split_core_substituent(ligand, mols=True, debug=debug)
    core_atoms, substituent_atoms             = split_core_substituent(ligand, mols=False, debug=debug)
    core_weights                              = weights[core_atoms]
    if metal_connected_ring_core_atoms(ligand_mol, get_metal_bound_ligand_atoms(ligand)):
        core_fp, core_contributions = atom_centered_weighted_morgan_fp(core_mol, core_weights, radius=radius, n_bits=n_bits)
    else:
        core_fp, core_contributions = atom_centered_weighted_morgan_fp(ligand_mol, weights, centers=core_atoms, radius=radius, n_bits=n_bits)
    substituent_fp, substituent_contributions = atom_centered_weighted_morgan_fp(ligand_mol, weights, centers=substituent_atoms, radius=radius, n_bits=n_bits)

    return {'core_fp': core_fp, 'substituent_fp': substituent_fp, 'core_contributions': core_contributions, 'substituent_contributions': substituent_contributions}

#################################
#### Descriptor Construction ####
#################################
## Main
def collect_ligands_data(systems, descriptor_type: int=2, alpha=0.25, radius=2, n_bits=1024):
    records  = []
    failures = []

    for system in systems:
        refcode = getattr(system, 'refcode_wo_digits', None) or getattr(system, 'refcode', None)
        try:
            hs_ref_mol             = get_hs_ref_mol(system)
            lig                    = hs_ref_mol.ligands[0]
            lig_mol, smiles        = get_rdkit_info(lig)
            core_mol, subs_mol     = split_core_substituent(lig, mols=True)
            core_atoms, subs_atoms = split_core_substituent(lig, mols=False)
            family                 = core_family_label(core_mol)
            is_aromatic            = int(any(atom.GetIsAromatic() for atom in lig_mol.GetAtoms()))
            denticity              = lig.get_denticity()

            if descriptor_type == 1:
                desc = get_ligand_descriptor_type1(core_mol, subs_mol, radius=radius, n_bits=n_bits)
            elif descriptor_type == 2:
                desc = get_ligand_descriptor_type2(lig, lig_mol, alpha=alpha, radius=radius, n_bits=n_bits)
            else:
                raise ValueError('descriptor_type must be 1 or 2')

            records.append({
                'id'            : refcode,
                'smiles'        : smiles,
                'denticity'     : denticity,
                'is_aromatic'   : is_aromatic,
                'family'        : family,
                'core_mol'      : core_mol,
                'subs_mol'      : subs_mol,
                'core_atoms'    : core_atoms,
                'subs_atoms'    : subs_atoms,
                'lig_mol'       : lig_mol,
                'core_fp'       : desc['core_fp'],
                'substituent_fp': desc['substituent_fp'],
            })
        except Exception as exc:
            failures.append((refcode, repr(exc)))

    return records, failures

def display_svg_mol(mol, legend=None):
    svg = Draw.MolsToGridImage([mol], molsPerRow=1, subImgSize=(520, 360), legends=[legend], useSVG=True)
    display(SVG(svg.data if hasattr(svg, 'data') else str(svg)))

def color_atom_labels(svg_text, atom_indices, color='#D62728'):
    for atom_idx in atom_indices:
        pattern = rf'(<path class="atom-{int(atom_idx)}"[^>]*fill=")#[0-9A-Fa-f]{{6,8}}'
        svg_text = re.sub(pattern, rf'\1{color}', svg_text)
    return svg_text

def display_svg_mol_extra(record, size=(620, 420)):
    highlight_atoms = sorted(set(record['core_atoms']) | set(record['subs_atoms']))
    highlight_colors = {idx: (0.65, 0.90, 0.70) for idx in record['core_atoms']}
    highlight_colors.update({idx: (0.70, 0.82, 1.00) for idx in record['subs_atoms']})
    donor_atoms = {
        int(record['core_atoms'][atom.GetIdx()])
        for atom in record['core_mol'].GetAtoms()
        if atom.HasProp('_MetalDonor') and atom.GetBoolProp('_MetalDonor')
    }
    svg = Draw.MolsToGridImage(
        [record['lig_mol']],
        molsPerRow=1,
        subImgSize=size,
        legends=[f"{record['id']} | green=core, blue=substituent, red=coordinating N labels"],
        highlightAtomLists=[highlight_atoms],
        highlightAtomColors=[highlight_colors],
        useSVG=True,
    )
    svg_text = svg.data if hasattr(svg, 'data') else str(svg)
    display(SVG(color_atom_labels(svg_text, donor_atoms)))

def get_hs_ref_mol(system):
    found, molecule = system.find_source('ref_hs_mol')
    if not found: raise ValueError(f'System {system.refcode} is missing ref_hs_mol')
    return molecule

def get_rdkit_info(ligand):
    if hasattr(ligand, 'rdkit_obj') and ligand.rdkit_obj is not None:
        rdkit_mol = Chem.Mol(ligand.rdkit_obj)
    elif hasattr(ligand, 'smiles') and ligand.smiles:
        rdkit_mol = Chem.MolFromSmiles(ligand.smiles)
    else:
        raise ValueError('Could not build an RDKit object for ligand')

    if hasattr(ligand, 'smiles') and ligand.smiles:
        smiles = ligand.smiles
    else:
        smiles = Chem.MolToSmiles(rdkit_mol)
    return rdkit_mol, smiles

def morgan_fp_array(mol, radius=2, n_bits=1024):
    from rdkit import DataStructs
    fp = rdMolDescriptors.GetMorganFingerprintAsBitVect(mol, radius, nBits=n_bits)
    arr = np.zeros((n_bits,), dtype=int)
    DataStructs.ConvertToNumpyArray(fp, arr)
    return arr

def get_ligand_descriptor_type1(core_mol, subs_mol, radius=2, n_bits=1024):
    core_fp = morgan_fp_array(core_mol, radius=radius, n_bits=n_bits)
    substituent_fp = morgan_fp_array(subs_mol, radius=0, n_bits=n_bits)
    return {'core_fp': core_fp, 'substituent_fp': substituent_fp}

####################
#### Map Layout ####
####################
def pairwise_distances(matrix, metric='euclidean'):
    matrix = np.asarray(matrix, dtype=float)
    if metric == 'euclidean':
        differences = matrix[:, None, :] - matrix[None, :, :]
        return np.linalg.norm(differences, axis=2)
    if metric == 'cosine':
        norms = np.linalg.norm(matrix, axis=1)
        denominator = np.outer(norms, norms)
        similarities = np.divide(
            matrix @ matrix.T,
            denominator,
            out=np.zeros_like(denominator),
            where=denominator != 0,
        )
        distances = 1.0 - np.clip(similarities, -1.0, 1.0)
        np.fill_diagonal(distances, 0.0)
        return distances
    raise ValueError("metric must be 'euclidean' or 'cosine'")

def metric_mds(dissimilarities, n_components=2, random_state=7, max_iter=300, eps=1e-6):
    dissimilarities = np.asarray(dissimilarities, dtype=float)
    if dissimilarities.ndim != 2 or dissimilarities.shape[0] != dissimilarities.shape[1]:
        raise ValueError('dissimilarities must be a square matrix')
    if not np.allclose(dissimilarities, dissimilarities.T):
        raise ValueError('dissimilarities must be symmetric')

    n_samples = len(dissimilarities)
    rng = np.random.RandomState(random_state)
    coords = rng.uniform(size=(n_samples, n_components))
    distances = pairwise_distances(coords)
    old_stress = None

    for iteration in range(max_iter):
        distances[distances == 0] = 1e-5
        ratio = dissimilarities / distances
        matrix_b = -ratio
        diagonal = np.arange(n_samples)
        matrix_b[diagonal, diagonal] += ratio.sum(axis=1)
        coords = matrix_b @ coords / n_samples

        distances = pairwise_distances(coords)
        stress = np.sum((distances - dissimilarities) ** 2) / 2
        if old_stress is not None:
            squared_distances = np.sum(distances ** 2) / 2
            if squared_distances and (old_stress - stress) / squared_distances < eps:
                break
        old_stress = stress

    return coords, stress

def pca_projection(matrix, n_components=2):
    matrix = np.asarray(matrix, dtype=float)
    centered = matrix - np.mean(matrix, axis=0)
    _, singular_values, components = np.linalg.svd(centered, full_matrices=False)
    n_available = min(n_components, len(components))
    coords = centered @ components[:n_available].T
    if n_available < n_components:
        coords = np.pad(coords, ((0, 0), (0, n_components - n_available)))
    model = {'components': components[:n_available], 'singular_values': singular_values[:n_available]}
    return coords, model

def project_map(matrix, method='mds', metric='cosine'):
    matrix = np.asarray(matrix, dtype=float)
    if len(matrix) == 0:
        return np.empty((0, 2)), None
    if len(matrix) == 1 or np.all(np.var(matrix, axis=0) == 0):
        return np.zeros((len(matrix), 2)), None
    if len(matrix) == 2:
        distance = pairwise_distances(matrix, metric=metric)[0, 1]
        distance = 1.0 if distance == 0 else distance
        return np.array([[-0.5 * distance, 0.0], [0.5 * distance, 0.0]]), None
    if method == 'mds':
        distances = pairwise_distances(matrix, metric=metric)
        coords, stress = metric_mds(distances, n_components=2, random_state=7)
        return coords, {'stress': stress}
    if method == 'pca':
        return pca_projection(matrix, n_components=2)
    raise ValueError("embedding_method must be 'mds' or 'pca'")

def spread_points(coords, min_distance=None, strength=0.35, iterations=450, anchor_strength=0.025):
    coords = np.asarray(coords, dtype=float)
    if len(coords) <= 1:
        return coords.copy()

    spread = coords.copy()
    anchor = coords.copy()
    span = max(np.ptp(coords[:, 0]), np.ptp(coords[:, 1]), 1.0)
    if min_distance is None:
        min_distance = 0.12 * span

    rng = np.random.default_rng(7)
    spread += rng.normal(scale=1e-4 * span, size=spread.shape)

    for _ in range(iterations):
        forces = np.zeros_like(spread)
        for i in range(len(spread)):
            delta = spread[i] - spread
            distances = np.linalg.norm(delta, axis=1)
            crowded = (distances > 1e-12) & (distances < min_distance)
            if np.any(crowded):
                unit = delta[crowded] / distances[crowded, None]
                overlap = (min_distance - distances[crowded])[:, None]
                forces[i] += np.sum(unit * overlap, axis=0)
        forces -= anchor_strength * (spread - anchor)
        spread += strength * forces / max(len(spread), 1)

    return spread

def assign_grid_centers(family_rows, raw_centers, grid_cell_spacing=1.0):
    for row, center in zip(family_rows, raw_centers):
        row['raw_family_x'] = center[0]
        row['raw_family_y'] = center[1]

    ordered_families = sorted(family_rows, key=lambda row: (row['raw_family_y'], row['raw_family_x']))
    n_families = len(ordered_families)
    n_cols = int(np.ceil(np.sqrt(n_families * 1.35)))
    n_rows = int(np.ceil(n_families / n_cols))

    grid_sites = []
    for row_idx in range(n_rows):
        for col_idx in range(n_cols):
            x = (col_idx - (n_cols - 1) / 2) * grid_cell_spacing
            y = - (row_idx - (n_rows - 1) / 2) * grid_cell_spacing
            grid_sites.append((x, y, row_idx + 1, col_idx + 1))

    for family_row, (x, y, row_idx, col_idx) in zip(ordered_families, grid_sites):
        family_row['family_x'] = x
        family_row['family_y'] = y
        family_row['grid_row'] = row_idx
        family_row['grid_column'] = col_idx

    return family_rows

def assign_position_family_labels(family_rows, prefix='F'):
    ordered_rows = sorted(
        family_rows,
        key=lambda row: (-float(row['family_y']), float(row['family_x']), row['Core-ring family']),
    )
    for idx, row in enumerate(ordered_rows, start=1):
        row['family_position_index'] = idx
        if row.get('grid_row') is not None and row.get('grid_column') is not None:
            row['family_id'] = f"{prefix}({int(row['grid_row'])},{int(row['grid_column'])})"
        else:
            row['family_id'] = f'{prefix}{idx}'
    return family_rows

def build_core_ring_map_data(
    descriptor_records,
    target_denticity,
    data=None,
    embedding_method='mds',
    grid=False,
    spread_family_centers=True,
    family_min_distance=0.12,
    grid_cell_spacing=1.0,
    grid_member_radius=0.28,
    local_member_min_distance=0.08,
):
    data = data or {}
    selected_records = [record for record in descriptor_records if str(record['denticity']) == str(target_denticity)]
    if not selected_records:
        raise ValueError(f'No homoleptic systems found for denticity {target_denticity}')

    family_to_records = {}
    for record in selected_records:
        family_to_records.setdefault(record['family'], []).append(record)

    family_rows = []
    for family, records in family_to_records.items():
        branch_vectors = np.asarray([record['core_fp'] for record in records], dtype=float)
        if not np.any(branch_vectors):
            branch_vectors = np.asarray([record['substituent_fp'] for record in records], dtype=float)
        family_rows.append({
            'Core-ring family': family,
            'family_count': len(records),
            'branch_vector': np.mean(branch_vectors, axis=0),
        })

    family_rows = sorted(family_rows, key=lambda row: (-row['family_count'], row['Core-ring family']))
    raw_centers, embedding = project_map([row['branch_vector'] for row in family_rows], method=embedding_method)

    if grid:
        family_rows = assign_grid_centers(family_rows, raw_centers, grid_cell_spacing=grid_cell_spacing)
        local_scale = grid_member_radius
    else:
        family_centers = spread_points(
            raw_centers,
            min_distance=family_min_distance * max(np.ptp(raw_centers[:, 0]), np.ptp(raw_centers[:, 1]), 1.0),
        ) if spread_family_centers else raw_centers
        for row, raw_center, center in zip(family_rows, raw_centers, family_centers):
            row['raw_family_x'] = raw_center[0]
            row['raw_family_y'] = raw_center[1]
            row['family_x'] = center[0]
            row['family_y'] = center[1]
            row['grid_row'] = None
            row['grid_column'] = None

        if len(family_centers) > 1:
            center_distances = pairwise_distances(family_centers)
            nonzero_distances = center_distances[center_distances > 1e-12]
            local_scale = 0.18 * np.median(nonzero_distances) if len(nonzero_distances) else 0.35
        else:
            local_scale = 0.35

    family_rows = assign_position_family_labels(family_rows)
    family_rows = sorted(family_rows, key=lambda row: row['family_position_index'])

    map_rows = []
    for family_row in family_rows:
        family = family_row['Core-ring family']
        records = family_to_records[family]
        local_coords, _ = project_map([record['substituent_fp'] for record in records], method=embedding_method)
        max_radius = np.max(np.linalg.norm(local_coords, axis=1)) if len(local_coords) else 0.0
        if max_radius > 0:
            scale = local_scale if grid else local_scale * min(1.4, 0.75 + 0.07 * np.sqrt(len(records)))
            local_coords = local_coords / max_radius * scale
        else:
            scale = local_scale if grid else local_scale

        if local_member_min_distance and len(local_coords) > 1:
            local_coords = spread_points(
                local_coords,
                min_distance=local_member_min_distance * local_scale,
                strength=0.25,
                iterations=250,
                anchor_strength=0.08,
            )
            max_radius = np.max(np.linalg.norm(local_coords, axis=1)) if len(local_coords) else 0.0
            max_allowed_radius = max(scale * 1.15, local_member_min_distance * local_scale)
            if max_radius > max_allowed_radius:
                local_coords = local_coords / max_radius * max_allowed_radius

        for record, local_coord in zip(records, local_coords):
            row = {
                'Refcode': record['id'],
                'family_id': family_row['family_id'],
                'map_x': family_row['family_x'] + local_coord[0],
                'map_y': family_row['family_y'] + local_coord[1],
                'family_x': family_row['family_x'],
                'family_y': family_row['family_y'],
                'raw_family_x': family_row['raw_family_x'],
                'raw_family_y': family_row['raw_family_y'],
                'local_x': local_coord[0],
                'local_y': local_coord[1],
                'denticity': record['denticity'],
                'is_aromatic': bool(record['is_aromatic']),
                'family_position_index': family_row['family_position_index'],
                'grid_row': family_row.get('grid_row'),
                'grid_column': family_row.get('grid_column'),
                'Core-ring family': family,
                'family_count': family_row['family_count'],
                'core_atom_count': len(record['core_atoms']),
                'substituent_atom_count': len(record['subs_atoms']),
                'SMILES': record['smiles'],
            }
            row.update(data.get(record['id'], {}))
            map_rows.append(row)

    map_rows.sort(key=lambda row: (row['family_position_index'], row['Refcode']))
    return map_rows, family_rows, selected_records, embedding

##########################
#### Interactive HTML ####
##########################
def format_value(value, digits=2, unit=None):
    if value in (None, '', 'None') or (isinstance(value, (float, np.floating)) and np.isnan(value)):
        return '-'
    try:
        formatted = f'{float(value):.{digits}f}'
    except (TypeError, ValueError):
        formatted = str(value)
    return formatted if unit in (None, '') else f'{formatted} {unit}'

def mol_to_svg(mol, size=(260, 180)):
    if mol is None:
        return '<div class="missing-mol">Molecule unavailable</div>'
    svg = Draw.MolsToGridImage([mol], molsPerRow=1, subImgSize=size, useSVG=True)
    return svg.data if hasattr(svg, 'data') else str(svg)

def normalize_property_specs(properties_to_show):
    specs = []
    for label, spec in properties_to_show.items():
        if isinstance(spec, str):
            specs.append({'label': label, 'key': spec, 'unit': None, 'digits': 2})
        else:
            specs.append({
                'label': label,
                'key': spec.get('key', label),
                'unit': spec.get('unit'),
                'digits': spec.get('digits', 2),
            })
    return specs


def resolve_property_key(property_name, property_specs, available_columns):
    available_columns = list(available_columns)
    if property_name in available_columns:
        return property_name

    candidates = [property_name]
    for spec in property_specs:
        if property_name in (spec['label'], spec['key']):
            candidates.extend([spec['key'], spec['label']])
            if spec.get('unit'):
                candidates.extend([
                    f"{spec['key']} [{spec['unit']}]",
                    f"{spec['label']} [{spec['unit']}]",
                ])

    for candidate in candidates:
        if candidate in available_columns:
            return candidate

    def compact(text):
        return ''.join(char.lower() for char in str(text) if char.isalnum())

    compact_candidates = [compact(candidate) for candidate in candidates]
    for column in available_columns:
        compact_column = compact(column)
        if compact_column in compact_candidates or any(compact_column.startswith(candidate) for candidate in compact_candidates):
            return column

    return property_name

def resolve_property_specs(property_specs, available_columns):
    resolved_specs = []
    for spec in property_specs:
        resolved = dict(spec)
        resolved['key'] = resolve_property_key(spec['key'], property_specs, available_columns)
        resolved_specs.append(resolved)
    return resolved_specs

def is_numeric_series(values):
    present = [value for value in values if value not in (None, '', 'None')]
    if not present:
        return False
    try:
        [float(value) for value in present]
    except (TypeError, ValueError):
        return False
    return True


def numeric_colorbar_ticks(values, n_ticks=5, include_value=None, use_unique_when_few=False, calendar_years=False):
    values = np.asarray([
        float(value) for value in values
        if value not in (None, '', 'None') and not (isinstance(value, (float, np.floating)) and np.isnan(value))
    ], dtype=float)
    if not len(values):
        return []

    n_ticks = max(2, int(n_ticks or 5))
    unique_values = np.unique(values)
    integer_like = bool(np.all(np.isclose(unique_values, np.round(unique_values))))

    if calendar_years and unique_values.max() > unique_values.min():
        raw_step = (unique_values.max() - unique_values.min()) / max(n_ticks - 1, 1)
        magnitude = 10 ** np.floor(np.log10(max(raw_step, 1)))
        candidate_steps = magnitude * np.array([1, 2, 5, 10], dtype=float)
        step = candidate_steps[np.argmin(np.abs(candidate_steps - raw_step))]
        first_tick = np.floor(unique_values.min() / step) * step
        last_tick = np.ceil(unique_values.max() / step) * step
        tick_values = np.arange(first_tick, last_tick + 0.5 * step, step)
    elif use_unique_when_few and len(unique_values) <= n_ticks:
        tick_values = unique_values
    else:
        locator = MaxNLocator(nbins=max(n_ticks - 1, 1), integer=integer_like, min_n_ticks=2)
        tick_values = np.asarray(locator.tick_values(float(unique_values.min()), float(unique_values.max())), dtype=float)

    vmin = float(unique_values.min())
    vmax = float(unique_values.max())
    if include_value is not None and vmin <= float(include_value) <= vmax:
        tick_values = np.unique(np.append(tick_values, float(include_value)))
        tick_values = np.array(sorted(tick_values), dtype=float)

    def tick_label(value):
        if np.isclose(value, round(value)):
            return str(int(round(value)))
        return f'{value:.2f}'.rstrip('0').rstrip('.')

    return [
        {
            'value': tick_label(value),
            'numeric_value': float(value),
            'position': float((float(value) - vmin) / (vmax - vmin)) if vmax != vmin else 0.5,
        }
        for value in tick_values
    ]

def compact_property_name(color_property):
    return ''.join(char.lower() for char in str(color_property) if char.isalnum())

def is_t12_property(color_property):
    compact_name = compact_property_name(color_property)
    return compact_name in {'t12', 't12k', 't1over2', 't1slash2'}

def is_dhelec_property(color_property):
    return compact_property_name(color_property).startswith('dhelec')

def is_year_property(color_property):
    return compact_property_name(color_property) in {'year', 'publicationyear'}

def numeric_color_limits(color_property, vmin, vmax):
    if is_dhelec_property(color_property):
        return vmin, min(vmax, 80.0)
    return vmin, vmax

def numeric_colormap_name(color_property):
    return 'blue_white_red' if is_t12_property(color_property) or is_dhelec_property(color_property) else 'turbo'

def blue_white_red_colormap(center_position=0.5):
    center_position = float(np.clip(center_position, 0.0, 1.0))
    transition = min(0.04, center_position, 1.0 - center_position)
    return LinearSegmentedColormap.from_list(
        'blue_white_red',
        [
            (0.00, '#08306b'),  # dark blue, low T1/2
            (center_position - transition, '#deebf7'),
            (center_position, '#f7f7f7'),
            (center_position + transition, '#fee0d2'),
            (1.00, '#67000d'),  # dark red, high T1/2
        ],
    )

def numeric_colormap_and_norm(color_property, vmin, vmax):
    if is_t12_property(color_property) or is_dhelec_property(color_property):
        center = 298.0 if is_t12_property(color_property) else 0.0
        if vmin < center < vmax:
            center_position = (center - vmin) / (vmax - vmin)
            norm = Normalize(vmin=vmin, vmax=vmax)
            return blue_white_red_colormap(center_position), norm
        if vmin >= center:
            cmap = LinearSegmentedColormap.from_list('white_red', ['#f7f7f7', '#67000d'])
        else:
            cmap = LinearSegmentedColormap.from_list('blue_white', ['#08306b', '#f7f7f7'])
        norm = Normalize(vmin=vmin, vmax=vmax) if vmax != vmin else Normalize(vmin=0, vmax=1)
        return cmap, norm

    cmap = colormaps.get_cmap('turbo')
    norm = Normalize(vmin=vmin, vmax=vmax) if vmax != vmin else Normalize(vmin=0, vmax=1)
    return cmap, norm

def add_color_column(map_rows, color_property, top_n_categories=20, cmap_name='tab20'):
    map_rows = [dict(row) for row in map_rows]
    values = [row.get(color_property, 'unknown') for row in map_rows]

    if is_numeric_series(values):
        numeric = np.asarray([
            float(value) if value not in (None, '', 'None') else np.nan
            for value in values
        ], dtype=float)
        valid_numeric = numeric[~np.isnan(numeric)]
        cmap_name_numeric = numeric_colormap_name(color_property)
        if not len(valid_numeric):
            for row in map_rows:
                row['_color_group'] = 'missing'
                row['_point_color'] = '#d9d9d9'
            legend_items = [{'kind': 'category', 'label': 'missing', 'color': '#d9d9d9'}]
            return map_rows, legend_items

        vmin = float(valid_numeric.min())
        vmax = float(valid_numeric.max())
        color_vmin, color_vmax = numeric_color_limits(color_property, vmin, vmax)
        tick_payload = numeric_colorbar_ticks(
            [color_vmin, color_vmax],
            n_ticks=top_n_categories,
            calendar_years=is_year_property(color_property),
        )
        if is_dhelec_property(color_property):
            tick_payload = [tick for tick in tick_payload if tick['numeric_value'] <= color_vmax]
        if tick_payload:
            color_vmin = min(color_vmin, tick_payload[0]['numeric_value'])
            color_vmax = max(color_vmax, tick_payload[-1]['numeric_value'])
        cmap, norm = numeric_colormap_and_norm(color_property, color_vmin, color_vmax)
        missing_color = '#000000' if is_t12_property(color_property) else '#d9d9d9'
        for row, value, numeric_value in zip(map_rows, values, numeric):
            row['_color_group'] = str(value)
            row['_point_color'] = (
                to_hex(cmap(norm(np.clip(numeric_value, color_vmin, color_vmax))))
                if not np.isnan(numeric_value) else missing_color
            )
        for tick in tick_payload:
            tick['position'] = float(norm(tick['numeric_value'])) if color_vmax != color_vmin else 0.5
        legend_items = [{
            'kind': 'colorbar',
            'label': color_property,
            'cmap_name': cmap_name_numeric,
            'vmin': format_value(color_vmin),
            'vmax': format_value(color_vmax),
            'color_vmin': color_vmin,
            'color_vmax': color_vmax,
            'ticks': tick_payload,
            'colors': [to_hex(cmap(x)) for x in np.linspace(0, 1, 24)],
        }]
        return map_rows, legend_items

    clean_values = ['unknown' if value in (None, '', 'None') else str(value) for value in values]
    counts = Counter(clean_values)
    top_values = [value for value, _ in counts.most_common(top_n_categories)]
    groups = [value if value in top_values else 'Other' for value in clean_values]
    palette = colormaps.get_cmap(cmap_name)(np.linspace(0, 1, max(len(top_values), 1)))
    colors = {value: to_hex(palette[i % len(palette)]) for i, value in enumerate(top_values)}
    colors['Other'] = '#d9d9d9'
    for row, group in zip(map_rows, groups):
        row['_color_group'] = group
        row['_point_color'] = colors[group]
    order = top_values + (['Other'] if 'Other' in set(groups) else [])
    legend_items = [
        {'kind': 'category', 'label': f'{value} (n={groups.count(value)})', 'color': colors[value]}
        for value in order
    ]
    return map_rows, legend_items

def write_interactive_html(
    map_rows,
    family_rows,
    descriptor_records,
    output_html_path,
    properties_to_show,
    color_property='last_author',
    color_properties=None,
    title=None,
    caption=None,
    top_n_categories=20,
    plot_width=980,
    plot_height=680,
    iframe_height=760,
):
    output_html_path = Path(output_html_path)
    output_html_path.parent.mkdir(parents=True, exist_ok=True)
    available_columns = set().union(*(row.keys() for row in map_rows))
    property_specs = resolve_property_specs(normalize_property_specs(properties_to_show), available_columns)
    color_properties = list(dict.fromkeys(color_properties or [color_property]))
    if color_property not in color_properties:
        color_properties.insert(0, color_property)

    point_colors = {}
    color_configs = {}
    for property_name in color_properties:
        color_key = resolve_property_key(property_name, property_specs, available_columns)
        colored_rows, legend_items = add_color_column(
            map_rows,
            color_key,
            top_n_categories=top_n_categories,
        )
        point_colors[property_name] = [row['_point_color'] for row in colored_rows]

        property_spec = next(
            (spec for spec in property_specs if spec['key'] == color_key),
            {'label': property_name, 'key': color_key, 'unit': None, 'digits': 2},
        )
        legend_payload = []
        for idx, item in enumerate(legend_items):
            payload = dict(item)
            payload['x'] = plot_width - 230 + 30
            payload['y'] = 18 + 25 + idx * 22
            legend_payload.append(payload)

        color_configs[property_name] = {
            'color_key': color_key,
            'plain_label': plain_property_label(property_name),
            'label': html_property_label_payload(property_name),
            'legend': legend_payload,
            'unit': property_spec.get('unit'),
            'digits': property_spec.get('digits', 2),
        }

    margin = {'left': 42, 'right': 230, 'top': 18, 'bottom': 38}
    x_values = np.asarray([row['map_x'] for row in map_rows] + [row['family_x'] for row in family_rows])
    y_values = np.asarray([row['map_y'] for row in map_rows] + [row['family_y'] for row in family_rows])
    x_pad = 0.025 * (x_values.max() - x_values.min() or 1.0)
    y_pad = 0.025 * (y_values.max() - y_values.min() or 1.0)
    x_min, x_max = x_values.min() - x_pad, x_values.max() + x_pad
    y_min, y_max = y_values.min() - y_pad, y_values.max() + y_pad

    def scale_x(x):
        return margin['left'] + (x - x_min) / (x_max - x_min) * (plot_width - margin['left'] - margin['right'])

    def scale_y(y):
        return plot_height - margin['bottom'] - (y - y_min) / (y_max - y_min) * (plot_height - margin['top'] - margin['bottom'])

    ligand_mols_by_refcode = {record['id']: record.get('lig_mol', record.get('mol')) for record in descriptor_records}

    family_payload = []
    for row in family_rows:
        family_payload.append({
            'family_id': row['family_id'],
            'family': row['Core-ring family'],
            'family_count': int(row['family_count']),
            'x': float(scale_x(row['family_x'])),
            'y': float(scale_y(row['family_y'])),
        })

    point_payload = []
    for row_idx, row in enumerate(map_rows):
        mol_svg = mol_to_svg(ligand_mols_by_refcode.get(row['Refcode']))
        property_rows = []
        for spec in property_specs:
            value = row.get(spec['key'], None)
            property_rows.append(
                f"<tr><td>{escape(spec['label'])}</td><td>{escape(format_value(value, digits=spec['digits'], unit=spec['unit']))}</td></tr>"
            )
        property_html = '\n'.join(property_rows)
        tooltip_html = f"""
            <div class="tooltip-title">{escape(row['Refcode'])}</div>
            <div class="tooltip-subtitle" data-color-summary></div>
            <div class="mol-svg">{mol_svg}</div>
            <table>
                <tr><td>Family</td><td>{escape(row['family_id'])} ({escape(row['Core-ring family'])})</td></tr>
                <tr><td>Core atoms</td><td>{int(row['core_atom_count'])}</td></tr>
                {property_html}
            </table>
        """
        point_payload.append({
            'refcode': row['Refcode'],
            'family_id': row['family_id'],
            'x': float(scale_x(row['map_x'])),
            'y': float(scale_y(row['map_y'])),
            'family_x': float(scale_x(row['family_x'])),
            'family_y': float(scale_y(row['family_y'])),
            'colors': {
                property_name: point_colors[property_name][row_idx]
                for property_name in color_properties
            },
            'color_values': {
                property_name: format_value(
                    row.get(config['color_key']),
                    digits=config['digits'],
                    unit=config['unit'],
                )
                for property_name, config in color_configs.items()
            },
            'tooltip_html': tooltip_html,
        })

    title = title or 'Interactive core-ring ligand map'
    caption = caption or 'Grey radial lines connect each compound to its family center. Hover points to inspect ligand and metadata.'
    option_html = '\n'.join(
        f'<option value="{escape(property_name)}"{" selected" if property_name == color_property else ""}>{escape(config["plain_label"])}</option>'
        for property_name, config in color_configs.items()
    )

    html_doc = f"""<!DOCTYPE html>
<html>
<head>
<meta charset="utf-8">
<title>{escape(title)}</title>
<style>
body {{ margin: 0; font-family: Helvetica, Arial, sans-serif; background: #ffffff; color: #111827; }}
.wrapper {{ position: relative; width: {plot_width}px; margin: 0 auto; }}
.title {{ font-size: 18px; font-weight: 700; margin: 14px 0 2px 0; }}
.caption {{ font-size: 12px; color: #4b5563; margin-bottom: 8px; max-width: 820px; }}
.controls {{ display: flex; align-items: center; gap: 8px; margin: 8px 0; font-size: 12px; font-weight: 700; }}
.controls select {{ font: inherit; font-weight: 400; padding: 3px 24px 3px 7px; border: 1px solid #9ca3af; border-radius: 4px; background: white; }}
svg {{ border: 1.4px solid black; background: white; }}
.branch-line {{ stroke: #c7c7c7; stroke-width: 0.8; pointer-events: none; }}
.family-center {{ fill: #111827; stroke: white; stroke-width: 1.2; opacity: 0.85; }}
.point {{ cursor: pointer; stroke: black; stroke-width: 0.45; opacity: 0.92; }}
.point:hover {{ stroke-width: 1.8; opacity: 1; }}
.axis-label {{ font-size: 13px; fill: #111827; }}
.legend-title {{ font-size: 12px; font-weight: 700; }}
.legend-unit {{ font-weight: 400; }}
.legend-text {{ font-size: 11px; dominant-baseline: middle; }}
.tooltip {{ position: absolute; display: none; width: 330px; max-height: 560px; overflow: auto; background: rgba(255,255,255,0.98); border: 1px solid #111827; box-shadow: 0 10px 30px rgba(0,0,0,0.18); border-radius: 8px; padding: 10px 12px; pointer-events: none; z-index: 10; }}
.tooltip-title {{ font-size: 15px; font-weight: 800; margin-bottom: 2px; }}
.tooltip-subtitle {{ font-size: 12px; color: #374151; margin-bottom: 6px; }}
.mol-svg svg {{ border: 0; display: block; margin: 0 auto 4px auto; }}
.missing-mol {{ height: 80px; display: flex; align-items: center; justify-content: center; color: #6b7280; }}
table {{ border-collapse: collapse; width: 100%; font-size: 11px; }}
td {{ padding: 2px 4px; border-top: 1px solid #e5e7eb; vertical-align: top; }}
td:first-child {{ color: #4b5563; width: 35%; }}
</style>
</head>
<body>
<div class="wrapper">
  <div class="title">{escape(title)}</div>
  <div class="caption">{escape(caption)}</div>
  <div class="controls"><label for="color-property">Color by</label><select id="color-property">{option_html}</select></div>
  <svg id="plot" width="{plot_width}" height="{plot_height}" viewBox="0 0 {plot_width} {plot_height}"></svg>
  <div id="tooltip" class="tooltip"></div>
</div>
<script>
const points = {json.dumps(point_payload)};
const families = {json.dumps(family_payload)};
const colorConfigs = {json.dumps(color_configs)};
const svg = document.getElementById('plot');
const tooltip = document.getElementById('tooltip');
const wrapper = document.querySelector('.wrapper');
const colorSelector = document.getElementById('color-property');
const ns = 'http://www.w3.org/2000/svg';
function add(tag, attrs, text) {{
  const element = document.createElementNS(ns, tag);
  for (const [key, value] of Object.entries(attrs || {{}})) element.setAttribute(key, value);
  if (text !== undefined) element.textContent = text;
  svg.appendChild(element);
  return element;
}}
for (const point of points) {{ add('line', {{x1: point.family_x, y1: point.family_y, x2: point.x, y2: point.y, class: 'branch-line'}}); }}
for (const family of families) {{
  const center = add('circle', {{cx: family.x, cy: family.y, r: 4.4, class: 'family-center'}});
  const centerTitle = document.createElementNS(ns, 'title');
  centerTitle.textContent = family.family_id + ' | n=' + family.family_count + ' | ' + family.family;
  center.appendChild(centerTitle);
}}
let activeProperty = colorSelector.value;
const pointElements = [];
function positionTooltip(event) {{
  const wrapperRect = wrapper.getBoundingClientRect();
  const localX = event.clientX - wrapperRect.left;
  const localY = event.clientY - wrapperRect.top;
  const tooltipWidth = tooltip.offsetWidth || 330;
  const tooltipHeight = tooltip.offsetHeight || 360;
  let left = localX + 18;
  let top = localY + 18;
  if (left + tooltipWidth > wrapperRect.width - 8) left = localX - tooltipWidth - 18;
  if (top + tooltipHeight > wrapperRect.height - 8) top = localY - tooltipHeight - 18;
  tooltip.style.left = Math.max(8, left) + 'px';
  tooltip.style.top = Math.max(8, top) + 'px';
}}
for (const point of points) {{
  const circle = add('circle', {{cx: point.x, cy: point.y, r: 5.2, fill: point.colors[activeProperty], class: 'point'}});
  pointElements.push(circle);
  circle.addEventListener('mousemove', (event) => {{
    tooltip.style.display = 'block';
    tooltip.innerHTML = point.tooltip_html;
    const summary = tooltip.querySelector('[data-color-summary]');
    summary.textContent = colorConfigs[activeProperty].plain_label + ': ' + point.color_values[activeProperty];
    requestAnimationFrame(() => positionTooltip(event));
  }});
  circle.addEventListener('mouseleave', () => {{ tooltip.style.display = 'none'; }});
}}
<!-- add('text', {{x: {margin['left'] + (plot_width - margin['left'] - margin['right']) / 2}, y: {plot_height - 16}, class: 'axis-label', 'text-anchor': 'middle'}}, 'Map coordinate 1') -->;
<!-- add('text', {{x: 18, y: {margin['top'] + (plot_height - margin['top'] - margin['bottom']) / 2}, class: 'axis-label', transform: 'rotate(-90 18 {margin['top'] + (plot_height - margin['top'] - margin['bottom']) / 2})', 'text-anchor': 'middle'}}, 'Map coordinate 2') -->;
function renderLegendTitle(legendTitle) {{
  const legendTitleElement = add('text', {{x: {plot_width - margin['right'] + 30}, y: {margin['top']}, class: 'legend-title dynamic-legend'}});
  if (legendTitle.kind === 't12') {{
    legendTitleElement.appendChild(document.createTextNode('T'));
    const sub = document.createElementNS(ns, 'tspan');
    sub.setAttribute('baseline-shift', 'sub');
    sub.setAttribute('font-size', '75%');
    sub.textContent = '1/2';
    legendTitleElement.appendChild(sub);
    const unit = document.createElementNS(ns, 'tspan');
    unit.setAttribute('class', 'legend-unit');
    unit.textContent = ' [K]';
    legendTitleElement.appendChild(unit);
  }} else if (legendTitle.kind === 'dhelec') {{
    legendTitleElement.appendChild(document.createTextNode('ΔH'));
    const sub = document.createElementNS(ns, 'tspan');
    sub.setAttribute('baseline-shift', 'sub');
    sub.setAttribute('font-size', '75%');
    sub.textContent = 'elec';
    legendTitleElement.appendChild(sub);
    const unit = document.createElementNS(ns, 'tspan');
    unit.setAttribute('class', 'legend-unit');
    unit.textContent = ' [kJ/mol]';
    legendTitleElement.appendChild(unit);
  }} else {{
    legendTitleElement.textContent = legendTitle.text;
  }}
}}
function renderLegend(propertyName) {{
  svg.querySelectorAll('.dynamic-legend').forEach((element) => element.remove());
  const config = colorConfigs[propertyName];
  renderLegendTitle(config.label);
  for (const item of config.legend) {{
  if (item.kind === 'colorbar') {{
    const defs = add('defs', {{class: 'dynamic-legend'}});
    const gradientId = 'colorbar-gradient-' + Math.random().toString(36).slice(2);
    const gradient = document.createElementNS(ns, 'linearGradient');
    const barWidth = 22;
    const barHeight = 260;
    const barX = item.x;
    const barY = item.y + 8;
    gradient.setAttribute('id', gradientId);
    gradient.setAttribute('x1', '0%');
    gradient.setAttribute('x2', '0%');
    gradient.setAttribute('y1', '100%');
    gradient.setAttribute('y2', '0%');
    item.colors.forEach((color, idx) => {{
      const stop = document.createElementNS(ns, 'stop');
      stop.setAttribute('offset', (100 * idx / Math.max(item.colors.length - 1, 1)) + '%');
      stop.setAttribute('stop-color', color);
      gradient.appendChild(stop);
    }});
    defs.appendChild(gradient);
    add('rect', {{x: barX, y: barY, width: barWidth, height: barHeight, fill: 'url(#' + gradientId + ')', stroke: 'black', 'stroke-width': 0.6, class: 'dynamic-legend'}});
    for (const tick of item.ticks) {{
      const tickY = barY + barHeight * (1 - tick.position);
      add('line', {{x1: barX + barWidth, y1: tickY, x2: barX + barWidth + 7, y2: tickY, stroke: 'black', 'stroke-width': 0.7, class: 'dynamic-legend'}});
      add('text', {{x: barX + barWidth + 11, y: tickY + 3.5, class: 'legend-text dynamic-legend'}}, tick.value);
    }}
  }} else {{
    add('circle', {{cx: item.x, cy: item.y, r: 5, fill: item.color, stroke: 'black', 'stroke-width': 0.4, class: 'dynamic-legend'}});
    add('text', {{x: item.x + 11, y: item.y, class: 'legend-text dynamic-legend'}}, item.label);
  }}
  }}
}}
function updateColorProperty() {{
  activeProperty = colorSelector.value;
  pointElements.forEach((pointElement, idx) => pointElement.setAttribute('fill', points[idx].colors[activeProperty]));
  tooltip.style.display = 'none';
  renderLegend(activeProperty);
}}
colorSelector.addEventListener('change', updateColorProperty);
updateColorProperty();
</script>
</body>
</html>"""
    output_html_path.write_text(html_doc)
    return output_html_path




def pretty_property_label(label):
    if is_t12_property(label):
        return r'$\mathbf{T}_{\mathbf{1/2}}$ [K]'
    if is_dhelec_property(label):
        return r'$\Delta H_{\mathrm{elec}}$ [kJ/mol]'
    return str(label).replace('_', ' ').strip().title()

def plain_property_label(label):
    if is_t12_property(label):
        return 'T1/2'
    if is_dhelec_property(label):
        return 'ΔHelec'
    return str(label).replace('_', ' ').strip().title()

def html_property_label_payload(label):
    if is_t12_property(label):
        return {'kind': 't12'}
    if is_dhelec_property(label):
        return {'kind': 'dhelec'}
    return {'kind': 'text', 'text': pretty_property_label(label)}

def write_static_png(
    map_rows,
    family_rows,
    output_png_path,
    properties_to_show,
    color_property='last_author',
    top_n_categories=20,
    figsize=(7.9, 6.6),
    dpi=400,
):
    output_png_path = Path(output_png_path)
    output_png_path.parent.mkdir(parents=True, exist_ok=True)

    available_columns = set().union(*(row.keys() for row in map_rows))
    property_specs = resolve_property_specs(normalize_property_specs(properties_to_show), available_columns)
    color_key = resolve_property_key(color_property, property_specs, available_columns)
    plot_rows, legend_items = add_color_column(map_rows, color_key, top_n_categories=top_n_categories)

    fig, ax = plt.subplots(figsize=figsize)
    fig.subplots_adjust(left=0.025, right=0.825, top=0.995, bottom=0.005)

    for row in plot_rows:
        ax.plot(
            [row['family_x'], row['map_x']],
            [row['family_y'], row['map_y']],
            color='#c7c7c7',
            linewidth=0.55,
            zorder=1,
        )

    if family_rows:
        ax.scatter(
            [row['family_x'] for row in family_rows],
            [row['family_y'] for row in family_rows],
            s=20,
            color='#111827',
            edgecolor='white',
            linewidth=0.5,
            zorder=2,
        )

    ax.scatter(
        [row['map_x'] for row in plot_rows],
        [row['map_y'] for row in plot_rows],
        s=38,
        c=[row['_point_color'] for row in plot_rows],
        edgecolor='black',
        linewidth=0.35,
        alpha=0.94,
        zorder=3,
    )

    x_values = np.asarray([row['map_x'] for row in plot_rows] + [row['family_x'] for row in family_rows])
    y_values = np.asarray([row['map_y'] for row in plot_rows] + [row['family_y'] for row in family_rows])
    x_pad = 0.025 * (x_values.max() - x_values.min() or 1.0)
    y_pad = 0.10 * (y_values.max() - y_values.min() or 1.0)
    ax.set_xlim(x_values.min() - x_pad, x_values.max() + x_pad)
    ax.set_ylim(y_values.min() - y_pad, y_values.max() + y_pad)
    ax.set_aspect('equal', adjustable='box')

    grid_families = [row for row in family_rows if row.get('grid_row') is not None and row.get('grid_column') is not None]
    if grid_families:
        col_positions = sorted({int(row['grid_column']): float(row['family_x']) for row in grid_families}.items())
        row_positions = sorted({int(row['grid_row']): float(row['family_y']) for row in grid_families}.items())
        ax.set_xticks([pos for _, pos in col_positions])
        ax.set_xticklabels([str(idx) for idx, _ in col_positions], fontsize=10)
        ax.set_yticks([pos for _, pos in row_positions])
        ax.set_yticklabels([str(idx) for idx, _ in row_positions], fontsize=10)
        ax.set_xlabel('Column', fontsize=11)
        ax.set_ylabel('Row', fontsize=11)
        ax.tick_params(axis='both', length=3.5, width=0.9, direction='out')
    else:
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_xlabel('')
        ax.set_ylabel('')

    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_linewidth(1.3)
        spine.set_color('black')

    if legend_items and legend_items[0].get('kind') == 'colorbar':
        numeric = np.asarray([
            float(row[color_key])
            for row in plot_rows
            if row.get(color_key) not in (None, '', 'None')
        ], dtype=float)
        numeric = numeric[~np.isnan(numeric)]
        if len(numeric):
            vmin = float(legend_items[0].get('color_vmin', numeric.min()))
            vmax = float(legend_items[0].get('color_vmax', numeric.max()))
            cmap_name_numeric = legend_items[0].get('cmap_name', numeric_colormap_name(color_key))
            cmap, norm = numeric_colormap_and_norm(color_key, vmin, vmax)
            scalar_mappable = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
            scalar_mappable.set_array([])
            fig.canvas.draw_idle()
            ax_box = ax.get_position()
            cbar_height = 0.84 * ax_box.height
            cbar_y = ax_box.y0 + 0.08 * ax_box.height
            cbar_ax = fig.add_axes([ax_box.x1 + 0.025, cbar_y, 0.026, cbar_height])
            cbar = fig.colorbar(scalar_mappable, cax=cbar_ax)
            cbar.set_label(pretty_property_label(color_property), fontsize=12)
            tick_payload = legend_items[0].get('ticks', [])
            if tick_payload:
                cbar.set_ticks([tick['numeric_value'] for tick in tick_payload])
                cbar.set_ticklabels([tick['value'] for tick in tick_payload])
            else:
                cbar.set_ticks(np.linspace(vmin, vmax, min(max(int(top_n_categories), 2), 8)) if vmax != vmin else [vmin])
            cbar.ax.tick_params(labelsize=11, length=4, width=0.9)
            cbar.outline.set_linewidth(1.0)
    elif legend_items:
        handles = [
            plt.Line2D(
                [0], [0],
                marker='o',
                linestyle='none',
                markerfacecolor=item['color'],
                markeredgecolor='black',
                markeredgewidth=0.4,
                markersize=7.0,
                label=item['label'],
            )
            for item in legend_items
        ]
        legend = ax.legend(
            handles=handles,
            title=pretty_property_label(color_property),
            loc='center left',
            bbox_to_anchor=(1.03, 0.5),
            frameon=False,
            fontsize=11,
            title_fontsize=12,
            borderaxespad=0.0,
        )
        legend.get_title().set_fontweight('bold')

    fig.savefig(output_png_path, dpi=dpi, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    return output_png_path

def generate_core_ring_interactive_map(
    descriptor_records,        # Precomputed ligand descriptor dictionaries from collect_ligands_data().
    data,                      # Metadata keyed by refcode; values are added to the map table/tooltips.
    properties_to_show,        # Tooltip specification: display label -> metadata key/unit/digits.
    color_property,            # Column/property used to color the points and build the legend.
    target_denticity,          # Only ligands with this denticity are included in the map.
    output_path,               # Shared output basename for the generated .html and .csv files.
    embedding_method='mds',    # Method for placing families/members by fingerprint similarity: 'mds' or 'pca'.
    grid=False,                # If True, snap family centers to a sparse grid; otherwise keep continuous layout.
    top_n_categories=20,       # Number of categorical color groups kept before assigning the rest to 'Other'.
    grid_cell_spacing=1.0,     # Distance between family bins when grid=True.
    grid_member_radius=0.50,   # Maximum radius occupied by compounds around each family bin when grid=True.
    local_member_min_distance=0.08, # Minimum visual separation between compounds inside the same family, relative to local scale.
    spread_family_centers=True,# If True, gently separate continuous family centers to reduce overlap.
    family_min_distance=0.12,  # Minimum relative separation enforced between family centers in continuous maps.
    display_iframe=False,      # If True, display the generated HTML in the notebook immediately.
    png_figsize=(7.9, 6.6),    # Manuscript PNG size in inches; legend is placed outside the map frame.
    png_dpi=400,               # Resolution of the automatically generated PNG.
    color_properties=None,     # Optional properties exposed by the HTML color selector.
    create_png=True,           # If True, also write a static PNG for color_property.
):
    map_rows, family_rows, selected_records, embedding = build_core_ring_map_data(
        descriptor_records,
        target_denticity=target_denticity,
        data=data,
        embedding_method=embedding_method,
        grid=grid,
        spread_family_centers=spread_family_centers,
        family_min_distance=family_min_distance,
        grid_cell_spacing=grid_cell_spacing,
        grid_member_radius=grid_member_radius,
        local_member_min_distance=local_member_min_distance,
    )

    html_path = str(output_path)+'.html'
    png_path = str(output_path)+'.png' if create_png else None
    csv_path = str(output_path)+'.csv'

    # Creates .csv
    csv_columns = [
        'Refcode',
        'family_id',
        'family_count',
        'grid_row',
        'grid_column',
        'Core-ring family',
    ]
    available_columns = set().union(*(row.keys() for row in map_rows))
    property_specs = resolve_property_specs(normalize_property_specs(properties_to_show), available_columns)
    property_columns = [spec['key'] for spec in property_specs if spec['key'] in available_columns]
    csv_columns = [col for col in csv_columns + property_columns if col in available_columns]
    csv_columns = list(dict.fromkeys(csv_columns))
    with open(csv_path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=csv_columns, extrasaction='ignore')
        writer.writeheader()
        writer.writerows(map_rows)

    layout_name = 'grid' if grid else 'tree'
    title = f'Denticity {target_denticity} homoleptic ligands: {layout_name} map'
    if grid: 
        caption = (
        f'Families are ordered using the {embedding_method.upper()} embedding method, then snapped onto a sparse grid.'
        ' This means that neighboring families should be more similar than far-away families.'
        ' Grey radial lines connect each compound to its family bin center; local positions use fingerprint similarity.'
        )
    else:
        caption = (f'Coordination-core ring families are positioned by fingerprint similarity using {embedding_method.upper()}, then spread for readability.'
        ' Grey radial lines connect each compound to its family center; local positions use fingerprint similarity.'
        )

    # Creates .html
    write_interactive_html(
        map_rows,
        family_rows,
        selected_records,
        output_html_path=html_path,
        properties_to_show=properties_to_show,
        color_property=color_property,
        color_properties=color_properties,
        title=title,
        caption=caption,
        top_n_categories=top_n_categories,
    )

    # Optionally creates a manuscript-style .png without interactive elements.
    if create_png:
        write_static_png(
            map_rows,
            family_rows,
            output_png_path=png_path,
            properties_to_show=properties_to_show,
            color_property=color_property,
            top_n_categories=top_n_categories,
            figsize=png_figsize,
            dpi=png_dpi,
        )


    if display_iframe:
        display(IFrame(str(html_path), width='100%', height=760))

    return {
        'map_df': map_rows,
        'family_rows': family_rows,
        'descriptor_records': descriptor_records,
        'selected_records': selected_records,
        'html_path': html_path,
        'csv_path': csv_path,
        'png_path': png_path,
        'embedding': embedding,
    }

#################################
#### Ligand-Family Statistics ###
#################################
def ligand_family_dhelec_statistics(
    map_data,
    minimum_main_members=6,
    minimum_main_fraction=0.60,
    maximum_gap=15.0,
    outlier_iqr_factor=1.5,
    value_column='dHelec [kJ/mol]',
):
    statistics_rows = []
    membership_rows = []

    for denticity, denticity_rows in map_data.items():
        families = {}
        for row_index, row in enumerate(denticity_rows):
            value = as_float(row.get(value_column))
            if value is not None:
                families.setdefault(row['family_id'], []).append((row_index, row, value))

        for family_id, family_rows in families.items():
            if len(family_rows) < minimum_main_members:
                continue

            values = np.asarray([item[2] for item in family_rows], dtype=float)
            q1, q3 = np.quantile(values, [0.25, 0.75])
            iqr = q3 - q1
            lower_fence = q1 - outlier_iqr_factor * iqr
            upper_fence = q3 + outlier_iqr_factor * iqr
            preliminary = sorted(
                [item for item in family_rows if lower_fence <= item[2] <= upper_fence],
                key=lambda item: item[2],
            )
            if not preliminary:
                continue

            clusters = []
            cluster_start = 0
            preliminary_values = np.asarray([item[2] for item in preliminary], dtype=float)
            for split_point in np.flatnonzero(np.diff(preliminary_values) > maximum_gap) + 1:
                clusters.append(preliminary[cluster_start:split_point])
                cluster_start = split_point
            clusters.append(preliminary[cluster_start:])

            largest_size = max(map(len, clusters))
            largest_clusters = [cluster for cluster in clusters if len(cluster) == largest_size]
            if len(largest_clusters) != 1:
                continue

            main_rows = largest_clusters[0]
            if len(main_rows) < minimum_main_members or len(main_rows) / len(family_rows) < minimum_main_fraction:
                continue

            main_indices = {item[0] for item in main_rows}
            main_values = np.asarray([item[2] for item in main_rows], dtype=float)
            for row_index, row, value in sorted(family_rows, key=lambda item: item[2]):
                is_main = row_index in main_indices
                membership_rows.append({
                    'Denticity': denticity,
                    'Family': family_id,
                    'Refcode': row['Refcode'],
                    value_column: value,
                    'Main Distribution': is_main,
                    'Classification': 'main' if is_main else 'outside',
                })

            statistics_rows.append({
                'Denticity': denticity,
                'Family': family_id,
                'Num Entries': len(family_rows),
                'Num Analysed Entries': len(family_rows),
                'Num Central Entries': len(main_rows),
                'Num Outside Central Distribution': len(family_rows) - len(main_rows),
                'Central Window Min dHelec [kJ/mol]': float(main_values.min()),
                'Central Window Max dHelec [kJ/mol]': float(main_values.max()),
                'dHelec Base [kJ/mol]': float(np.median(main_values)),
                'dHelec Width [kJ/mol]': float(main_values.max() - main_values.min()),
            })

    statistics_rows.sort(key=lambda row: (int(row['Denticity']), row['Family']))
    membership_rows.sort(key=lambda row: (int(row['Denticity']), row['Family'], row[value_column]))
    return statistics_rows, membership_rows

def ligand_family_variance_decomposition(
    map_data,
    minimum_family_size=1,
    value_column='dHelec [kJ/mol]',
    n_permutations=2000,
    random_seed=7,
):
    """Partition within-denticity variance into between- and within-family terms."""
    records = []
    for denticity, rows in map_data.items():
        for row in rows:
            value = as_float(row.get(value_column))
            if value is not None:
                records.append((str(denticity), row['family_id'], value))

    family_sizes = Counter((denticity, family) for denticity, family, _ in records)
    records = [
        record for record in records
        if family_sizes[(record[0], record[1])] >= minimum_family_size
    ]
    if not records:
        raise ValueError('No ligand families satisfy the requested minimum size')

    denticities = np.asarray([record[0] for record in records])
    families = np.asarray([f'{record[0]}::{record[1]}' for record in records])
    values = np.asarray([record[2] for record in records], dtype=float)
    unique_denticities, denticity_codes = np.unique(denticities, return_inverse=True)
    unique_families, family_codes = np.unique(families, return_inverse=True)
    denticity_counts = np.bincount(denticity_codes)
    family_counts = np.bincount(family_codes)
    family_denticity_codes = np.asarray([
        denticity_codes[np.flatnonzero(family_codes == family_index)[0]]
        for family_index in range(len(unique_families))
    ])

    def sums_of_squares(candidate_values):
        denticity_means = np.bincount(denticity_codes, weights=candidate_values) / denticity_counts
        family_means = np.bincount(family_codes, weights=candidate_values) / family_counts
        ss_family = np.sum(
            family_counts * (family_means - denticity_means[family_denticity_codes]) ** 2
        )
        ss_within = np.sum((candidate_values - family_means[family_codes]) ** 2)
        return float(ss_family), float(ss_within)

    ss_family, ss_within = sums_of_squares(values)
    family_df = len(unique_families) - len(unique_denticities)
    within_df = len(values) - len(unique_families)
    if family_df <= 0 or within_df <= 0:
        raise ValueError('Insufficient repeated families for variance decomposition')

    family_ms = ss_family / family_df
    within_ms = ss_within / within_df
    effective_family_size = (
        len(values)
        - sum(np.sum(family_counts[family_denticity_codes == denticity_index] ** 2) / denticity_counts[denticity_index]
              for denticity_index in range(len(unique_denticities)))
    ) / family_df
    family_variance = max((family_ms - within_ms) / effective_family_size, 0.0)
    within_variance = within_ms
    intraclass_correlation = family_variance / (family_variance + within_variance)

    rng = np.random.default_rng(random_seed)
    permutation_ss = np.empty(n_permutations, dtype=float)
    for permutation in range(n_permutations):
        shuffled = values.copy()
        for denticity_index in range(len(unique_denticities)):
            indices = np.flatnonzero(denticity_codes == denticity_index)
            shuffled[indices] = rng.permutation(shuffled[indices])
        permutation_ss[permutation] = sums_of_squares(shuffled)[0]
    permutation_p = (1 + np.sum(permutation_ss >= ss_family)) / (n_permutations + 1)

    return {
        'Minimum Family Size': minimum_family_size,
        'Num Entries': len(values),
        'Num Families': len(unique_families),
        'Between-Family SS Fraction': ss_family / (ss_family + ss_within),
        'Between-Family Variance [kJ2/mol2]': family_variance,
        'Within-Family Variance [kJ2/mol2]': within_variance,
        'Family ICC': intraclass_correlation,
        'Permutation p-value': permutation_p,
        'Num Permutations': n_permutations,
    }

def rounded_rows(rows, digits_by_column):
    rounded = []
    for row in rows:
        result = dict(row)
        for column, digits in digits_by_column.items():
            if column in result and result[column] not in (None, ''):
                result[column] = round(float(result[column]), digits)
        rounded.append(result)
    return rounded

def draw_family_distribution(
    ax,
    y,
    members,
    statistics,
    value_column='dHelec [kJ/mol]',
    baseline_column='dHelec Base [kJ/mol]',
    minimum_column='Central Window Min dHelec [kJ/mol]',
    maximum_column='Central Window Max dHelec [kJ/mol]',
    colors=None,
):
    colors = colors or {
        'central_window': '#9bc9c2',
        'central_point': '#0b4f57',
        'outside': '#dc7b35',
    }
    values = np.asarray([float(row[value_column]) for row in members], dtype=float)
    central_min = float(statistics[minimum_column])
    central_max = float(statistics[maximum_column])

    if values.min() < central_min:
        ax.plot([values.min(), central_min], [y, y], color=colors['outside'], lw=16,
                alpha=0.35, solid_capstyle='butt', zorder=1)
    if values.max() > central_max:
        ax.plot([central_max, values.max()], [y, y], color=colors['outside'], lw=16,
                alpha=0.35, solid_capstyle='butt', zorder=1)

    ax.plot([central_min, central_max], [y, y], color=colors['central_window'], lw=16,
            solid_capstyle='butt', zorder=1)
    main_mask = (values >= central_min) & (values <= central_max)
    for accepted, color in [(True, colors['central_point']), (False, colors['outside'])]:
        selected = np.sort(values[main_mask == accepted])
        if not len(selected):
            continue
        offsets = np.linspace(-0.09, 0.09, len(selected)) if len(selected) > 1 else [0.0]
        ax.scatter(selected, y + np.asarray(offsets), s=36, color=color,
                   edgecolor='white', linewidth=0.55, zorder=3)

    ax.scatter(float(statistics[baseline_column]), y, marker='D', s=60, color='red',
               edgecolor='white', linewidth=0.6, zorder=4)

def add_family_legends(fig, colors=None, fontsize=9):
    colors = colors or {
        'central_window': '#9bc9c2',
        'central_point': '#0b4f57',
        'outside': '#dc7b35',
    }
    area_handles = [
        plt.Line2D([0], [0], color=colors['central_window'], lw=10,
                   label='Main distribution'),
        plt.Line2D([0], [0], color=colors['outside'], lw=10, alpha=0.35,
                   label='Outliers'),
    ]
    point_handles = [
        plt.Line2D([0], [0], marker='o', color='none', markerfacecolor=colors['central_point'],
                   markeredgecolor='white', markersize=7, label='Accepted'),
        plt.Line2D([0], [0], marker='o', color='none', markerfacecolor=colors['outside'],
                   markeredgecolor='white', markersize=7, label='Discarded (outlier)'),
        plt.Line2D([0], [0], marker='D', color='none', markerfacecolor='red',
                   markeredgecolor='white', markersize=7,
                   label=r'Baseline, $\Delta H_{\mathrm{elec}}^{base}$'),
    ]
    fig.legend(
        handles=area_handles + point_handles,
        loc='upper center',
        bbox_to_anchor=(0.5, 0.98),
        ncol=5,
        frameon=False,
        fontsize=fontsize,
        columnspacing=1.6,
        handlelength=1.8,
    )

###################
#### For Debug ####
###################
def atom_debug_table(record):
    mol = record['lig_mol']
    core_atoms = set(record['core_atoms'])
    substituent_atoms = set(record['subs_atoms'])
    ring_membership = {idx: [] for idx in range(mol.GetNumAtoms())}
    for ring_id, ring in enumerate(mol.GetRingInfo().AtomRings()):
        for atom_idx in ring:
            ring_membership[atom_idx].append(ring_id)

    rows = []
    for atom in mol.GetAtoms():
        idx = atom.GetIdx()
        rows.append({
            'atom_idx': idx,
            'symbol': atom.GetSymbol(),
            'is_aromatic': atom.GetIsAromatic(),
            'ring_ids': ring_membership[idx],
            'core': idx in core_atoms,
            'substituent': idx in substituent_atoms,
            'neighbors': [nbr.GetIdx() for nbr in atom.GetNeighbors()],
        })
    return rows

def draw_core_assignment(record):
    mol = record['lig_mol']
    core_atoms = sorted(record['core_atoms'])
    substituent_ring_atoms = sorted(
        atom_idx for ring in mol.GetRingInfo().AtomRings() for atom_idx in ring
        if atom_idx not in set(core_atoms)
    )

    highlight_atoms = sorted(set(core_atoms) | set(substituent_ring_atoms))
    highlight_colors = {}
    for idx in substituent_ring_atoms:
        highlight_colors[idx] = (0.30, 0.55, 1.00)
    for idx in core_atoms:
        highlight_colors[idx] = (0.20, 0.75, 0.35)

    svg = Draw.MolsToGridImage(
        [mol],
        molsPerRow=1,
        subImgSize=(520, 320),
        legends=[f"{record['id']} | green=core atoms, blue=substituent ring atoms"],
        highlightAtomLists=[highlight_atoms],
        highlightAtomColors=[highlight_colors],
        useSVG=True,
    )
    display(SVG(svg.data if hasattr(svg, 'data') else str(svg)))

def summarize_core_families(records):
    return [
        {
            'Refcode': record['id'],
            'denticity': record['denticity'],
            'SMILES': record['smiles'],
            'Core-ring family': record['family'],
            'core_atom_count': len(record['core_atoms']),
            'substituent_atom_count': len(record['subs_atoms']),
        }
        for record in records
    ]

###############
#### Other ####
###############
def as_float(value):
    if value is None or str(value).strip().lower() in ('', 'none', 'nan', '-'):
        return None
    converted = float(value)
    return None if np.isnan(converted) else converted

def get_last_author(authors):
    author_parts = [part.strip() for part in str(authors).split(';') if part.strip()]
    return author_parts[-1] if author_parts else 'unknown'
