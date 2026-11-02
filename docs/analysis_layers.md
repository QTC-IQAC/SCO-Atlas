# Analysis Layers

The SCO ATLAS analyses are organized as chronological, self-contained layers.
Each layer has a descriptive numbered name and owns its notebooks, derived
outputs, figures, and any analysis-specific support files.

```text
SCO-Atlas/
└── 1-TMC_Thermodynamics/
    ├── Bibliography/
    ├── Notebooks/
    ├── Outputs/
    └── Figures/
        └── HTML/
```

Future analyses should be added as sibling directories, for example
`2-<Analysis_Name>/`, rather than mixed into the existing layer. Shared code can
be introduced later if two or more layers genuinely require it.

Static figures are reproducible notebook products and remain untracked. HTML
maps and their associated CSV tables are versioned under each layer's
`Figures/HTML/` directory because they are direct interactive research outputs.
