import csv
import json
import urllib.parse
import urllib.request
from pathlib import Path

HERE = Path(__file__).resolve().parent
JOURNAL_ABBREVIATIONS_PATH = HERE / "journal_abbreviations_chemistry_merged.csv"

def clean_text(value, debug: int = 0):
    bad_values = [None, "", "-", "None", "nan", "NaN"]
    if value in bad_values:
        if debug:
            print(f"clean_text: rejected empty value {value!r}")
        return ""
    cleaned = str(value).strip()
    if debug and cleaned != value:
        print(f"clean_text: {value!r} -> {cleaned!r}")
    return cleaned


def normalize_text(value, debug: int = 0):
    normalized = "".join(ch.lower() for ch in clean_text(value) if ch.isalnum())
    if debug:
        print(f"normalize_text: {value!r} -> {normalized!r}")
    return normalized


def load_journal_abbreviations(path=JOURNAL_ABBREVIATIONS_PATH, debug: int = 0):
    path = Path(path)
    if not path.exists():
        if debug:
            print(f"load_journal_abbreviations: file not found: {path}")
        return {}

    abbreviations = {}
    with path.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            full = clean_text(row.get("full_journal"))
            abbr = clean_text(row.get("abbreviation"))
            if full and abbr:
                abbreviations[normalize_text(abbr)] = full

    if debug:
        print(f"load_journal_abbreviations: loaded {len(abbreviations)} abbreviations")
    return abbreviations


journal_abbreviations = load_journal_abbreviations()


def crossref_journal_name(journal_name, debug: int = 0):
    aliases = {
        "Chem.-Eur.J.": "Chemistry - A European Journal",
        "Eur.J.Inorg.Chem.": "European Journal of Inorganic Chemistry",
        "Int.J.Inorg.Chem": "International Journal of Inorganic Chemistry",
        "Int.J.Inorg.Chem.": "International Journal of Inorganic Chemistry",
    }
    journal_name = clean_text(journal_name)
    mapped = aliases.get(
        journal_name,
        journal_abbreviations.get(normalize_text(journal_name), journal_name),
    )
    if debug and mapped != clean_text(journal_name):
        print(f"crossref_journal_name: {journal_name!r} -> {mapped!r}")
    return mapped


def journal_names_for_matching(journal_name):
    names = [clean_text(journal_name), crossref_journal_name(journal_name)]
    return [name for i, name in enumerate(names) if name and name not in names[:i]]


def first_cif_author(authors, debug: int = 0):
    if not authors:
        if debug:
            print("first_cif_author: no authors available")
        return ""
    first = str(authors[0]).replace(",", " ").strip()
    if debug:
        print(f"first_cif_author: {first!r}")
    return first


def candidate_year(candidate, debug: int = 0):
    date_parts = candidate.get("issued", {}).get("date-parts", [])
    year = date_parts[0][0] if date_parts and date_parts[0] else None
    if debug:
        print(f"candidate_year: {year}")
    return year


def page_matches(cif_page, crossref_page, debug: int = 0):
    original_cif_page = cif_page
    original_crossref_page = crossref_page
    cif_page = clean_text(cif_page).lower()
    crossref_page = clean_text(crossref_page).lower()
    if not cif_page or not crossref_page:
        if debug:
            print(
                "page_matches: False; incomplete page data "
                f"cif={original_cif_page!r}, crossref={original_crossref_page!r}"
            )
        return False
    if cif_page in crossref_page:
        if debug:
            print(f"page_matches: True; {cif_page!r} found in {crossref_page!r}")
        return True
    # Handles pages like m212-m213 or 11042-11050 by comparing first page tokens.
    match = crossref_page.replace(" ", "").startswith(cif_page.replace(" ", ""))
    if debug:
        print(f"page_matches: {match}; cif={cif_page!r}, crossref={crossref_page!r}")
    return match


def author_matches(cif_authors, candidate, debug: int = 0):
    if not cif_authors:
        if debug:
            print("author_matches: False; no CIF authors")
        return False
    cif_first = normalize_text(first_cif_author(cif_authors), debug=debug)
    if not cif_first:
        if debug:
            print("author_matches: False; empty first CIF author after normalization")
        return False
    if debug:
        print(f"author_matches: looking for CIF first author {cif_first!r}")
    for author in candidate.get("author", []):
        family = normalize_text(author.get("family", ""))
        given = normalize_text(author.get("given", ""))
        full = given + family
        if debug:
            print(
                f"  candidate author: given={given!r}, family={family!r}, "
                f"full={full!r}"
            )
        if family and family in cif_first:
            if debug:
                print("  author_matches: True by family-name match")
            return True
        if full and full in cif_first:
            if debug:
                print("  author_matches: True by full-name match")
            return True
    if debug:
        print("author_matches: False; no candidate author matched")
    return False


def candidate_score(
    candidate,
    journal,
    year,
    volume,
    page,
    authors=None,
    debug: int = 0,
):
    score = 0
    bad_titles = [
        "contents list",
        "front cover",
        "back cover",
        "inside front cover",
        "inside back cover",
        "graphical abstract",
        "erratum",
        "corrigendum",
    ]
    bad_years = [None, "", "-", "None", "nan", "NaN"]

    if debug:
        print("Scoring candidate:")
        print(f"  DOI: {candidate.get('DOI', '')}")
        print(f"  title: {' '.join(candidate.get('title', []))}")
        print(f"  journal: {'; '.join(candidate.get('container-title', []))}")
        print(
            "  volume/page/year: "
            f"{candidate.get('volume', '')} / {candidate.get('page', '')} / "
            f"{candidate_year(candidate)}"
        )

    title = " ".join(candidate.get("title", []))
    if any(bad_title in title.lower() for bad_title in bad_titles):
        score -= 5
        if debug:
            print("  -5 title looks like cover/contents/erratum")

    if candidate.get("type") == "journal-article":
        score += 1
        if debug:
            print("  +1 type is journal-article")
    elif debug:
        print(f"  +0 type is {candidate.get('type')!r}")

    cand_year = candidate_year(candidate)
    if cand_year in bad_years:
        score -= 1
        if debug:
            print("  -1 missing candidate year")
    else:
        try:
            if cand_year == int(year):
                score += 2
                if debug:
                    print(f"  +2 year match: {cand_year}")
            elif debug:
                print(f"  +0 year mismatch: candidate={cand_year}, cif={year}")
        except (TypeError, ValueError):
            if debug:
                print(f"  +0 CIF year not comparable: {year!r}")

    cand_volume = clean_text(candidate.get("volume"))
    if clean_text(volume) and cand_volume.lower() == clean_text(volume).lower():
        score += 2
        if debug:
            print(f"  +2 volume match: {cand_volume}")
    elif debug:
        print(
            "  +0 volume mismatch/incomplete: "
            f"candidate={cand_volume!r}, cif={volume!r}"
        )

    if page_matches(page, candidate.get("page", ""), debug=debug):
        score += 3
        if debug:
            print("  +3 page match")
    elif debug:
        print("  +0 page mismatch/incomplete")

    cand_journals = (
        candidate.get("container-title", [])
        + candidate.get("short-container-title", [])
    )
    journal_norms = [
        normalize_text(name) for name in journal_names_for_matching(journal)
    ]
    journal_match = any(
        journal_norm
        and (
            journal_norm in normalize_text(candidate_journal)
            or normalize_text(candidate_journal) in journal_norm
        )
        for journal_norm in journal_norms
        for candidate_journal in cand_journals
    )
    if journal_match:
        score += 1
        if debug:
            print("  +1 journal match")
    elif debug:
        print(
            f"  +0 journal mismatch: cif={journal!r}, "
            f"aliases={journal_names_for_matching(journal)!r}, "
            f"candidates={cand_journals!r}"
        )

    if author_matches(authors, candidate, debug=debug):
        score += 1
        if debug:
            print("  +1 author match")
    elif debug:
        print("  +0 author mismatch/not found")

    if debug:
        print(f"  final score: {score}")
    return score


def search_doi(
    journal_name,
    journal_year,
    journal_volume,
    journal_page,
    authors=None,
    rows=5,
    debug: int = 0,
):
    query_journal = crossref_journal_name(journal_name, debug=debug)
    values = [
        query_journal,
        journal_year,
        journal_volume,
        journal_page,
        first_cif_author(authors, debug=debug),
    ]
    query = " ".join(clean_text(value) for value in values if clean_text(value))
    if debug:
        print("Crossref DOI search")
        print(f"  journal_name: {journal_name!r}")
        print(f"  query_journal: {query_journal!r}")
        print(f"  journal_year: {journal_year!r}")
        print(f"  journal_volume: {journal_volume!r}")
        print(f"  journal_page: {journal_page!r}")
        print(f"  first_author: {first_cif_author(authors)!r}")
        print(f"  query: {query!r}")
    if not query:
        if debug:
            print("  search_doi: aborting; query is empty")
        return None

    url = "https://api.crossref.org/works?" + urllib.parse.urlencode(
        {"query.bibliographic": query, "rows": rows}
    )
    if debug:
        print(f"  url: {url}")
    request = urllib.request.Request(
        url,
        headers={"User-Agent": "SCO-ATLAS bibliographic search"},
    )
    with urllib.request.urlopen(request, timeout=20) as response:
        data = json.load(response)

    candidates = data.get("message", {}).get("items", [])
    if debug:
        print(f"  candidates retrieved: {len(candidates)}")
    if not candidates:
        if debug:
            print("  search_doi: no candidates returned")
        return None

    scored_candidates = [
        (
            candidate_score(
                candidate,
                journal_name,
                journal_year,
                journal_volume,
                journal_page,
                authors,
                debug=debug,
            ),
            candidate,
        )
        for candidate in candidates
    ]
    score, best = max(scored_candidates, key=lambda item: item[0])

    if debug:
        print("Best Crossref candidate")
        print(f"  score: {score}")
        print(f"  DOI: {best.get('DOI', '')}")
        print(f"  URL: {best.get('URL', '')}")
        print(f"  title: {' '.join(best.get('title', []))}")

    return {
        "query": query,
        "score": score,
        "doi": best.get("DOI", ""),
        "url": best.get("URL", ""),
        "title": " ".join(best.get("title", [])),
    }


def parse_cif_authors(cif_authors):
    authors = []
    for entry in cif_authors or []:
        clean_entry = str(entry).strip().strip("[").strip("]").replace("'", "")
        for author in clean_entry.split(","):
            author = author.strip()
            if author:
                authors.append(author)
    return authors


def get_biblio_info(
    source,
    score_threshold: int = 8,
    crossref_rows: int = 3,
    debug: int = 0,
):
    journal = source.cif.journal_name
    doi, url = "", ""
    match = None
    if journal != "CSD Communication(Private Communication)":
        try:
            match = search_doi(
                source.cif.journal_name,
                source.cif.journal_year,
                source.cif.journal_volume,
                source.cif.journal_page,
                source.cif.authors,
                rows=crossref_rows,
                debug=debug,
            )
            if match and match["score"] >= score_threshold:
                doi = match["doi"]
                url = match["url"]
        except Exception as exc:
            print(f"DOI search failed: {exc}")

    if journal == "CSD Communication(Private Communication)":
        journal = "CSD Communication"

    return {
        "match": match,
        "query": match["query"] if match else "",
        "score": match["score"] if match else "",
        "authors": "; ".join(parse_cif_authors(source.cif.authors)),
        "year": source.cif.journal_year,
        "journal": journal,
        "doi": doi,
        "url": url,
    }


def get_ligands_formula(source):
    unique_ligands = []
    for molecule in source.molecules:
        if molecule.iscomplex:
            for ligand in molecule.ligands:
                if ligand.formula not in unique_ligands:
                    unique_ligands.append(ligand.formula)
    return unique_ligands


def get_complex_formula(source):
    unique_molecules = []
    for molecule in source.molecules:
        if molecule.formula not in unique_molecules:
            unique_molecules.append(molecule.formula)
    return unique_molecules
