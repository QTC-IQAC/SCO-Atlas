import csv
import re
from pathlib import Path


HERE = Path(__file__).resolve().parent
OUTPUT = HERE / "journal_abbreviations_chemistry_merged.csv"


INPUT_FILES = {
    "journal_abbreviations_acs.csv": "acs",
    "journal_abbreviations_general.csv": "general",
    "journal_abbreviations_lifescience.csv": "lifescience",
}


CHEMISTRY_TERMS = (
    "chem",
    "inorg",
    "organic",
    "organometal",
    "coordination",
    "molec",
    "mol.",
    "mater",
    "material",
    "cryst",
    "crystall",
    "catal",
    "polymer",
    "spectro",
    "physical",
    "physico",
    "biochem",
    "biophys",
    "pharm",
    "nano",
    "surface",
    "solid state",
    "electrochem",
    "magnet",
    "supramol",
    "analytical",
    "dalton",
    "tetrahedron",
    "acta chem",
    "angew",
    "langmuir",
    "zeitschrift",
)


LIFESCIENCE_EXCLUDE_TERMS = (
    "medical",
    "medicine",
    "clinical",
    "clin.",
    "nursing",
    "nurse",
    "hospital",
    "health",
    "disease",
    "cancer",
    "oncology",
    "tumor",
    "aids",
    "hiv",
    "virus",
    "viral",
    "bacteria",
    "microbiology",
    "immunology",
    "neuroscience",
    "neurology",
    "brain",
    "psychology",
    "psychiatry",
    "behavior",
    "behaviour",
    "surgery",
    "dental",
    "dentistry",
    "veterinary",
    "nutrition",
    "diabetes",
    "cardiology",
    "radiology",
    "pathology",
    "toxicology",
    "epidemiology",
    "pediatrics",
    "paediatrics",
    "therapy",
    "therapeutics",
    "drug",
    "genetics",
    "genomics",
    "proteomics",
    "cell",
    "botany",
    "zoology",
    "ecology",
    "agriculture",
    "forestry",
    "plant",
)


NON_ENGLISH_HINTS = (
    "revue ",
    "rev ",
    "rivista",
    "annales",
    "annali",
    "zeitschrift",
    "z. ",
    " fuer ",
    " für ",
    "berichte",
    "ber. ",
    "beiheft",
    "bulletin de",
    "boletin",
    "boletín",
    "gazzetta",
    "journal de",
    "journal für",
    "journal fuer",
    "societe",
    "société",
    "sociedad",
    "deutsche",
    "deutscher",
    "dtsch",
    "mitteilungen",
    "verhandlungen",
    "wissenschaftliche",
    "monatshefte",
    "chemie der",
    "chemie ingenieur",
    "pharmazeutische",
    "naturwissenschaften",
    "francaise",
    "française",
    "espanola",
    "española",
    "italiana",
    "japanese journal",
    "chinese journal",
    "russian journal",
)


NON_ENGLISH_EXCEPTIONS = (
    "angewandte chemie, international edition",
    "angewandte chemie international edition",
)


def normalize_key(value):
    value = value.lower().replace("&", "and")
    value = re.sub(r"[^a-z0-9]+", "", value)
    return value


def clean_value(value):
    return " ".join(str(value).strip().split())


def is_chemistry_related(full_name, abbreviation):
    text = f"{full_name} {abbreviation}".lower()
    return any(term in text for term in CHEMISTRY_TERMS)


def looks_non_english(full_name):
    text = full_name.lower()
    if any(exception in text for exception in NON_ENGLISH_EXCEPTIONS):
        return False
    if any(ord(char) > 127 for char in text):
        return True
    return any(hint in text for hint in NON_ENGLISH_HINTS)


def keep_row(full_name, abbreviation, source):
    if not full_name or not abbreviation:
        return False

    if looks_non_english(full_name):
        return False

    if source == "acs":
        return True

    if source == "general":
        return is_chemistry_related(full_name, abbreviation)

    if source == "lifescience":
        text = f"{full_name} {abbreviation}".lower()
        if not is_chemistry_related(full_name, abbreviation):
            return False
        if any(term in text for term in LIFESCIENCE_EXCLUDE_TERMS):
            return False
        return True

    return False


def read_rows():
    for filename, source in INPUT_FILES.items():
        path = HERE / filename
        with path.open(newline="", encoding="utf-8-sig") as handle:
            for row in csv.reader(handle):
                if len(row) < 2:
                    continue
                full_name = clean_value(row[0])
                abbreviation = clean_value(row[1])
                if keep_row(full_name, abbreviation, source):
                    yield {
                        "full_journal": full_name,
                        "abbreviation": abbreviation,
                        "source": source,
                    }


def merge_rows(rows):
    merged = {}
    for row in rows:
        key = (normalize_key(row["full_journal"]), normalize_key(row["abbreviation"]))
        if key not in merged:
            row["sources"] = row.pop("source")
            merged[key] = row
        else:
            sources = set(merged[key]["sources"].split(";"))
            sources.add(row["source"])
            merged[key]["sources"] = ";".join(sorted(sources))
    return sorted(
        merged.values(),
        key=lambda row: (row["full_journal"].lower(), row["abbreviation"].lower()),
    )


def main():
    rows = merge_rows(read_rows())
    for row in rows:
        row.pop("sources", None)
    fieldnames = ["full_journal", "abbreviation"]
    with OUTPUT.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    print(f"Wrote {len(rows)} rows to {OUTPUT}")


if __name__ == "__main__":
    main()
