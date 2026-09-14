import csv
import html
import json
import re
import urllib.parse
import urllib.request
from pathlib import Path


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
INPUT_CSV = ROOT / "Outputs" / "unique_bibliographic_data.csv"
OUTPUT_CSV = ROOT / "Outputs" / "experimental_transition_candidates_first10.csv"


def clean_html(text):
    text = re.sub(r"<script.*?</script>|<style.*?</style>", " ", text, flags=re.I | re.S)
    text = re.sub(r"<[^>]+>", " ", text)
    return re.sub(r"\s+", " ", html.unescape(text)).strip()


def fetch_url_text(url, max_bytes=400_000, timeout=25):
    request = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
    with urllib.request.urlopen(request, timeout=timeout) as response:
        final_url = response.geturl()
        raw = response.read(max_bytes).decode("utf-8", "ignore")
    return final_url, clean_html(raw)


def fetch_crossref_text(doi):
    if not doi:
        return "", "", ""
    url = "https://api.crossref.org/works/" + urllib.parse.quote(doi, safe="")
    with urllib.request.urlopen(url, timeout=20) as response:
        item = json.load(response).get("message", {})
    title = " ".join(item.get("title", []))
    abstract = clean_html(item.get("abstract", ""))
    return url, title, abstract


def fetch_available_manuscript_text(doi, url):
    texts = []

    if doi:
        try:
            crossref_url, title, abstract = fetch_crossref_text(doi)
            if title:
                texts.append(("crossref_title", crossref_url, title))
            if abstract:
                texts.append(("crossref_abstract", crossref_url, abstract))
        except Exception as exc:
            texts.append(("crossref_error", "", str(exc)))

    target_url = url or (f"https://doi.org/{doi}" if doi else "")
    if target_url:
        try:
            final_url, text = fetch_url_text(target_url)
            if text:
                texts.append(("landing_page", final_url, text))
        except Exception as exc:
            texts.append(("landing_page_error", target_url, str(exc)))

    return texts


def extract_transition_mentions(text):
    patterns = [
        r"(?:T\s*1\s*/\s*2|T\s*1\/2|T1/2|T<sub>1/2</sub>|T\s*\(?LIESST\)?|spin transition|transition)[^.]{0,220}?\b\d{2,3}\s*K[^.]{0,120}",
        r"\b\d{2,3}\s*K[^.]{0,160}?(?:T\s*1\s*/\s*2|T\s*1\/2|T1/2|spin transition|transition)",
        r"centered at[^.]{0,160}?\b\d{2,3}\s*(?:and\s*)?\d{0,3}\s*K",
    ]
    mentions = []
    for pattern in patterns:
        for match in re.finditer(pattern, text, flags=re.I):
            snippet = re.sub(r"\s+", " ", match.group(0)).strip()
            if snippet and snippet not in mentions:
                mentions.append(snippet)
    return mentions


def extract_temperatures(mentions):
    values = []
    for mention in mentions:
        for match in re.finditer(r"\b(\d{2,3})(?:\s*(?:±|\+/-)\s*\d{1,3})?\s*K\b", mention):
            value = match.group(1)
            if value not in values:
                values.append(value)
        for before_k in re.findall(r"((?:\d{2,3}\s*(?:,|and)\s*)+\d{2,3})\s*K\b", mention):
            for value in re.findall(r"\d{2,3}", before_k):
                if value not in values:
                    values.append(value)
    return values


def grouped_refcodes(rows):
    grouped = {}
    for row in rows:
        grouped.setdefault(row["refcode"], []).append(row)
    return grouped


def compact_unique(values):
    compact = []
    for value in values:
        if value and value not in compact:
            compact.append(value)
    return compact


def main():
    with INPUT_CSV.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))

    output_rows = []
    for refcode, ref_rows in grouped_refcodes(rows).items():
        doi = next((row["doi"] for row in ref_rows if row.get("doi")), "")
        url = next((row["url"] for row in ref_rows if row.get("url")), "")
        texts = fetch_available_manuscript_text(doi, url)

        mentions = []
        evidence_sources = []
        errors = []
        for source_type, source_url, text in texts:
            if source_type.endswith("_error"):
                errors.append(f"{source_type}: {text}")
                continue
            source_mentions = extract_transition_mentions(text)
            if source_mentions:
                mentions.extend(source_mentions)
                evidence_sources.append(source_type)

        mentions = compact_unique(mentions)
        output_rows.append({
            "refcode": refcode,
            "doi": doi,
            "url": url,
            "mol_formulae": "; ".join(compact_unique(row.get("mol_formulae", "") for row in ref_rows)),
            "lig_formulae": "; ".join(compact_unique(row.get("lig_formulae", "") for row in ref_rows)),
            "sym_groups": "; ".join(compact_unique(row.get("sym_group", "") for row in ref_rows)),
            "diff_temps[K]": "; ".join(compact_unique(row.get("diff_temp", "") for row in ref_rows)),
            "candidate_transition_temperatures[K]": "; ".join(extract_temperatures(mentions)),
            "candidate_transition_mentions": " || ".join(mentions),
            "evidence_sources": "; ".join(compact_unique(evidence_sources)),
            "retrieval_errors": " || ".join(errors),
            "curation_status": "candidate_mentions_found" if mentions else "needs_manual_or_llm_reading",
        })

    fieldnames = [
        "refcode",
        "doi",
        "url",
        "mol_formulae",
        "lig_formulae",
        "sym_groups",
        "diff_temps[K]",
        "candidate_transition_temperatures[K]",
        "candidate_transition_mentions",
        "evidence_sources",
        "retrieval_errors",
        "curation_status",
    ]
    with OUTPUT_CSV.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(output_rows)

    print(f"Wrote {len(output_rows)} rows to {OUTPUT_CSV}")


if __name__ == "__main__":
    main()
