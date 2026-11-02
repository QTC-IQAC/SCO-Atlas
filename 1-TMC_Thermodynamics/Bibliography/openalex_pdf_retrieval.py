import argparse
import csv
import json
import os
import re
import shutil
import urllib.parse
import urllib.request
from pathlib import Path

HERE = Path(__file__).resolve().parent
PROJECT_ROOT = HERE.parents[2]
OPENALEX_API = "https://api.openalex.org"
DEFAULT_INPUT_CSV = PROJECT_ROOT / "Outputs" / "unique_bibliographic_data.csv"
DEFAULT_INPUT_BIB = PROJECT_ROOT / "Manuscript" / "references.bib"
DEFAULT_OUTPUT_DIR = PROJECT_ROOT / "Literature" / "PDFs"
DEFAULT_REPORT_CSV = PROJECT_ROOT / "Outputs" / "openalex_pdf_retrieval.csv"
DEFAULT_BIB_REPORT_CSV = PROJECT_ROOT / "Outputs" / "openalex_references_pdf_retrieval.csv"
REPORT_FIELDS = [
    "identifiers",
    "doi",
    "input_url",
    "input_file",
    "title",
    "openalex_id",
    "status",
    "source",
    "version",
    "license",
    "pdf_url",
    "local_path",
    "error",
]

def normalize_doi(doi):
    doi = str(doi or "").strip()
    doi = re.sub(r"^https?://(?:dx\.)?doi\.org/", "", doi, flags=re.I)
    doi = re.sub(r"^doi:\s*", "", doi, flags=re.I)
    return doi.strip()

def safe_filename(value):
    value = re.sub(r"[^A-Za-z0-9._-]+", "_", str(value or "").strip())
    return value.strip("._-") or "article"


def parse_biblatex(bib_path):
    """Read citation keys and DOI/URL/file fields from a BibLaTeX export."""
    records = []
    current = None
    with Path(bib_path).open(encoding="utf-8") as handle:
        for line in handle:
            entry_match = re.match(r"^\s*@\w+\s*[\{(]\s*([^,]+),", line)
            if entry_match:
                if current:
                    records.append(current)
                current = {
                    "citation_key": entry_match.group(1).strip(),
                    "doi": "",
                    "url": "",
                    "file": "",
                }
                continue
            if current is None:
                continue

            field_match = re.match(r"^\s*(doi|url|file)\s*=\s*(.*?)\s*,?\s*$", line, re.I)
            if not field_match:
                continue
            value = field_match.group(2).strip()
            if len(value) >= 2 and (value[0], value[-1]) in [("{", "}"), ('"', '"')]:
                value = value[1:-1]
            current[field_match.group(1).lower()] = value.strip()

    if current:
        records.append(current)
    return records


def existing_pdf_attachment(file_value):
    for value in str(file_value or "").split(";"):
        value = re.sub(r":application/pdf$", "", value.strip(), flags=re.I)
        if value.startswith(":"):
            value = value[1:]
        path = Path(value).expanduser()
        if path.suffix.lower() == ".pdf" and path.is_file():
            return path
    return None


def write_report(results, report_csv):
    report_csv = Path(report_csv)
    report_csv.parent.mkdir(parents=True, exist_ok=True)
    with report_csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=REPORT_FIELDS, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(results)
    print(f"Wrote {len(results)} DOI retrieval results to {report_csv}")

def openalex_request(
    endpoint,
    params=None,
    email="",
    api_key="",
    timeout=30,
    debug: int = 0,
):
    params = dict(params or {})
    if email:
        params["mailto"] = email
    if api_key:
        params["api_key"] = api_key

    url = f"{OPENALEX_API}/{endpoint.lstrip('/')}"
    if params:
        url += "?" + urllib.parse.urlencode(params)
    if debug:
        print(f"OpenAlex request: {url}")

    agent = "SCO-ATLAS literature retrieval"
    if email:
        agent += f" (mailto:{email})"
    request = urllib.request.Request(url, headers={"User-Agent": agent})
    with urllib.request.urlopen(request, timeout=timeout) as response:
        return json.load(response)

def get_work_by_doi(doi, email="", api_key="", timeout=30, debug: int = 0):
    doi = normalize_doi(doi)
    if not doi:
        raise ValueError("A DOI is required for exact OpenAlex matching")
    endpoint = "works/doi:" + urllib.parse.quote(doi, safe="")
    return openalex_request(
        endpoint,
        email=email,
        api_key=api_key,
        timeout=timeout,
        debug=debug,
    )

def search_works(
    query,
    per_page=25,
    filters="",
    email="",
    api_key="",
    timeout=30,
    debug: int = 0,
):
    """Search OpenAlex by topic without treating the results as verified matches."""
    params = {"search": query, "per-page": per_page}
    if filters:
        params["filter"] = filters
    data = openalex_request(
        "works",
        params=params,
        email=email,
        api_key=api_key,
        timeout=timeout,
        debug=debug,
    )
    return data.get("results", [])

def get_pdf_candidates(work):
    """Return unique, OA PDF locations in OpenAlex preference order."""
    locations = [
        ("best_oa_location", work.get("best_oa_location")),
        ("primary_location", work.get("primary_location")),
    ]
    locations.extend(("location", location) for location in work.get("locations", []))

    candidates = []
    seen = set()
    for source, location in locations:
        if not location:
            continue
        pdf_url = location.get("pdf_url")
        if not pdf_url or pdf_url in seen:
            continue
        if source != "best_oa_location" and location.get("is_oa") is not True:
            continue
        seen.add(pdf_url)
        candidates.append(
            {
                "source": source,
                "pdf_url": pdf_url,
                "landing_page_url": location.get("landing_page_url", ""),
                "version": location.get("version", ""),
                "license": location.get("license", ""),
            }
        )
    return candidates

def download_pdf(
    url,
    output_path,
    email="",
    timeout=60,
    overwrite=False,
    debug: int = 0,
):
    """Download a URL only when its response starts with a PDF signature."""
    output_path = Path(output_path)
    if output_path.exists() and not overwrite:
        return "existing"

    output_path.parent.mkdir(parents=True, exist_ok=True)
    partial_path = output_path.with_suffix(output_path.suffix + ".part")
    agent = "SCO-ATLAS literature retrieval"
    if email:
        agent += f" (mailto:{email})"
    request = urllib.request.Request(
        url,
        headers={
            "User-Agent": agent,
            "Accept": "application/pdf",
        },
    )

    if debug:
        print(f"Downloading: {url}")
    try:
        with urllib.request.urlopen(request, timeout=timeout) as response:
            signature = response.read(5)
            if signature != b"%PDF-":
                content_type = response.headers.get("Content-Type", "")
                raise ValueError(
                    f"Response is not a PDF (Content-Type: {content_type or 'unknown'})"
                )
            with partial_path.open("wb") as handle:
                handle.write(signature)
                shutil.copyfileobj(response, handle)
        partial_path.replace(output_path)
    finally:
        if partial_path.exists():
            partial_path.unlink()
    return "downloaded"


def retrieve_pdf(
    doi,
    output_path,
    input_url="",
    email="",
    api_key="",
    timeout=60,
    overwrite=False,
    debug: int = 0,
):
    """Resolve one DOI through OpenAlex and try each verified OA PDF location."""
    result = {
        "doi": normalize_doi(doi),
        "input_url": input_url,
        "title": "",
        "openalex_id": "",
        "status": "",
        "source": "",
        "version": "",
        "license": "",
        "pdf_url": "",
        "local_path": "",
        "error": "",
    }
    if not result["doi"]:
        result["status"] = "missing_doi"
        return result

    try:
        work = get_work_by_doi(
            result["doi"],
            email=email,
            api_key=api_key,
            timeout=timeout,
            debug=debug,
        )
    except Exception as exc:
        result["status"] = "openalex_error"
        result["error"] = str(exc)
        return result

    result["title"] = work.get("display_name", "")
    result["openalex_id"] = work.get("id", "")
    candidates = get_pdf_candidates(work)
    if not candidates:
        result["status"] = "no_open_pdf"
        return result

    errors = []
    for candidate in candidates:
        try:
            status = download_pdf(
                candidate["pdf_url"],
                output_path,
                email=email,
                timeout=timeout,
                overwrite=overwrite,
                debug=debug,
            )
            result.update(candidate)
            result["status"] = status
            result["local_path"] = str(Path(output_path).resolve())
            return result
        except Exception as exc:
            errors.append(f"{candidate['pdf_url']}: {exc}")

    result["status"] = "download_failed"
    result["error"] = " | ".join(errors)
    return result


def retrieve_pdfs_from_csv(
    input_csv,
    output_dir,
    report_csv,
    doi_column="doi",
    url_column="url",
    id_column="refcode",
    email="",
    api_key="",
    timeout=60,
    overwrite=False,
    selected_dois=None,
    selected_ids=None,
    debug: int = 0,
):
    """Retrieve one PDF per unique DOI and write an auditable CSV report."""
    input_csv = Path(input_csv)
    output_dir = Path(output_dir)
    report_csv = Path(report_csv)
    with input_csv.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))

    selected_dois = {
        normalize_doi(doi).lower(): normalize_doi(doi)
        for doi in selected_dois or []
        if normalize_doi(doi)
    }
    selected_ids = {
        str(identifier).strip().upper()
        for identifier in selected_ids or []
        if str(identifier).strip()
    }
    if selected_dois or selected_ids:
        rows = [
            row
            for row in rows
            if normalize_doi(row.get(doi_column, "")).lower() in selected_dois
            or str(row.get(id_column, "")).strip().upper() in selected_ids
        ]

        # A DOI can be retrieved directly even when it is absent from the input CSV.
        present_dois = {
            normalize_doi(row.get(doi_column, "")).lower() for row in rows
        }
        for doi_key, doi in selected_dois.items():
            if doi_key not in present_dois:
                rows.append({doi_column: doi, id_column: doi})

    grouped = {}
    missing_doi_rows = []
    for row in rows:
        doi = normalize_doi(row.get(doi_column, ""))
        if not doi:
            missing_doi_rows.append(row)
            continue
        grouped.setdefault(doi.lower(), []).append(row)

    results = []
    for doi_key, doi_rows in grouped.items():
        identifiers = []
        for row in doi_rows:
            identifier = str(row.get(id_column, "")).strip()
            if identifier and identifier not in identifiers:
                identifiers.append(identifier)
        identifier = identifiers[0] if identifiers else doi_key
        output_path = output_dir / f"{safe_filename(identifier)}.pdf"
        result = retrieve_pdf(
            doi_key,
            output_path,
            input_url=next(
                (row.get(url_column, "") for row in doi_rows if row.get(url_column)),
                "",
            ),
            email=email,
            api_key=api_key,
            timeout=timeout,
            overwrite=overwrite,
            debug=debug,
        )
        result["identifiers"] = "; ".join(identifiers)
        results.append(result)
        if debug:
            print(f"{result['identifiers'] or doi_key}: {result['status']}")

    for row in missing_doi_rows:
        results.append(
            {
                "identifiers": row.get(id_column, ""),
                "doi": "",
                "input_url": row.get(url_column, ""),
                "title": "",
                "openalex_id": "",
                "status": "missing_doi",
                "source": "",
                "version": "",
                "license": "",
                "pdf_url": "",
                "local_path": "",
                "error": "",
            }
        )

    write_report(results, report_csv)
    return results


def retrieve_pdfs_from_bib(
    input_bib,
    output_dir,
    report_csv,
    email="",
    api_key="",
    timeout=60,
    overwrite=False,
    selected_dois=None,
    selected_ids=None,
    debug: int = 0,
):
    """Retrieve PDFs for a BibLaTeX library, preserving existing attachments."""
    records = parse_biblatex(input_bib)
    output_dir = Path(output_dir)

    selected_dois = {
        normalize_doi(doi).lower(): normalize_doi(doi)
        for doi in selected_dois or []
        if normalize_doi(doi)
    }
    selected_ids = {
        str(identifier).strip().lower()
        for identifier in selected_ids or []
        if str(identifier).strip()
    }
    if selected_dois or selected_ids:
        records = [
            record
            for record in records
            if normalize_doi(record.get("doi", "")).lower() in selected_dois
            or record.get("citation_key", "").lower() in selected_ids
        ]
        present_dois = {
            normalize_doi(record.get("doi", "")).lower() for record in records
        }
        for doi_key, doi in selected_dois.items():
            if doi_key not in present_dois:
                records.append(
                    {
                        "citation_key": doi,
                        "doi": doi,
                        "url": "",
                        "file": "",
                    }
                )

    grouped = {}
    missing_doi_records = []
    for record in records:
        doi = normalize_doi(record.get("doi", ""))
        if not doi:
            missing_doi_records.append(record)
            continue
        grouped.setdefault(doi.lower(), []).append(record)

    results = []
    for doi_key, doi_records in grouped.items():
        identifiers = [
            record["citation_key"]
            for record in doi_records
            if record.get("citation_key")
        ]
        identifier = identifiers[0] if identifiers else doi_key
        attachment = None
        for record in doi_records:
            attachment = existing_pdf_attachment(record.get("file", ""))
            if attachment:
                break
        input_url = next(
            (record.get("url", "") for record in doi_records if record.get("url")),
            "",
        )
        input_file = next(
            (record.get("file", "") for record in doi_records if record.get("file")),
            "",
        )

        if attachment:
            result = {
                "doi": doi_key,
                "input_url": input_url,
                "input_file": input_file,
                "title": "",
                "openalex_id": "",
                "status": "existing_attachment",
                "source": "zotero_attachment",
                "version": "",
                "license": "",
                "pdf_url": "",
                "local_path": str(attachment.resolve()),
                "error": "",
            }
        else:
            output_path = output_dir / f"{safe_filename(identifier)}.pdf"
            result = retrieve_pdf(
                doi_key,
                output_path,
                input_url=input_url,
                email=email,
                api_key=api_key,
                timeout=timeout,
                overwrite=overwrite,
                debug=debug,
            )
            result["input_file"] = input_file

        result["identifiers"] = "; ".join(identifiers)
        results.append(result)
        if debug:
            print(f"{result['identifiers'] or doi_key}: {result['status']}")

    for record in missing_doi_records:
        attachment = existing_pdf_attachment(record.get("file", ""))
        results.append(
            {
                "identifiers": record.get("citation_key", ""),
                "doi": "",
                "input_url": record.get("url", ""),
                "input_file": record.get("file", ""),
                "title": "",
                "openalex_id": "",
                "status": "existing_attachment" if attachment else "missing_doi",
                "source": "zotero_attachment" if attachment else "",
                "version": "",
                "license": "",
                "pdf_url": "",
                "local_path": str(attachment.resolve()) if attachment else "",
                "error": "",
            }
        )

    write_report(results, report_csv)
    return results


def main():
    parser = argparse.ArgumentParser(
        description="Retrieve open-access PDFs for bibliography records through OpenAlex."
    )
    inputs = parser.add_mutually_exclusive_group()
    inputs.add_argument(
        "--input",
        type=Path,
        help=f"Bibliographic CSV input (default: {DEFAULT_INPUT_CSV}).",
    )
    inputs.add_argument(
        "--bib",
        type=Path,
        help=f"BibLaTeX input, for example {DEFAULT_INPUT_BIB}.",
    )
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--report", type=Path)
    parser.add_argument(
        "--doi",
        action="append",
        default=[],
        help="Retrieve this DOI; repeat the option to request several DOIs.",
    )
    parser.add_argument(
        "--id",
        action="append",
        default=[],
        help="Retrieve this CSV refcode or BibLaTeX citation key; repeat as needed.",
    )
    parser.add_argument(
        "--email",
        default=os.environ.get("OPENALEX_MAIL_ADDRESS", ""),
        help="Optional contact email sent to OpenAlex; it is not authentication.",
    )
    parser.add_argument(
        "--api-key",
        default=os.environ.get("OPENALEX_API_KEY", ""),
    )
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--debug", action="store_true")
    args = parser.parse_args()
    report = args.report or (
        DEFAULT_BIB_REPORT_CSV if args.bib else DEFAULT_REPORT_CSV
    )

    common = {
        "email": args.email,
        "api_key": args.api_key,
        "overwrite": args.overwrite,
        "selected_dois": args.doi,
        "selected_ids": args.id,
        "debug": int(args.debug),
    }
    if args.bib:
        retrieve_pdfs_from_bib(
            args.bib,
            args.output_dir,
            report,
            **common,
        )
    else:
        retrieve_pdfs_from_csv(
            args.input or DEFAULT_INPUT_CSV,
            args.output_dir,
            report,
            **common,
        )


if __name__ == "__main__":
    main()
