# -*- coding: utf-8 -*-
"""Cross-reference selected bacterial proteins using KEGG and UniProtKB.

The script evaluates selected bacterial species for annotated MiaA, MnmE/TrmE
and MnmG/GidA homologues using two complementary annotation resources:

1. KEGG Orthology assignments among KEGG genome representatives matching each
   species.
2. UniProtKB searches based on gene names, protein names and functional/domain
   annotations.

No genome sequence files are downloaded and no sequence-similarity search is
performed. The analysis therefore evaluates database annotations rather than
establishing homology de novo.

Requirements
------------
Python 3.9 or later and the following packages:
    requests, pandas, openpyxl

Outputs
-------
A user-selected output directory containing a timestamped results folder with:
    actinomycetota_MnmEG_validation_results.xlsx
    actinomycetota_MnmEG_validation_summary.csv
    actinomycetota_MnmEG_validation_all_uniprot_hits.csv

The species, targets, scoring thresholds and network settings can be modified
in the configuration section below.
"""

# =============================================================================
# USER SETTINGS
# =============================================================================

OUTPUT_BASE_DIR = None  # Selected interactively at runtime.
OUTPUT_FOLDER_PREFIX = "actinomycetota_MnmEG_validation"
CREATE_TIMESTAMPED_OUTPUT_FOLDER = True

# Keep the same stable cache folder as the KEGG conservation script.
CACHE_FOLDER_NAME = "_cache"

# Selected Actinomycetota species. Taxonomy IDs are NCBI species IDs.
# Species included in the validation analysis.
SPECIES_TO_CHECK = [
    {"species": "Mycobacterium tuberculosis", "taxid": "1773"},
    {"species": "Mycobacterium smegmatis", "taxid": "1772"},
    {"species": "Corynebacterium glutamicum", "taxid": "1718"},
    {"species": "Streptomyces coelicolor", "taxid": "1902"},
    {"species": "Streptomyces lividans", "taxid": "1916"},
    {"species": "Nocardia farcinica", "taxid": "37329"},
    {"species": "Bifidobacterium longum", "taxid": "216816"},
]

# Include MiaA as a useful positive-control target, while focusing interpretation
# mainly on MnmE/MnmG.
TARGETS = {
    "MiaA": {
        "ko": "K00791",
        "gene_synonyms": ["miaA"],
        "search_terms": [
            "miaA",
            "tRNA dimethylallyltransferase",
            "dimethylallyltransferase",
            "tRNA isopentenyltransferase",
        ],
        "protein_keywords": [
            "dimethylallyltransferase",
            "isopentenyltransferase",
            "miaa",
        ],
        "domain_keywords": [
            "trna dimethylallyltransferase",
            "trna isopentenyltransferase",
            "ippt",
            "miaa",
        ],
    },
    "MnmE": {
        "ko": "K03650",
        "gene_synonyms": ["mnmE", "trmE", "thdF"],
        "search_terms": [
            "mnmE",
            "trmE",
            "thdF",
            "tRNA modification GTPase",
            "GTPase MnmE",
            "GTPase TrmE",
        ],
        "protein_keywords": [
            "trna modification gtpase",
            "gtpase mnme",
            "gtpase trme",
            "mnme",
            "trme",
            "thdf",
        ],
        "domain_keywords": [
            "trme",
            "mnme",
            "gtpase",
            "mss1",
            "trna modification",
        ],
    },
    "MnmG": {
        "ko": "K03495",
        "gene_synonyms": ["mnmG", "gidA", "mto1"],
        "search_terms": [
            "mnmG",
            "gidA",
            "mto1",
            "glucose-inhibited division protein A",
            "tRNA uridine 5-carboxymethylaminomethyl modification enzyme",
            "carboxymethylaminomethyl",
        ],
        "protein_keywords": [
            "trna uridine 5-carboxymethylaminomethyl modification enzyme",
            "glucose-inhibited division protein a",
            "carboxymethylaminomethyl",
            "mnmg",
            "gida",
            "mto1",
        ],
        "domain_keywords": [
            "mnmg",
            "gida",
            "mto1",
            "fad",
            "nad-binding",
            "carboxymethylaminomethyl",
            "trna modification",
        ],
    },
}

# UniProt settings.
UNIPROT_BASE_URL = "https://rest.uniprot.org/uniprotkb/search"
UNIPROT_FIELDS_FULL = (
    "accession,id,reviewed,protein_name,gene_names,organism_name,organism_id,"
    "length,xref_interpro,xref_pfam,xref_kegg,cc_function"
)
UNIPROT_FIELDS_MINIMAL = "accession,id,reviewed,protein_name,gene_names,organism_name,organism_id,length"
UNIPROT_SIZE_PER_QUERY = 100
REVIEWED_ONLY = False

# Score thresholds for calling hits.
STRONG_HIT_SCORE = 5
POSSIBLE_HIT_SCORE = 2

# Network/cache settings.
REQUEST_TIMEOUT_SECONDS = 120
POLITE_DELAY_SECONDS = 0.5
USE_CACHED_DOWNLOADS = True

# =============================================================================
# SCRIPT
# =============================================================================

import io
import tkinter as tk
from tkinter import filedialog
import json
import re
import sys
import time
from datetime import datetime
from pathlib import Path

import pandas as pd
import requests


KEGG_ENDPOINTS = {
    "genomes": "list/genome/bacteria",
}

HEADERS = {
    "User-Agent": "tRNA-modifier-annotation-validation/1.0",
    "Accept": "text/plain,*/*;q=0.8",
}


def select_output_base_folder() -> Path:
    """Prompt the user to select the directory in which results are created."""
    root = tk.Tk()
    root.withdraw()
    root.update()
    selected = filedialog.askdirectory(title="Select output directory")
    root.destroy()
    if not selected:
        raise RuntimeError("No output directory was selected.")
    return Path(selected).expanduser().resolve()


def output_base_folder() -> Path:
    """Return the output base directory selected during program startup."""
    if OUTPUT_BASE_DIR is None:
        raise RuntimeError("The output directory has not been initialized.")
    base = Path(OUTPUT_BASE_DIR)
    base.mkdir(parents=True, exist_ok=True)
    return base


def cache_folder() -> Path:
    cache = output_base_folder() / CACHE_FOLDER_NAME
    cache.mkdir(parents=True, exist_ok=True)
    return cache


def make_output_folder() -> Path:
    base = output_base_folder()
    if CREATE_TIMESTAMPED_OUTPUT_FOLDER:
        stamp = datetime.now().strftime("%Y_%m_%d_%Hh%Mm%Ss")
        out = base / f"{OUTPUT_FOLDER_PREFIX}_{stamp}"
    else:
        out = base / OUTPUT_FOLDER_PREFIX
    out.mkdir(parents=True, exist_ok=True)
    return out


def kegg_url(endpoint: str) -> str:
    return "https://rest.kegg.jp/" + endpoint.lstrip("/")


def endpoint_to_file_name(endpoint: str) -> str:
    return endpoint.replace("/", "_") + ".tsv"


def download_text_cached(url: str, cache_name: str) -> str:
    path = cache_folder() / cache_name
    if USE_CACHED_DOWNLOADS and path.exists() and path.stat().st_size > 0:
        print(f"Using cached file: {path}")
        return path.read_text(encoding="utf-8")

    print(f"Downloading: {url}")
    response = requests.get(url, headers=HEADERS, timeout=REQUEST_TIMEOUT_SECONDS)
    if response.status_code != 200 or not response.text.strip():
        preview = response.text[:500].replace("\n", " ").replace("\r", " ")
        raise RuntimeError(
            f"Download failed: {url}\n"
            f"HTTP status: {response.status_code}\n"
            f"Response preview: {preview!r}"
        )
    path.write_text(response.text, encoding="utf-8")
    time.sleep(POLITE_DELAY_SECONDS)
    return response.text


def get_kegg_text(endpoint: str) -> str:
    return download_text_cached(kegg_url(endpoint), endpoint_to_file_name(endpoint))


def parse_kegg_bacterial_genome_list(text: str) -> pd.DataFrame:
    rows = []
    for line in text.splitlines():
        if not line.strip() or "\t" not in line:
            continue
        left, right = line.split("\t", 1)
        t_match = re.search(r"T\d{5}", left)
        if not t_match or ";" not in right:
            continue
        org_code, organism_name = right.split(";", 1)
        rows.append({
            "t_number": t_match.group(0),
            "org_code": org_code.strip(),
            "organism_name": organism_name.strip(),
        })
    df = pd.DataFrame(rows)
    if df.empty:
        raise RuntimeError("Could not parse KEGG list/genome/bacteria output.")
    return df.drop_duplicates("t_number").reset_index(drop=True)


def parse_kegg_ko_gene_links(text: str) -> pd.DataFrame:
    rows = []
    for line in text.splitlines():
        parts = line.strip().split("\t") if "\t" in line else line.strip().split()
        if len(parts) < 2 or ":" not in parts[1]:
            continue
        gene_prefix, gene_id = parts[1].split(":", 1)
        rows.append({
            "ko_entry": parts[0].strip(),
            "gene_entry": parts[1].strip(),
            "gene_prefix": gene_prefix,
            "gene_id": gene_id,
        })
    return pd.DataFrame(rows)


def normalize_text(x) -> str:
    return str(x if x is not None else "").lower()


def normalize_gene_tokens(x) -> set:
    txt = normalize_text(x)
    txt = re.sub(r"[;,|()/\[\]{}]", " ", txt)
    return {t.strip().lower() for t in txt.split() if t.strip()}


def quote_term_for_uniprot(term: str) -> str:
    term = str(term).strip()
    if " " in term or "-" in term:
        return '"' + term.replace('"', "") + '"'
    return term


def build_uniprot_query(taxid: str, target: dict, broad: bool = False) -> str:
    terms = []

    if not broad:
        # Prioritize exact gene names while retaining broader annotation terms.
        terms.extend([f"gene_exact:{g}" for g in target["gene_synonyms"]])

    terms.extend([quote_term_for_uniprot(t) for t in target["search_terms"]])

    or_block = " OR ".join(terms)
    query = f"(taxonomy_id:{taxid}) AND ({or_block})"
    if REVIEWED_ONLY:
        query += " AND (reviewed:true)"
    return query


def parse_uniprot_tsv(text: str) -> pd.DataFrame:
    if not text.strip():
        return pd.DataFrame()
    return pd.read_csv(io.StringIO(text), sep="\t")


def run_uniprot_search(query: str, fields: str) -> pd.DataFrame:
    params = {
        "query": query,
        "format": "tsv",
        "fields": fields,
        "size": str(UNIPROT_SIZE_PER_QUERY),
    }
    response = requests.get(
        UNIPROT_BASE_URL,
        params=params,
        headers=HEADERS,
        timeout=REQUEST_TIMEOUT_SECONDS,
    )
    if response.status_code != 200:
        preview = response.text[:500].replace("\n", " ").replace("\r", " ")
        raise RuntimeError(f"UniProt query failed with HTTP {response.status_code}: {preview}")
    time.sleep(POLITE_DELAY_SECONDS)
    return parse_uniprot_tsv(response.text)


def uniprot_search_with_fallback(query: str) -> pd.DataFrame:
    """Retry with a minimal field set if UniProt rejects an optional field."""
    try:
        return run_uniprot_search(query, UNIPROT_FIELDS_FULL)
    except Exception as exc:
        print(f"  Full UniProt field set failed; retrying minimal fields. Reason: {exc}")
        return run_uniprot_search(query, UNIPROT_FIELDS_MINIMAL)


def get_column(row, possible_names):
    for name in possible_names:
        if name in row.index:
            return row.get(name, "")
    return ""


def standardize_uniprot_hits(df: pd.DataFrame) -> pd.DataFrame:
    if df.empty:
        return pd.DataFrame(columns=[
            "accession", "entry_name", "reviewed", "protein_name", "gene_names",
            "organism_name", "organism_id", "length", "interpro", "pfam",
            "kegg", "function_cc",
        ])

    rows = []
    for _, row in df.iterrows():
        rows.append({
            "accession": get_column(row, ["Entry", "accession"]),
            "entry_name": get_column(row, ["Entry Name", "id"]),
            "reviewed": get_column(row, ["Reviewed", "reviewed"]),
            "protein_name": get_column(row, ["Protein names", "protein_name"]),
            "gene_names": get_column(row, ["Gene Names", "gene_names"]),
            "organism_name": get_column(row, ["Organism", "organism_name"]),
            "organism_id": get_column(row, ["Organism (ID)", "organism_id"]),
            "length": get_column(row, ["Length", "length"]),
            "interpro": get_column(row, ["InterPro", "xref_interpro"]),
            "pfam": get_column(row, ["Pfam", "xref_pfam"]),
            "kegg": get_column(row, ["KEGG", "xref_kegg"]),
            "function_cc": get_column(row, ["Function [CC]", "cc_function"]),
        })
    return pd.DataFrame(rows)


def score_uniprot_hit(row: pd.Series, target: dict) -> tuple[int, str]:
    score = 0
    flags = []

    gene_tokens = normalize_gene_tokens(row.get("gene_names", ""))
    gene_syns = {g.lower() for g in target["gene_synonyms"]}
    exact_gene_matches = sorted(gene_tokens.intersection(gene_syns))
    if exact_gene_matches:
        score += 5
        flags.append("gene=" + "/".join(exact_gene_matches))

    protein_txt = normalize_text(row.get("protein_name", ""))
    protein_hits = [k for k in target["protein_keywords"] if k.lower() in protein_txt]
    if protein_hits:
        score += 4
        flags.append("protein_name=" + "/".join(protein_hits[:3]))

    domain_txt = " ".join([
        normalize_text(row.get("interpro", "")),
        normalize_text(row.get("pfam", "")),
        normalize_text(row.get("function_cc", "")),
    ])
    domain_hits = [k for k in target["domain_keywords"] if k.lower() in domain_txt]
    if domain_hits:
        score += 2
        flags.append("domain_or_function=" + "/".join(domain_hits[:3]))

    reviewed_txt = normalize_text(row.get("reviewed", ""))
    if "reviewed" in reviewed_txt or "swiss-prot" in reviewed_txt:
        score += 1
        flags.append("reviewed")

    return score, "; ".join(flags)


def classify_status(score: int) -> str:
    if score >= STRONG_HIT_SCORE:
        return "strong UniProt hit"
    if score >= POSSIBLE_HIT_SCORE:
        return "possible UniProt hit"
    return "no convincing UniProt hit"


def search_uniprot_for_species_target(species: dict, target_name: str, target: dict) -> pd.DataFrame:
    taxid = species["taxid"]
    species_name = species["species"]

    all_hits = []

    queries = [
        ("gene_plus_terms", build_uniprot_query(taxid, target, broad=False)),
        ("broad_terms", build_uniprot_query(taxid, target, broad=True)),
    ]

    for label, query in queries:
        print(f"  UniProt {species_name} / {target_name} / {label}")
        try:
            df = uniprot_search_with_fallback(query)
        except Exception as exc:
            print(f"    Query failed: {exc}")
            df = pd.DataFrame()
        std = standardize_uniprot_hits(df)
        if not std.empty:
            std["query_label"] = label
            all_hits.append(std)

    if not all_hits:
        return pd.DataFrame()

    hits = pd.concat(all_hits, ignore_index=True)
    hits = hits.drop_duplicates("accession").copy()

    scores = hits.apply(lambda row: score_uniprot_hit(row, target), axis=1)
    hits["score"] = [x[0] for x in scores]
    hits["evidence_flags"] = [x[1] for x in scores]
    hits["status"] = hits["score"].apply(classify_status)
    hits["target"] = target_name
    hits["species"] = species_name
    hits["species_taxid"] = taxid

    hits["reviewed_rank"] = hits["reviewed"].astype(str).str.contains(
        "reviewed", case=False, na=False
    ).astype(int)

    hits = hits.sort_values(
        ["score", "reviewed_rank", "accession"],
        ascending=[False, False, True],
    )

    return hits.drop(columns=["reviewed_rank"])


def build_kegg_species_status() -> pd.DataFrame:
    genome_text = get_kegg_text(KEGG_ENDPOINTS["genomes"])
    genomes = parse_kegg_bacterial_genome_list(genome_text)

    ko_positive_org_codes = {}
    for target_name, target in TARGETS.items():
        ko = target["ko"]
        link_text = get_kegg_text(f"link/genes/{ko}")
        links = parse_kegg_ko_gene_links(link_text)
        if links.empty:
            ko_positive_org_codes[target_name] = set()
        else:
            ko_positive_org_codes[target_name] = set(links["gene_prefix"])

    rows = []
    for species in SPECIES_TO_CHECK:
        species_name = species["species"]
        species_lower = species_name.lower()
        matched = genomes[
            genomes["organism_name"]
            .str.lower()
            .str.contains(re.escape(species_lower), na=False)
        ].copy()

        matched_codes = set(matched["org_code"])

        for target_name, positives in ko_positive_org_codes.items():
            pos_codes = matched_codes.intersection(positives)
            rows.append({
                "species": species_name,
                "species_taxid": species["taxid"],
                "target": target_name,
                "target_ko": TARGETS[target_name]["ko"],
                "kegg_matching_genomes_n": len(matched),
                "kegg_matching_org_codes": ";".join(sorted(matched_codes)),
                "kegg_positive_genomes_n": len(pos_codes),
                "kegg_positive_org_codes": ";".join(sorted(pos_codes)),
                "kegg_positive_percent": (100 * len(pos_codes) / len(matched)) if len(matched) else None,
            })

    return pd.DataFrame(rows)


def build_validation_tables() -> tuple[pd.DataFrame, pd.DataFrame]:
    kegg_status = build_kegg_species_status()

    all_hits = []
    summary_rows = []

    for species in SPECIES_TO_CHECK:
        for target_name, target in TARGETS.items():
            hits = search_uniprot_for_species_target(species, target_name, target)

            if not hits.empty:
                all_hits.append(hits)
                top = hits.iloc[0]
                best_score = int(top["score"])
                best_status = classify_status(best_score)
                n_strong = int((hits["score"] >= STRONG_HIT_SCORE).sum())
                n_possible = int((hits["score"] >= POSSIBLE_HIT_SCORE).sum())

                summary_rows.append({
                    "species": species["species"],
                    "species_taxid": species["taxid"],
                    "target": target_name,
                    "target_ko": target["ko"],
                    "uniprot_status": best_status,
                    "uniprot_best_score": best_score,
                    "uniprot_hits_total": len(hits),
                    "uniprot_hits_possible_or_better": n_possible,
                    "uniprot_hits_strong": n_strong,
                    "best_accession": top.get("accession", ""),
                    "best_entry_name": top.get("entry_name", ""),
                    "best_reviewed": top.get("reviewed", ""),
                    "best_protein_name": top.get("protein_name", ""),
                    "best_gene_names": top.get("gene_names", ""),
                    "best_organism_name": top.get("organism_name", ""),
                    "best_length": top.get("length", ""),
                    "best_interpro": top.get("interpro", ""),
                    "best_pfam": top.get("pfam", ""),
                    "best_kegg_xref": top.get("kegg", ""),
                    "best_evidence_flags": top.get("evidence_flags", ""),
                })

            else:
                summary_rows.append({
                    "species": species["species"],
                    "species_taxid": species["taxid"],
                    "target": target_name,
                    "target_ko": target["ko"],
                    "uniprot_status": "no UniProt hit returned",
                    "uniprot_best_score": 0,
                    "uniprot_hits_total": 0,
                    "uniprot_hits_possible_or_better": 0,
                    "uniprot_hits_strong": 0,
                    "best_accession": "",
                    "best_entry_name": "",
                    "best_reviewed": "",
                    "best_protein_name": "",
                    "best_gene_names": "",
                    "best_organism_name": "",
                    "best_length": "",
                    "best_interpro": "",
                    "best_pfam": "",
                    "best_kegg_xref": "",
                    "best_evidence_flags": "",
                })

    summary = pd.DataFrame(summary_rows)
    summary = summary.merge(
        kegg_status,
        on=["species", "species_taxid", "target", "target_ko"],
        how="left",
    )

    def interpret(row):
        kegg_pos = row.get("kegg_positive_genomes_n", 0) or 0
        score = row.get("uniprot_best_score", 0) or 0

        if kegg_pos > 0 and score >= STRONG_HIT_SCORE:
            return "KEGG and UniProt support canonical ortholog"
        if kegg_pos == 0 and score >= STRONG_HIT_SCORE:
            return "KEGG KO absent but UniProt suggests possible homolog; inspect manually"
        if kegg_pos == 0 and score < POSSIBLE_HIT_SCORE:
            return "No KEGG KO and no convincing UniProt hit; candidate true absence/divergence"
        if kegg_pos > 0 and score < POSSIBLE_HIT_SCORE:
            return "KEGG KO present but UniProt text search weak; inspect KEGG-linked protein"
        return "Ambiguous; inspect manually"

    summary["interpretation"] = summary.apply(interpret, axis=1)

    if all_hits:
        all_hits_df = pd.concat(all_hits, ignore_index=True)
    else:
        all_hits_df = pd.DataFrame()

    return summary, all_hits_df


def export_results(out: Path, summary: pd.DataFrame, all_hits: pd.DataFrame):
    summary_csv = out / "actinomycetota_MnmEG_validation_summary.csv"
    hits_csv = out / "actinomycetota_MnmEG_validation_all_uniprot_hits.csv"
    xlsx_path = out / "actinomycetota_MnmEG_validation_results.xlsx"

    summary.to_csv(summary_csv, index=False, encoding="utf-8-sig")
    all_hits.to_csv(hits_csv, index=False, encoding="utf-8-sig")

    settings = pd.DataFrame([
        ["output_base_dir", str(output_base_folder())],
        ["cache_folder", str(cache_folder())],
        ["uniprot_endpoint", UNIPROT_BASE_URL],
        ["kegg_genome_endpoint", kegg_url(KEGG_ENDPOINTS["genomes"])],
        ["targets", json.dumps({k: v["ko"] for k, v in TARGETS.items()})],
        ["species", json.dumps(SPECIES_TO_CHECK)],
        ["reviewed_only", REVIEWED_ONLY],
        ["strong_hit_score", STRONG_HIT_SCORE],
        ["possible_hit_score", POSSIBLE_HIT_SCORE],
    ], columns=["setting", "value"])

    with pd.ExcelWriter(xlsx_path, engine="openpyxl") as writer:
        summary.to_excel(writer, sheet_name="summary", index=False)
        all_hits.to_excel(writer, sheet_name="all_uniprot_hits", index=False)
        settings.to_excel(writer, sheet_name="settings", index=False)

    print(f"Saved: {summary_csv}")
    print(f"Saved: {hits_csv}")
    print(f"Saved: {xlsx_path}")


def print_console_summary(summary: pd.DataFrame):
    print("\nValidation summary")
    print("=" * 80)

    cols = [
        "species", "target", "kegg_positive_genomes_n",
        "kegg_matching_genomes_n", "uniprot_status", "best_accession",
        "best_gene_names", "interpretation",
    ]

    view = summary[cols].copy()

    for _, row in view.iterrows():
        print(
            f"{row['species']} / {row['target']}: "
            f"KEGG {row['kegg_positive_genomes_n']}/{row['kegg_matching_genomes_n']} genomes; "
            f"UniProt = {row['uniprot_status']}; "
            f"best = {row['best_accession']} {row['best_gene_names']}; "
            f"{row['interpretation']}"
        )


def main():
    global OUTPUT_BASE_DIR
    OUTPUT_BASE_DIR = select_output_base_folder()
    out = make_output_folder()

    print(f"Output folder: {out}")
    print(f"Cache folder: {cache_folder()}")
    print(f"Python: {sys.version.split()[0]}")

    summary, all_hits = build_validation_tables()

    export_results(out, summary, all_hits)
    print_console_summary(summary)

    print("\nDone.")
    print("Review the summary sheet and the complete UniProtKB hit table for detailed results.")


if __name__ == "__main__":
    main()