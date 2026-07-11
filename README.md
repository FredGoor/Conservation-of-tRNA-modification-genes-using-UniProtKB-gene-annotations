## tRNA modifier annotation validation using UniProtKB and KEGG

### Overview

`Conservation of tRNA modification genes using UniProtKB gene annotations.py` evaluates the presence of selected tRNA-modification enzymes in a defined set of bacterial species using complementary information from KEGG and UniProtKB.

The script was designed to investigate cases in which a KEGG Orthology assignment is absent or unexpectedly rare in a taxonomic group. For each selected species, it:

1. identifies KEGG genome representatives matching the species name;
2. determines whether these genomes contain the specified KEGG Orthology identifiers;
3. searches UniProtKB for proteins matching relevant gene names, protein names, synonyms and functional annotations;
4. assigns a simple evidence score to each UniProtKB hit; and
5. generates a combined summary of KEGG and UniProtKB evidence.

The analysis does not download genome sequences and does not perform sequence-based homology searches such as BLAST or HMMER. It should therefore be interpreted as an annotation-based validation rather than definitive evidence of gene presence or absence.

### Requirements

The script requires Python 3.9 or later and the following packages:

```bash
pip install pandas requests openpyxl
```

A graphical desktop environment is required for the output-directory selection dialog.

### Usage

Run the script from a Python environment or command line:

```bash
python tRNA_modifier_annotation_validation_UniProt_KEGG.py
```

At startup, a folder-selection dialog prompts the user to choose an output directory. The script then creates a timestamped analysis folder within the selected directory.

The list of species, target enzymes, KEGG Orthology identifiers and UniProtKB search terms can be edited in the configuration section near the beginning of the script.

### Default targets

The default analysis includes:

* MiaA — KEGG Orthology `K00791`
* MnmE/TrmE — KEGG Orthology `K03650`
* MnmG/GidA — KEGG Orthology `K03495`

MiaA may be used as a reference target when the primary objective is to assess MnmE and MnmG annotation coverage.

### UniProtKB evidence scoring

UniProtKB records are scored using evidence from:

* exact matches to known gene names or synonyms;
* matches to expected protein names;
* matches to domain or functional annotations; and
* reviewed UniProtKB status.

Hits are classified as strong, possible or unconvincing according to configurable score thresholds. These scores are intended to prioritize records for inspection and should not be treated as formal homology statistics.

### Output files

The script generates:

* `tRNA_modifier_annotation_validation_results.xlsx`
  Excel workbook containing the summary table, all UniProtKB hits and analysis settings.

* `tRNA_modifier_annotation_validation_summary.csv`
  Combined KEGG and UniProtKB summary for each species and target.

* `tRNA_modifier_annotation_validation_all_uniprot_hits.csv`
  Full list of retrieved and scored UniProtKB records.

* `_cache/`
  Cached KEGG responses used to reduce repeated downloads.

### Interpretation

Typical result patterns include:

* **KEGG positive and strong UniProtKB evidence:** annotation is supported by both resources.
* **KEGG negative and strong UniProtKB evidence:** the apparent absence may reflect incomplete KEGG Orthology assignment or annotation differences.
* **KEGG negative and weak UniProtKB evidence:** the target may be absent, highly divergent or poorly annotated.
* **KEGG positive and weak UniProtKB evidence:** the KEGG-linked record should be inspected directly.

For uncertain cases, sequence-based validation using BLAST, HMMER or profile-domain searches is recommended.
