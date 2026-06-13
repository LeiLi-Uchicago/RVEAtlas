# FLUExplorer Methods and Information

## Overview

FLUExplorer is an interactive Shiny application for exploring amino acid variation across Influenza genomes. It is designed to help users examine subtype-specific sequence diversity, compare genetic clades, inspect variation at individual positions, and visualize gene-wide patterns such as conservation and fixed mutations.

The current version of the app is loaded with curated human Influenza A subtype H1N1, H3N2, human Influenza B Yam and Vic lineage amino acid usage data organized by gene, genetic clades, year, and year-month.

FLUExplorer is intended for researchers, bioinformaticians, genomic epidemiologists, and other users who want to investigate FLU evolutionary patterns, mutation dynamics, and lineage-specific amino acid changes through an accessible visual interface.

---
## Nextclade References

Data were sourced from GISAID, with all sequences annotated via Nextclade 3. The reference datasets used for each lineage are listed below:

#### H1N1 pdm09

| Segment | Nextclade reference dataset |
|---|---|
| HA | `flu_h1n1pdm_ha` |
| NA | `flu_h1n1pdm_na` |
| MP | `nextstrain/flu/h1n1pdm/mp` |
| NP | `nextstrain/flu/h1n1pdm/np` |
| NS | `nextstrain/flu/h1n1pdm/ns` |
| PA | `nextstrain/flu/h1n1pdm/pa` |
| PB1 | `nextstrain/flu/h1n1pdm/pb1` |
| PB2 | `nextstrain/flu/h1n1pdm/pb2` |

#### H1N1 seasonal

| Segment | Nextclade reference dataset |
|---|---|
| HA | `flu_h1n1_ha` |
| NA | `flu_h1n1_na` |
| MP | `flu_h1n1_mp` |
| NP | `flu_h1n1_np` |
| NS | `flu_h1n1_ns` |
| PA | `flu_h1n1_pa` |
| PB1 | `flu_h1n1_pb1` |
| PB2 | `flu_h1n1_pb2` |

#### H3N2

| Segment | Nextclade reference dataset |
|---|---|
| HA | `nextstrain/flu/h3n2/ha/EPI1857216` |
| NA | `nextstrain/flu/h3n2/na/EPI1857215` |
| MP | `nextstrain/flu/h3n2/mp` |
| NP | `nextstrain/flu/h3n2/np` |
| NS | `nextstrain/flu/h3n2/ns` |
| PA | `nextstrain/flu/h3n2/pa` |
| PB1 | `nextstrain/flu/h3n2/pb1` |
| PB2 | `nextstrain/flu/h3n2/pb2` |

#### B/Victoria

| Segment | Nextclade reference dataset |
|---|---|
| HA | `nextstrain/flu/vic/ha/KX058884` |
| NA | `nextstrain/flu/vic/na/CY073894` |
| MP | `nextstrain/flu/vic/mp` |
| NP | `nextstrain/flu/vic/np` |
| NS | `nextstrain/flu/vic/ns` |
| PA | `nextstrain/flu/vic/pa` |
| PB1 | `nextstrain/flu/vic/pb1` |
| PB2 | `nextstrain/flu/vic/pb2` |

#### B/Yamagata

| Segment | Nextclade reference dataset |
|---|---|
| HA | `nextstrain/flu/yam/ha/JN993010` |
| NA | `nextstrain/flu/b/na/CY073894` |
| MP | `nextstrain/flu/b/mp` |
| NP | `nextstrain/flu/b/np` |
| NS | `nextstrain/flu/b/ns` |
| PA | `nextstrain/flu/b/pa` |
| PB1 | `nextstrain/flu/b/pb1` |
| PB2 | `nextstrain/flu/b/pb2` |


---

## Update Log

### 2026-05-25
- Fix a BUG in count for HA protein

### 2026-05-12

- Update clade annotation and codon.

### 2026-05-05

- Added the Genetic Clade tab. Users can select a subtype-specific clade annotation, search ranked clade/group choices, and review clade-level summary statistics.

### 2026-04-28

- Added DuckDB-backed usage table loading to reduce memory pressure from large amino acid and nucleotide tables. When DuckDB is available, FLUExplorer now queries only the selected subtype, gene, position, grouping, and time range instead of loading full tables into memory.
- Improved the Single Position Explorer time filter by replacing the draggable Year-Month range slider with explicit Start and End selectors.
- Fixed Year-Month grouping so plots and tables display real `YYYY-MM` groups instead of collapsing into `Unknown`.
- Fixed a Year-Month plotting error that could appear when a selected position had sparse or special time values.

### 2026-04-10

- Updated UI.

### 2026-03-06

- 1st version on-line. Using NextStrain data (Downloaded on 2025).
