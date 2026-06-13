# Extended Gene Coexpression Analysis Pipeline

> **Status: early scaffold** — structure and config contract only; no analysis logic implemented yet.

Takes any preprocessed 10x single-cell/single-nucleus dataset (currently Seurat objects)
and produces context-robust co-expression modules with biological interpretation,
plus a gene-of-interest (GOI) lookup resource.

## Install

```r
# install.packages("devtools")
devtools::install_github("noborilab/coexpression_arabidopsis")
```

## CLI usage

```bash
# Validate a config file (no pipeline run)
Rscript inst/scripts/run_pipeline.R --config config/my_config.yaml --validate-only

# Run the full pipeline
Rscript inst/scripts/run_pipeline.R --config config/my_config.yaml
```

## Intended usage

1. Copy `config/example_config.yaml` → `config/my_config.yaml`.
2. Fill in your dataset path, stratum variables, estimation mode, etc.
3. Run the CLI entrypoint above.
4. Results appear in the output directory specified in the config.

## First analysis targets

### Dev atlas — pseudobulk
**Lee, Illouz-Eliaz, Nobori et al. 2025** *Nature Plants*  
GEO [GSE226097](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE226097)  
~34 k genes, clusters nested in 10 organs/stages (seed, seedling, flower, rosette, stem, silique).
Pseudobulk available; cell-level objects archived off-machine.

### Pathogen multiome / PRIMER — SingleCellGGM
**Nobori et al. 2024** *Nature*  
GEO [GSE226826](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE226826) + [GSE248054](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE248054)  
Mock / DC3000 / AvrRpt2 / AvrRpm1, ~65 k nuclei, 15 samples.
SingleCellGGM already run casually; needs review before use as pipeline core (Phase 0).

## Design principles

- **Input adapter boundary** (`R/adapter_seurat.R`): the core never assumes Seurat forever.
  All Seurat-specific logic is isolated in one file. Future adapters (AnnData,
  raw-count re-normalization from GEO) implement the same interface without touching core logic.

- **Two first-class estimation modes**: pseudobulk (`R/estimate_pseudobulk.R`) and
  SingleCellGGM (`R/estimate_singlecellggm.R`).

- **Optional robustness layer** (`R/robustness.R`): cross-stratum R_score + cross-dataset
  replication. Enabled/disabled per config; skipped when data doesn't support it.

- **Shared output/interpretation layer** (`R/interpret.R`, `R/goi_lookup.R`): modules,
  preservation, hierarchy/hubs, GO + TF + curated-anchor enrichment, GOI lookup.

- **Swappable domain plugins** (`plugins/`): TF lists, curated gene sets, anchor sets —
  all referenced by path in config, never hardcoded.

## Package structure

```
R/
  adapter_seurat.R          # Input adapter: Seurat → core abstraction (only Seurat-aware file)
  estimate_pseudobulk.R     # Estimation mode 1: pseudobulk marginal correlation
  estimate_singlecellggm.R  # Estimation mode 2: cell-level graphical Gaussian model
  robustness.R              # Optional cross-context robustness (R_score, replication)
  interpret.R               # Modules, preservation, hubs, GO/TF/anchor enrichment
  goi_lookup.R              # Gene-of-interest lookup resource
inst/scripts/
  run_pipeline.R            # Config-driven CLI entrypoint
config/
  example_config.yaml       # Annotated example config — copy and edit
docs/
  BACKGROUND.md             # Scientific background and design substrate
  ARCHITECTURE.md           # Architecture design
plugins/                    # Domain-specific gene lists (swappable, not committed)
tests/testthat/             # Test harness
```

## Background

See [docs/BACKGROUND.md](docs/BACKGROUND.md) for scientific context and design rationale.
See [docs/ARCHITECTURE.md](docs/ARCHITECTURE.md) for the pipeline architecture.

## License

MIT © Nobori Lab
