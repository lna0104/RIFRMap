# RIFRMap

An interactive platform for exploring rifampicin resistance mutations in the bacterial *rpoB* gene.

## Summary

RIFRMap is an interactive Shiny application for exploring rifampicin (RIF) resistance mutations in the DNA-directed RNA polymerase β-subunit gene, *rpoB*, across bacterial species. It integrates seven public and curated resistance datasets — CARD, CRyPTIC, DRAGdb, Yang-RIF, RifReport, PointFinder, and the NCBI Reference Gene Catalog — into a unified framework by standardising all mutation positions onto a single canonical numbering scheme, the amino-acid position in *Escherichia coli* MG1655 *rpoB*. This standardisation enables direct comparison of resistance-associated variants reported across heterogeneous sources and bacterial taxa. RIFRMap provides dynamic visualisation, filtering, and cross-dataset comparison tools to support research of antimicrobial resistance.


## Installation & Usage

RIFRMap requires a working R installation and the R package dependencies listed in `code/install_deps.R`.

```r
# Install required dependencies
source("code/install_deps.R")

# Launch the application
shiny::runApp()
```

## Implementation

RIFRMap is implemented as an R Shiny application, organised into a data-processing pipeline and a Shiny front end. The pipeline scripts live under `code/`; the Shiny app's entrypoints (`ui.R`, `server.R`, `global.R`) and its front-end assets (`www/`) sit at the repo root, alongside `raw_data/`, `processed_data/`, `output/`, and `plots/`. Pipeline scripts are run with the working directory set to the repo root (e.g. `source("code/update_all.R")`); the Shiny app's working directory is also the repo root at runtime, set automatically by `shiny::runApp()`.

**Data.** Raw records retrieved from each source are stored under `raw_data/`; the standardised, analysis-ready datasets produced from them are stored under `processed_data/`.

**Data-processing pipeline.** A set of R scripts retrieve, validate, and standardise mutation records from each of the seven sources onto the canonical *E. coli* numbering scheme. Refreshing data is a two-step, validate-then-release process: `source("code/validate_all.R")` downloads/processes everything and writes a validation report without publishing anything, then `source("code/update_all.R")` re-validates and, only if the entire run is clean, publishes it.

- `code/install_deps.R` — installs required R package dependencies
- `code/update_data.R` — downloads and updates raw datasets from their original sources
- `code/process_refseqs.R` — generates and maintains the reference-sequence coordinate map used for cross-species standardisation
- `code/process_data.R` — cleans, filters, and standardises each raw dataset onto the canonical numbering scheme
- `code/pipeline_run.R` — shared download/process/validate core used by both `validate_all.R` and `update_all.R`
- `code/validate_all.R` — Stage 1: validates a full refresh without publishing anything
- `code/update_all.R` — Stage 2: re-validates, then publishes if and only if the whole run is clean
- `code/update_report.R` — validation checks, plus the validation/release reports and changelog
- `code/helpers.R` — shared utility functions used throughout the pipeline

**Shiny application.** The processed data are loaded once at application start-up and served through:

- `ui.R` — user interface
- `server.R` — server logic and reactive behaviour
- `global.R` — shared objects, themes, and data loading
- `www/` — front-end assets (CSS, HTML, JavaScript)

## Data Sources

### 1. Comprehensive Antibiotic Resistance Database (CARD)

The **CARD** database is a hand-curated resource covering the full range of antimicrobial resistance genes (ARGs) [1]. To be included, each AMR determinant must meet all of the following:

- Reported in a peer-reviewed publication
- Associated with a GenBank sequence
- Demonstrates increased minimum inhibitory concentration (MIC) compared to control strains through experimental evidence


For *Mycobacterium tuberculosis*, where experimental antimicrobial susceptibility testing is particularly challenging, CARD employed a likelihood‑based framework that infers resistance associations from single‑nucleotide polymorphisms (SNPs) and observed phenotypic resistance across large collections of sequenced genomes.

### 2. CRyPTIC (Comprehensive Resistance Prediction for Tuberculosis: an International Consortium)

The **CRyPTIC Consortium** profiled MICs for 13 anti-tuberculosis drugs — including rifampicin — in over 15,000 *M. tuberculosis* isolates worldwide [2]. Using multivariable linear mixed models and whole-genome sequencing, they identified:

- **492** variants associated with increased MICs (resistance)
- **91** variants linked to decreased MICs (hypersusceptibility)

RIFRMap allows interactive filtering of this dataset by effect size, variant count, and Benjamini–Hochberg adjusted p-value.

### 3. DRAGdb (Drug Resistance Associated Genes Database)

**DRAGdb** is a manually curated repository of mutational data for drug resistance–associated genes (DRAGs) across:

- **ESKAPE pathogens**: *Enterococcus faecium*, *Staphylococcus aureus*, *Klebsiella pneumoniae*, *Acinetobacter baumannii*, *Pseudomonas aeruginosa*, and *Enterobacter* spp.
- Other bacteria, with emphasis on *Mycobacterium tuberculosis*

The database includes **PROVEAN scores** to estimate the predicted functional effects of mutations [3].

### 4. Yang-RIF (2023)

The **Yang-RIF dataset** originates from a saturation mutagenesis study by Yang et al., targeting the rifampicin-binding region in *E. coli* RNA polymerase. Researchers constructed **760** single-residue *rpoB* mutants across positions **510–537** and **563–572** to systematically assess resistance phenotypes [4]. Within RIFRMap, users can filter the dataset based on mutation frequency changes under rifampicin selection, using thresholds on log₂(fold-change) to identify enriched variants.

### 5. RifReport

RifReport was compiled through a comprehensive review of the scientific literature [5]. Relevant studies were identified using Google Scholar with search terms including `"rifampicin resistance"`, `"RNA polymerase subunit beta"`, `"point mutation"`, and `"amino acid substitution"`. From these studies, single amino acid substitutions in *rpoB* reported to confer rifampicin resistance were extracted. Within RIFRMap, users can filter this dataset based on the number of independent supporting studies and the number of bacterial species in which a mutation has been reported.

### 6. PointFinder

**PointFinder** is a curated database of chromosomal point mutations associated with antimicrobial resistance, spanning multiple clinically relevant bacterial pathogens, and forms part of the broader **ResFinder** framework, which combines acquired resistance genes, chromosomal mutations, and analysis tools to support antimicrobial resistance detection [6, 7]. Within RIFRMap, this dataset is filtered to chromosomal point mutations in *rpoB* associated with rifampicin resistance.

### 7. Reference Gene Catalog

The **Reference Gene Catalog** is a curated resource developed as part of NCBI's Pathogen Detection project and underlies NCBI's AMRFinderPlus tool, providing a standardised set of antimicrobial resistance genes and proteins maintained in the Bacterial Antimicrobial Resistance Reference Gene Database (BioProject PRJNA313047) [8]. Within RIFRMap, this dataset is restricted to mutations in the *rpoB* gene family associated with the rifampin subclass.

## References

[1]: <https://doi.org/10.1093/nar/gkac920> Alcock, B. P., Raphenya, A. R., Tyrrell, G. J., & McArthur, A. G. (2023). CARD 2023: expanded curation, support for machine learning, and resistome prediction at the Comprehensive Antibiotic Resistance Database. *Nucleic Acids Research*, 51, D690–D699.

[2]: <https://doi.org/10.1038/s41467-023-44325-5> The CRyPTIC Consortium. (2024). Quantitative measurement of antibiotic resistance in *Mycobacterium tuberculosis* reveals genetic determinants of resistance and susceptibility in a target gene approach. *Nature Communications*, 15, 488.

[3]: <https://doi.org/10.1038/s41598-020-65766-8> Ghosh, A., N, S., & Saha, S. (2020). Survey of drug resistance associated gene mutations in *Mycobacterium tuberculosis*, ESKAPE and other bacterial species. *Scientific Reports*, 10(1), 8957.

[4]: <https://doi.org/10.1038/s41586-023-06495-6> Yang, K. B., Cameranesi, M., Gowder, M., et al. (2023). High-resolution landscape of an antibiotic binding site. *Nature*, 622, 180–187.

[5]: <https://doi.org/10.1073/pnas.2424307122> Bolourchi, N., Brown, C. R. P., Letten, A. D., & Engelstädter, J. (2025). Evolution and evolvability of rifampicin resistance across the bacterial tree of life. *Proceedings of the National Academy of Sciences of the United States of America*, 122(31), e2424307122.

[6]: <https://doi.org/10.1093/jac/dkx217> Zankari, E., Allesøe, R., Joensen, K. G., Cavaco, L. M., Lund, O., & Aarestrup, F. M. (2017). PointFinder: a novel web tool for WGS-based detection of antimicrobial resistance associated with chromosomal point mutations in bacterial pathogens. *Journal of Antimicrobial Chemotherapy*, 72(10), 2764–2768.

[7]: <https://doi.org/10.1099/mgen.0.000748> Florensa, A. F., Kaas, R. S., Clausen, P. T. L. C., Aytan-Aktug, D., & Aarestrup, F. M. (2022). ResFinder – an open online resource for identification of antimicrobial resistance genes in next-generation sequencing data and prediction of phenotypes from genotypes. *Microbial Genomics*, 8(1), 000748.

[8]: <https://doi.org/10.1038/s41598-021-91456-0> Feldgarden, M., Brover, V., Gonzalez-Escalona, N., et al. (2021). AMRFinderPlus and the Reference Gene Catalog facilitate examination of the genomic links among antimicrobial resistance, stress response, and virulence. *Scientific Reports*, 11, 12728.

## Data Availability & Reproducibility

Raw, unprocessed data as retrieved from each source are stored under `raw_data/`; cleaned and standardised datasets used by the application are stored under `processed_data/`. Retrieval and processing of all seven datasets are performed by versioned, auditable R scripts, allowing the underlying data to be refreshed and reprocessed from original sources at any time. The deployed application itself makes no external network requests at runtime: all data are precomputed and bundled with the application.

## How to Cite

RIFRMap does not yet have an associated peer-reviewed publication or archival DOI. In the interim, please cite this software by reference to its repository, [github.com/lna0104/RIFRMap](https://github.com/lna0104/RIFRMap), including the version or commit hash used. A formal citation (publication and/or Zenodo DOI) will be added here upon publication.

## License

No license has yet been assigned to this repository. Licensing terms will be specified by the maintainer prior to publication/release.