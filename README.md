# RIFRMap

**RIFRMap** is an interactive Shiny-based platform for exploring rifampicin (RIF) resistance mutations across bacterial species. It integrates multiple public and curated datasets into a unified framework, enabling systematic analysis and comparison of resistance-associated mutations.

The platform standardises mutation nomenclature across heterogeneous sources and provides tools to examine:
- mutation distribution across species  
- frequency and prevalence of variants  
- overlap and consistency between data sources  

RIFRMap provides dynamic visualisations and flexible filtering capabilities, enabling users to explore mutations by species, supporting evidence, and dataset-specific criteria informed by the literature.
This tool is designed to support research into antimicrobial resistance by providing a comprehensive and accessible view of rifampicin resistance mutations.

## Project Structure

- **raw_data/**  
  Contains all original downloaded data files.

- **processed_data/**  
  Stores cleaned and processed datasets used for analysis and the app.

### Data Processing Scripts
- `install_deps.R` – Install required dependencies  
- `update_data.R` – Download and update raw datasets  
- `process_data.R` – Clean and process datasets  
- `process_refseqs.R` – Generate reference sequence coordinates (`coordinates.RData`)  
- `helpers.R` – Utility/helper functions  

### Shiny App
- `ui.R` – User interface  
- `server.R` – Server logic  
- `global.R` – Shared objects and data loading  

### Frontend
- **www/**  
  Custom CSS, HTML, and JavaScript for the Shiny app  
  
## Run the App

```r
shiny::runApp()  
```

## Data Sources

### Comprehensive Antibiotic Resistance Database (CARD)

The **CARD** database is a hand-curated resource covering the full range of antimicrobial resistance genes (ARGs) [1](https://doi.org/10.1093/nar/gkac920).
To be included, each AMR determinant must meet all of the following:

-   Reported in a peer-reviewed publication
-   Associated with a GenBank sequence
-   Demonstrates increased minimum inhibitory concentration (MIC) compared to control strains through experimental evidence

For *Mycobacterium tuberculosis*, additional curated mutation data were integrated from the **CRyPTIC**, **ReSeqTB**, and **WHO 2023** resistance catalogs.


### CRyPTIC (Comprehensive Resistance Prediction for Tuberculosis: an International Consortium)

The **CRyPTIC Consortium** profiled MICs for 13 anti-tuberculosis drugs—including rifampicin—in over 15,000 *M. tuberculosis* isolates worldwide [2](https://doi.org/10.1038/s41467-023-44325-5).
Using multivariable linear mixed models and whole-genome sequencing, they identified:

-   **492** variants associated with increased MICs (resistance)
-   **91** variants linked to decreased MICs (hypersusceptibility)

**RIFRMap** allows interactive filtering of this dataset by:

-   Effect size
-   Variant count
-   Benjamini–Hochberg adjusted p-value

### DRAGdb (Drug Resistance Associated Genes Database)

**DRAGdb** is a manually curated repository of mutational data for drug resistance–associated genes (DRAGs) across:

-   **ESKAPE pathogens**:
    *Enterococcus faecium*, *Staphylococcus aureus*, *Klebsiella pneumoniae*, *Acinetobacter baumannii*, *Pseudomonas aeruginosa*, and *Enterobacter* spp.
-   Other bacteria, with emphasis on *Mycobacterium tuberculosis*

The database includes **PROVEAN scores** to estimate the predicted functional effects of mutations [3](https://doi.org/10.1038/s41598-020-65766-8).


### Yang-RIF (2023)

The **Yang-RIF dataset** originates from a saturation mutagenesis study by Yang et al., targeting the rifampicin-binding region in *E. coli* RNA polymerase.
Researchers constructed **760** single-residue *rpoB* mutants across positions **510–537** and **563–572** to systematically assess resistance phenotypes [4](https://doi.org/10.1038/s41586-023-06495-6)

Within **RIFRMap**, users can filter the dataset based on mutation frequency changes under rifampicin selection, using thresholds on log₂(fold-change) to identify enriched variants.


### Curated Literature Data

This dataset was manually compiled by our lab through a comprehensive review of the scientific literature [5](https://doi.org/10.1073/pnas.2424307122).
Relevant studies were identified using Google Scholar with search terms such as:

-   `"rifampicin resistance"`
-   `"RNA polymerase subunit beta"`
-   `"point mutation"`
-   `"amino acid substitution"`

From these studies, they extracted single amino acid substitutions in *rpoB* reported to confer rifampicin resistance.

In **RIFRMap**, users can filter this dataset based on:

-   Number of independent supporting studies
-   Number of bacterial species in which the mutation has been reported


## Citation

[1]: <https://doi.org/10.1093/nar/gkac920> Alcock, B. P., Raphenya, A. R., Tyrrell, G. J., & McArthur, A. G. (2023). CARD 2023: expanded curation, support for machine learning, and resistome prediction at the Comprehensive Antibiotic Resistance Database. Nucleic Acids Research, 51, D690–D699.

[2]: <https://doi.org/10.1038/s41467-023-44325-5> The CRyPTIC Consortium. (2024, January). Quantitative measurement of antibiotic resistance in Mycobacterium tuberculosis reveals genetic determinants of resistance and susceptibility in a target gene approach. Nature Communications, 15, 488.

[3]: <https://doi.org/10.1038/s41598-020-65766-8> Ghosh, A., N, S., & Saha, S. (2020). Survey of drug resistance associated gene mutations in *Mycobacterium tuberculosis*, ESKAPE and other bacterial species. *Scientific Reports*, 10(1), 8957.

[4]: <https://doi.org/10.1038/s41586-023-06495-6> Yang, K.B., Cameranesi, M., Gowder, M., et al. (2023). High-resolution landscape of an antibiotic binding site. *Nature*, 622, 180–187.

[5]: <https://doi.org/10.1073/pnas.2424307122> Bolourchi, N., Brown, C.R.P., Letten, A.D., & Engelstädter, J. (2025). Evolution and evolvability of rifampicin resistance across the bacterial tree of life. *Proceedings of the National Academy of Sciences of the United States of America*, 122(31), e2424307122.
