# RIF Resistance Database

This repository contains an end-to-end R workflow for compiling, processing, analysing, and visualising rifampicin (RIF) resistance mutations from multiple public and curated sources. The project integrates heterogeneous datasets, standardises mutation nomenclature, builds a harmonised RIF mutation database, and provides an interactive Shiny application for exploration.

## Structure

- data/raw_data - Contains all original downloaded data files.
- data/RIFR_data - Stores cleaned datasets used in the analysis.
- code - The main analysis pipeline scripts.

## Data Sources

1.  CARD (Comprehensive Antibiotic Resistance Database)

    -   CARD includes only experimentally validated AMR determinants: genes must be reported in a peer-reviewed publication, have a sequence in GenBank, and show clear experimental evidence of increased MIC over controls [1](https://pubmed.ncbi.nlm.nih.gov/36263822/) (Curated-R).

    -   For *M. tuberculosis*, additional curated mutation data were incorporated from several major catalogs:

        -   CRyPTIC (CRyPTIC-R): mutation data acquired from the [CRyPTIC catalog](https://card.mcmaster.ca/ontology/47971).

        -   ReSeqTB (ReSeqTB-High, ReSeqTB-Minimal, ReSeqTB-Moderate): mutation data acquired from the [ReSeqTB catalog](https://card.mcmaster.ca/ontology/47970), evaluated as conferring resistance (Minimal, Moderate, High), not conferring resistance (None), or Indeterminate.

        -   WHO (WHO-R): mutation data acquired from the [WHO 2023 catalog](https://card.mcmaster.ca/ontology/47972).

    -   Raw CARD datasets used in this project were downloaded from: <https://card.mcmaster.ca/download>.


2.  CRyPTIC

    -   The Comprehensive Resistance Prediction for Tuberculosis: an International Consortium (CRyPTIC) is a Mycobacterium tuberculosis AMR mutation catalog published in 2022.
    -   The study performed WGS and determined the MICs of 13 drugs for 15211 Mtb isolates selected from patient samples gathered from 23 countries over five continents using a previously validated microtiter plate [2](https://doi.org/10.1038/s41467-023-44325-5). Using a censored interval regression model and correcting for multiple testing, the authors identified mutations significantly associated with elevated MICs.
    -   The raw dataset used in this project was downloaded from Supplementary Data 3 of the publication.
    -   For this pipeline, the processed dataset (CRyPTIC_processed.csv) includes only rifampicin-associated mutations that passed the Benjamini–Hochberg false discovery rate threshold of \< 0.05.

3.  DRAGdb

    -   DRAGdb is a database of mutational information of DRAGs across MTB clinical strains, ESKAPE bacteria and other pathogenic and non-pathogenic bacterial species with special reference to MTB H37Rv.
    
    
4. Yang's data

    -   A high-throughput screen using multiplex automated genome engineering to generate all possible single–amino acid substitutions within the RIF-binding region of rpoB.


## Processed Data

- All processed datasets retain only single-variant rifampicin-associated mutations.


## Citation

[1]: <https://pubmed.ncbi.nlm.nih.gov/36263822/> Alcock, B. P., Raphenya, A. R., Tyrrell, G. J., & McArthur, A. G. (2023). CARD 2023: expanded curation, support for machine learning, and resistome prediction at the Comprehensive Antibiotic Resistance Database. Nucleic Acids Research, 51, D690–D699.

[2]: <https://doi.org/10.1038/s41467-023-44325-5> The CRyPTIC Consortium. (2024, January). Quantitative measurement of antibiotic resistance in Mycobacterium tuberculosis reveals genetic determinants of resistance and susceptibility in a target gene approach. Nature Communications, 15, 488.
