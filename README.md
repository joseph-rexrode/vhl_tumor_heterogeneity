# Tumor Heterogeneity in VHL-specific ccRCCs - Repository

This repository contains all code and relevant functions used in this project, which was undergone to determine whether highly variable genes exhibit different expression variability levels across distinct VHL ccRCC populations.

## Data

All data was obtained from a publicly available dataset from GEO, namely <a href="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE254444">GSE254444</a>. 

Data was sourced from 6 samples via a nephrectomy from a single VHL patient, 5 of which are
distinct tumors, and the last of which is a normal epithelial sample. In total, around 20,000 genes and 11,000 cells were present in this dataset.


## Distribution Modeling

A Negative Binomial distribution was used to model this data, as is customary for scRNA-seq data, due to the large number of zeros present. 

$P(k) = \int_{0}^{\infty}{Poisson(k | \lambda) \cdot Gamma(\lambda | \alpha, \theta) d\lambda}$

This distribution was modeled for each unique tumor sample, as well as the normal epithelial sample, and visually compared. Dispersion parameters were also estimated across all models and compared to determine significant heterogeneity between populations.

## Sources

[1] Find variable features — FindVariableFeatures. (2026). Satijalab.org.
      <a>https://satijalab.org/seurat/reference/findvariablefeatures</a>

[2] Gerlinger, M. et al. (2012). Intratumor Heterogeneity and Branched Evolution Revealed by
Multiregion Sequencing. New England Journal of Medicine, 366(10), 883–892.
<a>https://doi.org/10.1056/nejmoa1113205</a>

[3] Maria Scotti, G. (2024 ). GEO Accession viewer. Nih.gov.
<a>https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE254444</a>

[4] Rowe, I. et al. (2026). Distinct genomic, microenvironmental, and nephron signatures in VHL
kidney cysts and tumors. Scientific Reports, 16(1).
<a>https://doi.org/10.1038/s41598-025-34381-w</a>

[5] Sarkar, A., & Stephens, M. (2021). Separating measurement and expression models clarifies
confusion in single-cell RNA sequencing analysis. Nature Genetics, 53( 6), 770– 777.
<a>https://doi.org/10.1038/s41588-021-00873-4</a>

[6] Sato, Y. et al. (2013). Integrated molecular analysis of clear-cell renal cell carcinoma. Nature Genetics, 45(8), 860–867. <a>https://doi.org/10.1038/ng.2699</a>

[7] Wu, H., & Li, Z. (2021). Analysis of single-cell RNA-seq data.
<a>https://www.haowulab.org/teaching/ENAR2021/ENAR2021_scRNA-seq_lecture2.pdf</a>