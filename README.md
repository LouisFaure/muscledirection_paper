[![DOI](https://img.shields.io/badge/DOI-10.1038/s41467--023--38647--7-blue)](https://doi.org/10.1038/s41467-023-38647-7)
[![GEO](https://img.shields.io/badge/scRNAseq-download-green)](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE160098)
[![GEO](https://img.shields.io/badge/bulkRNAseq-download-green)](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE199093)

# muscledirection_paper

This repository contains code for reproducing analysis of scRNAseq and bulkRNAseq from cultivate muscle cells, for the following study:

```
Sunadome, K., Erickson, A.G., Kah, D. et al. 
Directionality of developing skeletal muscles is set by mechanical forces.
Nat Commun 14, 3060 (2023). https://doi.org/10.1038/s41467-023-38647-7
```

## Prepare conda environment

```bash
mamba env create -f env.yml
mamba activate muscle
python -m ipykernel install --user --name muscle_sc --display-name "Muscle"
Rscript -e "IRkernel::installspec(name = 'muscle_sc', displayname = 'R Muscle )"
```
