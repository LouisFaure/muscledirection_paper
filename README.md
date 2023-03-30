# muscledirection_paper

This repository contains code for reproducing analysis of scRNAseq and bulkRNAseq from cultivate muscle cells, for the following study:

```
paper in publication process, citation coming soon!
```

## Prepare conda environment

```bash
mamba env create -f env.yml
mamba activate muscle
python -m ipykernel install --user --name muscle_sc --display-name "Muscle"
Rscript -e "IRkernel::installspec(name = 'muscle_sc', displayname = 'R Muscle )"
```
