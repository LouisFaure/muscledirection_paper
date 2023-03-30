# scRNAseq analysis

## Downloading BAM files and converting back to fastq

```bash
wget -O bam2fastq https://github.com/10XGenomics/bamtofastq/releases/download/v1.4.1/bamtofastq_linux
chmod +x bam2fastq
while read srr; do
        ffq --ftp $srr | jq -r '.[] | .url' | head -n1 | xargs curl -O
	name=$(basename $name .1)
	./bam2fastq --nthreads=20 $name $srr
	mv $srr/*/*.fastq.gz $srr/
	rm -rf $srr/*/
done <SRR_Acc_List.txt
```

## Preparing cellranger reference

```bash
wget http://cf.10xgenomics.com/supp/cell-exp/refdata-cellranger-mm10-1.2.0.tar.gz
tar xvf refdata-cellranger-mm10-1.2.0.tar.gz
```

## Run cellranger

```bash
while read srr; do
        cellranger count --id=$srr --fastqs=$srr --transcriptome=refdata-cellranger-mm10-1.2.0 --sample=$srr
done <SRR_Acc_List.txt
```

## CellBender

*Disclaimer*: this part of the analysis is not fully reproducible, possibly due to different package versions,  slight shifts in training losses will lead to differing count matrices, it is recommended to reuse the cleaned matrices from GEO for reproducing analysis.

```bash
mamba create -n cellbender pytorch==1.5.1 python=3.8 cudatoolkit=10.1 -c pytorch
mamba activate cellbender
git clone https://github.com/broadinstitute/CellBender
cd CellBender
git checkout v0.1.0
pip install -e .
cd ..

cellbender remove-background \
                 --epochs 150 \
                 --input raw_feature_bc_matrix.h5 \
                 --output cellbended.h5 \
                 --expected-cells 2000 \
                 --total-droplets-included 10000 \
                 --cuda
```

