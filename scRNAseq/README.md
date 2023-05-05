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

```bash
mamba create -n cellbender python=3.6 pytorch==1.4.0 -c pytorch -y
mamba activate cellbender

pip install -r cellbender_reqs_2020-02-25.txt
git clone https://github.com/broadinstitute/CellBender
cd CellBender && git checkout v0.1.0
pip install -e .

cd ..
wget https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM4859nnn/GSM4859876/suppl/GSM4859876_100pka_raw_feature_bc_matrix.h5
cellbender remove-background \
                 --epochs 150 \
                 --input GSM4859876_100pka_raw_feature_bc_matrix.h5 \
                 --output 100pka_cellbender.h5 \
                 --expected-cells 2100 \
                 --total-droplets-included 10000 \
                 --cuda

wget https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM4859nnn/GSM4859875/suppl/GSM4859875_05pka_raw_feature_bc_matrix.h5
cellbender remove-background \
                 --epochs 150 \
                 --input GSM4859875_05pka_raw_feature_bc_matrix.h5 \
                 --output 05pka_cellbender.h5 \
                 --expected-cells 1900 \
                 --total-droplets-included 10000 \
                 --cuda
```
