# scRNAseq analysis

## Downloading BAM files and converting back to fastq

```bash
wget -O bam2fastq https://github.com/10XGenomics/bamtofastq/releases/download/v1.4.1/bamtofastq_linux
chmod +x bam2fastq
while read srr; do
        ffq --ftp $srr | jq -r '.[] | .url' | head -n1 | xargs curl -O
	name=$(basename $name .1)
	./bam2fastq --nthreads=20 $name $srr
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
