# scRNAseq analysis

## Downloading BAM files and converting back to fastq

```bash

while read srr; do
        mkdir -p $srr
        ffq --ftp $srr | jq -r '.[] | .url' | head -n1 | xargs curl -O
	name=$(basename $name .1)
	bam2fastq --nthreads=20 $name fastq/
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
        $pathtocellranger3.1/cellranger-cs/3.1.0/bin/count --id=$srr --fastqs=$srr --transcriptome=refdata-cellranger-mm10-1.2.0 --sample=$srr
done <SRR_Acc_List.txt
```
