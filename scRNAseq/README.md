# scRNAseq analysis

## Setup

```bash
mkdir -p fastq
while read srr; do
	ffq --ftp $srr | jq -r '.[] | .url' | head -n1 | xargs curl -O
	name=$(ffq --gcp $srr | jq -r '.[] | .filename' | head -n 1)
	name=$(basename $name .1)
	bam2fastq --nthreads=20 $name fastq/
done <SRR_Acc_List.txt
```
