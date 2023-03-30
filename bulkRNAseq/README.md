# bulkRNAseq analysis

## Preparations

It is possible to run these preparation steps in parallel to gain time.

### Download fastq files

```bash
while read srr; do
	ffq --ftp $srr | jq -r '.[] | .url' | xargs curl -O
	oldname=$(ffq --ftp $srr | jq -r '.[] | .filename' | head -n 1)
	name=$(ffq --gcp $srr | jq -r '.[] | .filename' | head -n 1)
	name=$(basename $name .1)
	mv $oldname fastq/$name
done <SRR_Acc_List.txt
```

### Make STAR index

```bash
wget -O - https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M27/gencode.vM27.annotation.gtf.gz | gunzip -c > gencode.vM27.annotation.gtf
wget -O - https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M32/GRCm39.primary_assembly.genome.fa.gz | gunzip -c > GRCm39.primary_assembly.genome.fa

STAR --runMode genomeGenerate --runThreadN 40 --genomeDir GRCm39_vM27_Muscle --genomeFastaFiles GRCm39.primary_assembly.genome.fa --sjdbGTFfile gencode.vM27.annotation.gtf --sjdbOverhang 75
```

## Generate transcriptome counts

This needs to be run on a powerful system:

```bash
bt -f fastq -s GRCm39_vM27_Muscle -g gencode.vM27.annotation.gtf
```