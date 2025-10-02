# kDNA-Assembly

This repository contains a comprehensive pipeline for assembling the complete **kinetoplast DNA (kDNA)** of *Leishmania*, including both **maxicircle** and **minicircle** genomes. The pipeline supports both long-read only (Oxford Nanopore Technologies, ONT) and hybrid assemblies (ONT + Illumina), and includes downstream steps for gRNA annotation.

> 🚨 This is a complete pipeline that enables **direct minicircle assembly from long ONT reads**

---

## 📦 Requirements

- `Minimap2`, `BWA`, `SAMtools`, `BEDtools`
- `MetaFlye`, `Canu`, `Pilon`, `Unicycler`
- `Bandage`, `BLAST+`, `CD-HIT`, `MEME`
- `T-Aligner`, `Clustal Omega`
- Python 3, Biopython, NumPy, etc.

---

## 📘 Maxicircle Assembly Pipeline

### 🔹 Input:
- ONT (DNAseq) reads
- Illumina (DNAseq) reads

### 🔹 Steps:

1. **Map reads to nuclear assembly to extract kDNA-unmapped reads:**
```bash
minimap2 -ax map-ont nuclear.fasta ONT_reads.fastq > aln.sam
samtools view -b -f 4 aln.sam > unmapped.bam
bedtools bamtofastq -i unmapped.bam -fq unmapped_ONT.fastq
```
2. **Assemble unmapped ONT reads (MetaFlye):**
```bash
flye --nano-raw unmapped_ONT.fastq --out-dir flye_output --genome-size 30k
```
3. **Polish contig using Illumina reads (Pilon):**
```bash
bwa mem maxicircle.fasta illumina_R1.fastq illumina_R2.fastq > aln.sam
pilon --genome maxicircle.fasta --frags aln.bam --output maxicircle_polished
```
4.	**Circularity Check:**

	•	Load contigs in Bandage GUI to verify circularity.
5. **Gene Annotation**
```bash
blastn -query maxicircle.fasta -db nt -outfmt 6 -evalue 1e-5 -out annotations.tsv
```
## 📘 Minicircle Assembly Pipeline
This pipeline performs **de novo assembly and annotation of minicircle DNA** from *Leishmania* kinetoplast genomes. It supports both:

- **Long-read-only assembly** (Oxford Nanopore)
- **Hybrid assembly** (ONT + Illumina)
- **gRNA annotation** using RNA-seq and transcript alignment

## 📁 Folder Structure
minicircle_pipeline/
├── long_read_only/         # Long-read-only assembly workflow
├── hybrid_assembly/        # Hybrid ONT + Illumina workflow
├── scripts/                # Custom Python scripts (CSB1 rotation, clustering filters, etc.)
├── gRNA_annotation/        # gRNA annotation using T-Aligner
└── results/                # Output FASTA, clusters, trees, annotations


---

## 🔧 Prerequisites

Install the following tools (via conda or system packages):

- `minimap2`, `bwa`, `samtools`, `bedtools`, `seqkit`
- `CD-HIT`, `Canu`, `Unicycler`, `MEME Suite`
- `Clustal Omega`, `T-Aligner`, `Python 3`

---

## 🧬 Minicircle Assembly Pipeline

### 📌 Option 1: Long-Read Only Assembly

```bash
# Step 1: Filter unmapped ONT reads (500–1500 bp)
seqkit seq -m 500 -M 1500 unmapped_ONT.fastq > filtered_ONT.fastq

# Step 2: Cluster at 95% identity
cd-hit-est -i filtered_ONT.fastq -o clustered_reads.fasta -c 0.95 -n 10

# Step 3: Filter clusters with ≥10 reads (custom script)
python scripts/filter_clusters.py --input clustered_reads.fasta --clstr clustered_reads.clstr --min_reads 10 --output cluster_filtered/

# Step 4: Assemble each cluster individually
for cluster in cluster_filtered/*.fastq; do
    cluster_id=$(basename "$cluster" .fastq)
    canu -p $cluster_id -d canu_$cluster_id genomeSize=1k -nanopore-raw $cluster
done

# Step 5: Circularize assembled contigs (custom script)
for asm in canu_*/$cluster_id.contigs.fasta; do
    python scripts/simple_circularize.py --input $asm --output ${asm%.fasta}_circular.fasta
done

# Step 6: Final clustering of circularized minicircles
cat *_circular.fasta > all_circular_minicircles.fasta
cd-hit-est -i all_circular_minicircles.fasta -o minicircle_MSCs.fasta -c 0.95 -n 10

### 📌 Option 2: Hybrid Assembly (ONT + Illumina)
# Step 1: Extract unmapped reads from nuclear genome
# ONT reads
minimap2 -ax map-ont nuclear.fasta ONT.fastq | samtools view -b -f 4 - > unmapped_ONT.bam
bedtools bamtofastq -i unmapped_ONT.bam -fq unmapped_ONT.fastq

# Illumina reads
bwa mem nuclear.fasta R1.fastq R2.fastq | samtools view -b -f 4 - > unmapped_Illumina.bam
bedtools bamtofastq -i unmapped_Illumina.bam -fq unmapped_R1.fastq -fq2 unmapped_R2.fastq

# Step 2: Hybrid assembly using Unicycler
unicycler -1 unmapped_R1.fastq -2 unmapped_R2.fastq -l unmapped_ONT.fastq -o unicycler_output/

# Step 3: Identify CSB motifs (MEME)
meme unicycler_output/assembly.fasta -oc meme_out -mod zoops -nmotifs 3 -minw 15 -maxw 30 -dna

# Step 4: Rotate contigs to start at CSB1 (custom script)
python scripts/restart_with_csb1.py --input unicycler_output/assembly.fasta --csb1 ATGGT...

# Step 5: Cluster final minicircles
cd-hit-est -i restarted_minicircles.fasta -o minicircle_MSCs.fasta -c 0.95 -n 10

## 🧬 gRNA Annotation Pipeline
# Step 1: Transcript assembly from RNA-seq (T-Aligner)
taligner -i RNAseq.bam -o transcripts.fasta

# Step 2: Annotate gRNAs using minicircle sequences
taligner annotate -r transcripts.fasta -q minicircle_MSCs.fasta -o gRNA_annotations.tsv
