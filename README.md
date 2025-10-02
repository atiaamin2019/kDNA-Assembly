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

