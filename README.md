# Translocator: local realignment and global remapping enabling accurate translocation detection using third generation sequencing data
[![License: AGPL v3](https://img.shields.io/badge/License-AGPL%20v3-blue.svg)](https://www.gnu.org/licenses/agpl-3.0)

Contact: Ye Wu	
Email: ywu@cs.hku.hk  

## Introduction
Translocation is an important class of structural variants known to be associated with cancer formation and treatment. The recent development in single-molecule sequencing technologies that produce long reads has promised an advance in detecting translocations accurately. However, existing tools struggled with the high base error-rate of the long reads. Figuring out the correct translocation breakpoints is especially challenging due to suboptimally aligned reads. To address the problem, we developed Translocator, a robust and accurate translocation detection method that implements an effective realignment algorithm to recover the correct alignments. For benchmarking, we analyzed using NA12878 long reads against a modified GRCh38 reference genome embedded with translocations at known locations. Our results show that Translocator significantly outperformed other state-of-the-art methods, including Sniffles and PBSV. On Oxford Nanopore data, the recall improved from 48.2% to 87.5% and the precision from 88.7% to 92.7%.


---

## Contents
- [Installation](#installation)
- [Quick Start](#quick-demo)
- [Usage](#usage)

---

## Installation

### Option 1. Build Translocator using cmake

```bash
git clone https://github.com/HKU-BAL/Translocator.git
cd Translocator/
mkdir build/
cd build/
cmake ..
make

cd ../bin/translocator*
./translocator
```

### Option 2. Bioconda (in progress)


## Quick Start
You need to have a sorted bam file preferably aligned using NGMLR or minimap2 and the reference file used for the alignment. Translocator's algorithm is not dependent on specific aligners, but we haven't tested out other aligners yet.

```bash
./translocator -m sorted.bam -a ref.fa -v output.vcf
```

---

## Usage

### For PacBio data
* Align the PacBio reads to a reference genome (using NGMLR)
```bash
ngmlr -t threads -r ref.fa -q PacBio.fq | samtools view -Sb mapped.bam
samtools sort -@ threads -O bam -o sorted.bam mapped.bam
```
* Call translocations and other SVs
```bash
translocator -m sorted.bam -a ref.fa -v out.vcf
```
### For Oxford Nanopore data
* Align the ONT reads to a reference genome (using minimap2)
```bash
minimap2 -t threads -ax map-ont ref.fa ont.fq.gz --MD | samtools view -Sb > mapped.bam
samtools sort -@ threads -O bam -o sorted.bam mapped.bam
```
* Call translocations and other SVs
```bash
translocator -m sorted.bam -a ref.fa -v out.vcf --global_remap
```

### Parameters

| ﻿Parameter           | Default | Description                                                              |
|---------------------|---------|--------------------------------------------------------------------------|
| -m/--mapped_reads   | NA      | Sorted .bam file either from NGMLR or minimap2                           |
| -a/--reference      | NA      | Reference used for realignment, should be consistent with the one used in the mapped file                                          |
| -v/--vcf            | NA      | Name of the vcf file to be reported                                      |
| -l/--min_length     | 100     | Minimum length of the SVs to be reported                                 |
| -g/--realign_length | 100     | Minimum length of the sequence to be realigned                           |
| --global_remap      | false   | if global remapping is enabled. Recommended to set as true for ONT data. |
| --realign_clipped   | true    | if clipped reads are realigned in local realignment.                     |
