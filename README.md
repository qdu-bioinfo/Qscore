# Introduction
Qscore is a comprehensive method to evaluate the performance of amplicons by integrating the amplification rate, multi-tier taxonomic annotation, sequence type and length. It takes taxonomic annotation and its abundance at any taxonomic level as input, and predicts the optimal configuration of the ecological habitat through a large number of simulated amplification data that we analyze in advance.

# System Requirement and dependency
## Hardware Requirements
Qscore only requires a standard computer with sufficient RAM to support the operations defined by a user. For typical users, this would be a computer with about 2 GB of RAM. For optimal performance, we recommend a computer with the following specs:

  RAM: 8+ GB  
  CPU: 4+ cores

## Software Requirements
The package depends C++ (>= 4.8.2) we recommend using the Homebrew package manager:
```
brew install gcc
```
# Online service
We recommend it directly at [Qscore Online](http://qscore.single-cell.cn/)

# Installation guide
## Automatic Installation (recommended)
#### **a. Download the package**
```
git clone https://github.com/qdu-bioinfo/qscore.git	
```
#### **b. Install by installer**
```
cd Qscore
source install.sh
```

# Usage
### a. Calculate the optimal configuration of specific habitat
You can run this following command to get detailed parameters.
```
Qscore -h
```
Calculate Qscore needs input a habitat with taxonomy annotation and abundance (required)，such as
```
#Taxonomy	Count
p__Dictyoglomi;	337
p__Elusimicrobia;	420
p__Euryarchaeota;	130
p__Fibrobacteres;	480
p__Firmicutes;	751
```

```
Qscore -D R -i taxonomy.abd -W 1 1 1 -o optimal_configure.txt
```
Here taxonomy.abd is a taxonomic annotation file at any level with abundance, and optimal_configure.txt is the output file. Then choose the reference database (R: NCBI RefSeq, G: GreenGenes-13-8, C: GreenGenes-13-8-99, S: Silva), and the weight of different evaluation indicators (Senstivity: Wprecision: Cost, default is 1:1:1) according to the user's preferences.

### b. Simulate the next-generation 16S rRNA sequencing data or shotgun sequencing
Here we provide the method of generating 16S rRNA sequencing data and shotgun sequencing data based on the full-length metagenome sequence
### 16S rRNA sequencing：
You can run this following command to get detailed parameters.
```
Extract_16S_rRNA -h
```
And this following command provide a demo to get all 16S rRNA (default 300bp) of primer extract from a full-length metagenome sequence file. 
```
Extract_16S_rRNA -i genome.fa -o 16S_rRNA
```
### Shotgun sequencing:
You can run this following command to get detailed parameters.
```
Extract_WGS -h
```
And this following command provide a demo to get a simulated genome sequencing specimen (default 300bp) of use silider window extract from a full-length metagenome sequence file. 
```
Extract_WGS -i genome.fa -o WGS
```
Here genome.fa is a full-length metagenome sequence. These simulated sequences can also be used to test the accuracy of Qscore prediction through [NCBI-opal](https://pubmed.ncbi.nlm.nih.gov/34837585/)
