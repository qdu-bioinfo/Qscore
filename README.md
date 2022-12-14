# Introduction
Qscore is a comprehensive method to evaluate the performance of amplicons by integrating the amplification rate, multi-tier taxonomic annotation, sequence type and length. It takes taxonomic annotation and its abundance at any taxonomic level as input, and predicts the optimal configuration of the ecological habitat through a large number of simulated amplification data that we analyze in advance.

# Software Requirements
The package depends C++ (>= 4.8.2) we recommend using the Homebrew package manager:
```
brew install gcc
```
# Online service
We recommend it directly at [Qscore Online](http://qscore.single-cell.cn/)

# Usage
### a. Calculate the optimal configuration of specific habitat
Calculate Qscore needs input a habitat with taxonomy annotation and abundance (required)，such as
```
#Taxonomy	Count
p__Dictyoglomi;	337
p__Elusimicrobia;	420
p__Euryarchaeota;	130
p__Fibrobacteres;	480
p__Firmicutes;	751
```
Then choose the reference database (R: NCBI RefSeq, G: Greengenes cut-off 0.97, C: Greengenes cut-off 0.99, S: Silva), and the weight of different evaluation indicators (Senstivity: Wprecision: Cost, default is 1:1:1) according to the user's preferences
```
Qscore -i taxonomy.abd -o optimal_configure.txt
```
### b. Simulate the second generation 16S rRNA sequencing data or shotgun sequencing
Here we provide the method of generating 16S rRNA sequencing data and shotgun sequencing data based on the full-length metagenome sequence
16S rRNA sequencing：
```
Generate_16S_rRNA -i genome.fa -o 16S_rRNA
```
shotgun sequencing:
```
Generate_WGS -i genome.fa -o WGS
```

