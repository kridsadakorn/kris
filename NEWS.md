---
title: "NEWS"
author: "Kridsadakorn Chaichoompu"
date: "27/04/2018"
output:
  html_document:
    keep_md: yes
  pdf_document: default
---



# kris 1.0.0.9000

---

# kris 1.0.0

## Initial functions

* ```cal.pc.linear``` A function for linear principal component analysis (PCA)
* ```fst.each.snp.hudson``` A function for fixation index (Fst) calculation for 
all SNPs between two groups.
* ```fst.hudson``` A function for average fixation index (Fst) calculation 
between two groups.
* ```plot3views``` A function to create scatter plots in three views.
* ```read.bed``` Read the binary PLINK format (BED, BIM, and FAM)
* ```rubikclust``` A function for unsupervised clustering to detect rough 
structures and outliers.
* ```write.bed``` Write an list of SNP object to the binary PLINK format (BED, 
BIM, and FAM)
* ```xxt``` A function for calculating matrix multipication between a matrix and 
its transpose for large data.

## Initial R data 

* ```simsnp``` Synthetic dataset containing single nucleotide polymorphisms 
(SNP)
* ```sample_labels``` Synthetic dataset containing population labels for the 
dataset simsnp.

## Initial example files

* ```example_SNP.bed``` Synthetic dataset containing single nucleotide polymorphisms 
(SNP) in binary format
* ```example_SNP.bim``` Simulated SNP information
* ```example_SNP.fam``` Simulated sample information

## Updates

From the initial idea, some functions were changed their names:

* The name of function ```cal.PC.linear``` was changed to ```cal.pc.linear```.
* The name of function ```plot.3views``` was changed to ```plot3views```.
* The name of function ```rubikClust``` was changed to ```rubikclust```.
* The name of function ```XXt``` was changed to ```xxt```.