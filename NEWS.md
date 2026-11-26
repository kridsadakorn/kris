# KRIS 1.1.11

* `read.bed()` no longer returns wrong genotypes when the size of the BED file does not split evenly into its internal chunks (for example 10 individuals and 6 SNPs). It now also reads individual-major BED files and stops with a message if the BED file is too short (ISSUE-008).
* `cal.pc.linear(XXT = FALSE)` now returns the same eigenvalues and PCs as `XXT = TRUE`; before, it returned singular values instead of eigenvalues (ISSUE-008).
* `cal.pc.linear()` and `cal.pc.projection()` work with `no.pc = 1` (ISSUE-008).
* `xxt()` works on a matrix with one row (ISSUE-008).
* Fixed spelling, grammar and wrong facts in the documentation: the `simsnp` and `sample_labels` help pages now give the real size of the data (103 individuals, 2,000 SNPs), `read.bed()` names the FAM file correctly, `cal.pc.projection()` documents `$id`, and the `plot3views()` examples use as many colors as there are labels (ISSUE-007).
* Fixed typos in messages printed by `xxt()` and `write.bed()` (ISSUE-007).

# KRIS 1.1.10

* CITATION uses `bibentry()` and `c()` instead of the old-style `citEntry()` and `personList()` (ISSUE-001).
* Added `.Rbuildignore` so repository-only files are left out of the package (ISSUE-002).
* Updated GitHub Actions to r-lib/actions v2 (ISSUE-003).
* Fixed a broken PLINK link, dead README links and the DOI format (ISSUE-004).
* `xxt()` no longer fails when `window.size` equals the number of columns (ISSUE-005).
* `plot3views()` restores graphical parameters and accepts matrix input for `labels` (ISSUE-006).

# KRIS 1.1.9

* Updated README.Rmd after the pull request from fossabot; the pull request itself was not merged.

# KRIS 1.1.8

* Updated the manual.

# KRIS 1.1.7

* Added the DOI of the reference to the DESCRIPTION file.
* Updated Markdown files.
* GitLab only: built the documentation site with `pkgdown::build_site()`.
* GitLab only: added GitLab CI with `use_gitlab_ci()`.

# KRIS 1.1.6

* Updated the package for R 4.0.0.

# KRIS 1.1.5

* Added a function to calculate p-values between 2 groups with a regression model, with optional correction for covariates.
* Changed how `plot3views` is called and fixed bugs in it.

# KRIS 1.1.4

* Fixed some incorrect descriptions in the manual.

# KRIS 1.1.3

* Fixed mismatched colors and patterns in `plot3views`.

# KRIS 1.1.2

* Fixed bugs in `read.bed`.

# KRIS 1.1.1

* Updated and added examples in `plot3views`.

# KRIS 1.1.0

* Added the function ```cal.pc.projection``` 
* Moved ```replace.missing``` from principal.component.R to for_parallel.R
* Fixed bugs in test_plot3views.R

# KRIS 1.0.4
* Changed the URL in the CITATION file from GitLab to CRAN.
* Changed the email in the DESCRIPTION file from Gmail to the official email.
* Fixed a misspelling in the data manual.

# KRIS 1.0.3

* Updated URLs and added details of package in the DESCRIPTION file
* Updated the manual of all functions following CRAN's comments.

# KRIS 1.0.2

* Fixed the help file for ```read.bed```, ```write.bed```, ```fst.each.snp.hudson``` and ```fst.hudson```.

# KRIS 1.0.1

* Fixed the help file for ```cal.pc.linear``` and ```rubikclust```.

# KRIS 1.0.0

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
* ```write.bed``` Write a list of SNP objects to the binary PLINK format (BED, 
BIM, and FAM)
* ```xxt``` A function for calculating matrix multiplication between a matrix and 
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

Some functions were renamed from the initial version:

* The name of function ```cal.PC.linear``` was changed to ```cal.pc.linear```.
* The name of function ```plot.3views``` was changed to ```plot3views```.
* The name of function ```rubikClust``` was changed to ```rubikclust```.
* The name of function ```XXt``` was changed to ```xxt```.
