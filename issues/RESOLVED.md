# Resolved issues (newest first)

### ISSUE-008: Wrong results in read.bed(), cal.pc.linear() and xxt()
- Reported: 2026-11-26
- Type: bug
- Planned bump: patch
- Status: resolved
- Description: Found by probing the code after a local test run.
  1. read.bed() split the raw bytes into 10 equal parts with round(); when the byte count did not divide evenly, parts overlapped or ran backwards, so genotypes came back wrong (write.bed/read.bed round trip failed for 10x6, 5x3, 7x1).
  2. cal.pc.linear(XXT = FALSE) returned singular values as eigenvalues, so evalue and PC scores differed from XXT = TRUE.
  3. no.pc = 1 failed with "non-conformable arguments": diag() of a single value builds an identity matrix.
  4. xxt() on a one-row matrix returned the wrong shape because column subsetting dropped dimensions.
- Analysis: see description.
- Plan: Decode BED data in whole SNP blocks (also handle individual-major files, check file length); square singular values; diag(x, nrow = length(x)); drop = FALSE.
- Resolved: 2026-11-26
- Version: 1.1.11
- Files: R/plinkformat.R, R/principal.component.R, R/matrix.multiply.R, tests/testthat/test_ISSUE-008.R, tests/testthat/test_calpcprojection.R (set.seed)
- Commit: the commit "... (v1.1.11)"
- Verification: test_ISSUE-008.R fails 9 times on the old code and passes on the new code (R 4.1.3). Example data reads back identical to simsnp.

### ISSUE-007: Spelling, grammar and factual errors in the documentation
- Reported: 2026-11-26
- Type: docs
- Planned bump: patch
- Status: resolved
- Description: User asked for a check of all documents. Found typos (patially, pricipal, charactors, multipication, devided, possitive, calcualte), grammar errors, and wrong facts.
- Analysis: Wrong facts: simsnp/sample_labels docs said 753 individuals and 3,000 SNPs (data has 103 and 2,000) and listed wrong column names and labels; read.bed() said "BAM file"; cal.pc.projection() documented `$individual_id` but returns `$id`; fst.each.snp.hudson() said it returns a matrix (it returns a vector); plot3views() examples gave 4 colors for 3 labels, so the custom colors were silently ignored.
- Plan: Fix roxygen, README.Rmd/README.md, DESCRIPTION, NEWS.md; fix message typos in xxt() and write.bed(). Follow-up: rewrite all help pages and the README for readability; add the missing @return to plot3views().
- Resolved: 2026-11-26
- Version: 1.1.11
- Files: DESCRIPTION, NEWS.md, README.Rmd, README.md, R/*.R, man/*.Rd
- Commit: the commit "ISSUE-007: ... (v1.1.11)"
- Verification: local R 4.1.3 `R CMD check --as-cran`: 0 errors, 0 warnings, same 3 environment NOTEs as 1.1.10; CI not yet run on 1.1.11

### ISSUE-006: plot3views() leaves par() changed and breaks on matrix input
- Reported: 2026-11-26
- Type: cran
- Planned bump: patch
- Status: resolved
- Description: CRAN requires par() to be restored; class(x) == "..." fails for matrices on R >= 4.2.
- Analysis: No on.exit(); class comparisons.
- Plan: Restore par with on.exit(); use is.*().
- Resolved: 2026-11-26
- Version: 1.1.10
- Files: R/plot.3views.R, tests/testthat/test_ISSUE-006.R
- Commit: 06b2c8b
- Verification: regression test passes on R 4.1.3 and in CI on R 4.5.3, 4.6.2 and devel

### ISSUE-005: xxt() fails when window.size equals ncol(X)
- Reported: 2026-11-26
- Type: bug
- Planned bump: patch
- Status: resolved
- Description: Error "subscript out of bounds".
- Analysis: The sum loop ran over 2:length() even with a single block.
- Plan: Sum blocks with Reduce().
- Resolved: 2026-11-26
- Version: 1.1.10
- Files: R/matrix.multiply.R, tests/testthat/test_ISSUE-005.R
- Commit: 06b2c8b
- Verification: regression test passes on R 4.1.3

### ISSUE-004: Broken links and DOI format
- Reported: 2026-11-26
- Type: docs
- Planned bump: patch
- Status: resolved
- Description: Local check flags the PLINK URL; README has dead links; DOI has a space.
- Analysis: Old URLs.
- Plan: Update links, write <doi:...>.
- Resolved: 2026-11-26
- Version: 1.1.10
- Files: DESCRIPTION, R/plinkformat.R, man/read.bed.Rd, man/write.bed.Rd, README.Rmd, README.md
- Commit: 06b2c8b
- Verification: links checked with curl; urlchecker not run (no pandoc)

### ISSUE-003: GitHub Actions workflow cannot run
- Reported: 2026-11-26
- Type: ci
- Planned bump: patch
- Status: resolved
- Description: Workflow uses ubuntu-16.04 and r-lib/actions@v1.
- Analysis: Those runners no longer exist.
- Plan: Replace with r-lib/actions v2 matrix.
- Resolved: 2026-11-26
- Version: 1.1.10
- Files: .github/workflows/R-CMD-check.yaml
- Commit: 06b2c8b
- Verification: GitHub Actions run 45777507299 passed all 6 jobs (R 4.6.2 release, 4.5.3 oldrel-1, devel)

### ISSUE-002: No .Rbuildignore in the repository
- Reported: 2026-11-26
- Type: cran
- Planned bump: patch
- Status: resolved
- Description: Local check: NOTEs for hidden files and non-standard top-level files.
- Analysis: .gitignore listed .Rbuildignore, so it was never committed.
- Plan: Track .Rbuildignore; remove it from .gitignore.
- Resolved: 2026-11-26
- Version: 1.1.10
- Files: .Rbuildignore, .gitignore
- Commit: 06b2c8b
- Verification: local R 4.1.3 check; hidden-files and non-standard-files NOTEs are gone

### ISSUE-001: Old-style CITATION
- Reported: 2026-11-26
- Type: cran
- Planned bump: patch
- Status: resolved
- Description: CRAN r-devel Debian NOTE on 1.1.6: personList() and citEntry() are old-style.
- Analysis: inst/CITATION still used both on GitHub.
- Plan: Rewrite with bibentry() and c(person()).
- Resolved: 2026-11-26
- Version: 1.1.10
- Files: inst/CITATION
- Commit: 06b2c8b
- Verification: file rewritten; the CRAN NOTE can only be confirmed on R-devel
