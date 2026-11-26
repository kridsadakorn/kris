# Resolved issues (newest first)

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
- Commit: the commit "ISSUE-001..006 (v1.1.10)"
- Verification: regression test passes on R 4.1.3; R >= 4.2 behaviour needs CI

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
- Commit: the commit "ISSUE-001..006 (v1.1.10)"
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
- Commit: the commit "ISSUE-001..006 (v1.1.10)"
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
- Commit: the commit "ISSUE-001..006 (v1.1.10)"
- Verification: NOT verified; needs a push and a CI run

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
- Commit: the commit "ISSUE-001..006 (v1.1.10)"
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
- Commit: the commit "ISSUE-001..006 (v1.1.10)"
- Verification: file rewritten; the CRAN NOTE can only be confirmed on R-devel
