## Submission

Update of KRIS 1.1.6 (on CRAN) to 1.1.11. It fixes the NOTE about old-style citEntry() and
personList() in CITATION, bugs in read.bed(), cal.pc.linear() and xxt(), and errors in the documentation. See NEWS.md.

## Test environments

* local macOS (arm64), R 4.1.3: 0 errors, 0 warnings, 3 NOTEs (DOI check returned 403; unable to verify current time; no pandoc). Older than the current release, so not sufficient alone.
* GitHub Actions (1.1.10): macOS release, Windows release and devel, Ubuntu devel, release (4.6.2) and oldrel-1 (4.5.3): all passed. TODO: re-run on 1.1.11.
* TODO win-builder devel
* TODO rhub

## R CMD check results

TODO: fill in from the environments above.

## Reverse dependencies

TODO: FILEST (imports KRIS).
