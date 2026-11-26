# KRIS — project rules for Claude

R package **KRIS** (Keen and Reliable Interface Subroutines for Bioinformatic Analysis).
Goal: keep it passing `R CMD check --as-cran` on every platform and publish updates to CRAN.

## Project facts (checked 2026-11-26)

- GitHub (source of truth): https://github.com/kridsadakorn/kris, branch `master`
- CRAN: https://cran.r-project.org/package=KRIS. **On CRAN** (not archived), version 1.1.6, published 2021-01-21.
- CRAN has seen 1.0.3, 1.1.1, 1.1.6. GitHub reached 1.1.9 without a CRAN release.
- Reverse dependency: **FILEST** imports KRIS. Check it before changing any export.
- Maintainer email: `kridsadakorn@biostatgen.org`. CRAN sends every notice there, so it must be an inbox the maintainer reads.
- Local R is 4.1.3 (old). Current R release on 2026-11-26: 4.6.2. No `gh` CLI, no pandoc, no `spelling` package.
- CRAN submission winter break: 2026-12-22 to 2027-01-07.

## Layout

| Path | What it is | Who edits it |
|---|---|---|
| `R/` | package source (roxygen2 comments); no `src/` | r-coder |
| `man/`, `NAMESPACE` | **generated** by `roxygen2::roxygenise()`, never hand-edit (`man/figures/` holds README images) | regenerate only |
| `tests/testthat/` | unit tests | tester |
| `vignettes/` | none at present | docs-qc |
| `inst/CITATION` | citation, `bibentry()` + `c(person())` only | docs-qc |
| `DESCRIPTION`, `NEWS.md` | metadata, version, changelog | main session |
| `README.Rmd` → `README.md` | edit the .Rmd, knit to .md | docs-qc |
| `docs/` | pkgdown site, excluded from build | docs-qc |
| `cran-comments.md` | notes for the CRAN submission form | release-architect |
| `.github/workflows/` | CI | ci-setup |
| `issues/` | `ISSUES.md` (open), `RESOLVED.md` (closed, newest first) | main session |

Claude-only files (`CLAUDE.md`, `.claude/`, `issues/`, `cran-comments.md`) must stay in `.Rbuildignore`. `.Rbuildignore` is tracked in git; do not add it to `.gitignore`.

## Agents (`.claude/agents/`)

| Agent | Model | Use for |
|---|---|---|
| `r-coder` | sonnet | writing/fixing R code and roxygen |
| `tester` | sonnet | tests, local `R CMD check` |
| `ci-fixer` | sonnet | failing GitHub Actions run |
| `docs-qc` | haiku | all documentation QC |
| `ci-setup` | haiku | workflow YAML |
| `release-architect` | opus | go/no-go, `cran-comments.md`; use sparingly |

Rules: use the cheapest agent that can do the job. Use opus only for release decisions or a problem sonnet could not solve. The main session does small tasks itself.

## Ask the user first

Deleting files (move to Trash, never `rm`), pushing to GitHub, opening or sending email, submitting to CRAN.

## Issue workflow

On "report issue ..." or "มี issue ...":

1. Record it in `issues/ISSUES.md` as `ISSUE-NNN`.
2. Analyse the cause, choose agents, fix.
3. Test.
4. Bump the version and write `NEWS.md`, citing the issue ID.
5. Move the entry to the top of `issues/RESOLVED.md` with Resolved / Version / Files / Commit / Verification.
6. Commit as `ISSUE-NNN: <title> (vX.Y.Z)`. Push only when the user says so.

## Automatic version bump (SemVer)

Every resolved issue bumps `Version:` in `DESCRIPTION` without asking.

- PATCH: bug fix, CRAN NOTE/WARNING, docs, CI, internal refactor
- MINOR: backward-compatible function, argument or feature
- MAJOR: removed/renamed export, changed default or output format, raised minimum R
- Several issues at once: one bump at the highest level.
- Never reuse a version CRAN has seen; a new version must be higher than every submitted one.
- `NEWS.md` heading `# KRIS X.Y.Z` must match DESCRIPTION. Update the `inst/CITATION` year if needed.
- Tag `vX.Y.Z` only after CRAN accepts.

## Testing against the latest R

- Never submit unless the package passes on the current R release **and** R-devel. CRAN checks with R-devel first, then release/oldrel on Linux, macOS and Windows.
- A local R older than the current release is not evidence. Before each submission look up the current release at https://cran.r-project.org/ and confirm CI ran on it and on devel.
- Tools: GitHub Actions (all OS), `devtools::check_win_devel()`, rhub, `R CMD check --as-cran`.

## CRAN release checklist

1. `roxygen2::roxygenise()` leaves no diff in `man/` or `NAMESPACE`.
2. `R CMD build` + `R CMD check --as-cran` outside the repo: 0 errors, 0 warnings, only acceptable NOTEs.
3. CI green on all OS including R-devel; win-builder devel and rhub pass.
4. `urlchecker::url_check()` and spell check pass.
5. Reverse dependency check of FILEST.
6. release-architect writes the go/no-go list and updates `cran-comments.md`.
7. Commit everything; ask before pushing.
8. The user submits at https://cran.r-project.org/submit.html and confirms the email. Claude never does.

## Known history

- Never archived as of 2026-11-26. Open CRAN finding on 1.1.6: old-style `personList()`/`citEntry()` in CITATION (r-devel Debian NOTE).
- GitHub CI used `ubuntu-16.04` and `r-lib/actions@v1` until 1.1.10.
- This checkout (`kris_github/kris`) is the only working copy. Old copies, do not use: `../../kris`, `../../kris_gitlab/kris`, `../../kris_gitlab/kris_doc` (old website repo), `../../kris_rel_backup/kris`.
- Local branch `backup/gitlab-6b92a87` keeps an unpushed GitLab commit (LICENSE wording only). Do not push it.
