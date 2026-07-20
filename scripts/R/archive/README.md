# Archived R analyses

Files in this directory are retained for historical reference and are not part
of the numbered analysis pipeline.

- `04_baseline_characteristics.Rmd` is an obsolete exploratory notebook. It
  contains ineffective TableOne label arguments and a reference to an object
  removed earlier in the document, so it is not expected to knit.
- `05_QALYs.Rmd` is an obsolete exploratory notebook. It conditionally creates
  an object that is printed unconditionally and relies on chunk-local working
  directory changes, so it is not expected to knit.

The reusable Danish and UK EQ-5D-5L scoring functions formerly defined in
`05_QALYs.Rmd` now live in `scripts/R/functions/eq5d5l_utility.R` and have a
standalone regression test in `scripts/R/tests/test_eq5d5l_utility.R`.
