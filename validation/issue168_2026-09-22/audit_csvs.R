# Compare every tracked CSV with the freshly rendered object-verified artifact.
# The writer checks exact CSV serialization of each generating object; this
# audit additionally records differences from the pre-fix Git tree.
root <- here::here()
manifest <- read.csv(file.path(root, "outputs/manifest.csv"), stringsAsFactors = FALSE)
files <- system2("git", c("ls-files", "outputs/tables/*.csv"), stdout = TRUE)
rows <- lapply(files, function(file) {
  entry <- manifest[manifest$file == file, , drop = FALSE]
  stopifnot(nrow(entry) == 1L, entry$status == "csv_verified",
    entry$artifact_md5 == unname(tools::md5sum(file.path(root, file))))
  original <- system2("git", c("show", paste0("b8e1b02:", file)), stdout = TRUE)
  stopifnot(is.null(attr(original, "status")))
  before <- read.csv(text = paste(original, collapse = "\n"), check.names = FALSE,
    colClasses = "character", na.strings = "NA")
  after <- read.csv(file.path(root, file), check.names = FALSE,
    colClasses = "character", na.strings = "NA")
  same_schema <- identical(names(before), names(after)) && identical(dim(before), dim(after))
  changed <- if (same_schema) {
    a <- as.matrix(before); b <- as.matrix(after)
    sum(xor(is.na(a), is.na(b)) | (!is.na(a) & !is.na(b) & a != b))
  } else NA_integer_
  data.frame(file = file, generating_vignette = entry$generating_vignette,
    verified_against_generating_object = TRUE, before_rows = nrow(before),
    after_rows = nrow(after), schema_unchanged = same_schema,
    changed_cells = changed, values_unchanged = identical(before, after),
    artifact_md5 = entry$artifact_md5)
})
audit <- do.call(rbind, rows)
write.csv(audit, file.path(root, "validation/issue168_2026-09-22/csv_audit.csv"), row.names = FALSE)
print(audit[, c("file", "values_unchanged", "changed_cells")], row.names = FALSE)
cat("PASS:", nrow(audit), "tracked publication CSVs verified against generating objects.\n")
