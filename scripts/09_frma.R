# Single-array frozen RMA (McCall et al. 2010) for every CEL file of the retained series.
# fRMA uses only frozen, externally estimated vectors => no information shared across samples,
# so it can be run once for all arrays without leaking between folds/test set.
# Output: data/frma/<GSE>__<chiptype>.tsv.gz  (probeset x GSM, log2)
suppressMessages({library(affy); library(frma); library(affyio)})
args <- commandArgs(TRUE)
raw <- "data/raw"; out <- "data/frma"; dir.create(out, showWarnings = FALSE)
chunk <- 40
for (gse in args) {
  d <- file.path(raw, gse)
  if (!dir.exists(d)) { dir.create(d); untar(file.path(raw, paste0(gse, "_RAW.tar")), exdir = d) }
  cels <- list.files(d, pattern = "\\.cel(\\.gz)?$", ignore.case = TRUE, full.names = TRUE)
  ct <- vapply(cels, function(f) read.celfile.header(f)$cdfName, "")
  cat(gse, "CELs:", length(cels), "chip types:", paste(names(table(ct)), table(ct), collapse = "; "), "\n")
  for (chip in unique(ct)) {
    of <- file.path(out, paste0(gse, "__", chip, ".tsv.gz"))
    if (file.exists(of)) { cat("  exists", of, "\n"); next }
    fs <- cels[ct == chip]; res <- list()
    for (i in seq(1, length(fs), by = chunk)) {
      sub <- fs[i:min(i + chunk - 1, length(fs))]
      if (grepl("HuGene", chip, ignore.case = TRUE)) {
        suppressMessages(library(oligo))
        raw_fs <- oligo::read.celfiles(sub)
        e <- exprs(frma(raw_fs, target = "core"))
      } else {
        ab <- ReadAffy(filenames = sub)
        e <- exprs(frma(ab))
      }
      res[[length(res) + 1]] <- e
      cat("  ", chip, i, "-", i + length(sub) - 1, "\n")
    }
    m <- do.call(cbind, res)
    colnames(m) <- sub("^(GSM\\d+).*$", "\\1", basename(colnames(m)))
    con <- gzfile(of, "w"); write.table(round(m, 5), con, sep = "\t", quote = FALSE, col.names = NA); close(con)
    cat("  wrote", of, dim(m), "\n")
  }
}
