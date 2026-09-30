# Build probe -> Entrez/symbol maps from Bioconductor annotation packages.
# Probes mapping to >1 Entrez gene are dropped. Output: data/annot/<platform>_map.tsv
suppressMessages({library(AnnotationDbi); library(jsonlite)})
out <- "data/annot"
mk <- function(pkg, tag) {
  suppressMessages(library(pkg, character.only = TRUE))
  db <- get(pkg)
  k <- keys(db, keytype = "PROBEID")
  m <- suppressMessages(AnnotationDbi::select(db, keys = k, columns = c("ENTREZID", "SYMBOL", "CHR"), keytype = "PROBEID"))
  m <- m[!is.na(m$ENTREZID), ]
  multi <- names(which(tapply(m$ENTREZID, m$PROBEID, function(x) length(unique(x))) > 1))
  m <- m[!(m$PROBEID %in% multi), ]
  m <- m[!duplicated(m$PROBEID), ]
  write.table(m, file.path(out, paste0(tag, "_map.tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
  cat(tag, pkg, as.character(packageVersion(pkg)), "probes:", nrow(m), "genes:", length(unique(m$ENTREZID)), "dropped multi:", length(multi), "\n")
}
mk("hgu133a.db", "GPL96")
mk("hgu133plus2.db", "GPL570")
mk("hugene10sttranscriptcluster.db", "GPL6244")
# Brainarray ENTREZG custom CDFs (GPL10526/GPL17027): IDs are "<entrez>_at"; symbols via org.Hs.eg.db
suppressMessages(library(org.Hs.eg.db))
cat("org.Hs.eg.db", as.character(packageVersion("org.Hs.eg.db")), "\n")
eg <- keys(org.Hs.eg.db, keytype = "ENTREZID")
s <- suppressMessages(AnnotationDbi::select(org.Hs.eg.db, keys = eg, columns = c("SYMBOL", "CHR"), keytype = "ENTREZID"))
s <- s[!duplicated(s$ENTREZID), ]
write.table(s, file.path(out, "entrez_symbol.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
for (p in c("frma","oligo","affy","limma","GEOquery","hgu133afrmavecs","hgu133plus2frmavecs","pd.hugene.1.0.st.v1","hgu133a.db","hgu133plus2.db","hugene10sttranscriptcluster.db","org.Hs.eg.db"))
  cat(p, tryCatch(as.character(packageVersion(p)), error=function(e) "NOT INSTALLED"), "\n")
