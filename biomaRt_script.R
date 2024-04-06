#SUBSTITUIR O CHUNK DOS LOGCOUNTS POR ESTe

## distributions - log transform
logcounts = cpm(dgeObj,log=TRUE)

# Set up the connection to Ensembl
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
ensembl_dataset <- useDataset('hsapiens_gene_ensembl', mart = ensembl)

# Remove version numbers from Ensembl gene IDs in logcounts
ensembl_ids <- gsub("\\..*", "", rownames(logcounts))
rownames(logcounts) <- ensembl_ids

# Retrieve gene symbols for the Ensembl IDs
gene_symbols <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"),
                      filters = "ensembl_gene_id",
                      values = ensembl_ids,
                      mart = ensembl)

# Loop through each Ensembl gene ID in seqdata
for (i in seq_along(ensembl_ids)) {
  # Find the index of the Ensembl gene ID in gene_symbols
  idx <- match(ensembl_ids[i], gene_symbols$ensembl_gene_id)
  # If a corresponding gene symbol is found and it's not an empty string, replace the row name
  if (!is.na(idx) && gene_symbols$external_gene_name[idx] != "") {
    rownames(logcounts)[i] <- gene_symbols$external_gene_name[idx]
  }
}