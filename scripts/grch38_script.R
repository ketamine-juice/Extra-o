#Variável

my_genes = resultados_mod$row_names

# Filter 'grch38' based on 'row_names' from 'resultados_mod'
resultados_id <- grch38 %>%
  filter(ensgene %in% resultados_mod$row_names) %>%
  # Select only the 'symbol' column
  select(symbol)

resultados_ensembl_id <- grch38 %>%
  filter(ensgene %in% resultados_mod$row_names) %>%
  # Select both 'ensgene' and 'symbol' columns
  select(ensgene, symbol)

# Print the filtered table
print(resultados_id)

resultados_id = as.data.frame(resultados_id)

resultados_ensembl_id

# Loop through each gene in my_genes
for (i in seq_along(my_genes)) {
  # Check if the gene is not present in resultados_ensembl_id
  if (!my_genes[i] %in% resultados_ensembl_id$ensgene) {
    # Create a new row in resultados_ensembl_id with NA for ensgene and the gene name
    resultados_ensembl_id <- rbind(resultados_ensembl_id, c(NA, my_genes[i]))
  }
}
