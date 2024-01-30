### Assign Taxonomy
Now we assign taxonomy to our sequences using the 16S database. You can download the fasta and txt file [here](https://ucedna.com/reference-databases-for-metabarcoding).
We're gonna use local BLASTN. 

```sh

```

For inspecting the classification, removing sequence rownames for display only.
```sh
taxa_print <- taxa  
rownames(taxa_print) <- NULL
head(taxa_print)
```

Formatting data and read tables, transpose reads table
```
seq_table_nochim_transpose <- t(seq_table_nochim) 
```

Combine taxa and read tables
```
combined_table <- cbind(taxa, seq_table_nochim_transpose) 
```

Write into csv file
```
write.csv(taxa, file="taxa.csv") 
write.csv(seq_table_nochim, file="reads.csv") 
write.csv(seq_table_nochim_transpose, file="reads_t.csv") 
write.csv(combined_table, file="Combined_raw_ASVs_table.csv") 
```
