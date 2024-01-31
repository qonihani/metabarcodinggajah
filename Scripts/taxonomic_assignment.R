#BLAST on Windows di command prompt/powershell
#BLASTN on Windows

#--BLASTN--

makeblastdb -in CO1_Blast.fasta -dbtype nucl -out CO1_Blast.fasta -title "COI_Database"

blastdbcmd -db CO1_Blast.fasta -info

blastn -query ASV_Sequence.fasta -db CO1_Blast.fasta -out taxonomic_assignment.txt -max_target_seqs 1 -perc_identity 97 -outfmt "6 qseqid sseqid"

# Merging Taxonomic result with dada2 result in Rstudio
taxonomic_assignment <- read.table("taxonomic_assignment.txt")


Taxonomy = right_join(seq_table_nochim_trans, taxonomic_assignment, by = c("ASVNumber" = "V1"))

head(Taxonomy)

X16S_taxo <- read_delim("16S_taxo.txt", delim = ";", 
                        +     escape_double = FALSE, trim_ws = TRUE)
taxo_16S <- X16S_taxo

Taxonomy = merge(Taxonomy, taxo_16S, by.x = "V2", by.y = "X1", all.x=TRUE)


colnames(Taxonomy) <- c("Accession","ASVNumber","Sequence","Anak-F","Anak-M","Induk-F","Induk-M","Kingdom","Phylum","Class","Order","Family","Genus","Species")
head(Taxonomy)

write.csv(Taxonomy, file="Results.csv")
write_xlsx(Taxonomy,"Results.xlsx")
