# About
# This script reads in raw data > preprocesses data > to form a data matrix
# Data matrix is used for further stats analysis

# ============== Clears Memory ======
# clears all objects includes hidden objects
rm(list = ls(all.names = TRUE)) 

# frees up memory and reports the memory usage.
gc() 

# ============== Loads Packages =======
library(readxl)
library(tidyverse)
library(data.table)
library(destiny)
library(tidyr)
library(stringr)
library(janitor)
library(IMIFA)

# ============== 1. Read & Selects from Excel File ======
chick_raw <- read_excel('Chick_Analysis_Dataset.xlsx', 
                        sheet = '20220708_105819_15Jun22_Chick_B') %>%
  select(-c(`PG.ProteinDescriptions`,
            `PG.ProteinNames`,
            `PG.MolecularFunction`))
  
# splits accession by ";" delimiter (ie "Q9JHU4-1; Q9JHU4" --> "Q9JHU4-1")
chick_raw$PG.ProteinAccessions <- sapply(strsplit(chick_raw$PG.ProteinAccessions,";"), `[`, 1)

# renames column to "Accession"
chick_export <- chick_raw %>%
  rename("Accession" = "PG.ProteinAccessions")

# exports accession numbers to upload to Uniprot
fwrite(data.frame(chick_export$`Accession`), "Chick_Accession.csv", sep = ",")

  
# ============== 2. Combines Uniprot Data To Combined Matrix =====
# reads in Gene Symbol table downloaded from Uniprot
gene_symbol_map <- fread('Chick_Uniprot.tsv')

colnames(gene_symbol_map) <- c("Accession", "Gene Symbol") 

# merges gene symbol column to main df
ratio_combined <- chick_export %>%
  
  # merges gene symbol column to main df
  left_join(gene_symbol_map,
            by = "Accession") %>%
  
  # relocates columns and removes NAs
  relocate(c(`Accession`,
             `Gene Symbol`)) %>%
  na.omit() %>%
  
  # adds number to the end of duplicate gene symbols (ie Sptbn1-2)
  group_by(`Gene Symbol`) %>%
  mutate(`GS_count` = 1:n()) %>%
  mutate(`Gene Symbol` = ifelse(`GS_count` == 1, 
                                `Gene Symbol`, 
                                paste0(`Gene Symbol`, "-", `GS_count`))) %>%
  
  # converts all NaN into NAs, removes NAs
  na_if("NaN") %>%
  na.omit() #3626 rows left

# exports data matrix for further analysis (using MetaboAnalyst)  
fwrite(data.frame(ratio_combined), "Output/Chick_Combined.csv", sep = ",")
