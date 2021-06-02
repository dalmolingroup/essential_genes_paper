# Caenorhabditis elegans --------------------------------------------------

library(data.table)
library(tidyverse)
library(biomaRt)
library(geneplast)

# Read table containing only lethal genes
table_lethal <- as.data.frame(fread("datasets/c_elegans_filtered.csv"))

# Get all genes
ensembl <- useMart("ensembl")
ensembl <- useDataset("celegans_gene_ensembl", ensembl)
all_genes <- getBM(attributes = c("ensembl_peptide_id", "wormbase_gene"),
                   mart = ensembl)

# Get lethal genes 
lethal_genes <- semi_join(all_genes, table_lethal, by = c("wormbase_gene" = "wb_gene_id"))

# Get non-lethal genes
nonlethal_genes <- anti_join(all_genes, table_lethal, by = c("wormbase_gene" = "wb_gene_id"))

# Retrieve COG ids
load("scripts_geneplast/geneplast_objects/geneplast_objects/geneplastData.RData")
cog_celegans <- unique(cogdata[cogdata$ssp_id == 6239, c("protein_id", "cog_id")])
colnames(cog_celegans)[1] <- "ensembl_peptide_id"

# Save tables
save(lethal_genes, nonlethal_genes, cog_celegans, 
     file = "scripts_geneplast/c_elegans/tables_cogs_genes_celegans.RData")

######################

# GENEPLAST ANALYSIS ----

# Prepare data for geneplast functions
cog_ids <- data.frame(cog_id = unique(cog_celegans$cog_id), stringsAsFactors = F)
rownames(cog_ids) <- cog_ids$cog_id

# Create object ogp
ogp_celegans <- gplast.preprocess(cogdata = cogdata, 
                                  sspids = ssp, 
                                  cogids = cog_ids, verbose = FALSE)

# GplastTest
ogp_celegans <- gplast(ogp_celegans, verbose = FALSE)

# GplastRes
res_plast_celegans <- gplast.get(ogp_celegans, what = "results")

# Greate ogr object
ogr_celegans <- groot.preprocess(cogdata = cogdata, phyloTree = phyloTree,
                                 spid = '6239', cogids = cog_ids, verbose=FALSE)

# GrootTest
set.seed(1)
library(snow)
options(cluster = makeCluster(3, "SOCK"))
ogr_celegans <- groot(ogr_celegans, nPermutations = 10000, verbose = FALSE)
stopCluster(getOption("cluster"))

# grootRes
res_root_celegans <- groot.get(ogr_celegans, what = "results")

# Phylotree
groot.plot(ogr_celegans, plot.lcas = T, fname = "tree_celegans.pdf")

# MERGE TABLES OF ROOT INFERENCE AND GENES
res_root_celegans$cog_id <- rownames(res_root_celegans)
res_plast_celegans$cog_id <- rownames(res_plast_celegans)
table_cog_root <- inner_join(cog_celegans, res_root_celegans, by = "cog_id")
table_cog_root <- inner_join(table_cog_root, res_plast_celegans,
                             by = "cog_id")

table_cog_root_lethal <- table_cog_root %>% 
  mutate(lethal_nonlethal = case_when(
    ensembl_peptide_id %in% lethal_genes$ensembl_peptide_id ~ "lethal"
  ), organism = "celegans") %>% 
  filter(!is.na(lethal_nonlethal))

table_cog_root_nonlethal <- table_cog_root %>% 
  mutate(lethal_nonlethal = case_when(
    ensembl_peptide_id %in% nonlethal_genes$ensembl_peptide_id ~ "nonlethal"
  ), organism = "celegans") %>% 
  filter(!is.na(lethal_nonlethal))

table_cog_root <- rbind(table_cog_root_lethal, table_cog_root_nonlethal)
table_plot_celegans <- table_cog_root
save(table_plot_celegans, file = "scripts_geneplast/c_elegans/table_plot_celegans.RData")















