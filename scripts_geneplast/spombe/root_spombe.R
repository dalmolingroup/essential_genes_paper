# S. POMBE ANALYSIS 

library(tidyverse)
library(biomaRt)

spombe_table <- read.delim("datasets/FYPOviability.tsv", 
                           col.names = c("gene", "phenotype"))

# find inviables genes
idx <- grep("inviable", spombe_table$phenotype)
spombe_inviables <- spombe_table[idx,]

# find nonlethal genes
spombe_nonlethal <- anti_join(spombe_table, spombe_inviables)

# MAP COGS
load("scripts_geneplast/geneplast_objects/geneplastData.RData")
spombe_cogs <- unique(cogdata[cogdata$ssp_id == 4896, c("protein_id", "cog_id")])
colnames(spombe_cogs)[1] <- "ensembl_peptide_id"
#spombe_cogs$ensembl_peptide_id <- gsub("\\.1$", "", spombe_cogs$ensembl_peptide_id)

# save cog tables 
save(spombe_inviables, spombe_nonlethal, spombe_cogs, 
     file = "scripts_geneplast/spombe/cog_spombe_tables.RData")

# Prepare data for geneplast functions
cog_ids <- data.frame(cog_id = unique(spombe_cogs$cog_id), stringsAsFactors = F)
rownames(cog_ids) <- cog_ids$cog_id

#### GENEPLAST ANALYSIS ----
library(geneplast)

# create object ogp
ogp_spombe <- gplast.preprocess(cogdata = cogdata, sspids = ssp, cogids = cog_ids)

# gplastTest
ogp_spombe <- gplast(ogp_spombe)

# gplastRes
res_plast_spombe <- gplast.get(ogp_spombe, what = "results")

# create ogr object
ogr_spombe <- groot.preprocess(cogdata = cogdata, phyloTree = phyloTree,
                               spid = '4896', cogids = cog_ids)

# grootTest
set.seed(1)
library(snow)
options(cluster=makeCluster(3, "SOCK"))
ogr_spombe <- groot(ogr_spombe, nPermutations = 10000, verbose = FALSE)
stopCluster(getOption("cluster"))

# grootResggsignif
res_root_spombe <- groot.get(ogr_spombe, what = "results")

# Phylotree
groot.plot(ogr_spombe, plot.lcas = T, fname = "scripts_geneplast/spombe/tree_spombe.pdf", width = 10, height = 20)


# MERGE TABLES OF ROOT INFERENCE AND GENES
res_root_spombe$cog_id <- rownames(res_root_spombe)
res_plast_spombe$cog_id <- rownames(res_plast_spombe)
table_cog_root <- inner_join(spombe_cogs, res_root_spombe, by = "cog_id")
# table_cog_root <- inner_join(table_cog_root, res_plast_spombe,
#                              by = "cog_id")

spombe_inviables$gene <- paste0(spombe_inviables$gene, ".1")
spombe_nonlethal$gene <- paste0(spombe_nonlethal$gene, ".1")

table_cog_root_lethal <- table_cog_root %>% 
  mutate(lethal_nonlethal = case_when(
    ensembl_peptide_id %in% spombe_inviables$gene ~ "lethal"
  ), organism = "spombe") %>% 
  filter(!is.na(lethal_nonlethal))

table_cog_root_nonlethal <- table_cog_root %>% 
  mutate(lethal_nonlethal = case_when(
    ensembl_peptide_id %in% spombe_nonlethal$gene ~ "nonlethal"
  ), organism = "spombe") %>% 
  filter(!is.na(lethal_nonlethal))

table_cog_root <- rbind(table_cog_root_lethal, table_cog_root_nonlethal)

table_plot_spombe <- table_cog_root
save(table_plot_spombe, file = "scripts_geneplast/spombe/table_plot_spombe.RData")
