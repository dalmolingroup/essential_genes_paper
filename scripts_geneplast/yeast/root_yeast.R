# YEAST ANALYSIS (S. scerevisae)

library(tidyverse)
library(biomaRt)

cols <- c('feature_name','feature_type','gene_name','sgd_id','reference',
          'experiment_type','mutant_type','allele','strain_bg','phenotype',
          'chemical','conditions','details','reporter')
scerevisae_table <- read.delim(file = "datasets/phenotype_data.tab", 
                               col.names = cols, stringsAsFactors = F)

# find lethal genes
idx <- scerevisae_table$mutant_type == 'null' & scerevisae_table$phenotype == 'inviable'
scerevisae_inviables <- scerevisae_table[idx, ]

# find nonlethal genes 
scerevisae_nonlethal <- anti_join(scerevisae_table, scerevisae_inviables, 
                                  by = "sgd_id")

ensembl <- useMart('ensembl')
ensembl <- useDataset('scerevisiae_gene_ensembl', ensembl)
all_genes <- getBM(attributes = c('ensembl_peptide_id', "external_gene_name"), 
                   mart = ensembl)

all_genes$ensembl_peptide_id <- gsub("_mRNA", "", all_genes$ensembl_peptide_id)
all_genes[all_genes$ensembl_peptide_id == "",] <- NA
all_genes <- na.omit(all_genes)

scerevisae_nonlethal <- anti_join(all_genes, scerevisae_inviables, 
                                  by = c("ensembl_peptide_id" = "feature_name"))

# MAP COGS
load("scripts_geneplast/geneplast_objects/geneplastData.RData")
yeast_cogs <- unique(cogdata[cogdata$ssp_id == 4932, c("protein_id", "cog_id")])
colnames(yeast_cogs)[1] <- "ensembl_peptide_id"

# save cog tables 
save(scerevisae_inviables, scerevisae_nonlethal, yeast_cogs, 
     file = "scripts_geneplast/yeast/cog_yeast_tables.RData")

# Prepare data for geneplast functions
cog_ids <- data.frame(cog_id = unique(yeast_cogs$cog_id), stringsAsFactors = F)
rownames(cog_ids) <- cog_ids$cog_id

#### GENEPLAST ANALYSIS ----
library(geneplast)

# create object ogp
ogp_yeast <- gplast.preprocess(cogdata = cogdata, 
                               sspids = ssp, 
                               cogids = cog_ids, verbose = FALSE)

# gplastTest
ogp_yeast <- gplast(ogp_yeast, verbose = FALSE)

# gplastRes
res_plast_yeast <- gplast.get(ogp_yeast, what = "results")

# create ogr object
ogr_yeast <- groot.preprocess(cogdata = cogdata, phyloTree = phyloTree,
                              spid = '4932', cogids = cog_ids, verbose=FALSE)

# grootTest
set.seed(1)
library(snow)
options(cluster=makeCluster(3, "SOCK"))
ogr_yeast <- groot(ogr_yeast, nPermutations = 10000, verbose = FALSE)
stopCluster(getOption("cluster"))

# groot
res_root_yeast <- groot.get(ogr_yeast, what = "results")

# Phylotree
groot.plot(ogr_yeast, plot.lcas = T, fname = "scripts_geneplast/yeast/tree_yeast.pdf", height = 20, width = 10)

# MERGE TABLES OF ROOT INFERENCE AND GENES
res_root_yeast$cog_id <- rownames(res_root_yeast)
res_plast_yeast$cog_id <- rownames(res_plast_yeast)
table_cog_root <- inner_join(yeast_cogs, res_root_yeast, by = "cog_id")
table_cog_root <- inner_join(table_cog_root, res_plast_yeast,
                             by = "cog_id")

table_cog_root_lethal <- table_cog_root %>% 
  mutate(lethal_nonlethal = case_when(
    ensembl_peptide_id %in% scerevisae_inviables$feature_name ~ "lethal"
  ), organism = "yeast") %>% 
  filter(!is.na(lethal_nonlethal))

table_cog_root_nonlethal <- table_cog_root %>% 
  mutate(lethal_nonlethal = case_when(
    ensembl_peptide_id %in% scerevisae_nonlethal$ensembl_peptide_id ~ "nonlethal"
  ), organism = "yeast") %>% 
  filter(!is.na(lethal_nonlethal))

table_cog_root <- rbind(table_cog_root_lethal, table_cog_root_nonlethal)

table_plot_yeast <- table_cog_root
save(table_plot_yeast, file = "scripts_geneplast/yeast/table_plot_yeast.RData")




















