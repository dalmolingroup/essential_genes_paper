### IDENTIFICATION OF LETHAL GENES IN DROSOPHILA

library(tidyverse)
library(biomaRt)

path <- "datasets/allele_phenotypic_data_fb_2019_01.tsv"
table_drosp <- read_delim(file = path, delim = "\t", col_names = T)

# import table with gene ids
table_genes <- read_tsv("scripts_geneplast/drosophila/fbal_to_fbgn_fb_2019_01.tsv")

# search for lethal genes
idx <- str_detect(table_drosp$phenotype, "lethal")
table_lethal_drosp1 <- table_drosp[idx, ]
idx <- str_detect(table_lethal_drosp1$phenotype, "partially")
table_lethal_drosp2 <- table_lethal_drosp1[!idx, ]
table_lethal_drosp2$allele_symbol <- str_replace(table_lethal_drosp2$allele_symbol, pattern = "\\[.*\\]", replacement = "")

# map the allele ids into gene ids used by FlyBase
table_lethal_mapped <- inner_join(table_lethal_drosp2, table_genes, by = c("allele_FBal" = "AlleleID"))

# get a list of all lethal genes
lethal_genes_drosp <- unique(table_lethal_mapped$GeneID)

# get a list of non-lethal genes
table_nonlethal_mapped <- anti_join(table_genes, table_lethal_mapped, by = "GeneID")
nonlethal_genes_drosp <- unique(table_nonlethal_mapped$GeneID)

# confirm if all the genes were mapped
n_distinct(nonlethal_genes_drosp) + n_distinct(lethal_genes_drosp) == n_distinct(table_genes$GeneID)

# map genes with biomaRt
ensembl <- useMart("ensembl")
ensembl <- useDataset("dmelanogaster_gene_ensembl", ensembl)
lethal_genes_map <- getBM(attributes = c("ensembl_peptide_id", "flybase_gene_id"),
                          filters = "flybase_gene_id",
                          values = lethal_genes_drosp, 
                          mart = ensembl)

# some ids have non-characters symbols and must be removed
nonlethal_genes_drosp <- nonlethal_genes_drosp[grep("^[A-Za-z0-9]+$", nonlethal_genes_drosp, perl = TRUE)]
non_lethal_map <- getBM(attributes = c("ensembl_peptide_id", "flybase_gene_id"),
                        filters = "flybase_gene_id",
                        values = nonlethal_genes_drosp,
                        mart = ensembl)

# map cogs
load("scripts_geneplast/geneplast_objects/geneplastData.RData")
cog_dros <- unique(cogdata[cogdata$ssp_id == 7227, c("protein_id", "cog_id")])
colnames(cog_dros)[1] <- "ensembl_peptide_id"

# save files
save(lethal_genes_map, non_lethal_map, cog_dros, 
     file = "scripts_geneplast/drosophila/tables_cogs_genes_drosp.RData")

# GENEPLAST ANALYSIS
# Prepare data for geneplast functions
cog_ids <- data.frame(cog_id = unique(cog_dros$cog_id), stringsAsFactors = F)
rownames(cog_ids) <- cog_ids$cog_id

#### GENEPLAST ANALYSIS ----
library(geneplast)

# create object ogp
ogp_drosp <- gplast.preprocess(cogdata = cogdata, 
                               sspids = ssp, 
                               cogids = cog_ids, verbose = FALSE)

# gplastTest
ogp_drosp <- gplast(ogp_drosp, verbose = FALSE)

# gplastRes
res_plast_drosp <- gplast.get(ogp_drosp, what = "results")

# create ogr object
ogr_drosp <- groot.preprocess(cogdata = cogdata, phyloTree = phyloTree,
                              spid = '7227', cogids = cog_ids, verbose=FALSE)

# grootTest
set.seed(1)
library(snow)
options(cluster=makeCluster(3, "SOCK"))
ogr_drosp <- groot(ogr_drosp, nPermutations = 10000, verbose = FALSE)
stopCluster(getOption("cluster"))

# grootResggsignif
res_root_drosp <- groot.get(ogr_drosp, what = "results")

# Phylotree
groot.plot(ogr_drosp, plot.lcas = T, fname = "scripts_geneplast/drosophila/tree_drosophila.pdf", height = 20, width = 10)

# MERGE TABLES OF ROOT INFERENCE AND GENES
res_root_drosp$cog_id <- rownames(res_root_drosp)
res_plast_drosp$cog_id <- rownames(res_plast_drosp)
table_cog_root <- inner_join(cog_dros, res_root_drosp, by = "cog_id")
table_cog_root <- inner_join(table_cog_root, res_plast_drosp,
                             by = "cog_id")

table_cog_root_lethal <- table_cog_root %>% 
  mutate(lethal_nonlethal = case_when(
    ensembl_peptide_id %in% lethal_genes_map$ensembl_peptide_id ~ "lethal"
  ), organism = "drosophila") %>% 
  filter(!is.na(lethal_nonlethal))

table_cog_root_nonlethal <- table_cog_root %>% 
  mutate(lethal_nonlethal = case_when(
    ensembl_peptide_id %in% non_lethal_map$ensembl_peptide_id ~ "nonlethal"
  ), organism = "drosophila") %>% 
  filter(!is.na(lethal_nonlethal))

table_cog_root <- rbind(table_cog_root_lethal, table_cog_root_nonlethal)

table_plot_drosophila <- table_cog_root
save(table_plot_celegans, file = "scripts_geneplast/drosophila/table_plot_drosophila.RData")













