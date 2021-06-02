### MUS MUSCULUS ROOT INFERENCE FOR LETHAL AND NON LETHAL GENES

# read table 
library(readr)
ontology <- read_delim("scripts_geneplast/mouse/mammalian_ontology.txt", delim = "\t", 
                       col_names = c("mp_id", "description", "details"),
                       col_types = "ccc")
mouse_geno <- read_delim("datasets/MGI_PhenoGenoMP.rpt",
                         delim = "\t", col_names = c("allelic_composition", "allele_symbol",
                                                     "genetic_backgrond", "mp_id", "pubmed_id",
                                                     "mgi_id"))

library(dplyr)
table_merged <- inner_join(ontology, mouse_geno, by = "mp_id")

# get one gene per row
library(tidyr)
table_merged <- mutate(table_merged, 
                       mgi_id = strsplit(as.character(mgi_id), ","))
table_merged <- unnest(table_merged, mgi_id )

# filter by lethality
idx <- grep('lethal', table_merged$description)
table_merged_lethal <- table_merged[idx, ]

# read table of defined categories 
# in this list, I excluded the ontology MP:0010831 because it was very unnespecific.
ontologies_lethal_categories <- read_delim(
  file = "scripts_geneplast/mouse/ontologies_lethal_categories.txt", delim = "\t") 

# I decided to eliminate the ontologies related with 'incomplete penetrance' and 
# the ones not well defined
idx <- grep('incomplete penetrance', ontologies_lethal_categories$description)
ontologies_lethal_categories <- ontologies_lethal_categories[-idx, ]
idx <- grep(',', ontologies_lethal_categories$categories)
ontologies_lethal_categories <- ontologies_lethal_categories[-idx, ]
write.csv(ontologies_lethal_categories, 
          file = "scripts_geneplast/mouse/ontologies_lethal_categories_final.csv", 
          quote = F, row.names = F)

# merge tables of lethal genes and ontologies
lethal_genes_table <- inner_join(table_merged_lethal, ontologies_lethal_categories,
                                 by = c("mp_id" = "mp.id"))
# final list of lethal genes
lethal_genes <- unique(lethal_genes_table$mgi_id)

# table of non-lethal genes
nonlethal_genes_table <- anti_join(table_merged, lethal_genes_table, by = "mgi_id")

# final list of non-lethal genes
nonlethal_genes <- unique(nonlethal_genes_table$mgi_id)

# get ensembl peptide ids for the lethal and non-lethal genes
library(biomaRt)
ensembl <- useMart('ensembl')
ensembl <- useDataset('mmusculus_gene_ensembl', ensembl)
lethal_genes_biomart <- getBM(attributes = c('mgi_symbol', 'mgi_id', 
                                             'mgi_description', 'ensembl_peptide_id'), 
                              filters = 'mgi_id', 
                              mart = ensembl, 
                              values = lethal_genes)
lethal_genes_biomart[lethal_genes_biomart$ensembl_peptide_id == "",] <- NA
lethal_genes_biomart <- na.omit(lethal_genes_biomart)

nonlethal_genes_biomart <- getBM(attributes = c('mgi_symbol', 'mgi_id', 
                                                'mgi_description', 'ensembl_peptide_id'), 
                                 mart = ensembl)
nonlethal_genes_biomart <- anti_join(nonlethal_genes_biomart, lethal_genes_biomart, by = "mgi_id")
nonlethal_genes_biomart[nonlethal_genes_biomart$ensembl_peptide_id == "",] <- NA
nonlethal_genes_biomart <- na.omit(nonlethal_genes_biomart)

all_genes <- getBM(attributes = c('mgi_symbol', 'mgi_id', 
                                  'mgi_description', 'ensembl_peptide_id'), mart = ensembl)
none_genes <- anti_join(all_genes, table_merged, by = "mgi_id")

# there are about 2500 mgi ids that do not match with all mgi ids retrieved from biomaRt. 
# is biomaRt correct?
sum(!unique(table_merged$mgi_id) %in% unique(all_genes$mgi_id))

# MAP COGS 
load("scripts_geneplast/geneplast_objects/geneplastData.RData")
mouse_cogs <- unique(cogdata[cogdata$ssp_id == 10090, c("protein_id", "cog_id")])
colnames(mouse_cogs)[1] <- "ensembl_peptide_id"

# save tables 
save(lethal_genes_biomart, nonlethal_genes_biomart, mouse_cogs,
     file = "scripts_geneplast/mouse/table_cogs_mouse.RData")

# prepare data for geneplast functions
cog_ids <- data.frame(cog_id = unique(mouse_cogs$cog_id), stringsAsFactors = F)
rownames(cog_ids) <- cog_ids$cog_id

#### GENEPLAST ANALYSIS FOR ROOT INFERING OF LETHAL GENES
library(geneplast)

# create object ogp
ogp_mouse <- gplast.preprocess(cogdata = cogdata, 
                               sspids = ssp, 
                               cogids = mouse_cogs, verbose = FALSE)

# gplastTest
ogp_mouse <- gplast(ogp_mouse, verbose = FALSE)

# gplastRes
res_plast_mouse <- gplast.get(ogp_mouse, what = "results")

# create ogr object
ogr_mouse <- groot.preprocess(cogdata = cogdata, phyloTree = phyloTree,
                              spid = '10090', cogids = mouse_cogs, verbose = FALSE)

# grootTest
set.seed(1)
library(snow)
options(cluster = makeCluster(4, "SOCK"))
ogr_mouse <- groot(ogr_mouse, nPermutations = 10000, verbose = FALSE)
stopCluster(getOption("cluster"))

# grootRes
res_root_mouse <- groot.get(ogr_mouse, what = "results")

# Phylotree
groot.plot(ogr_mouse, plot.lcas = T, fname = "scripts_geneplast/mouse/tree_mouse.pdf")

# MERGE TABLES OF ROOT INFERENCE AND GENES
res_root_mouse$cog_id <- rownames(res_root_mouse)
res_plast_mouse$cog_id <- rownames(res_plast_mouse)
table_cog_root <- inner_join(mouse_cogs, res_root_mouse, by = "cog_id")
table_cog_root <- inner_join(table_cog_root, res_plast_mouse,
                             by = "cog_id")

table_cog_root_lethal <- table_cog_root %>% 
  mutate(lethal_nonlethal = case_when(
    ensembl_peptide_id %in% lethal_genes_biomart$ensembl_peptide_id ~ "lethal"
  ), organism = "mouse") %>% 
  filter(!is.na(lethal_nonlethal))

table_cog_root_nonlethal <- table_cog_root %>% 
  mutate(lethal_nonlethal = case_when(
    ensembl_peptide_id %in% nonlethal_genes_biomart$ensembl_peptide_id ~ "nonlethal"
  ), organism = "mouse") %>% 
  filter(!is.na(lethal_nonlethal))

table_cog_root <- rbind(table_cog_root_lethal, table_cog_root_nonlethal)
table_plot_mouse <- table_cog_root

save(table_plot_mouse, file = "scripts_geneplast/mouse/table_plot_mouse.RData")
