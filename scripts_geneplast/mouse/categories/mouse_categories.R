# MOUSE LETHAL GENES

# Read tables
library(readr)
colnames_ont <- c("mp_id", "description","details")
ontology <-
  read_delim(
    file = "scripts_geneplast/mouse/mammalian_ontology.txt",
    delim = "\t",
    col_names = colnames_ont)

colnames_geno <- c("allelic_composition", "allele_symbol", "genetic_background", 
                   "mp_id", "pubmed_id", "mgi_id")
mouse_geno <-
  read_delim(file = "datasets/MGI_PhenoGenoMP.rpt",
             delim = "\t",
             col_names = colnames_geno)

# merge tables 
library(dplyr)
table_merged <- inner_join(ontology, mouse_geno, by = "mp_id")

# filter by lethality
idx <- grep("lethal", table_merged$description)
table_merged_filtered <- table_merged[idx, ]

# save
save(table_merged, table_merged_filtered, file = 'scripts_geneplast/mouse/mouse_table_ontology.RData')

# generate rows with single MGI ids
library(tidyr)
table_merged_filtered <- mutate(table_merged_filtered,
                                mgi_id = strsplit(as.character(mgi_id), ","))
table_merged_filtered <- unnest(table_merged_filtered, mgi_id)

# use the MGI id to retrieve ensembl peptide ids of lethal genes
library(biomaRt)
ensembl <- useMart("ensembl")
ensembl <- useDataset("mmusculus_gene_ensembl", ensembl)

lethal_genes <- table_merged_filtered$mgi_id
lethal_genes_table <- getBM(
  attributes = c('mgi_symbol', 'mgi_id', 'mgi_description', 'ensembl_peptide_id'),
  filters = "mgi_id",
  values = lethal_genes,
  mart = ensembl)
lethal_genes_table$ensembl_peptide_id[lethal_genes_table$ensembl_peptide_id == ""] <- NA
idx <- complete.cases(lethal_genes_table)
lethal_genes_table <- lethal_genes_table[idx,]

### DEFINE CATEGORIES
# Create a table with the ontology descriptions of lethal genes
ontologies <- data.frame(mp_id = unique(table_merged_filtered$mp_id),
                         stringsAsFactors = F)
ontologies_lethal <- inner_join(ontology, ontologies, by = "mp_id")

# After that, categories were defined manualy, i.e. we defined the categories and
# created a table relating them with MPID. The table with categories is named
# 'ontologies_lethal_categories.txt' saved 

# read table of defined categories
ontologies_lethal_categories <- read_tsv(
  file = "scripts_geneplast/mouse/ontologies_lethal_categories.txt",
  col_names = c("mp_id", "categories", "description", "details"),
  skip = 1)

# I decided to eliminate the ontologies related with 'incomplete penetrance'
idx <- grep("incomplete penetrance", ontologies_lethal_categories$description)
ontologies_lethal_categories <- ontologies_lethal_categories[-idx, ]

ontologies_lethal_categories <- inner_join(
  ontologies_lethal_categories, table_merged_filtered,
  by = "mp_id"
)

ontologies_lethal_categories <- inner_join(
  ontologies_lethal_categories, lethal_genes_table,
  by = "mgi_id"
)

# organize table 
working_lethal_mouse <- ontologies_lethal_categories[, c(11,12,14,1,2,3)]
colnames(working_lethal_mouse) <- 
  c('mgi_id', 'mgi_symbol', 'ensembl_peptide_id','mp_id','categories','description')

# eliminate redundancies of categories
idx <- grep(',', working_lethal_mouse$categories)
working_lethal_mouse <- working_lethal_mouse[-idx,]
working_lethal_mouse <- unique(working_lethal_mouse)

# MAP COGS
load("scripts_geneplast/geneplast_objects/geneplastData.RData")
mouse_cogs <- unique(cogdata[cogdata$ssp_id == 10090, c("protein_id", "cog_id")])
colnames(mouse_cogs)[1] <- "ensembl_peptide_id"

lethal_cogs <- inner_join(mouse_cogs, lethal_genes_table, by = "ensembl_peptide_id")
lethal_cogs <- inner_join(lethal_cogs, working_lethal_mouse, 
                          by = c("ensembl_peptide_id", "mgi_id", "mgi_symbol"))
lethal_cogs <- unique(lethal_cogs)

lethal_cog_ids_mouse <- data.frame(cog_id = unique(lethal_cogs$cog_id), stringsAsFactors = F)

save(lethal_cogs, working_lethal_mouse,
     file = "scripts_geneplast/mouse/categories/working_lethal_mouse.RData")

# Save lethal cogs by category
lethal_cog_ids_mouse <- data.frame(cog_id = unique(lethal_cogs$cog_id), stringsAsFactors = F)
temp <- split(lethal_cogs, f = lethal_cogs$categories)
early_mouse_cogs <- temp[[1]]
mild_mouse_cogs <- temp[[2]]
late_mouse_cogs <- temp[[3]]

save(early_mouse_cogs, mild_mouse_cogs, late_mouse_cogs, file = "scripts_geneplast/mouse/categories/lethal_cogs_categories_mouse.RData")

# GENEPLAST ANALYSIS
library(geneplast)

# create object ogp
ogp.mouse <- gplast.preprocess(cogdata = cogdata, sspids = ssp, cogids = lethal_cog_ids_mouse, 
                               verbose=FALSE)

# gplastTest
ogp.mouse <- gplast(ogp.mouse, verbose=FALSE)

# gplastRes
res.plast.mouse <- gplast.get(ogp.mouse, what = "results")

# create ogr object
ogr <- groot.preprocess(cogdata=cogdata, phyloTree=phyloTree, spid = '10090', 
                        cogids=lethal_cog_ids_mouse, verbose=FALSE)

# grootTest
set.seed(1)
ogr <- groot(ogr, nPermutations = 10000, verbose=FALSE)

# grootRes
res.root.mouse <- groot.get(ogr, what = "results")

# merge tables
res.root.mouse$orthologous_group <- rownames(res.root.mouse)
root_genes_mouse <- inner_join(res.root.mouse, lethal_cogs, by = c("orthologous_group" = "cog_id"))

save(res.plast.mouse, res.root.mouse,ogp.mouse, ogr, root_genes_mouse,
     file = "scripts_geneplast/mouse/categories/geneplast_results.RData")

groot.plot(ogr,plot.lcas = TRUE, fname = "scripts_geneplast/mouse/categories/tree_mouse_categories.pdf")


# Compare proportion of essential genes in the three categories -----------

# Load mouse data 
load("mouse/geneplast_results.RData")
load("table_plot_mouse.RData")

# Get roots for each category
root <- unique(table_plot_mouse[, c(2,3)])

# Get 90th root percentile 
old_root <- quantile(root$Root, probs = 0.9)

# Select the roots greater or equal than the 90th percentile
res <- root_genes_mouse %>% 
  group_by(categories) %>% 
  mutate(prop = ifelse(Root >= old_root, 1, 0))

# Count the number of older COGs and the total number of COGs by category
case <- tapply(res$prop[res$prop == 1], res$categories[res$prop == 1], sum)
total <- tapply(res$prop, res$categories, length)

# Perform proportion test
test <- prop.test(case, total)
prop <- as.vector(test$estimate)
names(prop) <- c("early", "late", "mild")
prop.trend.test(case, total)

dt <- data.frame(categories = c("Early", "Late", "Mild"), 
                 prop = prop, stringsAsFactors = F)
dt$categories <- factor(dt$categories, levels = c("Early", "Mild", "Late"))

# Plot
library(ggplot2)
library(ggpubr)
plot1 <- ggplot(dt, aes(x = categories, y = prop)) +
  geom_bar(stat = "identity") +
  theme_classic() +
  labs(title = "", x = "", y = "Proportion") + 
  theme(axis.text = element_text(face = "plain", size = 12, vjust = 0.5),
        axis.ticks.x = element_blank(),
        axis.line.x.top = element_blank(),
        panel.grid = element_blank(),
        axis.line.y.right = element_blank(),
        
        plot.title = element_text(face = "plain", hjust = 0.5)) 

plot1













