
# FUNCTIONAL ENRICHMENT ANALYSIS ------------------------------------------

library(tidyverse)
library(biomaRt)
library(clusterProfiler)
library(org.Ce.eg.db)
library(org.Dm.eg.db)
library(org.Mm.eg.db)
library(org.Sc.sgd.db)

# Read table
load("scripts_geneplast/table_all_org.RData")

# Select only the lethal genes
all_common_cogs <- purrr::reduce(split(table$cog_id[table$lethal_nonlethal == "lethal"], 
                                       table$organism[table$lethal_nonlethal == "lethal"]), dplyr::intersect)
table_cdm <- table[table$lethal_nonlethal == "lethal" & 
                 table$organism %in% c("celegans", "drosophila", "mouse"), ]

table_y <- table[table$lethal_nonlethal == "lethal" &
                   table$organism == "yeast",]


cogs_cdm <- setdiff(purrr::reduce(split(table_cdm$cog_id, table_cdm$organism),
                                  dplyr::intersect), unique(table_y$cog_id))

cogs_y <- setdiff(table_y$cog_id, table_cdm$cog_id)


# Table containing only the 196 lethal COGS in mouse, drosophila, and celegans
table_cdm <- table_cdm[table_cdm$cog_id %in% cogs_cdm,]

# Table containing only the 50 lethal COGS exclusive for yeast
table_y <- table_y[table_y$cog_id %in% cogs_y, ]


# MULTICELLULAR
# Enrichement analysis for common COGs for M. musculus, C. elegans and D. melanogaster --------

lethal_genes <- lapply(split(table_cdm$ensembl_peptide_id, table_cdm$organism), unique)

org_db <- c("org.Ce.eg.db", "org.Dm.eg.db", "org.Mm.eg.db")

lethal_genes_entrez <- lapply(seq_along(lethal_genes), function(i) {
  bitr(lethal_genes[[i]], fromType = "ENSEMBLPROT", toType = "ENTREZID", org_db[i])[[2]]
})

enrich <- lapply(seq_along(lethal_genes_entrez), function(i) {
  temp <- enrichGO(lethal_genes_entrez[[i]], OrgDb = org_db[i], 
           keyType = "ENTREZID", ont = "BP", pvalueCutoff = 0.001, minGSSize = 20)@result
  temp <- temp[temp$p.adjust < 0.05, ]
  return(temp)
})

go <- lapply(enrich, `[[`, 1)
go_cdm <- purrr::reduce(go, dplyr::intersect)

terms_cdm <- go2term(goid = go_cdm)

# UNICELLULAR
# Enrichment analysis for S. cerevisiae -------------------------------------------

lethal_genes_yeast <- paste0(table_y$ensembl_peptide_id, "_mRNA")
yeast_entrez <- bitr(lethal_genes_yeast, fromType = "ENSEMBLPROT", toType = "ENTREZID",
                     OrgDb = "org.Sc.sgd.db")[[2]]

go_yeast <- enrichGO(yeast_entrez, OrgDb = "org.Sc.sgd.db", 
                     keyType = "ENTREZID", ont = "BP", pvalueCutoff = 0.001)@result
go_yeast <- go_yeast[go_yeast$p.adjust < 0.05,]
go_yeast <- go_yeast[, 1:2]


# Enrichment analysis for all common COGs ---------------------------------

lethal_genes_all <- table[table$cog_id %in% all_common_cogs,]
lethal_genes_all <- split(lethal_genes_all$ensembl_peptide_id, lethal_genes_all$organism)
lethal_genes_all$yeast <- paste0(lethal_genes_all$yeast, "_mRNA")

org_db <- c("org.Ce.eg.db", "org.Dm.eg.db", "org.Mm.eg.db", "org.Sc.sgd.db")

lethal_genes_all_entrez <- lapply(seq_along(lethal_genes_all), function(i) {
  bitr(lethal_genes_all[[i]], fromType = "ENSEMBLPROT", toType = "ENTREZID", org_db[i])[[2]]
})

enrich_all <- lapply(seq_along(lethal_genes_entrez), function(i) {
  temp <- enrichGO(lethal_genes_all_entrez[[i]], OrgDb = org_db[i], 
                   keyType = "ENTREZID", ont = "BP", pvalueCutoff = 0.001)@result
  temp <- temp[temp$p.adjust < 0.05, ]
  return(temp)
})

go <- lapply(enrich_all, `[[`, 1)
go_all <- purrr::reduce(go, dplyr::intersect)

terms_all <- go2term(goid = go_all)

# ENRICHMENT ANALYSIS OF MOUSE CATEGORIES ---------------------------------

# Load all genes
load("scripts_geneplast/table_all_org.RData")

# Select only the yeast cogs
yeast_lethal_cogs <- unique(table$cog_id[table$organism == "yeast" &
                                           table$lethal_nonlethal == "lethal"])

# Load mouse categories
load("scripts_geneplast/mouse/categories/lethal_cogs_categories_mouse.RData")

# Intersect early mouse and yeast cogs
intersect(early_mouse_cogs$cog_id, yeast_lethal_cogs)

# Enrichment analysis of mouse categories
mouse_categories <- do.call(rbind, list(early_mouse_cogs, late_mouse_cogs, mild_mouse_cogs))
lethal_genes_mouse <- split(mouse_categories$ensembl_peptide_id, mouse_categories$categories)

lethal_genes_entrez <- lapply(seq_along(lethal_genes_mouse), function(i) {
  bitr(lethal_genes_mouse[[i]], fromType = "ENSEMBLPROT", toType = "ENTREZID", "org.Mm.eg.db")[[2]]
})

enrich <- lapply(seq_along(lethal_genes_entrez), function(i) {
  temp <- enrichGO(lethal_genes_entrez[[i]], OrgDb = "org.Mm.eg.db", 
                   keyType = "ENTREZID", ont = "BP", pvalueCutoff = 0.001, minGSSize = 20)@result
  temp <- temp[temp$p.adjust < 0.05, ]
  return(temp)
})

enrich_mild <- enrich[[2]]
enrich_late <- enrich[[3]]

save(enrich_late, enrich_mild, file = "GO_enrichment/enrich_mouse_categories.RData")

# Early genes did not enrich, because of the number of genes (entrez genes = 11)


# Enrichment analysis of percentiles ---------------------------------------

species <- c("celegans", "drosophila", "mouse", "yeast")

# Get the upper 30 pencentile
percent_30 <- sapply(split(table$Root, table$organism), function (i) {
  quantile(i, probs = 0.3)
})
names(percent_30) <- species

percent_70 <- sapply(split(table$Root, table$organism), function (i) {
  quantile(i, probs = 0.7)
})
names(percent_70) <- species

table$percent <- ""
lapply(species, function (i) {
  table$percent[table$organism == i] <<- ifelse(table$Root[table$organism == i] <= percent_30[i], "percent_30", 
                                          ifelse(table$Root[table$organism == i] >= percent_70[i], "percent_70", NA))
})
table_percent <- na.omit(table)

lethal_genes_30 <- split(table_percent$ensembl_peptide_id[table_percent$percent == "percent_30" & table_percent$lethal_nonlethal == "lethal"], 
                         table_percent$organism[table_percent$percent == "percent_30" & table_percent$lethal_nonlethal == "lethal"])


lethal_genes_70 <- split(table_percent$ensembl_peptide_id[table_percent$percent == "percent_70" & table_percent$lethal_nonlethal == "lethal"], 
                         table_percent$organism[table_percent$percent == "percent_70" & table_percent$lethal_nonlethal == "lethal"])


org_db <- c("org.Ce.eg.db", "org.Dm.eg.db", "org.Mm.eg.db", "org.Sc.sgd.db")

lethal_genes_entrez_30 <- lapply(seq_along(lethal_genes_30), function(i) {
  bitr(lethal_genes_30[[i]], fromType = "ENSEMBLPROT", toType = "ENTREZID", org_db[i])[[2]]
})

lethal_genes_entrez_70 <- lapply(seq_along(lethal_genes_70), function(i) {
  bitr(lethal_genes_70[[i]], fromType = "ENSEMBLPROT", toType = "ENTREZID", org_db[i])[[2]]
})

enrich_30 <- lapply(seq_along(lethal_genes_entrez_30), function(i) {
  temp <- enrichGO(lethal_genes_entrez_30[[i]], OrgDb = org_db[i], 
                   keyType = "ENTREZID", ont = "BP", pvalueCutoff = 0.001, minGSSize = 20)@result
  temp <- temp[temp$p.adjust < 0.05, ]
  temp <- temp[order(temp$p.adjust),]
  return(temp)
})

enrich_70 <- lapply(seq_along(lethal_genes_entrez_70), function(i) {
  temp <- enrichGO(lethal_genes_entrez_70[[i]], OrgDb = org_db[i], 
                   keyType = "ENTREZID", ont = "BP", pvalueCutoff = 0.001, minGSSize = 20)@result
  temp <- temp[temp$p.adjust < 0.05, ]
  temp <- temp[order(temp$p.adjust),]
  return(temp)
})

lapply(seq_along(enrich_30),  function (i) {
  write.table(enrich_30[[i]][, c("ID", "Description")], file = paste0("GO_enrichment/percentiles/percentiles_", species[i], "_enrich_30.csv"),
              sep = ";", quote = F, row.names = F)
})

lapply(seq_along(enrich_70),  function (i) {
  write.table(enrich_70[[i]][, c("ID", "Description")], file = paste0("GO_enrichment/percentiles/percentiles_", species[i], "_enrich_70.csv"),
              sep = ";", quote = F, row.names = F)
})


# Enrichment analysis for the exclusive genes only ------------------------

# Map the exclusive genes belonging to COGs exclusive to each species
species <- c("celegans", "drosophila", "mouse", "yeast")

filter_cogs <- function(table, lethality) {
  table <- table[table$lethal_nonlethal == lethality,]
  cogs <- lapply(split(table$cog_id, table$organism), unique)
  cogs
}
cogs_list <- filter_cogs(table, "lethal")
names(cogs_list) <- species
essential_exclusive <- lapply(species, function (i) {
  unique(setdiff(cogs_list[[i]], unique(unlist(unname(cogs_list[ !grepl(i, names(cogs_list)) ] )))))
})
names(essential_exclusive) <- species
table$exclusive <- NA

table <- table[table$lethal_nonlethal == "lethal",]
combined_cogs <- reduce(cogs_list, intersect)

invisible(lapply(species, function (i) {
  table$exclusive[table$organism == i] <<- ifelse(table$cog_id[table$organism == i & table$lethal_nonlethal == "lethal"] %in% essential_exclusive[[i]], "E", NA )
}))
table$exclusive[table$cog_id %in% combined_cogs] <- "C"
table <- na.omit(table)

n_distinct(table$cog_id[table$exclusive == "C"]) == length(combined_cogs)

# Enrichment analysis for exclusive essential genes
table_exclusive <- table[table$exclusive == "E",]
lethal_genes <- lapply(split(table_exclusive$ensembl_peptide_id, table_exclusive$organism), unique)
org_db <- c("org.Ce.eg.db", "org.Dm.eg.db", "org.Mm.eg.db", "org.Sc.sgd.db")
lethal_genes$celegans <- gsub("\\.(\\d|\\w)*", "", lethal_genes$celegans)
lethal_genes_entrez <- lapply(seq_along(lethal_genes), function(i) {
  bitr(lethal_genes[[i]], fromType = "ENSEMBLPROT", toType = "ENTREZID", org_db[i])[[2]]
})

enrich <- lapply(seq_along(lethal_genes_entrez), function(i) {
  temp <- enrichGO(lethal_genes_entrez[[i]], OrgDb = org_db[i], 
                   keyType = "ENTREZID", ont = "BP", pvalueCutoff = 0.001, minGSSize = 20)@result
  temp <- temp[temp$p.adjust < 0.05, ]
  return(temp)
})

lapply(seq_along(species), function(i) {
  enrich[[i]]$organism <<- species[i]
})

enrich_exclusive <- do.call(rbind, enrich_exclusive)

write.table(enrich_exclusive[,c("ID", "Description", "organism")], file = "GO_enrichment/enrich_exclusive.csv",
            quote = F, row.names = F, sep = ";")

save(terms_cdm, terms_all, go_yeast, enrich_30, enrich_30, enrich_exclusive, enrich_mild, enrich_late,
     file = "GO_enrichment/enrich_tables.RData")

