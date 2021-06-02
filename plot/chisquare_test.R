library(dplyr)
library(tidyr)
library(ggplot2)

# Load table
load("scripts_geneplast/table_all_org.RData")

# Define quantiles for old and new genes
species <- c("celegans", "drosophila", "mouse",  "spombe", "yeast")
percent30 <- sapply(split(table$Root, table$organism), function (i) {
  quantile(i, probs = 0.3)
})
percent70 <- sapply(split(table$Root, table$organism), function (i) {
  quantile(i, probs = 0.7)
})
names(percent30) <- species
names(percent70) <- species

# Select the older and newer genes
table30 <- lapply(species, function (i) {
  temp <- table[table$organism == i,]
  temp2 <- dplyr::slice(temp, which(temp$Root <= percent30[i]))
  temp2
})
table30 <- do.call(rbind, table30)
table30$percent <- "percent30"

table70 <- lapply(species, function (i) {
  temp <- table[table$organism == i,]
  temp2 <- dplyr::slice(temp, which(temp$Root >= percent70[i]))
  temp2
})
table70 <- do.call(rbind, table70)
table70$percent <- "percent70"

percentiles <- rbind(table30, table70)

# Apply Chi-squared test of independence
tables_list <- lapply(species, function (i) {
  temp <- percentiles[percentiles$organism == i,]
  table(temp$percent, temp$lethal_nonlethal)
})
names(tables_list) <- species

lapply(tables_list, function(i) {
  chisq.test(i, correct = T)
})

percentiles$organism <- factor(percentiles$organism, levels = c("yeast", "drosophila", "celegans", "mouse"), 
                               labels =  c("S. cerevisiae", "D. melanogaster", "C. elegans", "M. musculus") )
percentiles %>% 
  group_by(organism, percent) %>% 
  summarise(lethal_prop = length(lethal_nonlethal[lethal_nonlethal == "lethal"])/ n()) %>% 
  ggplot(aes(x = organism, y = lethal_prop, fill = percent)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9)) + 
  scale_fill_manual(values = c(percent30 = "#BFA454", percent70 = "#D93D1A"),
                    labels = c("Young", "Old"), guide = guide_legend("")) +
  scale_x_discrete(name = "") +
  scale_y_continuous("% of essential genes") +
  ggtitle("Proportion of essential genes in young and old roots") +
  theme_classic() + 
  theme(axis.text.x = element_text(face = "italic", size = 8), 
        axis.text.y = element_text(size = 10),
        axis.title = element_text(size = 10),
        title = element_text(size = 8),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank())

ggsave(file = "plot/proportion_essential_genes.pdf", width = 5, height = 3)









