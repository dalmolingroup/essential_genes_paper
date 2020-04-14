# PLOTS -------------------------------------------------------------------

library(ggplot2)
library(tidyverse)
library(ggpubr)
library(UpSetR)
library(grid)
library(VennDiagram)

load("scripts_geneplast/table_all_org.RData")

# Ancestry plot ----
table$organism <- factor(table$organism, levels = c("yeast", "drosophila", "celegans", "mouse"),
                         labels =  c("S. cerevisiae", "D. melanogaster", "C. elegans", "M. musculus"))
# Violin plot
pos.d <- position_dodge(width = 0.9)
plot2 <- ggplot(table, aes(x = organism, y = ancestry, fill = factor(lethal_nonlethal))) +
  geom_violin(position = pos.d) +
  stat_summary(fun.y = mean, geom = "point", position = pos.d, cex = 3, show.legend = F) + 
  scale_fill_manual(values = c(lethal = "#ff4a4aff", nonlethal = "#3939c0ff"), 
                    labels = c("Essential", "Others"), guide = guide_legend("Essentiality")) +
  scale_x_discrete("") +
  scale_y_continuous("Ancestry", limits = c(0, 1.1), breaks = seq(0, 1, 0.1)) +
  theme_classic() + 
  theme(axis.text.x = element_text(face = "italic", size = 10), 
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title = element_text(face = "plain", hjust = 0.5))

plot2 <- plot2 +
  stat_compare_means(method = "wilcox.test", label = "p.signif", label.x = 1.5)

plot2

ggsave(plot = plot2, filename = "plot/ancestry_plot_violin.pdf",
       width = 7, height = 5)
ggsave(plot = plot2, filename = "plot/ancestry_plot_violin.png",
       width = 7, height = 5)
ggsave(plot = plot2, filename = "plot/ancestry_plot_violin.svg",
       width = 7, height = 5)

# Ancestry for exclusive orthologs ----
load("scripts_geneplast/table_all_org.RData")

organisms <- c("celegans", "drosophila", "mouse", "yeast")
filter_cogs <- function(table, lethality) {
  table <- table[table$lethal_nonlethal == lethality,]
  cogs <- lapply(split(table$cog_id, table$organism), unique)
  cogs
}
cogs_list <- filter_cogs(table, "lethal")
names(cogs_list) <- organisms
essential_exclusive <- lapply(organisms, function (i) {
  unique(setdiff(cogs_list[[i]], unique(unlist(unname(cogs_list[ !grepl(i, names(cogs_list)) ] )))))
})
names(essential_exclusive) <- organisms
table$exclusive <- NA

table <- table[table$lethal_nonlethal == "lethal",]
combined_cogs <- reduce(cogs_list, intersect)

invisible(lapply(organisms, function (i) {
  table$exclusive[table$organism == i] <<- ifelse(table$cog_id[table$organism == i & table$lethal_nonlethal == "lethal"] %in% essential_exclusive[[i]], "E", NA )
}))
table$exclusive[table$cog_id %in% combined_cogs] <- "C"
table <- na.omit(table)

n_distinct(table$cog_id[table$exclusive == "C"]) == length(combined_cogs)

table$organism <- factor(table$organism, levels = c("yeast", "drosophila", "celegans", "mouse"), 
                         labels =  c("S. cerevisiae", "D. melanogaster", "C. elegans", "M. musculus") )

pos.d <- position_dodge(width = 0.5)
plot1 <- ggplot(table, aes(x = organism, y = ancestry, col = factor(exclusive))) +
  stat_summary(aes(x = organism, y = ancestry, col = factor(exclusive)),
               fun.data = mean_se, fun.args = list(mult = 1), geom = "errorbar", 
               lwd = 1.5, position = pos.d, width = 0.2) +
  stat_summary(aes(x = organism, y = ancestry, group = factor(exclusive)),
               fun.y = mean, geom = "point", position = pos.d, cex = 1.2, color = "black") + 
  scale_color_manual(values = c(E = "#ec3663ff", C = "#18df00ff"),
                     labels = c(E = "Exclusive", C = "Shared"), guide = guide_legend("Exclusivity")) +
  guides(colour = guide_legend(override.aes = list(alpha=1, size = 5))) +
  theme_classic() + 
  stat_compare_means(method = "wilcox.test", label = "p.signif", label.x = 1.5, label.y = 1.1) +
  labs("", x = "", y = "Ancestry", color = "Exclusivity") +
  theme(axis.text.x = element_text(face = "italic", size = 7),
        axis.title.y = element_text(size = 8),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title = element_text(face = "plain", hjust = 0.5))
plot1

ggsave(plot = plot1, filename = "plot/ancestry_exclusive_combined_meanerror.pdf",
       width = 5, height = 3)
ggsave(plot = plot1, filename = "plot/ancestry_exclusive_combined_meanerror.png",
       width = 5, height = 3)
ggsave(plot = plot1, filename = "plot/ancestry_exclusive_combined_meanerror.svg",
       width = 5, height = 3)


# Intersection lethal cogs ----
load("scripts_geneplast/table_all_org.RData")

# Get only the lethal genes
organisms <- c("celegans", "drosophila", "mouse", "yeast")
table_lethal <- unique(table[table$lethal_nonlethal == "lethal" &
                               table$organism %in% organisms,])
# Get all unique cogs
cogs <- unique(table_lethal$cog_id)

# Get COGs from all organisms
get_cog <- function(x) {
  unique(table_lethal$cog_id[table_lethal$organism == x])
}
organism_cogs <- lapply(organisms, get_cog)
names(organism_cogs) <- organisms

venn.diagram(x = organism_cogs, filename = "plot/intersect_cogs.pdf")

# Create upset table of intersections
upset_list <- fromList(organism_cogs)

pdf("plot/plot_intersect_all.pdf", width = 10, height = 5)
upset(upset_list, order.by = "freq", nsets = 6, point.size = 3.5,
      line.size = 2, text.scale = c(1.3, 1.3, 1, 1, 2, 1.2))
dev.off()

png("plot/plot_intersect_all.png", width = 20, height = 12, units = "cm", res = 100)
upset(upset_list, order.by = "freq", nsets = 6, point.size = 3.5,
      line.size = 2, text.scale = c(1.3, 1.3, 1, 1, 2, 1.2))
dev.off()


# Ancestry for mouse categories
load("scripts_geneplast/mouse/geneplast_results.RData")

root_genes_mouse$ancestry <- root_genes_mouse$Root / max(root_genes_mouse$Root)

plot4 <- ggplot(root_genes_mouse, aes(x = categories, y = ancestry, fill = categories)) +
  geom_violin(col = "black") + 
  stat_summary(fun.y = mean, fun.args = list(mult = 1), geom = "point", size = 3) +
  scale_fill_manual(values = c("#c17d11ff", "#4e9a06ff", "#3465a4ff")) +
  scale_y_continuous("Ancestry", limits = c(0, 1.3), breaks = seq(0, 1, 0.2)) +
  scale_x_discrete("", limits = c("Early lethality", "Mild lethality", "Late lethality")) +
  guides(fill = FALSE) +
  theme_bw() +
  theme(axis.text = element_text(face = "plain", size = 12, vjust = 0.5),
        axis.ticks.x = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.line.y.left = element_line(colour = "black"),
        axis.line.x = element_blank(),
        
        plot.title = element_text(face = "plain", hjust = 0.5)) 
plot4 <- plot4 + 
  stat_compare_means(method = "wilcox.test", comparisons = list(c("Early lethality", "Mild lethality"), 
                                                                c("Early lethality", "Late lethality"), 
                                                                c("Mild lethality", "Late lethality")), 
                     label = "p.signif")

plot4
ggsave(plot4, filename = "plot/plot_mouse2.pdf",
       width = 5, height = 5)
ggsave(plot4, filename = "plot/plot_mouse2.png",
       width = 5, height = 5)
ggsave(plot4, filename = "plot/plot_mouse2.svg",
       width = 5, height = 5)

