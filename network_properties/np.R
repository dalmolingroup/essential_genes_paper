library(igraph)
library(dplyr)
library(tidyr)
library(ggridges)
library(ggplot2)
library(purrr)
library(ggpubr)
library(scales)

# Species comprised in our study
species <- c("celegans", "drosophila", "mouse", "yeast", "spombe")

# Load interaction tables
invisible(lapply(species, function(i) {
  tryCatch( 
    assign(i, read.csv(paste0("network_properties/", i, "_int.csv"), stringsAsFactors = F, header = T), 
           envir = .GlobalEnv),
    error = function(e) NULL
  )
}))

# Load table with all inferred roots
load("scripts_geneplast/table_all_org.RData")

# Interaction networks ----------------------------------------------------

# Remove duplicates function
remove_duplicates <- function(df) {
  idx <- duplicated(t(apply(df, 1, sort)))
  df <- df[ idx ,]
  return(df)
}

lapply(species, function (x) {
  assign(paste0(x, "_filtered"), remove_duplicates(get(x)), envir = .GlobalEnv)
})

# Create networks
invisible(lapply(species, function (x) {
  assign(paste("g", x, sep = "_"), graph_from_edgelist( as.matrix(  
    get(paste(x, "filtered", sep = "_"))[, 1:2] ), directed = F)
    ,  envir = .GlobalEnv)
}))


# Calculate network properties --------------------------------------------

# Node degree ----
dg <- lapply(species, function(j) {
  protein <- names(degree(get(paste("g", j, sep = "_"))))
  degree <- degree(get(paste("g", j, sep = "_")))
  organism <- j
  return(data.frame(protein, degree, organism, stringsAsFactors = F))
})

dg <- do.call(rbind, dg)
dg <- inner_join(dg, table[, c("ensembl_peptide_id", "Root", "lethal_nonlethal")], 
                 by = c("protein" = "ensembl_peptide_id"))

# Filter outlier values for degree
dg <- dg %>% 
  group_by(organism, lethal_nonlethal) %>% 
  filter(degree < mean(degree) * 5)

density_plot_degree <- ggplot(dg, aes(x = degree,  y = organism, fill = lethal_nonlethal)) +
  geom_density_ridges(alpha = 0.5, scale = 0.95) +
  scale_x_log10() +
  theme_ridges()

scatter_degree <- ggplot(dg, aes(x = organism, y = degree, col = lethal_nonlethal, group = lethal_nonlethal)) + 
  geom_point(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.5), 
             alpha = 0.5) +
  # stat_summary(fun.data = mean_sdl, fun.args = list(mult = 1), geom = "errorbar", 
  #              color = "black", width = 0.1, position = position_dodge(width = 0.5)) +
  # stat_summary(fun.y = mean, geom = "point", color = "black", 
  #              position = position_dodge(width = 0.5)) +
  guides(col = guide_legend(title = "Category", override.aes = aes(label = ""))) +
  stat_compare_means(method = "wilcox.test", label = "p.signif", label.x = 1.5 )

boxplot_degree <- ggplot(dg, aes(x = organism, y = degree, fill = lethal_nonlethal)) +
  geom_boxplot() +
  stat_compare_means(method = "wilcox.test", label = "p.signif", label.x = 1.5 )

walk(c("density_plot_degree", "scatter_degree", "boxplot_degree"), function (x) {
  ggsave(paste0(x, ".pdf"), plot = get(x), width = 10, height = 7)
})

walk(c("density_plot_degree", "scatter_degree", "boxplot_degree"), function (x) {
  ggsave(paste0(x, ".png"), plot = get(x), width = 10, height = 7)
})

# Betweeness --------------------------------------------------------------

btw <- lapply(species, function (j) {
  protein <- names(betweenness(get(paste("g", j, sep = "_")), directed = F, normalized = T))
  btw <- betweenness(get(paste("g", j, sep = "_")), directed = F, normalized = T)
  organism <- j
  return(data.frame(protein, btw, organism, stringsAsFactors = F))
})

btw <- do.call(rbind, btw)
btw <- inner_join(btw, table[, c("ensembl_peptide_id", "Root", "lethal_nonlethal")], 
                  by = c("protein" = "ensembl_peptide_id"))

btw <- btw %>% 
  group_by(organism, lethal_nonlethal) %>% 
  filter(btw < mean(btw) * 5)

density_plot_btw <- ggplot(btw, aes(x = btw,  y = organism, fill = lethal_nonlethal)) +
  geom_density_ridges(alpha = 0.5, scale = 0.9) +
  scale_x_sqrt() +
  theme_ridges()

scatter_btw <- ggplot(btw, aes(x = organism, y = btw, col = lethal_nonlethal, group = lethal_nonlethal)) + 
  geom_point(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.5), 
             alpha = 0.5) +
  # stat_summary(fun.data = mean_sdl, fun.args = list(mult = 1), geom = "errorbar", 
  #              color = "black", width = 0.1, position = position_dodge(width = 0.5)) +
  # stat_summary(fun.y = mean, geom = "point", color = "black", 
  #              position = position_dodge(width = 0.5)) +
  guides(col = guide_legend(title = "Category", override.aes = aes(label = ""))) +
  stat_compare_means(method = "wilcox.test", label = "p.signif", label.x = 1.5 )

boxplot_btw <- ggplot(btw, aes(x = organism, y = btw, fill = lethal_nonlethal)) +
  geom_boxplot() +
  stat_compare_means(method = "wilcox.test", label = "p.signif", label.x = 1.5 )

walk(c("density_plot_btw", "scatter_btw", "boxplot_btw"), function (x) {
  ggsave(paste0(x, ".pdf"), plot = get(x), width = 10, height = 7)
})

walk(c("density_plot_btw", "scatter_btw", "boxplot_btw"), function (x) {
  ggsave(paste0( x, ".png"), plot = get(x), width = 10, height = 7)
})

# Dynamites ---------------------------------------------------------------
pos.d <- position_dodge(width = 0.9)
# dg$organism <- factor(dg$organism, levels = c("celegans", "drosophila", "mouse", "yeast", "spombe"),
#                       labels = c("Caenorhabditis \nelegans", "Drosophila \nmelanogaster", 
#                                  "Mus \nmusculus", "Saccharomyces \ncerevisiae", 
#                                  "Schizosaccharomyces \npombe"))
ggplot(dg, aes(x = organism, y = degree, fill = factor(lethal_nonlethal)) ) +
  geom_bar(stat = "summary", fun.y = "mean", position = pos.d) +
  stat_summary(fun.data = mean_se, fun.args = list(mult = 1), geom = "errorbar", 
               lwd = 0.6, position = pos.d, width = 0.2) +
  scale_fill_manual(values = c(lethal = "#ff4a4aff", nonlethal = "#3939c0ff"), 
                    labels = c("Essential", "Others"), guide = guide_legend("Essentiality")) +
  labs(x = "", y = "<k>") +
  theme_classic() + 
  theme(axis.text.x = element_text(face = "italic", size = 11), 
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 13),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title = element_text(face = "plain", hjust = 0.5)) +
  stat_compare_means(method = "wilcox.test", label = "p.signif", label.x = 1.5, label.y = 120)
ggsave("dyn_degree.pdf", width = 7, height = 4)
ggsave("dyn_degree.png", width = 7, height = 4)

pos.d <- position_dodge(width = 0.9)
# btw$organism <- factor(btw$organism, levels = c("celegans", "drosophila", "mouse", "yeast", "spombe"),
#                        labels = c("Caenorhabditis \nelegans", "Drosophila \nmelanogaster", 
#                                   "Mus \nmusculus", "Saccharomyces \ncerevisiae", "Schizosaccharomyces \npombe"))
ggplot(btw, aes(x = organism, y = btw, fill = factor(lethal_nonlethal)) ) +
  geom_bar(stat = "summary", fun.y = "mean", position = pos.d) +
  stat_summary(fun.data = mean_se, fun.args = list(mult = 1), geom = "errorbar", 
               lwd = 0.6, position = pos.d, width = 0.2) +
  scale_fill_manual(values = c(lethal = "#ff4a4aff", nonlethal = "#3939c0ff"), 
                    labels = c("Essential", "Others"), guide = guide_legend("Essentiality")) +
  labs(x = "", y = "Beetweeness") +
  theme_classic() + 
  theme(axis.text.x = element_text(face = "italic", size = 11), 
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 13),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title = element_text(face = "plain", hjust = 0.5)) +
  stat_compare_means(method = "wilcox.test", label = "p.signif", label.x = 1.5, label.y = 0.001)
ggsave("network_properties/dyn_btw.pdf", width = 7, height = 4)
ggsave("network_properties/dyn_btw.png", width = 7, height = 4)

# plot grid ----
ggarrange(density_plot_degree, 
          ggarrange(boxplot_degree, scatter_degree, ncol = 2, labels = c("B", "C")),
          nrow = 2, labels = "A")
ggsave(filename = "comb_degree.pdf", width = 10, height = 7)
ggsave(filename = "comb_degree.png", width = 10, height = 7)

ggarrange(density_plot_btw, 
          ggarrange(boxplot_btw, scatter_btw, ncol = 2, labels = c("B", "C")),
          nrow = 2, labels = "A")
ggsave(filename = "network_properties/comb_btw.pdf", width = 10, height = 7)
ggsave(filename = "network_properties/comb_btw.png", width = 10, height = 7)

# Mouse -------------------------------------------------------------------
load("scripts_geneplast/mouse/categories/working_lethal_mouse.RData")

dg_mouse <- inner_join(dg, working_lethal_mouse[, c("ensembl_peptide_id", "categories")], 
                       by = c("protein" = "ensembl_peptide_id"))

dg_mouse$categories <- factor(dg_mouse$categories, levels = c("Early lethality", "Mild lethality", "Late lethality"))

density_plot_degree_mouse <- ggplot(dg_mouse, aes(x = degree,  y = categories, fill = categories)) +
  geom_density_ridges(alpha = 0.5, scale = 0.95) +
  #scale_x_log10() +
  theme_ridges()+
  theme(legend.position = "none")

scatter_degree_mouse <- ggplot(dg_mouse, aes(x = categories, y = degree, col = categories, group = categories)) + 
  geom_point(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.5), 
             alpha = 0.5) +
  # stat_summary(fun.data = mean_sdl, fun.args = list(mult = 1), geom = "errorbar", 
  #              color = "black", width = 0.1, position = position_dodge(width = 0.5)) +
  # stat_summary(fun.y = mean, geom = "point", color = "black", 
  #              position = position_dodge(width = 0.5)) +
  #geom_text(aes(label=""), show.legend = F) +
  stat_compare_means(method = "wilcox.test", label = "p.signif", 
                     comparisons = list(c("Early lethality", "Mild lethality"), 
                                        c("Early lethality", "Late lethality"), c("Mild lethality", "Late lethality"))) +
  theme(legend.position = "none")

boxplot_degree_mouse <- ggplot(dg_mouse, aes(x = categories, y = degree, fill = categories)) +
  geom_boxplot() +
  stat_compare_means(method = "wilcox.test", label = "p.signif", 
                     comparisons = list(c("Early lethality", "Mild lethality"), 
                                        c("Early lethality", "Late lethality"), c("Mild lethality", "Late lethality"))) +
  theme(legend.position = "none")

walk(c("density_plot_degree_mouse", "scatter_degree_mouse", "boxplot_degree_mouse"), function (x) {
  ggsave(paste0(x, ".pdf"), plot = get(x), width = 10, height = 7)
})

walk(c("density_plot_degree_mouse", "scatter_degree_mouse", "boxplot_degree_mouse"), function (x) {
  ggsave(paste0(x, ".png"), plot = get(x), width = 10, height = 7)
})

# beetweeness -----
btw_mouse <- inner_join(btw, working_lethal_mouse[, c("ensembl_peptide_id", "categories")], 
                        by = c("protein" = "ensembl_peptide_id"))

btw_mouse$categories <- factor(btw_mouse$categories, levels = c("Early lethality", "Mild lethality", "Late lethality"))

density_plot_btw_mouse <- ggplot(btw_mouse, aes(x = btw,  y = categories, fill = categories)) +
  geom_density_ridges(alpha = 0.5, scale = 0.9) +
  #scale_x_sqrt() +
  theme_ridges()+
  theme(legend.position = "none")

scatter_btw_mouse <- ggplot(btw_mouse, aes(x = categories, y = btw, col = categories, group = categories)) + 
  geom_point(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.5), 
             alpha = 0.5) +
  # stat_summary(fun.data = mean_sdl, fun.args = list(mult = 1), geom = "errorbar", 
  #              color = "black", width = 0.1, position = position_dodge(width = 0.5)) +
  # stat_summary(fun.y = mean, geom = "point", color = "black", 
  #              position = position_dodge(width = 0.5)) +
  #guides(col = guide_legend(title = "Category", override.aes = aes(label = ""))) +
  stat_compare_means(method = "wilcox.test", label = "p.signif", 
                     comparisons = list(c("Early lethality", "Mild lethality"), 
                                        c("Early lethality", "Late lethality"), c("Mild lethality", "Late lethality"))) +
  theme(legend.position = "none")

boxplot_btw_mouse <- ggplot(btw_mouse, aes(x = categories, y = btw, fill = categories)) +
  geom_boxplot() +
  stat_compare_means(method = "wilcox.test", label = "p.signif", 
                     comparisons = list(c("Early lethality", "Mild lethality"), 
                                        c("Early lethality", "Late lethality"), c("Mild lethality", "Late lethality")),
  ) +
  theme(legend.position = "none")

walk(c("density_plot_btw_mouse", "scatter_btw_mouse", "boxplot_btw_mouse"), function (x) {
  ggsave(paste0("mouse/", x, ".pdf"), plot = get(x), width = 10, height = 7)
})

walk(c("density_plot_btw_mouse", "scatter_btw_mouse", "boxplot_btw_mouse"), function (x) {
  ggsave(paste0("mouse/", x, ".png"), plot = get(x), width = 10, height = 7)
})

# plot grids ----
ggarrange(density_plot_degree_mouse, 
          ggarrange(boxplot_degree_mouse, scatter_degree_mouse, ncol = 2, labels = c("B", "C")),
          nrow = 2, labels = "A")
ggsave(filename = "mouse_degree.pdf", width = 10, height = 7)
ggsave(filename = "mouse_degree.png", width = 10, height = 7)

ggarrange(density_plot_btw_mouse, 
          ggarrange(boxplot_btw_mouse, scatter_btw_mouse, ncol = 2, labels = c("B", "C")),
          nrow = 2, labels = "A")
ggsave(filename = "mouse_btw.pdf", width = 10, height = 7)
ggsave(filename = "mouse_btw.png", width = 10, height = 7)


# Network properties by percentiles ---------------------------------------

species <- c("celegans", "drosophila", "mouse", "spombe", "yeast")

# Get the upper 30 pencentile
percent_30 <- sapply(split(table$Root, table$organism), function (i) {
  quantile(i, probs = 0.3)
})
names(percent_30) <- species

percent_70 <- sapply(split(table$Root, table$organism), function (i) {
  quantile(i, probs = 0.7)
})
names(percent_70) <- species

dg$percent <- ""
lapply(species, function (i) {
  dg$percent[dg$organism == i] <<- ifelse(dg$Root[dg$organism == i] <= percent_30[i], "percent_30", 
                                          ifelse(dg$Root[dg$organism == i] >= percent_70[i], "percent_70", NA))
})
dg_percent <- na.omit(dg)

btw$percent <- ""
lapply(species, function (i) {
  btw$percent[btw$organism == i] <<- ifelse(btw$Root[btw$organism == i] <= percent_30[i], "percent_30", 
                                            ifelse(btw$Root[btw$organism == i] >= percent_70[i], "percent_70", NA))
})
btw_percent <- na.omit(btw)

dg_percent$organism <- factor(dg_percent$organism, levels = c("yeast","spombe", "drosophila", "celegans", "mouse"),
                              labels =  c("S. cerevisiae", "S. pombe", "D. melanogaster", "C. elegans", "M. musculus") )

btw_percent$organism <- factor(btw_percent$organism, levels = c("yeast", "spombe", "drosophila", "celegans", "mouse"),
                               labels =  c("S. cerevisiae", "S. pombe", "D. melanogaster", "C. elegans", "M. musculus") )

# BY PERCENT COMPARISON ----
# degree - percent 30
dg_percent30 <- ggplot(dg_percent[dg_percent$percent == "percent_30",], aes(x = organism, y = degree, fill = factor(lethal_nonlethal)) ) +
  geom_bar(stat = "summary", fun.y = "mean", position = pos.d) +
  stat_summary(fun.data = mean_se, fun.args = list(mult = 1), geom = "errorbar", 
               lwd = 0.6, position = pos.d, width = 0.2) +
  scale_fill_manual(values = c(lethal = "#ff4a4aff", nonlethal = "#3939c0ff"), 
                    labels = c("Essential", "Others"), guide = guide_legend("")) +
  labs(title = "degree - young", x = "", y = "<k>") +
  theme_classic() + 
  theme(axis.text.x = element_text(face = "italic", size = 8), 
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 13),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title = element_text(face = "plain", hjust = 0.5)) +
  stat_compare_means(method = "wilcox.test", label = "p.signif", label.x = 1.5, label.y = 200)

# degree - percent 70
dg_percent70 <- ggplot(dg_percent[dg_percent$percent == "percent_70",], aes(x = organism, y = degree, fill = factor(lethal_nonlethal)) ) +
  geom_bar(stat = "summary", fun.y = "mean", position = pos.d) +
  stat_summary(fun.data = mean_se, fun.args = list(mult = 1), geom = "errorbar", 
               lwd = 0.6, position = pos.d, width = 0.2) +
  scale_fill_manual(values = c(lethal = "#ff4a4aff", nonlethal = "#3939c0ff"), 
                    labels = c("Essential", "Others"), guide = guide_legend("")) +
  labs(title = "degree - old", x = "", y = "<k>") +
  theme_classic() + 
  theme(axis.text.x = element_text(face = "italic", size = 8), 
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 13),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title = element_text(face = "plain", hjust = 0.5)) +
  stat_compare_means(method = "wilcox.test", label = "p.signif", label.x = 1.5, label.y = 200)

# betweenness - percent 30
btw_percent30 <- ggplot(btw_percent[btw_percent$percent == "percent_30",], aes(x = organism, y = btw, fill = factor(lethal_nonlethal)) ) +
  geom_bar(stat = "summary", fun.y = "mean", position = pos.d) +
  stat_summary(fun.data = mean_se, fun.args = list(mult = 1), geom = "errorbar", 
               lwd = 0.6, position = pos.d, width = 0.2) +
  scale_fill_manual(values = c(lethal = "#ff4a4aff", nonlethal = "#3939c0ff"), 
                    labels = c("Essential", "Others"), guide = guide_legend("")) +
  scale_y_continuous(labels = comma) +
  labs(title = "btw - young", x = "", y = "Betweenness") +
  theme_classic() + 
  theme(axis.text.x = element_text(face = "italic", size = 8), 
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 13),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title = element_text(face = "plain", hjust = 0.5)) +
  stat_compare_means(method = "wilcox.test", label = "p.signif", label.x = 1.5, label.y = 1.2e-3)

# betweenness - percent 70
btw_percent70 <- ggplot(btw_percent[btw_percent$percent == "percent_70",], aes(x = organism, y = btw, fill = factor(lethal_nonlethal)) ) +
  geom_bar(stat = "summary", fun.y = "mean", position = pos.d) +
  stat_summary(fun.data = mean_se, fun.args = list(mult = 1), geom = "errorbar", 
               lwd = 0.6, position = pos.d, width = 0.2) +
  scale_fill_manual(values = c(lethal = "#ff4a4aff", nonlethal = "#3939c0ff"), 
                    labels = c("Essential", "Others"), guide = guide_legend("")) +
  scale_y_continuous(labels = comma) + 
  labs(title = "btw - old", x = "", y = "Betweenness") +
  theme_classic() + 
  theme(axis.text.x = element_text(face = "italic", size = 8), 
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 13),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title = element_text(face = "plain", hjust = 0.5)) +
  stat_compare_means(method = "wilcox.test", label = "p.signif", label.x = 1.5, label.y = 1.2e-3)


ggarrange(dg_percent30, dg_percent70, btw_percent30, btw_percent70, ncol = 2, nrow = 2, common.legend = T, legend = "right")
ggsave("network_properties/percentiles_percent_np.pdf", width = 10, height = 5)
ggsave("network_properties/percentiles_percent_np.svg", width = 10, height = 5)

# BY ESSENTIALITY COMPARISON ----

# degree - essential
dg_percent_essential <- ggplot(dg_percent[dg_percent$lethal_nonlethal == "lethal",], aes(x = organism, y = degree, fill = factor(percent)) ) +
  geom_bar(stat = "summary", fun.y = "mean", position = pos.d) +
  stat_summary(fun.data = mean_se, fun.args = list(mult = 1), geom = "errorbar", 
               lwd = 0.6, position = pos.d, width = 0.2) +
  scale_fill_manual(values = c(percent_30 = "#BFA454", percent_70 = "#D93D1A"), 
                    labels = c("Young", "Old"), guide = guide_legend("")) +
  labs(title = "degree - essential",x = "", y = "<k>") +
  theme_classic() + 
  theme(axis.text.x = element_text(face = "italic", size = 8), 
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 13),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title = element_text(face = "plain", hjust = 0.5)) +
  stat_compare_means(method = "wilcox.test", label = "p.signif", label.x = 1.5, label.y = 200)

# degree - nonessential
dg_percent_nonessential <- ggplot(dg_percent[dg_percent$lethal_nonlethal == "nonlethal",], aes(x = organism, y = degree, fill = factor(percent)) ) +
  geom_bar(stat = "summary", fun.y = "mean", position = pos.d) +
  stat_summary(fun.data = mean_se, fun.args = list(mult = 1), geom = "errorbar", 
               lwd = 0.6, position = pos.d, width = 0.2) +
  scale_fill_manual(values = c(percent_30 = "#BFA454", percent_70 = "#D93D1A"), 
                    labels = c("Young", "Old"), guide = guide_legend("")) +
  labs(title = "degree - others", x = "", y = "<k>") +
  theme_classic() + 
  theme(axis.text.x = element_text(face = "italic", size = 8), 
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 13),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title = element_text(face = "plain", hjust = 0.5)) +
  stat_compare_means(method = "wilcox.test", label = "p.signif", label.x = 1.5, label.y = 200)

# betweenness - essential
btw_percent_essential <- ggplot(btw_percent[btw_percent$lethal_nonlethal == "lethal",], aes(x = organism, y = btw, fill = factor(percent)) ) +
  geom_bar(stat = "summary", fun.y = "mean", position = pos.d) +
  stat_summary(fun.data = mean_se, fun.args = list(mult = 1), geom = "errorbar", 
               lwd = 0.6, position = pos.d, width = 0.2) +
  scale_fill_manual(values = c(percent_30 = "#BFA454", percent_70 = "#D93D1A"), 
                    labels = c("Young", "Old"), guide = guide_legend("")) +
  scale_y_continuous(labels = comma) +
  labs(title = "btw - essential", x = "", y = "Betweenness") +
  theme_classic() + 
  theme(axis.text.x = element_text(face = "italic", size = 8), 
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 13),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title = element_text(face = "plain", hjust = 0.5)) +
  stat_compare_means(method = "wilcox.test", label = "p.signif", label.x = 1.3, label.y = 1.2e-3)

# betweenness - nonessential
btw_percent_nonessential <- ggplot(btw_percent[btw_percent$lethal_nonlethal == "nonlethal",], aes(x = organism, y = btw, fill = factor(percent)) ) +
  geom_bar(stat = "summary", fun.y = "mean", position = pos.d) +
  stat_summary(fun.data = mean_se, fun.args = list(mult = 1), geom = "errorbar", 
               lwd = 0.6, position = pos.d, width = 0.2) +
  scale_fill_manual(values = c(percent_30 = "#BFA454", percent_70 = "#D93D1A"), 
                    labels = c("Young", "Old"), guide = guide_legend("")) +
  scale_y_continuous(labels = comma) +
  labs(title = "btw - others", x = "", y = "Betweenness") +
  theme_classic() + 
  theme(axis.text.x = element_text(face = "italic", size = 8), 
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 13),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title = element_text(face = "plain", hjust = 0.5)) +
  stat_compare_means(method = "wilcox.test", label = "p.signif", label.x = 1.5, label.y = 1.2e-3)

ggarrange(dg_percent_essential, dg_percent_nonessential, btw_percent_essential, btw_percent_nonessential, ncol = 2, nrow = 2, common.legend = T, legend = "right")
ggsave("network_properties/percentiles_essential_np.pdf", width = 12, height = 7)
ggsave("network_properties/percentiles_essential_np.svg", width = 12, height = 7)

