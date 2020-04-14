library(igraph)
library(GO.db)
library(GOSemSim)
library(TreeAndLeaf)


# Dendogram function
GOdendogram <- function(GOs, cutoff = 0, host = "127.0.0.1", port = 9091){
  xx <- as.list(GO.db::GOTERM)
  terms <- unique(GOs)
  desc <- sapply(terms, function(i){xx[[i]]@Term})
  if(length(terms)==0){
    stop("there are no enriched terms!")
  }
  cutoff <- ceiling(length(terms)*cutoff)
  message("computing similarity matrix... step 1 of 4")
  semData <- GOSemSim::godata(ont = "BP")
  mSim <- GOSemSim::mgoSim(terms, terms, semData = semData, measure = "Wang", combine = NULL)
  if(cutoff > 0){
    RS <- rowSums(mSim) - 1
    GOs <- names(sort(RS)[-1:-cutoff])
    mSim <- mSim[rownames(mSim) %in% GOs, colnames(mSim) %in% GOs]
    terms <- terms[terms %in% GOs]
    GOs <- terms
  }
  message("generating graph from matrix... step 2 of 4")
  hc <- hclust(dist(mSim))
  gg <- TreeAndLeaf::hclust2igraph(hc)
  graphSize <- length(terms)
  if(graphSize < 100){
    graphSize <- "S"
  }else if(graphSize < 500){
    graphSize <- "M"
  }else{
    graphSize <- "L"
  }
  graph <- suppressWarnings(TreeAndLeaf::formatTree(gg, as.data.frame(terms, stringsAsFactors = FALSE),
                                                    theme = paste0(graphSize, "SReds"),
                                                    nodeFontSize = 40,
                                                    nodeFontColor = "black", edgeColor = "grey80"))
  igraph::V(graph)$nodeAlias[match(terms,igraph::V(graph)$name)] <- desc[terms]
  igraph::V(graph)$nodeSize[match(terms,igraph::V(graph)$name)] <- 100
  rdp <- RedeR::RedPort(host = host, port = port)
  message("invoking RedeR... step 3 of 4")
  RedeR::calld(rdp)
  RedeR::resetd(rdp)
  if(graphSize == "S"){
    graphSize <- "small"
  }else if(graphSize == "M"){
    graphSize <- "medium"
  }else{
    graphSize <- "large"
  }
  igraph::V(graph)$nodeLineColor <- "grey80"
  message("adding graph into RedeR... step 4 of 4")
  suppressMessages(suppressWarnings(TreeAndLeaf::treeAndLeaf(rdp, graph, size = graphSize)))
  message("done!")
  return(rdp)
}

# IDS 
load("~/part_1/analysis/enrichment/enrich_tables.RData")

GO_ids_multicellular <- terms_cdm$go_id
GO_ids_yeast <- go_yeast$ID

rdp <- GOdendogram(GO_ids_multicellular)
rdp <- GOdendogram(GO_ids_yeast)
