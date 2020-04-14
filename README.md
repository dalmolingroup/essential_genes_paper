# Evolutionary analysis of essential genes in eukaryote model organisms

This repository has the scripts used in the evolutionary analysis of essential genes in *Saccharomyces cerevisiae*, *Drosophila melanogaster*, *Mus musculus* e *Caenorhabditis elegans*.

## Datasets selection
Essential genes datasets selection for each organism were conducted as follows (also depicted in Supplementary Methods):

- Saccharomyces cerevisiae
S. cerevisiae information was retrieved from Saccharomyces Genome Database (www.yeastgenome.org/, genome version R64-2-1). Phenotype-genotype information can be found at Downloads > curation > literature > phenotype_data.tab (https://downloads.yeastgenome.org/curation/literature/phenotype_data.tab). We selected the genes reported as producing the “inviable” phenotype. 

- Drosophila melanogaster
D. melanogaster data was retrieved from Flybase (flybase.org, FB2018_06 Dmel release 6.25). Phenotype-genotype information can be found at Downloads > releases > precomputed_files > alleles > allele_phenotypic_data_fb_2019_01.tsv.gz (ftp://ftp.flybase.net/releases/FB2019_01/precomputed_files/alleles/allele_phenotypic_data_fb_2019_01.tsv.gz). We selected only genes reported as producing the “lethal” phenotype. 
 
- Caenorhabditis elegans  
C. elegans data was retrieved from Wormbase (wormbase.org, release WS270). Steps for retrieving genes classified as lethal are the following:
 - On the WormBase site, go to “search directory” and change the search for “Phenotype”. Then we searched for “lethal”.
 - We selected the genes related to “embryonic lethal”, “larval lethal”, and “adult lethal”. This list of genes includes knock-out (marked as “Variants”)  and knock-down genes (marked as “RNAi”) and is saved as “c_elegans_filtered.csv”. 
  
- Mus musculus  
M. musculus information was retrieved from Mouse Genome Database (www.informatics.jax.org/, release 6.13).  In mouse, the phenotype produced by a mutation is classified under the Mammalian Phenotype Ontology (MPO). Phenotype data can be found at Downloads > Alleles and Phenotypes > MGI_PhenoGenoMP.rpt (http://www.informatics.jax.org/downloads/reports/MGI_GenePheno.rpt). We selected genes with MPO phenotype marked as lethal with complete penetrance. After identifying the lethal genes, we manually curated the mouse lethality categories based on MPO. At total, 16 MPO categories were considered and classified as following:

  - Early lethality: genes that cause lethality from fertilization to embryo implantation
  - Mild lethality: genes that cause lethality from embryo implantation to birth;
  - Late lethality: genes that cause lethality right after birth (a few weeks after pup birth).
  
  Mouse categories table can be found at `scripts_geneplast/mouse/ontologies_lethal_categories_final.csv`
  
  ## Directories in this repository
  This repository is organized in the following directories:
  - `scripts_geneplast`: holds geneplast scripts for all the considered organisms;
  - `GO_enrichment`: performs enrichment analysis using Gene Ontology;
  - `network_properties`: holds script for node degree and betweenness analysis of PPI networks;
  - `plot`: holds the code for the figures in the paper; 
  - `datasets`: stores the files acquires in the *Datasets selection* section. 
  
  
  
