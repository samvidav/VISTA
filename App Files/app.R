library(ggplot2)
theme_set(theme_bw())
library(UniProt.ws)
library(topGO)
library(reticulate)
library(networkDynamic)
library(RColorBrewer)
library(ndtv)
library(tsna)
library(shiny)
library(zip)
set.seed(123)

##### GET AND CLEAN DATA ####

links_filename <- list.files(path=".", pattern= "edges*")[1]
nodes_filename <- "nodes.csv"

links <- read.csv(links_filename, header=T, as.is=T, sep=",", stringsAsFactors = F)
nodes <- read.csv(nodes_filename, header=T, as.is=T, sep=",", stringsAsFactors = F)

# ensure all values are filled
links[is.na(links)] <- 0

# get edge weight threshold from input filename (GUI)
# get numeric portion of filename
CUTOFF_THRESHOLD <- as.numeric(unlist(regmatches(links_filename,
                                                 gregexpr("[[:digit:]]+\\.*[[:digit:]]*", 
                                                          links_filename))))
# weight threshold is set lower than cutoff to prevent artificial on-off behavior
if (CUTOFF_THRESHOLD <= 0.1 ) { 
  WEIGHT_THRESHOLD <- CUTOFF_THRESHOLD 
} else {
  WEIGHT_THRESHOLD <- CUTOFF_THRESHOLD - 0.1
} 

# get timepoint information from column names
num_link_attr <- ncol(links)
timepoints <- unique(as.numeric(unlist(regmatches(colnames(links)[3:num_link_attr],
                                                  gregexpr("[[:digit:]]+\\.*[[:digit:]]*", 
                                                           colnames(links)[3:num_link_attr])))))
NUM_TPS <- length(timepoints)

# set interval length to smallest interval between timepoints
interval <- {}
for (i in 1:(NUM_TPS-1)) {
  interval[i] <- timepoints[i+1] - timepoints[i]
}
INTERVAL <- min(interval)

# rename node columns
LOC_PROV <- F
if (ncol(nodes) == 3) { 
  colnames(nodes) <- c("accession", "gene_name", "taxid") 
} else { 
  LOC_PROV <- T
  colnames(nodes) <- c("accession", "gene_name", "taxid", "localization") 
  }

# rename link columns 
colnames(links)[c(1,2)] <- c("bait_accession", "prey_accession")
ABUND_PROV <- F
if (num_link_attr == NUM_TPS+2) {
  colnames(links)[3:num_link_attr] <- paste("w", timepoints, sep="")
} else { 
  ABUND_PROV <- T
  cnames <- c()
  for (i in 1:NUM_TPS) {
    cnames <- c(cnames, paste("w", timepoints[i], sep=""))
    cnames <- c(cnames, paste("a", timepoints[i], sep=""))
  }
  colnames(links)[3:num_link_attr] <- cnames
}

# only keep links where at least one timepoint is over the cutoff threshold
keep_links <- apply(links[, paste("w", timepoints, sep="")], 1, 
              function(x) { sum(x >= CUTOFF_THRESHOLD) > 0 } )
links <- links[keep_links, ]

# remove duplicate nodes and links
nodes <- nodes[!duplicated(nodes$accession), ]
nodes <- nodes[nodes$accession %in% links$bait_accession | 
                 nodes$accession %in% links$prey_accession, ]
NUM_NODES <- dim(nodes)[1]

links <- links[!duplicated(links[, c("bait_accession", "prey_accession")]), ]
NUM_LINKS <- dim(links)[1]

# assign unique ids to nodes
nodes$id <- c(1:NUM_NODES)
rownames(nodes) <- nodes$id

# attach ids and gene names to links
links <- merge(links, nodes[, c("accession", "gene_name", "id")], by.x="bait_accession", 
               by.y="accession")
colnames(links)[c(ncol(links)-1, ncol(links))] <- c("bait_gene_name", "bait_id")
links <- merge(links, nodes[, c("accession", "gene_name", "id")], by.x="prey_accession", 
               by.y="accession")
colnames(links)[c(ncol(links)-1, ncol(links))] <- c("prey_gene_name", "prey_id")

BAIT_IDS <- unique(links$bait_id)
NUM_BAITS <- length(BAIT_IDS)

TAXIDS <- unique(nodes$taxid)

##### CALCULATE SHARED NODES ####

for (tp in timepoints) {
  shared_tp <- c(mode='numeric', length=NUM_NODES)
  shared_bait_tp <- c("", length=NUM_NODES)
  for (i in 1:NUM_NODES) {
    shared_tp[i] <- sum(links[links$prey_accession == nodes$accession[i], 
                              paste('w', tp, sep="")] >= WEIGHT_THRESHOLD)
    if (shared_tp[i] > 0) {
      shared_baits <- links[which(links$prey_accession == nodes$accession[i] & 
                                    links[, paste('w', tp, sep="")] >= WEIGHT_THRESHOLD), 
                            "bait_gene_name"]
      shared_bait_tp[i] <- paste0(shared_baits, collapse = "; ")
    } else { shared_bait_tp[i] <- "" }
  }
  nodes <- cbind(nodes, shared_tp)
  colnames(nodes)[ncol(nodes)] <- paste("shared_", tp, sep="")
  nodes[, ncol(nodes)] <- as.character(nodes[, ncol(nodes)])
  
  nodes <- cbind(nodes, shared_bait_tp)
  colnames(nodes)[ncol(nodes)] <- paste("shared_baits_", tp, sep="")
  nodes[, ncol(nodes)] <- as.character(nodes[, ncol(nodes)])
}

##### CORUM COMPLEX INFORMATION ####  

corum_complexes <- read.table("coreComplexes.txt", header = T, sep = "\t", quote = "", 
                              fill = T, stringsAsFactors = F)
corum_complexes <- corum_complexes[, c("ComplexName", "subunits.UniProt.IDs.")]

list_complexes <- strsplit(corum_complexes[, "subunits.UniProt.IDs."], ";")
names(list_complexes) <- corum_complexes$ComplexName

# only keep complexes with at least 3 members (no duplexes)
list_complexes <- Filter(function(x) length(x) > 2, list_complexes)

# only keep complexes with at least 40% of its members in the provided dataset
keep_function <- function(x) {
  if (sum(x %in% nodes$accession) / length(x) >= 0.4 ) { 
    return (TRUE) 
  } else { return (FALSE) }
}
list_complexes <- Filter(keep_function, list_complexes)

chr2string <- function(x) {
  gene_names <- nodes[nodes$accession %in% x, "gene_name"]
  return (paste0(gene_names, collapse="; "))
}

percent_detected <- function(x) {
  return ( sum(x %in% nodes$accession) / length(x) )
}

# output complexes table (for download)
detected_complexes <- cbind.data.frame(names(list_complexes), 
                                          sapply(list_complexes, percent_detected),
                                          sapply(list_complexes, chr2string))
colnames(detected_complexes) <- c("Complex Name", "Fraction of Complex Detected",
                                  "Detected Complex Members")
detected_complexes <- detected_complexes[order(-detected_complexes$`Fraction of Complex Detected`), ]

# annotate nodes with complex membership
complex_membership <- rep("", NUM_NODES)
for (i in 1:NUM_NODES) {
  complexes <- which(grepl(nodes$gene_name[i], detected_complexes$`Detected Complex Members`))
  if (length(complexes) > 0) {
    complex_membership[i] <- paste0(names(list_complexes)[complexes], collapse="; ")
  }
}
nodes$complexes <- complex_membership


##### CLUSTER PROTEINS BASED ON ABUNDANCES OVER TIME ####

# if abundances are provided
if (ABUND_PROV) {
  links$cluster = NA
  relative_abundances <- {}
  
  for (bait_id in BAIT_IDS) {
    bait_gene_name <- nodes[which(nodes$id == bait_id), "gene_name"]
    
    # get abundances for all nodes associated with this bait (except bait itself)
    abundances <- links[which(links$bait_id == bait_id & links$prey_id != bait_id), 
                        c("prey_gene_name", paste("a", timepoints, sep=""))]
    
    # scale values from 0 to 1
    zero2one <- function(x) {
      return ((x - min(x)) / (max(x) - min(x)) )
    }
    
    # scale values to max 1
    maxOne <- function(x) {
      return ( x / max(x) )
    }
    
    to_scale <- abundances[, c(2:ncol(abundances))]
    abundances[, c(2:ncol(abundances))] <- t(apply(to_scale, 1, maxOne))
    
    # calculate pca
    abundances[is.na(abundances)] <- 0
    # only calculate pca on columns with differential expresson
    comp_pca <- which(dim(table(abundances[, c(2:ncol(abundances))])) > 1) + 1
    abundances_pca <- prcomp(scale(abundances[, comp_pca]))
    # find number of clusters as number of PCs that explain over 90% of the variance
    eigs <- abundances_pca$sdev^2 
    pve <- eigs / sum(eigs) # get percent of variance explained by each PC
    total_pve <- 0
    for (i in 1:length(pve)) {
      total_pve <- total_pve + pve[i]
      if (total_pve >= 0.9) { break }
    }
    num_clusters <- i
    
    # k-means clustering
    
    clusters <- kmeans(abundances[, c(2:ncol(abundances))], num_clusters)
    abundances <- cbind(abundances, clusters$cluster)
    colnames(abundances)[ncol(abundances)] <- "cluster"
    
    # add cluster numbers to links 
    for (i in 1:dim(abundances)[1]) {
      links$cluster[which(links$prey_gene_name == abundances$prey_gene_name[i] & 
                            links$bait_gene_name == bait_gene_name)] <- 
        abundances$cluster[i]
      }
      
    # save relative abundances to plot
    abundances$bait_gene_name <- bait_gene_name
    relative_abundances <- rbind(relative_abundances, abundances)
  }
    
  clusterPlot <- function(bait_id) {
    # plot profiles of proteins split by cluster
    melt_abund <- {}
    bait_gene_name <- nodes$gene_name[which(nodes$id == bait_id)]
    abundances <- relative_abundances[which(relative_abundances$bait_gene_name == 
                                              bait_gene_name), ]
    for (i in rownames(abundances)) {
      new <- cbind.data.frame(timepoints, rep(abundances[i, "cluster"], NUM_TPS))
      new <- cbind.data.frame(new, rep(abundances[i, "prey_gene_name"], NUM_TPS))
      new <- cbind.data.frame(new, as.vector(t(abundances[i, paste("a", timepoints, sep="")])))
      melt_abund <- rbind(melt_abund, new)
    }
    colnames(melt_abund) <- c("time", "cluster", "prey_gene_name", "abundance")
    rownames(melt_abund) <- 1:dim(melt_abund)[1]
    
    print(ggplot(data=melt_abund, aes(x=time, y=abundance, group=prey_gene_name)) +
            geom_line(aes(color=as.factor(cluster))) +
            geom_point(aes(color=as.factor(cluster))) +
            scale_x_continuous(breaks = timepoints) +
            labs(title=paste("Bait", bait_gene_name), x="Time", y="Scaled Relative Abundance", 
                 color="Cluster") +
            facet_wrap(~as.factor(cluster)))
  }
}

##### NORMALIZE BY PROTEOME ABUNDANCE ####

proteome_file <- list.files(path=".", pattern= "proteome_abundance.*")[1]
NORM_PROTEOME <- F
if (!is.na(proteome_file)) { NORM_PROTEOME <- T }

if (ABUND_PROV & NORM_PROTEOME) {
  library(gridExtra)
  
  ## NORMALIZE BY PROTEOME ABUNDANCE
  
  proteome_abundance <- read.table(proteome_file, header=T, sep="\t", 
                                   stringsAsFactors = F, quote="", comment.char = "")
  colnames(proteome_abundance) <- c("accession", "gene_name", paste("prot_a", timepoints, sep=""))
  
  proteome_abundance <- proteome_abundance[proteome_abundance$accession %in% 
                                             nodes$accession, ]
  heatmap_abundance <- links[, c("bait_gene_name", "prey_gene_name", "cluster",
                                 paste("a", timepoints, sep=""))]
  for (tp in timepoints) {
    norm_abund <- vector(mode='numeric', length=NUM_LINKS)
    for (i in 1:NUM_LINKS) {
      prey <- links$prey_accession[i]
      if (prey %in% proteome_abundance$accession) {
        norm_abund[i] <- links[i, paste("a", tp, sep="")] / 
          proteome_abundance[proteome_abundance$accession == prey, paste("prot_a", tp, sep="")]
      } else {
        norm_abund[i] <- NA
      }
    }
    heatmap_abundance <- cbind.data.frame(heatmap_abundance, norm_abund, stringsAsFactors = F)
    colnames(heatmap_abundance)[ncol(heatmap_abundance)] <- paste("norm_a", tp, sep="")
    heatmap_abundance[, paste("norm_a", tp, sep="")] <- as.numeric(heatmap_abundance[, paste("norm_a", tp, sep="")])
  }
  
  scaled_reg <- t(apply(heatmap_abundance[, c(paste("a", timepoints, sep=""))], 1, 
                        zero2one))
  scaled_norm <- t(apply(heatmap_abundance[, c(paste("norm_a", timepoints, sep=""))], 1, 
                         zero2one))
  
  heatmap_abundance_scaled <- cbind.data.frame(heatmap_abundance[, c("bait_gene_name",
                                                                       "prey_gene_name",
                                                                       "cluster")],
                                                scaled_reg, scaled_norm, stringsAsFactors = F)
  
  # append data to edge attributes file
  links <- cbind.data.frame(links, heatmap_abundance[, c(paste("norm_a", timepoints, 
                                                                sep=""))])
  
  plotHeatmaps <- function(bait_id) {
    bait_gene_name <- nodes$gene_name[which(nodes$id == bait_id)]
    bait_data <- heatmap_abundance_scaled[heatmap_abundance_scaled$bait_gene_name == bait_gene_name, ]
    
    clusters <- sort(unique(bait_data$cluster))
    grobList <- list()
    
    for (cluster in clusters) {
      cluster_data <- bait_data[bait_data$cluster == cluster, c("prey_gene_name", 
                                                                paste("a", timepoints, 
                                                                      sep=""),
                                                                paste("norm_a", 
                                                                      timepoints, sep=""))]
      prey_gene_name <- rep(cluster_data$prey_gene_name, each=NUM_TPS*2)
      type <- rep(c("Not Norm.", "Prot. Norm."), each=NUM_TPS, times=dim(cluster_data)[1])
      time <- rep(timepoints, times=dim(cluster_data)[1]*2)
      values <- c(t(cluster_data[, -1]))
      cluster_data_melt <- cbind.data.frame(prey_gene_name, type, time, values)
      
      grobList[[cluster]] <- ggplot(cluster_data_melt, aes(time, prey_gene_name)) + 
        geom_tile(aes(fill = values), colour = "white") + 
        scale_fill_gradient(low = "navyblue", high = "yellow3", 
                            name = "Scaled Abund.") + 
        facet_wrap(aes(type)) +
        xlab("Condition") + 
        ylab("Prey Gene Name") +
        scale_x_continuous(breaks = timepoints, labels = timepoints) +
        ggtitle(paste("Cluster", cluster))
    }
    do.call("grid.arrange", c(grobList, ncol=floor(sqrt(length(clusters))), 
                              top = paste("Bait", bait_gene_name)))
  }
}

##### SIZE PROTEINS BY ABUNDANCE FOR SINGLE-BAIT ####

if (length(BAIT_IDS) == 1 & ABUND_PROV) {
  # convert abundance to size and size by quartile
  
  overall_avg <- vector(mode='numeric', length=NUM_TPS)
  for (i in 1:NUM_TPS) {
    overall_avg[i] <- mean(links[, paste('a', timepoints[i], sep="")])
  }
  
  # set 2.5 as the average node size
  mult_factor <- 2.5 / mean(overall_avg)
  
  returnSize <- function(x) {
    y <- max(1, mult_factor*x) # prevent nodes from being too small
    y <- min(4, y) # prevent nodes from being too big
    return (y)
  }
  
  for (i in 1:NUM_TPS) {
    size <- lapply(links[, paste('a', timepoints[i], sep="")], returnSize)
    size <- cbind(links$prey_id, unlist(size))
    size <- rbind(size, c(BAIT_IDS, 2.5))
    colnames(size) <- c("id", paste("size_", timepoints[i], sep=""))
    nodes <- merge(nodes, size, by="id")
  }
}

##### GET NODE ANNOTATIONS FROM UNIPROT ####

# get UniProt data for proteins in each species and background list for enrichment

background_list <- read.table("background_gene_list.txt", header=T,
                              sep = "\t", stringsAsFactors = F)
colnames(background_list) <- c("accession", "taxid")

node_annotations <- {}
for (tx in TAXIDS) {
  tx_nodes <- nodes[which(nodes$taxid==tx), "accession"]
  tx_nodes <- c(tx_nodes, background_list[which(background_list$taxid == tx), 
                                          "accession"])
  # get data from UniProt
  up <- UniProt.ws(taxId=tx)
  up_selected <- select(up, keys=tx_nodes, columns = c("ORGANISM", "GO-ID", 
                                                       "SUBCELLULAR-LOCATIONS", 
                                                       keytype = "UNIPROTKB"))
  node_annotations <- rbind(node_annotations, up_selected)
}
colnames(node_annotations) <- c("accession", "organism", "GOid", "location")

# merge organism name into nodes file
nodes <- merge(nodes, node_annotations[, c("accession", "organism")], by="accession")
nodes <- nodes[!duplicated(nodes), ]
# create mapping from taxid to organism
txidToOrg <- list()
for (tx in TAXIDS) {
  txidToOrg[[toString(tx)]] <- na.omit(nodes[which(nodes$taxid == tx), "organism"])[1]
}
# if organism name was NA, replace with known 
unspecified <- which(is.na(nodes$organism))
for (u in unspecified) {
  nodes[u, "organism"] <- txidToOrg[toString(nodes[u, "taxid"])]
}
ORGANISMS <- unique(nodes$organism)

# create functional GO annotations list to later pipe into topGO
accessionToGO <- strsplit(node_annotations$GOid, "; ")
names(accessionToGO) <- node_annotations$accession

##### ASSIGN LOCALIZATIONS #####

# localizations already provided
if (LOC_PROV) { 
  # parse all given localizations along separator ";"
  locs <- unique(unlist(strsplit(nodes$localization, "; ")))
  NUM_LOCS <- length(locs)
  
  for (loc in locs) {
    check <- vector(mode="logical", length=NUM_NODES)
    check <- grepl(loc, nodes$localization, ignore.case=T)
    nodes <- cbind(nodes, check)
    colnames(nodes)[ncol(nodes)] <- paste("loc_", loc, sep="")
  }
} else { 
  locs <- c("Cell Membrane", "Cytoplasm", "Endoplasmic Reticulum", "Golgi", 
            "Mitochondrion", "Nucleus", "Peroxisome")
  NUM_LOCS <- length(locs)
  all_locs <- node_annotations$accession
  # assign locs based on whether 'subcellular location' comment in UniProt mentions location
  for (loc in locs) {
    x <- as.numeric(grepl(loc, node_annotations$location, ignore.case=T))
    all_locs <- cbind.data.frame(all_locs, x)
    colnames(all_locs)[ncol(all_locs)] <- paste("loc_", loc, sep="")
  }
  colnames(all_locs)[1] <- "accession"
  nodes <- merge(nodes, all_locs, by="accession")
  nodes <- nodes[!duplicated(nodes), ]
}

NUM_NODES <- dim(nodes)[1]

# convert localization information into single string
localization_string <- vector(mode='character', length=NUM_NODES)
localization_string_output <- vector(mode='character', length=NUM_NODES)

for (loc in locs) {
  index <- which(nodes[, paste('loc_', loc, sep="")] == 1)
  if (length(index) > 0) {
    localization_string_output[index] <- 
      paste(localization_string_output[index], loc, sep="; ")
    for (j in index) {
      if (localization_string[j] == "" ) { localization_string[j] <- loc; }
      # multiple possiblities
      else { localization_string[j] <- "Multiple" }
    }
  }
}
nodes$localization_string <- localization_string
nodes$localization_string_output <- localization_string_output

# if node lacks localization information replace with organism name
unspecified <- which(nodes$localization_string == "")
for (u in unspecified) {
  nodes[u, "localization_string"] <- txidToOrg[toString(nodes[u, "taxid"])]
}

##### ASSIGN FUNCTIONAL ANNOTATIONS AFTER GO ENRICHMENT #####

gene_universe <- names(accessionToGO)
interesting_genes <- nodes$accession

# assign whether a gene is "interesting" - in the provided data
gene_list <- factor(as.integer(gene_universe %in% interesting_genes))
names(gene_list) <- gene_universe

# create topGO object
# Biological Process (BP) ontology
GOdata <- new("topGOdata", ontology = "BP", allGenes = gene_list, 
              annot = annFUN.gene2GO, gene2GO = accessionToGO)

# use Fisher's test to determine enrichment
test <- new("classicCount", testStatistic = GOFisherTest, name = "Fisher test")
resultFisher <- getSigGroups(GOdata, test)

# perform multiple test correction using fdr (BH) method
fisher_p_values <- score(resultFisher)
corr_p_values <- p.adjust(fisher_p_values, method = "BH")
corrected_p_values <- as.data.frame(corr_p_values)
corrected_p_values$GO.ID <- rownames(corrected_p_values)
  
# generate results table
allRes <- GenTable(GOdata, Fisher_p_value = resultFisher,
                   topNodes = dim(corrected_p_values)[1])
allRes <- merge(allRes, corrected_p_values, by="GO.ID")
# only keep terms with p-value < 0.05 (if not more than 5 from multiple correction, then just Fisher)
if (length(which(allRes$corr_p_values < 0.05)) > 5) {
  allRes <- allRes[which(allRes$corr_p_values < 0.05), ]
  # sort by corrected p_values
  allRes <- allRes[order(allRes$corr_p_values), ]
} else {
  allRes <- allRes[which(allRes$Fisher_p_value < 0.05), ]
  # sort by Fisher p_values
  allRes <- allRes[order(allRes$Fisher_p_value), ]
}

colnames(allRes)[ncol(allRes)] <- "Corrected p-value"

# get list of top (upto 10) GO terms
num_go_terms <- min(dim(allRes)[1], 10)
TOP_GO_TERMS <- allRes$Term[1:num_go_terms]

# list the genes for each of the significant GO terms
annotated_genes <- lapply(allRes$GO.ID, 
                          function(x) as.character(
                            unlist(genesInTerm(object = GOdata, whichGO = x)))) 
names(annotated_genes) <- allRes$Term
# get GO terms for the original genes in the provided list
significant_genes <- lapply(annotated_genes, function(x) intersect(x, interesting_genes))

# create columns for each GO term designating whether a protein is associated with the term
# and string with all of top 10 significant GO terms
nodes$GO_term_string <- vector(mode = "character", length = NUM_NODES)
for (got in TOP_GO_TERMS) {
  genes_in_got <- rep(0, NUM_NODES)
  i <- which(nodes$accession %in% unlist(significant_genes[got]))
  nodes$GO_term_string[i] <- paste(nodes$GO_term_string[i], got, sep="; ")
  genes_in_got[i] <- 1
  nodes <- cbind(nodes, genes_in_got)
  colnames(nodes)[ncol(nodes)] <- paste("got_", got, sep="")
}
  
# output GO Table (for download)
GOTable <- allRes
NUM_GO_TERMS <- dim(GOTable)[1]
GOTable$Proteins <- vector(mode = "character", length = NUM_GO_TERMS)
for (i in 1:NUM_GO_TERMS) {
  GOTable$Proteins[i] <- paste0(nodes$gene_name[which(nodes$accession %in% 
                                                        unlist(significant_genes[GOTable$Term[i]]))], 
                                collapse = "; ")
}

##### EXTRACT LINK INFORMATION ####

links$type <- "Provided"

# calculate link onset and terminus information based on threshold
# onset when weight >= threshold, terminus when weight <= threshold
# duplicate link if there are non-overlapping spells of activity (ex: tp1-2, tp4-5)

NUM_LINKS <- dim(links)[1]

f <- 0
freq <- vector(mode='numeric', length=NUM_LINKS) 
# count number of times edge must repeat: non-overlapping spells
onset <- {}
terminus <- {} 

temp_timepoints <- c(timepoints, timepoints[NUM_TPS]+INTERVAL)
for (n in (1:NUM_LINKS)) {
  j <- 1
  freq[n] <- 0
  while (j <= NUM_TPS) {
    count <- 0
    time <- temp_timepoints[j]
    w <- links[n, paste('w', time, sep="")]
    if (w >= WEIGHT_THRESHOLD) {
      f <- f+1
      freq[n] <- freq[n] + 1
      onset[f] <- time
      while (w >= WEIGHT_THRESHOLD) {
        count <- count+1
        j <- j+1
        time <- temp_timepoints[j]
        if (time %in% timepoints) {
          w <- links[n, paste('w', time, sep="")] 
        } else { w <- -1 } # set a negative weight to break out of while loop
      }
      terminus[f] <- time
    }
    j <- j+1
    time <- temp_timepoints[j]
  }
}

links <- links[rep(rownames(links), freq), ]
links$onset <- onset
links$terminus <- terminus
rownames(links) <- 1:nrow(links)

##### ASSIGN NODE SPELLS AND SHARED NODES #####

# assign unique node onset and terminus times from edge info such that node activity encompasses 
# incident edge activity - examples:
# 1. edge tps=(1-2, 2-5) -> node tps=(1-5)
# 2. edge tps=(1-5, 2-4) -> node tps=(1-5)
# repeat until as condense as possible, then remove any remaining duplicates (on-off periods of edges) by
# extending node activity for entire edge spell

nodes <- nodes[order(nodes$id), ]
nodes <- nodes[!duplicated(nodes$id), ]
rownames(nodes) <- nodes$id
links <- links[order(links$prey_id), ]

# get the number of times to repeat each node (if it appears in multiple baits)
freq <- as.data.frame(table(links$prey_id), stringsAsFactors = F)
keep_nodes <- nodes[rep(freq$Var1, freq$Freq), ]
keep_nodes$onset <- links$onset
keep_nodes$terminus <- links$terminus

current_dim <- dim(keep_nodes)[1]
row.names(keep_nodes) <- 1:current_dim
new_dim <- 0

count <- 1
while (new_dim < current_dim) {
  count <- count+1
  remove_ids <- numeric()
  i <- 1
  while (i < dim(keep_nodes)[1]) {
    currentNode <- keep_nodes[i,"id"]
    while (i < dim(keep_nodes)[1] & keep_nodes[i+1, "id"] == currentNode) {
      # case 1
      if (keep_nodes[i, "onset"] == keep_nodes[i+1, "terminus"] | 
          keep_nodes[i, "terminus"] == keep_nodes[i+1, "onset"]) {
        keep_nodes[i, "onset"] <- min(keep_nodes[i, "onset"], keep_nodes[i+1, "onset"])
        keep_nodes[i, "terminus"] <- max(keep_nodes[i, "terminus"], 
                                         keep_nodes[i+1, "terminus"])
        remove_ids <- append(remove_ids, i+1)
      }
      # case 2
      else if ((keep_nodes[i, "onset"] <= keep_nodes[i+1, "onset"] & 
                keep_nodes[i, "terminus"] >= keep_nodes[i+1, "terminus"]) |
               (keep_nodes[i, "onset"] >= keep_nodes[i+1, "onset"] & 
                keep_nodes[i, "terminus"] <= keep_nodes[i+1, "terminus"])) {
        keep_nodes[i, "onset"] <- min(keep_nodes[i, "onset"], keep_nodes[i+1, "onset"])
        keep_nodes[i, "terminus"] <- max(keep_nodes[i, "terminus"], 
                                         keep_nodes[i+1, "terminus"])
        remove_ids <- append(remove_ids, i+1)
      }
      i <- i+1
    }
    i <- i+1
  }
  current_dim <- dim(keep_nodes)[1]
  if (length(remove_ids) == 0) { break }
  keep_nodes <- keep_nodes[-remove_ids, ]
  new_dim <- dim(keep_nodes)[1]
  rownames(keep_nodes) <- 1:new_dim
}

# remove remaining duplicates by extending them for max possible time when edges active
remove_ids <- numeric()
duplicates <- names(which(table(keep_nodes$id) > 1))
if (length(duplicates) > 0) {
  for (i in 1:length(duplicates)) {
    check <- which(keep_nodes$id == duplicates[i])
    onset <- min(keep_nodes[check, "onset"])
    terminus <- max(keep_nodes[check, "terminus"])
    keep_nodes[check[1], "onset"] <- onset
    keep_nodes[check[1], "terminus"] <- terminus
    remove_ids <- append(remove_ids, check[2:length(check)])
  }
  keep_nodes <- keep_nodes[-remove_ids, ]
}

# ensure that bait nodes are always present
missing_baits <- BAIT_IDS[!(BAIT_IDS %in% keep_nodes$id)]
if (length(missing_baits) > 0) {
  add_baits <- nodes[missing_baits, ]
  add_baits$onset <- timepoints[1]
  add_baits$terminus <- timepoints[NUM_TPS]
  
  keep_nodes <- rbind(keep_nodes, add_baits)
}

nodes <- keep_nodes 
nodes$onset[which(nodes$id %in% BAIT_IDS)] <- timepoints[1]
nodes$terminus[which(nodes$id %in% BAIT_IDS)] <- timepoints[NUM_TPS] + INTERVAL

nodes$duration <- nodes$terminus - nodes$onset
nodes <- nodes[order(nodes$id), ]

NUM_NODES <- dim(nodes)[1]

# convert old node ids to new ids
old_to_new <- setNames(as.list(1:NUM_NODES), nodes$id)
nodes$id <- 1:NUM_NODES

links$bait_id <- as.numeric(old_to_new[as.character(links$bait_id)])
links$prey_id <- as.numeric(old_to_new[as.character(links$prey_id)])

BAIT_IDS <- as.numeric(old_to_new[as.character(BAIT_IDS)])

##### GET STRING-DB INFORMATION ####

source_python("getStringInteractors.py")
source_python("gsi.py")

string_links <- {}

for (species in TAXIDS) {
  values <- nodes[which(nodes$taxid == species), "accession"]
  keys <- nodes[which(nodes$taxid == species), "gene_name"]
  node_list <- setNames(as.list(values), keys)
  interactions <- run_functions(dict(node_list), species)
  string_links <- rbind(string_links, interactions)
}

colnames(string_links) <- c("bait_accession", "prey_accession",
                            "bait_gene_name", "prey_gene_name", "weight")

# only keep interactions that have been experimentally verified 
string_links <- string_links[which(string_links$weight != 0), ]

NUM_STRING_EDGES <- dim(string_links)[1]

bait_id <- vector(mode = 'numeric', length = NUM_STRING_EDGES)
prey_id <- vector(mode = 'numeric', length = NUM_STRING_EDGES)
bait_gene_name <- vector(mode = 'character', length = NUM_STRING_EDGES)
prey_gene_name <- vector(mode = 'character', length = NUM_STRING_EDGES)

for (i in 1:NUM_STRING_EDGES) {
  bait_id[i] <- nodes[which(nodes$accession == string_links[i, "bait_accession"]), "id"]
  prey_id[i] <- nodes[which(nodes$accession == string_links[i, "prey_accession"]), "id"]
  bait_gene_name[i] <- nodes[which(nodes$accession == string_links[i, "bait_accession"]), 
                             "gene_name"][1]
  prey_gene_name[i] <- nodes[which(nodes$accession == string_links[i, "prey_accession"]), 
                             "gene_name"][1]
}

string_links$bait_id <- bait_id
string_links$prey_id <- prey_id
string_links$bait_gene_name <- bait_gene_name
string_links$prey_gene_name <- prey_gene_name

# give edge same weight and NA abundance at all time points
for (tp in timepoints) {
  weight_var <- string_links$weight
  string_links <- cbind(string_links, weight_var)
  colnames(string_links)[length(string_links)] <- paste('w', tp, sep="")
  abundance_var <- rep(NA, dim(string_links)[1])
  string_links <- cbind(string_links, abundance_var)
  colnames(string_links)[length(string_links)] <- paste('a', tp, sep="")
  norm_abundance_var <- rep(NA, dim(string_links)[1])
  string_links <- cbind(string_links, norm_abundance_var)
  colnames(string_links)[length(string_links)] <- paste('norm_a', tp, sep="")
}

# make onset and terminus infinity
string_links$onset <- Inf
string_links$terminus <- Inf

# assign type to be STRING inferred or both
both <- merge(string_links, links, by=c("bait_id", "prey_id"))
if (dim(both)[1] > 0) {
  both$type <- "Both"
  string_links <- merge(string_links, both, by=c("bait_id", "prey_id"), all=T)
  string_links$type[is.na(string_links$type)] <- "STRING"
} else {
  string_links$type <- "STRING"
}

# assign clusters to be NA 
string_links$cluster <- NA

# bind STRING links to provided links file
string_links <- string_links[, colnames(links)]
links <- rbind.data.frame(links, string_links)

NUM_LINKS <- dim(links)[1]

##### NETWORK ####

MAX_DISPLAY <- 500

if (NUM_NODES < MAX_DISPLAY) {
  # rearrange link dataframes to network
  if (ABUND_PROV) {
    links <- links[, c("bait_id", "prey_id", "bait_accession", "prey_accession", "bait_gene_name",
                       "prey_gene_name", "cluster", "onset", "terminus", "type", 
                       paste("w", timepoints, sep=""),
                       paste("a", timepoints, sep=""))]
  } else {
    links <- links[, c("bait_id", "prey_id", "bait_accession", "prey_accession", "bait_gene_name",
                       "prey_gene_name", "onset", "terminus", "type", 
                       paste("w", timepoints, sep=""))]
  }
  
  colnames(links)[c(1,2)] <- c("from", "to")
  
  # give localization, duration and edge type a factor numeric value (to later assign color)
  nodes$locanum <- as.numeric(factor(nodes$localization_string))
  # give duration factor based on what quartile of expression it is in
  nodes$durfactor <- cut(nodes$duration, unique(quantile(nodes$duration, probs=0:4/4)),
                         include.lowest=T, labels=F)
  links$edge_color <- c("black", "navyblue", "gray75")[as.numeric(factor(links$type))]
  
  # create network object
  
  net3 <- network(links, vertex.attr=nodes, matrix.type="edgelist", 
                  loops=F, multiple=F, directed=F, ignore.eval = F)
  
  # assign colors to localization, duration factors
  newPalette <- colorRampPalette(brewer.pal(8, "Set2"))
  loc_color <- newPalette(n = nlevels(factor(nodes$localization_string))+1)
  
  dur_color <- brewer.pal(n = nlevels(factor(nodes$durfactor)), name = "Reds")
  
  net3 %v% "loc_color" <- loc_color[net3 %v% "locanum"]
  net3 %v% "dur_color" <- dur_color[net3 %v% "durfactor"]
  
  # dynamic properties
  vs <- data.frame(onset=nodes$onset, terminus=nodes$terminus, vertex.id=1:NUM_NODES)
  es <- data.frame(onset=links$onset, terminus=links$terminus, tail=links$to, 
                   head=links$from)
  
  # create dynamic network object
  ## NET3.DYN HAS ALL THE ORIGINAL DATA, NEVER MODIFY THIS BUT ONLY MAKE COPIES TO MODIFY
  
  net3.dyn <- networkDynamic(base.net=net3, vertex.spells=vs, edge.spells=es)
  
  net3.dyn <- reconcile.edge.activity(net3.dyn, mode="reduce.to.vertices")
  
  # set dynamic edge weights
  for (i in 1:NUM_TPS) {
    temp_timepoints <- c(timepoints, timepoints[NUM_TPS]+INTERVAL)
    tp <- temp_timepoints[i]
    next_tp <- temp_timepoints[i+1]
    net3.dyn <- activate.edge.attribute(net3.dyn, 'edge_weight', 
                                        (links[, paste('w', tp, sep="")])*2, 
                                        onset=tp, terminus=next_tp, e=1:NUM_LINKS)
    # set dynamic vertex sizes (if single bait)
    if (length(BAIT_IDS) == 1 & ABUND_PROV) {
      net3.dyn <- activate.vertex.attribute(net3.dyn, 'node_size',
                                            nodes[, paste('size_', tp, sep="")],
                                            onset=tp, terminus=next_tp, v=1:NUM_NODES)
    }
  }
}

# output protein information
saveNodeAttributes <- function(filename) {
  output_nodes <- nodes[, c("gene_name", "accession", "organism", 
                            "localization_string_output", "complexes", 
                            "GO_term_string")]
  output_nodes$GO_term_string <- sub("^; ", "", output_nodes$GO_term_string)
  output_nodes$localization_string_output <- sub("^; ", "", 
                                                 output_nodes$localization_string_output)
  colnames(output_nodes) <- c("gene_name", "accession", "species", "localizations",
                              "complexes", "GO_terms")
  if (NUM_BAITS > 1) {
    output_nodes <- cbind(output_nodes, nodes[, paste("shared_baits_", timepoints, sep="")])
  }
  write.table(output_nodes, file = filename, sep = "\t", 
              quote = F, row.names = F)
}

# output edge information
saveEdgeAttributes <- function(filename) {
  # get edge information to output
  output_edges <- links[, c("bait_gene_name", "prey_gene_name",
                            "bait_accession", "prey_accession",
                            "onset", "terminus", "type", 
                            paste("w", timepoints, sep=""))]
  output_edges_colnames <- c("#node1", "node2", "node1_accession", "node2_accession", 
                             "onset", "terminus", "type", 
                             paste("confidence_score_", timepoints, sep=""))
  
  if (ABUND_PROV) {
    output_edges <- cbind(output_edges, links[, c("cluster", 
                                                  paste("a", timepoints, sep=""))])
    output_edges_colnames <- c(output_edges_colnames, "cluster",
                               paste("abundance_", timepoints, sep=""))
    if (NORM_PROTEOME) {
      output_edges <- cbind(output_edges, links[, paste("norm_a", timepoints, sep="")])
      output_edges_colnames <- c(output_edges_colnames, 
                                      paste("normalized_to_proteome_abundance_", 
                                            timepoints, 
                                            sep=""))
    }
  }
  
  # rename abundance columns
  colnames(output_edges) <- output_edges_colnames
  
  write.table(output_edges, file = filename, sep = "\t", 
              quote = F, row.names = F)
}
  
##### DEACTIVATE NODES AND EDGES BASED ON INPUT ####

if (NUM_NODES < MAX_DISPLAY) {
  sharedCalc <- function(min_shared) {
    plot_net <- net3.dyn
    # deactivate nodes based on number of baits they share
    for (tp in timepoints) {
      keep_nodes <- nodes[nodes[, paste("shared_", tp, sep="")] >= min_shared, "id"]
      keep_nodes <- c(keep_nodes, BAIT_IDS)
      keep_nodes <- unique(keep_nodes)
      deactivate_nodes <- nodes[!(nodes$id %in% keep_nodes), "id"]
      deactivate.vertices(plot_net, onset=tp, terminus=tp+INTERVAL, v=deactivate_nodes, 
                          deactivate.edges=T)
    }
    return(plot_net)
  }
  
  deactivateNE <- function(plot_net, st, locns, species, go_terms) {
    
    # deactivate edges (based on STRING confidence threshold)
    deactivate_link_ids <- which((links$type == "STRING") &
                                   (links[, paste("w", timepoints[1], sep="")] <= st))
    
    if (!is.null(deactivate_link_ids)) {
      deactivate.edges(plot_net, e=deactivate_link_ids)
    }
    
    # deactivate nodes (based on localizations, species, functions, etc.)
    
    # localizations
    keep_loc_vertex_ids <- {}
    for (loc in locns) {
      keep_loc_vertex_ids <- append(keep_loc_vertex_ids, 
                                    nodes[which(nodes[, paste("loc_", loc, sep="")] == 1), "id"])
    }
    
    # species 
    keep_spec_vertex_ids <- {}
    for (spec in species) {
      keep_spec_vertex_ids <- append(keep_spec_vertex_ids, 
                                     nodes[which(nodes$organism == spec), "id"])
    }
    
    # GO_terms
    keep_got_vertex_ids <- {}
    if ("All" %in% go_terms) {
      keep_got_vertex_ids <- nodes$id
    } else {
      for (got in go_terms) {
        keep_got_vertex_ids <- append(keep_got_vertex_ids, 
                                      nodes[which(nodes[, paste("got_", got, sep="")] == 1), "id"])
      }
    }
    
    keep_vertex_ids <- Reduce(intersect, list(keep_loc_vertex_ids, keep_spec_vertex_ids,
                                              keep_got_vertex_ids))
    deactivate_vertex_ids <- setdiff(nodes$id, keep_vertex_ids)
    
    if (!is.null(deactivate_vertex_ids)) {
      # don't deactivate baits
      deactivate_vertex_ids <- deactivate_vertex_ids[!deactivate_vertex_ids %in% BAIT_IDS]
      
      deactivate.vertices(plot_net, onset=timepoints[1], terminus=(timepoints[NUM_TPS]+INTERVAL),
                          v=deactivate_vertex_ids, deactivate.edges=T)
    }
    
    # compute network animation
    compute.animation(plot_net, animation.mode = "kamadakawai", default.dist=2, 
                      slice.par=list(start=timepoints[1], end=timepoints[NUM_TPS], 
                                     interval=INTERVAL, 
                                     aggregate.dur=INTERVAL, rule='any'))
    return (plot_net)
  }
}

##### PLOT DYNAMIC NETWORK ####

if (NUM_NODES < MAX_DISPLAY) {
  if (length(BAIT_IDS) == 1 & ABUND_PROV) {
    plotCompNet <- function(net3_network, lab) {
      # plot network animation
      render.d3movie(net3_network,
                     render.par=list(tween.frames = 10, show.time = F),
                     plot.par=list(bg='white'),
                     d3.options = list(animationDuration=6000,enterExitAnimationFactor=0.3),
                     usearrows = F,
                     displaylabels = lab, label=function(slice){slice%v%'gene_name'},
                     vertex.border=function(slice){slice%v%'dur_color'},
                     vertex.cex = function(slice){slice%v%'node_size'},
                     vertex.lwd = 3,
                     vertex.col = function(slice){slice%v%'loc_color'},
                     edge.lwd = 'edge_weight', 
                     edge.col = function(slice){slice%e%'edge_color'},
                     vertex.tooltip = function(slice){paste("<b>Gene Name:</b>", slice%v%'gene_name', 
                                                            "<br>", "<b>Localization:</b>", 
                                                            slice%v%'localization_string',
                                                            "<br>", "<b>Complexes:</b>", 
                                                            slice%v%'complexes',
                                                            "<br>", "<b>Time Active:</b>", 
                                                            slice%v%'duration')},
                     edge.tooltip = function(slice){paste('From:',slice%e%'bait_gene_name', ' To: ', 
                                                          slice%e%'prey_gene_name', '<br>',
                                                          'Weight:', (slice%e%'edge_weight')/2)},
                     output.mode = 'htmlWidget')
    } 
    
  } else {
    plotCompNet <- function(net3_network, lab) {
      # plot network animation
      render.d3movie(net3_network,
                     render.par=list(tween.frames = 10, show.time = F),
                     plot.par=list(bg='white'),
                     d3.options = list(animationDuration=6000,enterExitAnimationFactor=0.3),
                     usearrows = F,
                     displaylabels = lab, label=function(slice){slice%v%'gene_name'},
                     vertex.border=function(slice){slice%v%'dur_color'},
                     vertex.cex = 1.5,
                     vertex.lwd = 3,
                     vertex.col = function(slice){slice%v%'loc_color'},
                     edge.lwd = 'edge_weight', 
                     edge.col = function(slice){slice%e%'edge_color'},
                     vertex.tooltip = function(slice){paste("<b>Gene Name:</b>", slice%v%'gene_name', 
                                                            "<br>", "<b>Localization:</b>", 
                                                            slice%v%'localization_string',
                                                            "<br>", "<b>Complexes:</b>", 
                                                            slice%v%'complexes',
                                                            "<br>", "<b>Time Active:</b>", 
                                                            slice%v%'duration')},
                     edge.tooltip = function(slice){paste('From:',slice%e%'bait_gene_name', ' To: ', 
                                                          slice%e%'prey_gene_name', '<br>',
                                                          'Weight:', (slice%e%'edge_weight')/2)},
                     output.mode = 'htmlWidget')
    } 
  }
}

##### ANALYTICAL TOOLS REPORT ####

analyticalTools <- function() {
  # timeline of node activity separated by localization
  for (loc in locs) {
    vs <- nodes[which(nodes[, paste("loc_", loc, sep="")] == 1), "id"]
    print(timeline(net3.dyn, v=vs, plot.edge.spells=F, displaylabels = F, 
                   xlab = "Time", ylab = paste(loc, "Proteins")))
  }
  
  # timeline of node activity separated by GO terms
  for (got in TOP_GO_TERMS) {
    vs <- nodes[which(nodes[, paste("got_", got, sep="")] == 1), "id"]
    print(timeline(net3.dyn, v=vs, plot.edge.spells=F, displaylabels = F, 
                   xlab = "Time", ylab = paste(got, "proteins")))
  }
  
  # edge formation
  data <- as.data.frame(tEdgeFormation(net3.dyn))
  
  print(ggplot(data) +
          geom_line(mapping = aes(x = 1:dim(data)[1], y = as.numeric(data$x)), size = 1.5) +
          scale_x_continuous(breaks = timepoints) +
          labs(x = "Time", y = "Number of Edges Formed", title = "Time of Edge Formation"))
  
  # edge duration
  provided_edges <- which(links$type == "Provided")
  data <- data.frame(x=edgeDuration(net3.dyn, e=provided_edges))
  print(ggplot(data) +
    geom_histogram(mapping = aes(x = data$x)) +
    scale_x_continuous(breaks = timepoints) +
    labs(x = "Duration", y = "Number of Edges", title = "Duration of Edge Activity"))
}

##### ABUNDANCE PLOTS #####

if (ABUND_PROV) {
  # re-orient data frame to make it suitable for ggplot2 plotting
  abund_info <- links[which(links$type == "Provided"), ]
  abund_info <- abund_info[, c("bait_gene_name", "prey_gene_name", 
                               paste("a", timepoints, sep=""))]
  
  # melt abundances dataframe down
  abundances <- {}
  for (i in 1:dim(abund_info)[1]) {
    new <- cbind.data.frame(timepoints, rep(abund_info[i, "bait_gene_name"], NUM_TPS))
    new <- cbind.data.frame(new, rep(abund_info[i, "prey_gene_name"], NUM_TPS))
    new <- cbind.data.frame(new, as.vector(t(abund_info[i, paste("a", timepoints, sep="")])))
    abundances <- rbind(abundances, new)
  }
  colnames(abundances) <- c("time", "bait_gene_name", "prey_gene_name", "abundance")
  rownames(abundances) <- 1:dim(abundances)[1]
  
  # create list of neighbors for every node
  calcNeighbors <- function(confidence, locns, species, go_terms) {
  
    # localizations
    keep_loc_nodes <- {}
    if (all(locns == locs)) { 
      keep_loc_nodes <- nodes$accession 
      } else {
      for (loc in locns) {
        keep_loc_nodes <- append(keep_loc_nodes, 
                                 nodes[which(nodes[, paste("loc_", loc, sep="")] == 1), 
                                       "accession"])
      }
    }
    
    # species 
    keep_spec_nodes <- {}
    for (spec in species) {
      keep_spec_nodes <- append(keep_spec_nodes, 
                                nodes[which(nodes$organism == spec), "accession"])
    }
    
    # GO_terms
    keep_got_nodes <- {}
    if ("All" %in% go_terms) {
      keep_got_nodes <- nodes$accession
    } else {
      for (got in go_terms) {
        keep_got_nodes <- append(keep_got_nodes, 
                                 nodes[which(nodes[, paste("got_", got, sep="")] == 1), 
                                       "accession"])
      }
    }
    
    keep_nodes <- Reduce(intersect, list(keep_loc_nodes, keep_spec_nodes,
                                         keep_got_nodes))
    
    keep_links <- links[which(links$bait_accession %in% keep_nodes &
                                links$prey_accession %in% keep_nodes &
                                links[, paste('w', timepoints[1], sep="")] >= confidence &
                                links$type == "STRING"), ]
    neighbors <- list()
    for (i in 1:NUM_NODES) {
      add <- keep_links[which(keep_links$bait_accession == nodes[i, "accession"]), 
                        "prey_accession"]
      add <- c(add, keep_links[which(keep_links$prey_accession == nodes[i, "accession"]), 
                               "bait_accession"])
      neighbors[[nodes[i, "accession"]]] <- add
    }
    return (neighbors)
  }
  
  # decimal places to show in label
  scalefunction <- function(x) sprintf("%.2f", x)
  
  abundPlot <- function(genesyms, nhbr=F, neighborList) {
    genes <- unlist(strsplit(genesyms, "; "))
    # display neighboring nodes abundances as well
    if (nhbr == T) {
      accessions <- {}
      for (gene in genes) {
        accessions <- c(accessions, 
                        nodes[grepl(gene, nodes$gene_name, ignore.case = T), "accession"])
      }
      neighbor_genes <- nodes[nodes$accession %in% unlist(neighborList[accessions]), 
                              "gene_name"]
      genes <- as.vector(c(genes, neighbor_genes))
    }
    current_genes <- {}
    for (gene in genes) {
      current_genes <- rbind(current_genes, 
                        abundances[grepl(gene, abundances$prey_gene_name, 
                                         ignore.case = T), ])
    }
    print(ggplot(data = current_genes) +
            geom_line(mapping = aes(x=time, y=abundance, color=prey_gene_name, 
                                    linetype=bait_gene_name), size=1.5) +
            scale_y_continuous(trans='log', labels=scalefunction) +
            scale_x_continuous(breaks = timepoints) + 
            labs(x="Time", y="Abundance", color="Search Genes", 
                 linetype="Bait"))
  }
  
  printClusterNumber <- function(genesyms) {
    genes <- unlist(strsplit(genesyms, ", "))
    for (gene in genes) {
      for (bait in BAIT_IDS) {
        bait_gene_name <- nodes$gene_name[which(nodes$id == bait)]
        clusterNumber <- links[which(links$from == bait & 
                                       grepl(gene, links$prey_gene_name, ignore.case = T)),
                               "cluster"]
        if (length(clusterNumber) != 0) {
          out <- cat(paste(gene, " is in cluster ", clusterNumber, " for bait ", 
                           bait_gene_name, "\n", collapse=""))
        }
      }
    }
  }
}

##### SHINY INTERACTIVE WEB APP OUTPUT UI ####

# output if there are too many proteins to display (> 500) but abundance provided (assumed)

if (NUM_NODES >= MAX_DISPLAY) {
  ui <- fluidPage(

    titlePanel("Interactome"),
    
    # layout - sidebar + main
    sidebarLayout(
      
      # what to include in sidebar
      sidebarPanel(
        # slider for STRING confidence
        sliderInput(inputId = "confidence", label="STRING edges confidence threshold:",
                    min=0, max=1, value=0.5),
        
        # radio buttons for shared proteins
        radioButtons("shared", 
                     label="Show proteins shared between at least 'n' baits:",
                     choices=1:NUM_BAITS, selected=1, inline=T),
        
        # checkboxes for which localizations to include
        checkboxGroupInput("checkLocalizations", label="Show proteins localized to:",
                           choices=locs, selected=locs),
        
        # checkboxes for which GO terms to include
        checkboxGroupInput("checkGOTerms", label="Show proteins associated with:",
                           choices=c(TOP_GO_TERMS, "All"), 
                           selected="All"),
        
        # checkboxes for which species to include
        checkboxGroupInput("checkSpecies", label="Show species:",
                           choices=ORGANISMS, selected=ORGANISMS),
        
        # break to separate quantitative plotting tools
        br(),
        br(),
        
        # textbox for genes to plot
        textInput("genesym", label = "See quantitative values for:", 
                  value = "Ex: IFT27; CCL19; TRAV35"),
        
        # checkbox to include node neighbors in plot
        checkboxInput(inputId = "neighbors", "Display interactor quants", value=F),
        
        # action button to calculate plot and switch tabs to plots
        actionButton("abundGo", "Render Quantitative Plot")
      ),
      
      # what to include in the main (output) panel
      mainPanel(
        
        # tabs for the various kinds of outputs (network, plots, downloads)
        tabsetPanel(id = "inTabset",
                    
                    tabPanel("Network", value = "network",
                             h4("Too many proteins ( > 500 ) for network display; all other tool functions available.")
                    ),
                    
                    tabPanel("Quantitative Plots", value = "plots",
                             h4("Plot protein and interactor quantitative values"),
                             plotOutput("abundPlot"),
                             downloadButton("dAbund", "Download Quantitative Plot"),
                             br(),
                             br(),
                             verbatimTextOutput("clusterNumber"),
                             br(),
                             h4("Prey Proteins per Bait Clustered by Temporal Quantitative Profiles"),
                             uiOutput("clusterPlots"),
                             downloadButton("dCluster", "Download Cluster Plots"),
                             br(),
                             br(),
                             h4("Heatmaps of Scaled Quantitative Profiles, Prior to and After Normalization to Proteome"),
                             uiOutput("heatmapPlots"),
                             downloadButton("dHeatmaps", "Download Heatmaps"),
                             br(),
                             br()),
                    
                    tabPanel("Gene Ontology", value = "GO",
                             br(),
                             h4("Significant gene ontology terms (Corrected p-value < 0.05)"),
                             h6("If this stringent corrected p-value cut-off found fewer than 5 GO terms, all terms with p-values < 0.05 are displayed"),
                             br(),
                             downloadButton("dGOTable", "Download Gene Ontology Table"),
                             br(),
                             br(),
                             dataTableOutput("GOTable")),
                    
                    tabPanel("Detected Protein Complexes", value = "complexes",
                             br(),
                             h4("Detected Protein Complexes (> 40% membership)"),
                             br(),
                             downloadButton("dComplexes", "Download Detected Complexes Table"),
                             br(),
                             br(),
                             dataTableOutput("complexesTable")
                    ),
                    
                    tabPanel("Downloads", value = "downloads", 
                             verbatimTextOutput("downloadNetworkInfo"),
                             downloadButton("downloadNodes", "Node Attributes"),
                             downloadButton("downloadEdges", "Edge Attributes"),
                             br())
        )
      )
    )
    )
} else {
  # output if quantitative plots must be displayed
  
  if (ABUND_PROV) {
    ui <- fluidPage(
      # make progress bar more visible
      tags$head(tags$style(HTML("
                                .progress-striped .bar {
                                background-color: #149bdf;
                                background-image: -webkit-gradient(linear, 0 100%, 100% 0, color-stop(0.25, rgba(255, 255, 255, 0.6)), color-stop(0.25, transparent), color-stop(0.5, transparent), color-stop(0.5, rgba(255, 255, 255, 0.6)), color-stop(0.75, rgba(255, 255, 255, 0.6)), color-stop(0.75, transparent), to(transparent));
                                background-image: -webkit-linear-gradient(45deg, rgba(255, 255, 255, 0.6) 25%, transparent 25%, transparent 50%, rgba(255, 255, 255, 0.6) 50%, rgba(255, 255, 255, 0.6) 75%, transparent 75%, transparent);
                                background-image: -moz-linear-gradient(45deg, rgba(255, 255, 255, 0.6) 25%, transparent 25%, transparent 50%, rgba(255, 255, 255, 0.6) 50%, rgba(255, 255, 255, 0.6) 75%, transparent 75%, transparent);
                                background-image: -o-linear-gradient(45deg, rgba(255, 255, 255, 0.6) 25%, transparent 25%, transparent 50%, rgba(255, 255, 255, 0.6) 50%, rgba(255, 255, 255, 0.6) 75%, transparent 75%, transparent);
                                background-image: linear-gradient(45deg, rgba(255, 255, 255, 0.6) 25%, transparent 25%, transparent 50%, rgba(255, 255, 255, 0.6) 50%, rgba(255, 255, 255, 0.6) 75%, transparent 75%, transparent);
                                -webkit-background-size: 40px 40px;
                                -moz-background-size: 40px 40px;
                                -o-background-size: 40px 40px;
                                background-size: 40px 40px;
                                }
                                "))),
      
      titlePanel("Interactome"),
      
      # layout - sidebar + main
      sidebarLayout(
        
        # what to include in sidebar
        sidebarPanel(
          # checkbox to display labels
          checkboxInput(inputId = "labels", "Display node labels", value=F),
          
          # break to space out network plotting tools
          br(),
          
          # slider for STRING confidence
          sliderInput(inputId = "confidence", label="STRING edges confidence threshold:",
                      min=0, max=1, value=0.5),
          
          # radio buttons for shared proteins
          radioButtons("shared", 
                       label="Show proteins shared between at least 'n' baits:",
                       choices=1:NUM_BAITS, selected=1, inline=T),
          
          # checkboxes for which localizations to include
          checkboxGroupInput("checkLocalizations", label="Show proteins localized to:",
                             choices=locs, selected=locs),
          
          # checkboxes for which GO terms to include
          checkboxGroupInput("checkGOTerms", label="Show proteins associated with:",
                             choices=c(TOP_GO_TERMS, "All"), 
                             selected="All"),
          
          # checkboxes for which species to include
          checkboxGroupInput("checkSpecies", label="Show species:",
                             choices=ORGANISMS, selected=ORGANISMS),
          
          # action button to calculate and render plot
          actionButton("netGo", "Render Network Plot"),
          
          # break to separate quantitative plotting tools
          br(),
          br(),
          
          # textbox for genes to plot
          textInput("genesym", label = "See quantitative values for:", 
                    value = "Ex: IFT27; CCL19; TRAV35"),
          
          # checkbox to include node neighbors in plot
          checkboxInput(inputId = "neighbors", "Display interactor quants", value=F),
          
          # action button to calculate plot and switch tabs to plots
          actionButton("abundGo", "Render Quantitative Plot")
        ),
        
        # what to include in the main (output) panel
        mainPanel(
          
          # tabs for the various kinds of outputs (network, plots, downloads)
          tabsetPanel(id = "inTabset",
                      
                      tabPanel("Network", value = "network",
                               # stop boostrap css from messing up the tooltip in the widget
                               tags$style(HTML(".tooltip {opacity: 1}")), 
                               ndtv:::ndtvAnimationWidgetOutput("netPlot", 
                                                                width="100%", height="800px")
                      ),
                      
                      tabPanel("Quantitative Plots", value = "plots",
                               h4("Plot protein and interactor quantitative values"),
                               plotOutput("abundPlot"),
                               downloadButton("dAbund", "Download Quantitative Plot"),
                               br(),
                               br(),
                               verbatimTextOutput("clusterNumber"),
                               br(),
                               h4("Prey Proteins per Bait Clustered by Temporal Quantitative Profiles"),
                               uiOutput("clusterPlots"),
                               downloadButton("dCluster", "Download Cluster Plots"),
                               br(),
                               br(),
                               h4("Heatmaps of Scaled Quantitative Profiles, Prior to and After Normalization to Proteome"),
                               uiOutput("heatmapPlots"),
                               downloadButton("dHeatmaps", "Download Heatmaps"),
                               br(),
                               br()),
                      
                      tabPanel("Gene Ontology", value = "GO",
                               br(),
                               h4("Significant gene ontology terms (Corrected p-value < 0.05)"),
                               h6("If this stringent corrected p-value cut-off found fewer than 5 GO terms, all terms with p-values < 0.05 are displayed"),
                               br(),
                               downloadButton("dGOTable", "Download Gene Ontology Table"),
                               br(),
                               br(),
                               dataTableOutput("GOTable")),
                      
                      tabPanel("Detected Protein Complexes", value = "complexes",
                               br(),
                               h4("Detected Protein Complexes (> 40% membership)"),
                               br(),
                               downloadButton("dComplexes", "Download Detected Complexes Table"),
                               br(),
                               br(),
                               dataTableOutput("complexesTable")
                      ),
                      
                      tabPanel("Downloads", value = "downloads", 
                               verbatimTextOutput("downloadNetworkInfo"),
                               downloadButton("downloadNodes", "Node Attributes"),
                               downloadButton("downloadEdges", "Edge Attributes"),
                               br(),
                               br(),
                               verbatimTextOutput("downloadReportInfo"),
                               downloadButton("downloadReport", "Report"))
          )
        )
      )
      )
} else {
  ui <- fluidPage(
    # make progress bar more visible
    tags$head(tags$style(HTML("
                              .progress-striped .bar {
                              background-color: #149bdf;
                              background-image: -webkit-gradient(linear, 0 100%, 100% 0, color-stop(0.25, rgba(255, 255, 255, 0.6)), color-stop(0.25, transparent), color-stop(0.5, transparent), color-stop(0.5, rgba(255, 255, 255, 0.6)), color-stop(0.75, rgba(255, 255, 255, 0.6)), color-stop(0.75, transparent), to(transparent));
                              background-image: -webkit-linear-gradient(45deg, rgba(255, 255, 255, 0.6) 25%, transparent 25%, transparent 50%, rgba(255, 255, 255, 0.6) 50%, rgba(255, 255, 255, 0.6) 75%, transparent 75%, transparent);
                              background-image: -moz-linear-gradient(45deg, rgba(255, 255, 255, 0.6) 25%, transparent 25%, transparent 50%, rgba(255, 255, 255, 0.6) 50%, rgba(255, 255, 255, 0.6) 75%, transparent 75%, transparent);
                              background-image: -o-linear-gradient(45deg, rgba(255, 255, 255, 0.6) 25%, transparent 25%, transparent 50%, rgba(255, 255, 255, 0.6) 50%, rgba(255, 255, 255, 0.6) 75%, transparent 75%, transparent);
                              background-image: linear-gradient(45deg, rgba(255, 255, 255, 0.6) 25%, transparent 25%, transparent 50%, rgba(255, 255, 255, 0.6) 50%, rgba(255, 255, 255, 0.6) 75%, transparent 75%, transparent);
                              -webkit-background-size: 40px 40px;
                              -moz-background-size: 40px 40px;
                              -o-background-size: 40px 40px;
                              background-size: 40px 40px;
                              }
                              "))),
    
    titlePanel("Interactome"),
    
    # layout - sidebar + main
    sidebarLayout(
      
      # what to include in sidebar
      sidebarPanel(
        # checkbox to display labels
        checkboxInput(inputId = "labels", "Display node labels", value=F),
        
        # break to space out network plotting tools
        br(),
        
        # slider for STRING confidence
        sliderInput(inputId = "confidence", label="STRING edges confidence threshold:",
                    min=0, max=1, value=0.5),
        
        # radio buttons for shared proteins
        radioButtons("shared", 
                     label="Show proteins shared between at least 'n' baits:",
                     choices=1:NUM_BAITS, selected=1, inline=T),
        
        # checkboxes for which localizations to include
        checkboxGroupInput("checkLocalizations", label="Show proteins localized to:",
                           choices=locs, selected=locs),
        
        # checkboxes for which GO terms to include
        checkboxGroupInput("checkGOTerms", label="Show proteins associated with:",
                           choices=c(TOP_GO_TERMS, "All"), 
                           selected="All"),
        
        # checkboxes for which species to include
        checkboxGroupInput("checkSpecies", label="Show species:",
                           choices=ORGANISMS, selected=ORGANISMS),
        
        # action button to calculate and render plot
        actionButton("netGo", "Render Network Plot")
      ),
      
      # what to include in the main (output) panel
      mainPanel(
        
        # tabs for the various kinds of outputs (network, plots, downloads)
        tabsetPanel(id = "inTabset",
                    
                    tabPanel("Network", value = "network",
                             # stop boostrap css from messing up the tooltip in the widget
                             tags$style(HTML(".tooltip {opacity: 1}")), 
                             ndtv:::ndtvAnimationWidgetOutput("netPlot", 
                                                              width="100%", height="800px")
                    ),
                    
                    tabPanel("Gene Ontology", value = "GO",
                             br(),
                             h4("Significant gene ontology terms (Corrected p-value < 0.05)"),
                             br(),
                             downloadButton("dGOTable", "Download Gene Ontology Table"),
                             br(),
                             br(),
                             dataTableOutput("GOTable")),
                    
                    tabPanel("Detected Protein Complexes", value = "complexes",
                             br(),
                             h4("Detected Protein Complexes (> 40% membership)"),
                             br(),
                             downloadButton("dComplexes", "Download Detected Complexes Table"),
                             br(),
                             br(),
                             dataTableOutput("complexesTable")
                    ),
                    
                    tabPanel("Downloads", value = "downloads", 
                             verbatimTextOutput("downloadNetworkInfo"),
                             downloadButton("downloadNodes", "Node Attributes"),
                             downloadButton("downloadEdges", "Edge Attributes"),
                             br(),
                             br(),
                             verbatimTextOutput("downloadReportInfo"),
                             downloadButton("downloadReport", "Report"))
        )
      )
    )
    )
  }
}

##### SHINY INTERACTIVE WEB APP OUTPUT SERVER ####

server <- function(input, output, session) {
  
  ## NETWORK
  
  if (NUM_NODES < MAX_DISPLAY) {
    # switch tabs
    observeEvent(input$netGo, {
      updateTabsetPanel(session, "inTabset", "network")
    })
    
    # calculate inside a reactive to prevent multiple computations
    inputNet <- eventReactive(input$netGo, {
      deactivateNE(sharedCalc(input$shared), input$confidence, input$checkLocalizations, 
                   input$checkSpecies, input$checkGOTerms)
    })
    
    # output network animation
    output$netPlot <- ndtv:::renderNdtvAnimationWidget({
      withProgress(message="Calculating network", detail="This may take a while...", 
                   value=0.8, {
                     plotCompNet(inputNet(), lab=input$labels)
                   })
    })
  } 
  
  ## ABUNDANCE PLOTS
  
  if (ABUND_PROV) {
    # switch tabs
    observeEvent(input$abundGo, {
      updateTabsetPanel(session, "inTabset", "plots")
    })
    
    # find neighbors for given parameters
    confNeighbors <- eventReactive(input$abundGo, {
      calcNeighbors(input$confidence, input$checkLocalizations, 
                    input$checkSpecies, input$checkGOTerms)
    })
    
    # plot abundances
    output$abundPlot <- renderPlot({
      abundPlot(input$genesym, input$neighbors, confNeighbors())
    })
    
    # download abundance plot
    output$dAbund <- downloadHandler(
      filename= function() { 
        "abundance_plot.pdf"}
      ,
      content = function(file) {
        pdf(file)
        abundPlot(input$genesym, input$neighbors, confNeighbors())
        dev.off()
      }
    )
    
    # output cluster number for selected protein
    output$clusterNumber <- renderPrint({
      printClusterNumber(input$genesym)
    })
    
    # plot clusters and heatmaps
    output$clusterPlots <- renderUI({
      plot_output_list <- lapply(1:NUM_BAITS, function(i) {
        plotname <- paste0("plot", i)
        plotOutput(plotname)
      })
      # convert the list to a tagList - this is necessary for the list of 
      # items to display properly
      do.call(tagList, plot_output_list)
    })
    for(i in 1:NUM_BAITS) {
      # Need local so that each item gets its own number. Without it, the value
      # of i in the renderUI() will be the same across all instances, because
      # of when the expression is evaluated
      local({
        my_i <- i
        plotname <- paste0("plot", my_i)
        output[[plotname]] <- renderPlot({
          clusterPlot(BAIT_IDS[my_i])
        })
      })
    }
    
    # download cluster plots
    output$dCluster <- downloadHandler(
      filename= function() { 
        "cluster_plots.pdf"}
      ,
      content = function(file) {
        pdf(file)
        for (bait in BAIT_IDS) {
          clusterPlot(bait)
        }
        dev.off()
      }
    )
    
    if (NORM_PROTEOME) {
      # plot heatmaps
      output$heatmapPlots <- renderUI({
        heatmap_output_list <- lapply(1:NUM_BAITS, function(n) {
          heatmapname <- paste0("heatmap", n)
          plotOutput(heatmapname)
        })
        # convert the list to a tagList - this is necessary for the list of
        # items to display properly
        do.call(tagList, heatmap_output_list)
      })
      for(n in 1:NUM_BAITS) {
        # Need local so that each item gets its own number. Without it, the value
        # of n in the renderUI() will be the same across all instances, because
        # of when the expression is evaluated
        local({
          my_n <- n
          heatmapname <- paste0("heatmap", my_n)
          output[[heatmapname]] <- renderPlot({
            plotHeatmaps(BAIT_IDS[my_n])
          })
        })
      }

      # download heatmap plots
      output$dHeatmaps <- downloadHandler(
        filename= function() {
          "abundance_heatmaps.pdf"}
        ,
        content = function(file) {
          pdf(file)
          for (bait in BAIT_IDS) {
            plotHeatmaps(bait)
          }
          dev.off()
        }
      )
    }
  }
  
  ## GENE ONTOLOGY
  
  # download table
  output$dGOTable <- downloadHandler(
    filename = function() {
      "GO_enrichment_table.txt"}
    ,
    content = function(file) {
      write.table(GOTable, file = file, sep = "\t", 
                  quote = F, row.names = F)
    }
  )
  
  # print table output
  output$GOTable <- renderDataTable(GOTable)
  
  ## COMPLEXES
  
  # download table
  output$dComplexes <- downloadHandler(
    filename = function() {
      "complex_membership.txt"}
    ,
    content = function(file) {
      write.table(detected_complexes, file = file, sep = "\t", 
                  quote = F, row.names = F)
    }
  )
  
  # print table output
  output$complexesTable <- renderDataTable(detected_complexes)
  
  ## DOWNLOADS
  
  # download network information
  output$downloadNetworkInfo <- renderPrint(
    cat("Download tab-separated files with inferred and provided protein attributes 
        (localization, top GO terms, species name, abundances at every time point, etc.), 
        and edge attributes (STRING and provided edges, edge onset and terminus, 
        edge weights at every time point, etc.).")
  )
  
  # output nodes
  output$downloadNodes <- downloadHandler(
    filename= function() { 
      "interactome_node_attributes.txt"}
    ,
    content = function(file) {
      saveNodeAttributes(file)
    }
  )
  
  # output edges
  output$downloadEdges <- downloadHandler(
    filename= function() { 
      "interactome_edge_attributes.txt"}
    ,
    content = function(file) {
      saveEdgeAttributes(file)
    }
  )
  
  if (NUM_NODES < MAX_DISPLAY) {
    # download report information
    output$downloadReportInfo <- renderPrint(
      cat("Download network analysis graphs - spells of protein activity separated by 
          localization and functional assignment, plot of edge formation, and 
          histogram of edge durations.")
      )
    
    # output report
    output$downloadReport <- downloadHandler(
      filename= function() { 
        "interactome_analytical_graphs.pdf"}
      ,
      content = function(file) {
        pdf(file)
        analyticalTools()
        dev.off()
      }
    ) 
  }
  
  # close the R session when Chrome closes
  session$onSessionEnded(function() { 
    stopApp()
    q() 
  })
}

##### DEPLOY SHINY APP ####

shinyApp(ui, server)