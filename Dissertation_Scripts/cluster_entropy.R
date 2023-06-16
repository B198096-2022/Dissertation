#######################################################################
#                                                                     #
#                       Rquires user input                            #
#               specify directory and merged obejects                 #
#                                                                     #
#######################################################################

outdir <- "path/to/outdir"
all_merged <- "location/of/merged_object"  #The full pilot object 
small_merged <- "location/of/small_merged_object"  #The pilot object without tamburini data sets 
Library(Seurat)


#######################################################################
#                          First Harmony                              #
#                       Rquires user input                            #
#                    specify all patient labels                       #
#                                                                     #
#######################################################################

all_counts <- list()
all_probs <- list()
entropy <- list()
cluster_sizes <- list()
test_cluster <- list()
for (i in 1:length(levels(all_merged$harmony_clusters[1]))){
  clust <- i - 1
  test_cluster[[i]] <- subset(all_merged, subset = (harmony_clusters == paste0(clust)))
  total <- dim(test_cluster[[i]])[2]
  cluster_sizes[i] <- total
  counts <- list()
  counts[1] <- sum(test_cluster[[i]]$patient == "bfc_p1")
  counts[2] <- sum(test_cluster[[i]]$patient == "bfc_p5")
  counts[3] <- sum(test_cluster[[i]]$patient == "bhc_p2")
  counts[4] <- sum(test_cluster[[i]]$patient == "bhc_p4")
  counts[5] <- sum(test_cluster[[i]]$patient == "gfc_p4") 
  counts[6] <- sum(test_cluster[[i]]$patient == "ghc_p3") 
  counts[7] <- sum(test_cluster[[i]]$patient == "rfc_p1")
  counts[8] <- sum(test_cluster[[i]]$patient == "rhc_p3")
  counts[9] <- sum(test_cluster[[i]]$patient == "tfc_p1")
  counts[10] <- sum(test_cluster[[i]]$patient == "tfc_p2")
  counts[11] <- sum(test_cluster[[i]]$patient == "thc_p1") 
  counts[12] <- sum(test_cluster[[i]]$patient == "thc_p2")
  counts[13] <- sum(test_cluster[[i]]$patient == "zfc_p1")
  counts[14] <- sum(test_cluster[[i]]$patient == "zfc_p2")
  counts[15] <- sum(test_cluster[[i]]$patient == "zhc_p1") 
  counts[16] <- sum(test_cluster[[i]]$patient == "zhc_p2")
  counts[17] <- total
  probs <- c()
  for (x in 1:16){
    if (counts[x] == 0){
      counts[x] <- 1
    }
    probs <- c(probs, (as.integer(counts[x])/total) * log2((as.integer(counts[x])/total)))
  }
  all_counts[i] <- counts
  all_probs[i] <- probs
  entropy[i] <- 4 + sum(probs)
}

outdir <- "~/Desktop/big_pilot/entropy/"
saveRDS(entropy, file = paste0(outdir, "harmony_entropy.rds"))
saveRDS(all_counts, file = paste0(outdir, "harmony_cluster_counts.rds"))
saveRDS(all_probs, file = paste0(outdir, "harmony_cluster_probs.rds"))
saveRDS(cluster_sizes, file = paste0(outdir, "harmony_cluster_sizes.rds"))

ent_vec <- c()
for (i in 1:length(entropy)){
  ent_vec <- c(ent_vec, as.numeric(entropy[i]))
}

prob_vec <- c()
for (i in 1:length(all_probs)){
  prob_vec <- c(prob_vec, as.numeric(all_probs[i]))
}

clust_size_vec <- c()
for (i in 1:length(cluster_sizes)){
  clust_size_vec <- c(clust_size_vec, as.numeric(cluster_sizes[i]))
}

weight_ent_vec <- c()
for (i in 1:length(entropy)){
  weight_ent_vec <- c(weight_ent_vec, (clust_size_vec[i]/sum(clust_size_vec))*ent_vec[i])
}

assign("harmony_ent_vec", ent_vec)
assign("harmony_clust_size_vec", clust_size_vec)
assign("harmony_weight_ent_vec", weight_ent_vec)

#saveRDS(entropy, file = paste0(outdir, "unint_entropy.rds"))
#saveRDS(entropy, file = paste0(outdir, "harmony_entropy.rds"))

#SCORE OF 0 MEANS EVENLY SPREAD
#SCORE OF 4 MEANS ONLY A SINGLE PATIENT IN CLUSTER 

ent_vec <- c()
for (i in 1:length(entropy)){
  ent_vec <- c(ent_vec, as.numeric(entropy[i]))
}

  
  
#

#####################################################################################

#######################################################################
#                                                                     #
#                            Now Unintegrated                         #
#                     Use patient labels from abvoe                   #
#                                                                     #
#######################################################################

all_counts <- list()
all_probs <- list()
entropy <- list()
cluster_sizes <- list()
test_cluster <- list()
for (i in 1:length(levels(all_merged$unintegrated_clusters[1]))){
  clust <- i - 1
  test_cluster[[i]] <- subset(all_merged, subset = (unintegrated_clusters == paste0(clust)))
  total <- dim(test_cluster[[i]])[2]
  cluster_sizes[i] <- total
  counts <- list()
  counts[1] <- sum(test_cluster[[i]]$patient == "bfc_p1")
  counts[2] <- sum(test_cluster[[i]]$patient == "bfc_p5")
  counts[3] <- sum(test_cluster[[i]]$patient == "bhc_p2")
  counts[4] <- sum(test_cluster[[i]]$patient == "bhc_p4")
  counts[5] <- sum(test_cluster[[i]]$patient == "gfc_p4") 
  counts[6] <- sum(test_cluster[[i]]$patient == "ghc_p3") 
  counts[7] <- sum(test_cluster[[i]]$patient == "rfc_p1")
  counts[8] <- sum(test_cluster[[i]]$patient == "rhc_p3")
  counts[9] <- sum(test_cluster[[i]]$patient == "tfc_p1")
  counts[10] <- sum(test_cluster[[i]]$patient == "tfc_p2")
  counts[11] <- sum(test_cluster[[i]]$patient == "thc_p1") 
  counts[12] <- sum(test_cluster[[i]]$patient == "thc_p2")
  counts[13] <- sum(test_cluster[[i]]$patient == "zfc_p1")
  counts[14] <- sum(test_cluster[[i]]$patient == "zfc_p2")
  counts[15] <- sum(test_cluster[[i]]$patient == "zhc_p1") 
  counts[16] <- sum(test_cluster[[i]]$patient == "zhc_p2")
  counts[17] <- total
  probs <- c()
  for (x in 1:16){
    if (counts[x] == 0){
      counts[x] <- 1
    }
    probs <- c(probs, (as.integer(counts[x])/total) * log2((as.integer(counts[x])/total)))
  }
  all_counts[i] <- counts
  all_probs[i] <- probs
  entropy[i] <- 4 + sum(probs)
}

outdir <- "~/Desktop/big_pilot/entropy/"
saveRDS(entropy, file = paste0(outdir, "unintegrated_entropy.rds"))
saveRDS(all_counts, file = paste0(outdir, "unintegrated_cluster_counts.rds"))
saveRDS(all_probs, file = paste0(outdir, "unintegrated_cluster_probs.rds"))
saveRDS(cluster_sizes, file = paste0(outdir, "unintegrated_cluster_sizes.rds"))

ent_vec <- c()
for (i in 1:length(entropy)){
  ent_vec <- c(ent_vec, as.numeric(entropy[i]))
}

prob_vec <- c()
for (i in 1:length(all_probs)){
  prob_vec <- c(prob_vec, as.numeric(all_probs[i]))
}

clust_size_vec <- c()
for (i in 1:length(cluster_sizes)){
  clust_size_vec <- c(clust_size_vec, as.numeric(cluster_sizes[i]))
}

weight_ent_vec <- c()
for (i in 1:length(entropy)){
  weight_ent_vec <- c(weight_ent_vec, (clust_size_vec[i]/sum(clust_size_vec))*ent_vec[i])
}

assign("unintegrated_ent_vec", ent_vec)
assign("unintegrated_clust_size_vec", clust_size_vec)
assign("unintegrated_weight_ent_vec", weight_ent_vec)



#####################################################################################

#######################################################################
#                                                                     #
#                               Now scVI                              #
#                     Use patient labels from abvoe                   #
#                                                                     #
#######################################################################

all_counts <- list()
all_probs <- list()
entropy <- list()
cluster_sizes <- list()
test_cluster <- list()
for (i in 1:length(levels(all_merged$scvi_clusters[1]))){
  clust <- i - 1
  test_cluster[[i]] <- subset(all_merged, subset = (scvi_clusters == paste0(clust)))
  total <- dim(test_cluster[[i]])[2]
  cluster_sizes[i] <- total
  counts <- list()
  counts[1] <- sum(test_cluster[[i]]$patient == "bfc_p1")
  counts[2] <- sum(test_cluster[[i]]$patient == "bfc_p5")
  counts[3] <- sum(test_cluster[[i]]$patient == "bhc_p2")
  counts[4] <- sum(test_cluster[[i]]$patient == "bhc_p4")
  counts[5] <- sum(test_cluster[[i]]$patient == "gfc_p4") 
  counts[6] <- sum(test_cluster[[i]]$patient == "ghc_p3") 
  counts[7] <- sum(test_cluster[[i]]$patient == "rfc_p1")
  counts[8] <- sum(test_cluster[[i]]$patient == "rhc_p3")
  counts[9] <- sum(test_cluster[[i]]$patient == "tfc_p1")
  counts[10] <- sum(test_cluster[[i]]$patient == "tfc_p2")
  counts[11] <- sum(test_cluster[[i]]$patient == "thc_p1") 
  counts[12] <- sum(test_cluster[[i]]$patient == "thc_p2")
  counts[13] <- sum(test_cluster[[i]]$patient == "zfc_p1")
  counts[14] <- sum(test_cluster[[i]]$patient == "zfc_p2")
  counts[15] <- sum(test_cluster[[i]]$patient == "zhc_p1") 
  counts[16] <- sum(test_cluster[[i]]$patient == "zhc_p2")
  counts[17] <- total
  probs <- c()
  for (x in 1:16){
    if (counts[x] == 0){
      counts[x] <- 1
    }
    probs <- c(probs, (as.integer(counts[x])/total) * log2((as.integer(counts[x])/total)))
  }
  all_counts[i] <- counts
  all_probs[i] <- probs
  entropy[i] <- 4 + sum(probs)
}

outdir <- "~/Desktop/big_pilot/entropy/"
saveRDS(entropy, file = paste0(outdir, "scvi_entropy.rds"))
saveRDS(all_counts, file = paste0(outdir, "scvi_cluster_counts.rds"))
saveRDS(all_probs, file = paste0(outdir, "scvi_cluster_probs.rds"))
saveRDS(cluster_sizes, file = paste0(outdir, "scvi_cluster_sizes.rds"))

ent_vec <- c()
for (i in 1:length(entropy)){
  ent_vec <- c(ent_vec, as.numeric(entropy[i]))
}

prob_vec <- c()
for (i in 1:length(all_probs)){
  prob_vec <- c(prob_vec, as.numeric(all_probs[i]))
}

clust_size_vec <- c()
for (i in 1:length(cluster_sizes)){
  clust_size_vec <- c(clust_size_vec, as.numeric(cluster_sizes[i]))
}

weight_ent_vec <- c()
for (i in 1:length(entropy)){
  weight_ent_vec <- c(weight_ent_vec, (clust_size_vec[i]/sum(clust_size_vec))*ent_vec[i])
}

assign("scvi_ent_vec", ent_vec)
assign("scvi_clust_size_vec", clust_size_vec)
assign("scvi_weight_ent_vec", weight_ent_vec)


#####################################################################################

#######################################################################
#                                                                     #
#                              Now RPCA                               #
#                     Use patient labels from abvoe                   #
#                                                                     #
#######################################################################


all_counts <- list()
all_probs <- list()
entropy <- list()
cluster_sizes <- list()
test_cluster <- list()
for (i in 1:length(levels(small_merged$rpca_clusters[1]))){
  clust <- i - 1
  test_cluster[[i]] <- subset(small_merged, subset = (rpca_clusters == paste0(clust)))
  total <- dim(test_cluster[[i]])[2]
  cluster_sizes[i] <- total
  counts <- list()
  counts[1] <- sum(test_cluster[[i]]$patient == "bfc_p1")
  counts[2] <- sum(test_cluster[[i]]$patient == "bfc_p5")
  counts[3] <- sum(test_cluster[[i]]$patient == "bhc_p2")
  counts[4] <- sum(test_cluster[[i]]$patient == "bhc_p4")
  counts[5] <- sum(test_cluster[[i]]$patient == "gfc_p4") 
  counts[6] <- sum(test_cluster[[i]]$patient == "ghc_p3") 
  counts[7] <- sum(test_cluster[[i]]$patient == "rfc_p1")
  counts[8] <- sum(test_cluster[[i]]$patient == "rhc_p3")
  counts[9] <- sum(test_cluster[[i]]$patient == "tfc_p1")
  counts[10] <- sum(test_cluster[[i]]$patient == "tfc_p2")
  counts[11] <- sum(test_cluster[[i]]$patient == "thc_p1") 
  counts[12] <- sum(test_cluster[[i]]$patient == "thc_p2")
  counts[13] <- sum(test_cluster[[i]]$patient == "zfc_p1")
  counts[14] <- sum(test_cluster[[i]]$patient == "zfc_p2")
  counts[15] <- sum(test_cluster[[i]]$patient == "zhc_p1") 
  counts[16] <- sum(test_cluster[[i]]$patient == "zhc_p2")
  counts[17] <- total
  probs <- c()
  for (x in 1:16){
    if (counts[x] == 0){
      counts[x] <- 1
    }
    probs <- c(probs, (as.integer(counts[x])/total) * log2((as.integer(counts[x])/total)))
  }
  all_counts[i] <- counts
  all_probs[i] <- probs
  entropy[i] <- 4 + sum(probs)
}

outdir <- "~/Desktop/big_pilot/entropy/"
saveRDS(entropy, file = paste0(outdir, "rpca_entropy.rds"))
saveRDS(all_counts, file = paste0(outdir, "rpca_cluster_counts.rds"))
saveRDS(all_probs, file = paste0(outdir, "rpca_cluster_probs.rds"))
saveRDS(cluster_sizes, file = paste0(outdir, "rpca_cluster_sizes.rds"))

ent_vec <- c()
for (i in 1:length(entropy)){
  ent_vec <- c(ent_vec, as.numeric(entropy[i]))
}

prob_vec <- c()
for (i in 1:length(all_probs)){
  prob_vec <- c(prob_vec, as.numeric(all_probs[i]))
}

clust_size_vec <- c()
for (i in 1:length(cluster_sizes)){
  clust_size_vec <- c(clust_size_vec, as.numeric(cluster_sizes[i]))
}

weight_ent_vec <- c()
for (i in 1:length(entropy)){
  weight_ent_vec <- c(weight_ent_vec, (clust_size_vec[i]/sum(clust_size_vec))*ent_vec[i])
}

assign("rpca_ent_vec", ent_vec)
assign("rpca_clust_size_vec", clust_size_vec)
assign("rpca_weight_ent_vec", weight_ent_vec)

#####################################################################################

#######################################################################
#                                                                     #
#                              Now CCA                                #
#                     Use patient labels from abvoe                   #
#                                                                     #
#######################################################################

all_counts <- list()
all_probs <- list()
entropy <- list()
cluster_sizes <- list()
test_cluster <- list()
for (i in 1:length(levels(small_merged$cca_clusters[1]))){
  clust <- i - 1
  test_cluster[[i]] <- subset(small_merged, subset = (cca_clusters == paste0(clust)))
  total <- dim(test_cluster[[i]])[2]
  cluster_sizes[i] <- total
  counts <- list()
  counts[1] <- sum(test_cluster[[i]]$patient == "bfc_p1")
  counts[2] <- sum(test_cluster[[i]]$patient == "bfc_p5")
  counts[3] <- sum(test_cluster[[i]]$patient == "bhc_p2")
  counts[4] <- sum(test_cluster[[i]]$patient == "bhc_p4")
  counts[5] <- sum(test_cluster[[i]]$patient == "gfc_p4") 
  counts[6] <- sum(test_cluster[[i]]$patient == "ghc_p3") 
  counts[7] <- sum(test_cluster[[i]]$patient == "rfc_p1")
  counts[8] <- sum(test_cluster[[i]]$patient == "rhc_p3")
  counts[9] <- sum(test_cluster[[i]]$patient == "tfc_p1")
  counts[10] <- sum(test_cluster[[i]]$patient == "tfc_p2")
  counts[11] <- sum(test_cluster[[i]]$patient == "thc_p1") 
  counts[12] <- sum(test_cluster[[i]]$patient == "thc_p2")
  counts[13] <- sum(test_cluster[[i]]$patient == "zfc_p1")
  counts[14] <- sum(test_cluster[[i]]$patient == "zfc_p2")
  counts[15] <- sum(test_cluster[[i]]$patient == "zhc_p1") 
  counts[16] <- sum(test_cluster[[i]]$patient == "zhc_p2")
  counts[17] <- total
  probs <- c()
  for (x in 1:16){
    if (counts[x] == 0){
      counts[x] <- 1
    }
    probs <- c(probs, (as.integer(counts[x])/total) * log2((as.integer(counts[x])/total)))
  }
  all_counts[i] <- counts
  all_probs[i] <- probs
  entropy[i] <- 4 + sum(probs)
}

outdir <- "~/Desktop/big_pilot/entropy/"
saveRDS(entropy, file = paste0(outdir, "cca_entropy.rds"))
saveRDS(all_counts, file = paste0(outdir, "cca_cluster_counts.rds"))
saveRDS(all_probs, file = paste0(outdir, "cca_cluster_probs.rds"))
saveRDS(cluster_sizes, file = paste0(outdir, "cca_cluster_sizes.rds"))

ent_vec <- c()
for (i in 1:length(entropy)){
  ent_vec <- c(ent_vec, as.numeric(entropy[i]))
}

prob_vec <- c()
for (i in 1:length(all_probs)){
  prob_vec <- c(prob_vec, as.numeric(all_probs[i]))
}

clust_size_vec <- c()
for (i in 1:length(cluster_sizes)){
  clust_size_vec <- c(clust_size_vec, as.numeric(cluster_sizes[i]))
}

weight_ent_vec <- c()
for (i in 1:length(entropy)){
  weight_ent_vec <- c(weight_ent_vec, (clust_size_vec[i]/sum(clust_size_vec))*ent_vec[i])
}

assign("cca_ent_vec", ent_vec)
assign("cca_clust_size_vec", clust_size_vec)
assign("cca_weight_ent_vec", weight_ent_vec)


#####################################################################################

#######################################################################
#                                                                     #
#        Adding cluster entropy as meta data column                   #
#                                                                     #
#######################################################################

#Harmony
all_merged$harmony_entropy <- 0
for (i in 1:88341){
  cluster <- as.integer(as.character(all_merged$harmony_clusters[[i]]))
  all_merged$harmony_entropy[[i]] <-  harmony_ent_vec[(cluster + 1)]
}
harmony_ent_plot <- FeaturePlot(
  all_merged, 
  c("harmony_entropy"),
  reduction = "umap.harmony")
ggsave(filename = paste0(outdir, "pilot2_harmony_entropy_feature.pdf"), plot = harmony_ent_plot)

####RPCA
small_merged$rpca_entropy <- 0
for (i in 1:dim(small_merged)[2]){
  cluster <- as.integer(as.character(small_merged$rpca_clusters[[i]]))
  small_merged$rpca_entropy[[i]] <-  rpca_ent_vec[(cluster + 1)]
}
rpca_ent_plot <- FeaturePlot(
  small_merged, 
  c("rpca_entropy"),
  reduction = "umap.rpca")
ggsave(filename = paste0(outdir, "pilot2_rpca_entropy_feature.pdf"), plot = rpca_ent_plot)

####CCCA
small_merged$cca_entropy <- 0
for (i in 1:dim(small_merged)[2]){
  cluster <- as.integer(as.character(small_merged$cca_clusters[[i]]))
  small_merged$cca_entropy[[i]] <-  cca_ent_vec[(cluster + 1)]
}
cca_ent_plot <- FeaturePlot(
  small_merged, 
  c("cca_entropy"),
  reduction = "umap.cca")
ggsave(filename = paste0(outdir, "pilot2_cca_entropy_feature.pdf"), plot = cca_ent_plot)

####scVI
all_merged$scvi_entropy <- 0
for (i in 1:88341){
  cluster <- as.integer(as.character(all_merged$scvi_clusters[[i]]))
  all_merged$scvi_entropy[[i]] <-  scvi_ent_vec[(cluster + 1)]
}
scvi_ent_plot <- FeaturePlot(
  all_merged, 
  c("scvi_entropy"),
  reduction = "umap.scvi")
ggsave(filename = paste0(outdir, "pilot2_scvi_entropy_feature.pdf"), plot = scvi_ent_plot)

####Unintegrated
all_merged$unintegrated_entropy <- 0
for (i in 1:88341){
  cluster <- as.integer(as.character(all_merged$unintegrated_clusters[[i]]))
  all_merged$unintegrated_entropy[[i]] <-  unintegrated_ent_vec[(cluster + 1)]
}
unintegrated_ent_plot <- FeaturePlot(
  all_merged, 
  c("unintegrated_entropy"),
  reduction = "umap.unintegrated")
ggsave(filename = paste0(outdir, "pilot2_unintegrated_entropy_feature.pdf"), plot = unintegrated_ent_plot)








