library(SignallingSingleCell)
library(gplots)
# determine if virus+ cells in lowly infected donors are truly infected

####### load data in #######
data <- readRDS("~/nl/ycao/flu/NovaseqRuns/analysis/defineInfection/all_donors_pd.RData")
raw_counts <- readRDS("~/nl/ycao/flu/NovaseqRuns/analysis/raw_clean_donor_removed_102219.RData")
norm_clean <- readRDS("~/nl/ycao/flu/NovaseqRuns/analysis/norm_clean_donor_removed_102219.RData")

donors <- c("IAV1","IAV2","IAV3","IAV5","IAV6","IAV7", 
            "IBV1","IBV2","IBV3","IBV4","IBV5","IBV6")
pred_donors <- c("IAV1","IAV2","IAV5","IAV7","IBV1","IBV2")

pred_infected <- rownames(data)[which(data$pred_donors == "ZINB_applied" & data$pred_infection == "yes")]
unpred_infected <- rownames(data)[which(data$pred_donors == "ZINB_not_applied" & data$orig == "yes")]

norm_clean$pred_infection <- "ZINB_pred_uninfect"
norm_clean$pred_infection[norm_clean$cell_state == "2-Bystander"] <- "Bystander"
norm_clean$pred_infection[which(rownames(pData(norm_clean)) %in% pred_infected)] <- "ZINB_pred_infected"
norm_clean$pred_infection[which(rownames(pData(norm_clean)) %in% unpred_infected)] <- "flu_counts_infected"
norm_clean$pred_infection[norm_clean$cell_state == "1-Healthy" ] <- "Healthy"
pData(raw_counts)[colnames(norm_clean),"pred_infection"] <- norm_clean$pred_infection

norm_clean$metacelltype <- "Non-EPI"
norm_clean$metacelltype[norm_clean$Celltype_2 %in% c("CEP","GOB","MHCIIEPI","Squamous","BasalEPI")] <- "EPI"
pData(raw_counts)[colnames(norm_clean),"metacelltype"] <- norm_clean$metacelltype

saveRDS(pData(norm_clean), file = "~/nl/ycao/flu/NovaseqRuns/analysis/norm_clean_pData_102219.RData")
saveRDS(pData(raw_counts), file = "~/nl/ycao/flu/NovaseqRuns/analysis/raw_clean_pData_102219.RData")

####### run DE on pred infected vs bystander #######
# The input should be raw data, not normalized data
source("~/nl/ycao/code/SCanalysis/edgeR_DE.R")
raw_sub <- raw_counts[,which(raw_counts$pred_infection == "Bystander"|raw_counts$pred_infection =="ZINB_pred_infected")]
for (i in 1:2) {
  type <- unique(raw_counts$metacelltype)[i]
  print(type)
  # run DE
  z <- exprs(raw_sub[,which(raw_sub$metacelltype == type)])
  groupList <- rep(0, times=ncol(z))        # ZINB_pred_infected
  groupList[which(colnames(z) %in% rownames(pData(raw_counts)[which(raw_counts$pred_infection == "Bystander"),]))] <- 1 # Bystander
  group <- factor(groupList)
  tab <- edgeRDE(z, group, contrast = list(c(1,-1),c(-1,1)))
  filen <- paste("~/nl/ycao/flu/NovaseqRuns/analysis/DE/",type,"_predinfect_bystander_DE.txt", sep = "")
  print(filen)
  write.de(tab[['contrast_1']],file= filen)
}

de_epi <- read.table("~/nl/ycao/flu/NovaseqRuns/analysis/DE/EPI_predinfect_bystander_DE.txt", header = T, sep = "\t", row.names = 1)
de_nonepi <- read.table("~/nl/ycao/flu/NovaseqRuns/analysis/DE/Non-EPI_predinfect_bystander_DE.txt", header = T, sep = "\t", row.names = 1)
head(de_epi)
head(de_nonepi)

epi_inf_sig <- c(rownames(de_epi)[which(de_epi$logFC > 1 & de_epi$FDR <= 0.05)], 
                 rownames(de_epi)[which(de_epi$logFC < -1 & de_epi$FDR <= 0.05)])
epi_inf_sig <- epi_inf_sig[-grep("FluB", epi_inf_sig)]
epi_inf_sig <- epi_inf_sig[-grep("H3N2", epi_inf_sig)]
norm_h <- norm_clean[,which(norm_clean$cell_state == "1-Healthy" & norm_clean$Celltype_2 %in% c("CEP","GOB","MHCIIEPI","Squamous","BasalEPI"))]
norm_h <- id_markers(norm_h, id_by = "Celltype_2")
epi_markers <- unique(unname(unlist(return_markers(norm_h, return_by = "Celltype_2", num_markers = 50))))
epi_inf_sig <- epi_inf_sig[-which(epi_inf_sig %in% epi_markers)]

nonepi_inf_sig <- c(rownames(de_nonepi)[which(de_nonepi$logFC > 1 & de_nonepi$FDR <= 0.05)], 
                 rownames(de_nonepi)[which(de_nonepi$logFC < -1 & de_nonepi$FDR <= 0.05)])
nonepi_inf_sig <- nonepi_inf_sig[-grep("FluB", nonepi_inf_sig)]
nonepi_inf_sig <- nonepi_inf_sig[-grep("H3N2", nonepi_inf_sig)]

#### cluster pred infec, unpred infect and bystander cells for each meta-celltype using infection signature ####
# for each cell type, remove healthy cells
norm_flu <- norm_clean[,-which(norm_clean$cell_state == "1-Healthy")]
#norm_epi <- dim_reduce(norm_epi, genelist = epi_inf_sig, pre_reduce = "vPCA", iterations = 500, nVar = 0.85)
#ggplot(pData(norm_epi), aes(x,y, color= pred_infection)) +geom_point()
#ggplot(pData(norm_epi), aes(x,y, color= Celltype_2)) +geom_point()

###### NEU ######
MHCIIEPI <- norm_flu[,which(norm_flu$Celltype_2 == "MHCIIEPI")]
MHCIIEPI_bi <- MHCIIEPI[, which(MHCIIEPI$pred_infection == "Bystander" | MHCIIEPI$pred_infection == "ZINB_pred_infected")]
MHCIIEPI_bi <- id_markers(MHCIIEPI_bi, id_by = "pred_infection")
markerlist <- unname(unlist(return_markers(MHCIIEPI_bi, return_by = "pred_infection", num_markers =650)))
markerlist <- markerlist[-grep("FluB", markerlist)]
markerlist <- markerlist[-grep("H3N2", markerlist)]
markerlist <- unique(markerlist)

cSum <- colSums(exprs(MHCIIEPI)[markerlist,])
rSum <- rowSums(exprs(MHCIIEPI)[markerlist,])

cep_DE <- data.table::fread("~/nl/ycao/flu/NovaseqRuns/analysis/DE/test/CEPpred_bystander_DE.txt")
cep_DE <- cep_DE$V1[c(which(cep_DE$logFC> 1 & cep_DE$FDR <= 0.05),which(cep_DE$logFC < -1 & cep_DE$FDR <= 0.05))]
cep_DE <- cep_DE[-grep("FluB", cep_DE)]
cep_DE <- cep_DE[-grep("H3N2", cep_DE)]
cSum <- colSums(exprs(MHCIIEPI)[MHCIIEPI_DEG,])
unique(pData(MHCIIEPI)[which(cSum ==0),"pred_infection"])
MHCIIEPI_sub <- MHCIIEPI[,-which(cSum ==0)]

#### flow svm prediction ####
pData(MHCIIEPI)[which(MHCIIEPI$pred_infection == "Bystander"), "Pass_Gate"] <- "Bystander"
pData(MHCIIEPI)[which(MHCIIEPI$pred_infection == "ZINB_pred_infected"), "Pass_Gate"] <- "ZINB_pred_infected"
plot_tsne_metadata(MHCIIEPI, color_by = "Pass_Gate")
MHCIIEPI <- flow_svm(MHCIIEPI)
plot_tsne_metadata(MHCIIEPI, color_by = "SVM_Classify")

#### kNN louvain ####
gene_set = unique(MHCIIEPI_DEG)
MHCIIEPI_cor <- cell_cor(input = MHCIIEPI_sub, cor = "spearman", genelist = gene_set, n_blocks = 5)
MHCIIEPI_net <- cell_net(input_cor_matrix = MHCIIEPI_cor, number_edges = 10, weight = "standardized")
MHCIIEPI_sub <- reduce_cell_net(MHCIIEPI_sub, MHCIIEPI_net)
MHCIIEPI_sub <- cluster_cell_net(MHCIIEPI_sub, MHCIIEPI_net, cluster_method = "all")
colnames(pData(MHCIIEPI_sub))
plot_tsne_metadata(MHCIIEPI_sub,  color_by = "pred_infection", xcol = "louvain_x", ycol = "louvain_y")
plot_tsne_metadata(MHCIIEPI_sub,  color_by = "louvain_clusters", xcol = "louvain_x", ycol = "louvain_y")

pData(norm_clean)[colnames(MHCIIEPI_sub),"louvain_x"] <- MHCIIEPI_sub$louvain_x
pData(norm_clean)[colnames(MHCIIEPI_sub),"louvain_y"] <- MHCIIEPI_sub$louvain_y

plot_tsne_gene(norm_clean, gene = c("ITGAM", "FCGR3A","CSF3R"))


## use SVM classifier use pred_infected to predict bystander + unpred infect
norm_epi$Pass_Gate[norm_epi$pred_infection == "ZINB_pred_infected"] <- "ZINB_pred_infected"
norm_epi$Pass_Gate[norm_epi$pred_infection == "Bystander"] <- "Bystander"
norm_epi <- flow_svm(norm_epi)
ggplot(pData(norm_epi), aes(x,y, color= SVM_Classify)) +geom_point()

####### run DE to compare bystander and predicted infected for selected cell types #######
raw_flu <- raw_counts[,-which(raw_counts$cell_state == "1-Healthy")]
celltypes <- c("CEP","GOB","MAC","NEU","MHCIIEPI","Squamous")
for (i in 1:length(celltypes)) {
  type <- celltypes[i]
  print(type)
  # run DE
  z <- exprs(raw_flu[,which(raw_flu$Celltype_2 == type)])
  groupList <- rep(0, times=ncol(z))        # flu_counts_infected
  groupList[which(colnames(z) %in% rownames(pData(raw_flu)[which(raw_flu$pred_infection == "Bystander"),]))] <- 1 # Bystander
  group <- factor(groupList)
  tab <- edgeRDE(z, group, contrast = list(c(1,-1,0),c(1,0,-1)))
  filen <- paste("~/nl/ycao/flu/NovaseqRuns/analysis/DE/test/",type,"pred_bystander_DE.txt", sep = "")
  print(filen)
  write.de(tab[['contrast_1']],file= filen)
}

#raw_healthy <- raw_counts[,which(raw_counts$cell_state == "1-Healthy")]
#raw_healthy <- raw_healthy[,-which(raw_healthy$Celltype_2 == "MAST")]

findDEmarkers(raw_healthy, DEgroup = "Celltype_2", sizefactor = "size_factor", 
              outdir = "~/nl/ycao/flu/NovaseqRuns/analysis/DE/test", outsuffix = "healthy_DE.tsv",pVal = 1,
              batchID = "Patient", pd = pData(raw_healthy), lib_size = "UMI_sum")

###
####### run DE compare the unpred infect to pred infect and bystander #######
raw_sub <- raw_counts[,-which(raw_counts$pred_infection == "no")]
for (i in 1:2) {
  type <- unique(raw_counts$metacelltype)[i]
  print(type)
  # run DE
  z <- exprs(raw_sub[,which(raw_sub$metacelltype == type)])
  groupList <- rep(0, times=ncol(z))        # flu_counts_infected
  groupList[which(colnames(z) %in% rownames(pData(raw_counts)[which(raw_counts$pred_infection == "Bystander"),]))] <- 1 # Bystander
  groupList[which(colnames(z) %in% rownames(pData(raw_counts)[which(raw_counts$pred_infection == "ZINB_pred_infected"),]))] <- 2
  group <- factor(groupList)
  tab <- edgeRDE(z, group, contrast = list(c(1,-1,0),c(1,0,-1)))
  filen <- paste("~/nl/ycao/flu/NovaseqRuns/analysis/DE/",type, sep = "")
  print(filen)
  write.de(tab[['contrast_1']],file= paste(filen, "unpred_bystander_DE.txt", sep = ""))
  write.de(tab[['contrast_2']],file= paste(filen, "unpred_pred_DE.txt", sep = ""))
}

de_epi_test1 <-  read.table("~/nl/ycao/flu/NovaseqRuns/analysis/DE/EPI unpred_bystander_DE.txt", header = T, sep = "\t", row.names = 1)
de_epi_test2 <- read.table("~/nl/ycao/flu/NovaseqRuns/analysis/DE/EPI unpred_pred_DE.txt", header = T, sep = "\t", row.names = 1)
de_nonepi_test1 <-  read.table("~/nl/ycao/flu/NovaseqRuns/analysis/DE/Non-EPI unpred_bystander_DE.txt", header = T, sep = "\t", row.names = 1)
de_nonepi_test2 <- read.table("~/nl/ycao/flu/NovaseqRuns/analysis/DE/Non-EPI unpred_pred_DE.txt", header = T, sep = "\t", row.names = 1)

pred_upgenes <- rownames(de[which(de$logFC>1 & de$FDR < 0.05),])
length(which(pred_upgenes %in% epi_in_upgenes))
pred_upgenes[-which(pred_upgenes %in% epi_in_upgenes)]

unpred_upgenes <- rownames(de[which(de$logFC< -1 & de$FDR < 0.05),])
length(unpred_upgenes)

raw_counts$pred_infection[which(raw_counts$cell_state == "3-Infected" & raw_counts$pred_infection != "ZINB_pred_infected" &
                                  raw_counts$pred_infection != "flu_counts_infected")] <- "ZINB_pred_bystander"
raw_sub <- raw_counts[,-which(raw_counts$pred_infection == "no" | raw_counts$pred_infection == "flu_counts_infected")]
for (i in 1:2) {
  type <- unique(raw_counts$metacelltype)[i]
  print(type)
  # run DE
  z <- exprs(raw_sub[,which(raw_sub$metacelltype == type)])
  groupList <- rep(0, times=ncol(z))        # ZINB_pred_bystander
  groupList[which(colnames(z) %in% rownames(pData(raw_counts)[which(raw_counts$pred_infection == "Bystander"),]))] <- 1 # Bystander
  groupList[which(colnames(z) %in% rownames(pData(raw_counts)[which(raw_counts$pred_infection == "ZINB_pred_infected"),]))] <- 2
  group <- factor(groupList)
  tab <- edgeRDE(z, group, contrast = list(c(1,-1,0),c(1,0,-1)))
  filen <- paste("~/nl/ycao/flu/NovaseqRuns/analysis/DE/",type, sep = "")
  print(filen)
  write.de(tab[['contrast_1']],file= paste(filen, "predbystander_vs_bystander_DE.txt", sep = ""))
  write.de(tab[['contrast_2']],file= paste(filen, "predbystander_vs_infected_DE.txt", sep = ""))
}

### compare predicted low vload cells to unpredicted low vload cells us comparison to bystander cells as control
MHCIIEPI_DE <- data.table::fread("~/nl/ycao/flu/NovaseqRuns/analysis/DE/test/MHCIIEPIpred_bystander_DE.txt")
healthy_DEs <- Sys.glob("~/nl/ycao/flu/NovaseqRuns/analysis/DE/test/*healthy*tsv")
celltype_deg <- c()
for (f in 1:length(healthy_DEs)) {
  file <- data.table::fread(healthy_DEs[f], stringsAsFactors = F)
  celltype <- strsplit(healthy_DEs[f], split = "/")[[1]][11]
  celltype<- gsub("test","",strsplit(celltype, split = "_")[[1]][1])
  print(celltype)
  deg <- file$V1[which(file$logFC > 1 & file$FDR<=0.05)]
  celltype_deg <- c(celltype_deg, deg)
}
celltype_deg <- unique(celltype_deg)

dim(cep_cor)
saveRDS(cep_cor, file ="~/nl/ycao/flu/NovaseqRuns/analysis/DE/test/cep_corr_matrix.RData")
saveRDS(neu_cor, file ="~/nl/ycao/flu/NovaseqRuns/analysis/DE/test/cep_corr_matrix.RData")


###### 1. PCA ########
cep <- norm_clean[,norm_clean$Celltype_2 == "CEP"]
cep <- cep[,-which(cep$cell_state == "1-Healthy")]
cep_DE <- data.table::fread("~/nl/ycao/flu/NovaseqRuns/analysis/DE/test/CEPpred_bystander_DE.txt")
cep_DEG <- cep_DE$V1[c(which(cep_DE$logFC> 1 & cep_DE$FDR <= 0.05),which(cep_DE$logFC < -1 & cep_DE$FDR <= 0.05))]
cep_DEG <- cep_DEG[-grep("FluB", cep_DEG)]
cep_DEG <- cep_DEG[-grep("H3N2", cep_DEG)]

cep <- dim_reduce(cep, genelist = cep_DEG, pre_reduce = "iPCA", nComp = 50, iterations = 500)
plot_tsne_metadata(cep, color_by = "pred_infection", xcol = "iPC_Comp1", ycol="iPC_Comp2")

###### 2. graph based clustering ######
gene_set = unique(cep_DEG)
cSum <- colSums(exprs(cep)[cep_DEG,])
rSum <- rowSums(exprs(cep)[cep_DEG,])

cep_cor <- cell_cor(input = cep, cor = "spearman", genelist = gene_set, n_blocks = 5)
cep_net <- cell_net(input_cor_matrix = cep_cor, number_edges = 30, weight = "standardized")
cep <- reduce_cell_net(cep, cep_net)
cep <- cluster_cell_net(cep, cep_net, cluster_method = "all")
colnames(pData(cep))
plot_tsne_metadata(cep,  color_by = "pred_infection", xcol = "louvain_x", ycol = "louvain_y")
plot_tsne_metadata(cep,  color_by = "vload_state", xcol = "louvain_x", ycol = "louvain_y")
plot_tsne_metadata(cep,  color_by = "louvain_clusters", xcol = "louvain_x", ycol = "louvain_y")

###### 3. correlation #######
reorder_cormat <- function(cormat){
  # Use correlation between variables as distance
  dd <- as.dist((1-cormat)/2)
  hc <- hclust(dd)
  cormat <-cormat[hc$order, hc$order]
}
cep_cor <- reorder_cormat(cep_cor)
cep_cor_anno <- as.data.frame(cep_cor)
cep_cor_anno$cell <- rownames(cep_cor_anno)
cep_cor_anno$pred_infection_var1 <- pData(raw_counts)[cep_cor_anno$cell, "pred_infection"]
cep_cor_anno$louvain_cluster <- pData(cep)[cep_cor_anno$cell, "louvain_clusters"]
cep_cor_anno <- cep_cor_anno[,c(ncol(cep_cor_anno), ncol(cep_cor_anno)-1, ncol(cep_cor_anno)-2,1:(ncol(cep_cor_anno)-3))]
data.table::fwrite(cep_cor_anno, file = "~/Dropbox (UMass Medical School)/test/cep_correlation_matrix.csv", quote = F)
cep_cor_lng <- reshape2::melt(cep_cor_anno)
cep_cor_lng$pred_infection_var1 <- pData(raw_counts)[cep_cor_lng$Var1, "pred_infection"]
cep_cor_lng$pred_infection_var2 <- pData(raw_counts)[cep_cor_lng$Var2, "pred_infection"]
ggplot(cep_cor_lng) +
 geom_tile(aes(Var1, Var2, color=value))+
 # geom_bar(aes(x=100, y=1, fill = pred_infection_var1),stat="identity", width =0.1) +
  theme_classic() +
  scale_color_viridis() +
  theme(axis.text.x=element_blank(), axis.text.y = element_blank())
