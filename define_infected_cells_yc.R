require(pscl)
require(MASS)
require(boot)
require(ggplot2)
require(VGAM)
library(SignallingSingleCell)
library(fitdistrplus)
library(countreg)
###
### find viral load cutoff for each donor that marks true infection
# soup flu distribution and test if I see significantly more cells than expected 
# that comes from soup at a given UMI cutoff
# do first for dog spike in experiment
raw_filt <- readRDS("~/nl/ycao/flu/dogSpikein/hisat2_align/raw_cleaned_filt1000.RData")
pd <- pData(raw_filt)[which(raw_filt$donor == "flu31"),]
x=pd[pd$flu_gene_exp<10000,]
x=x[x$UMI_sum<50000,]

####### Zero inflated NB on human#######
model_zinb = hurdle((flu_gene_exp)~log(UMI_sum)|log(UMI_sum),data=x[x$species=="human",], dist="negbin")
summary(model_zinb)
saveRDS(model_zinb, file = "~/nl/ycao/flu/NovaseqRuns/analysis/defineInfection/flu29_zinb_model.RData")
pd$predict_infection_zinb <- predict(model_zinb, pd, type = "response")
pd$pzero <- round(predprob(model_zinb,newdata=pd)[,1])
phat = predprob(model_zinb,newdata=pd)[,1]
summary(phat)
pd$pvals = 1-pzanegbin(q=pd$flu_gene_exp,size=1/model_zinb$theta,munb=pd$predict_infection_zinb,pobs0=phat)
pd$padj_infection_zinb=p.adjust(pd$pvals,method = "fdr")

ggplot(pd, aes(padj_infection_zinb, color=species)) +
  geom_density()
table(pd$species)
table(pd$species[pd$padj_infection_zinb<=0.05])
countreg::rootogram(model_zinb,style="standing", max = 100, main = "HC4 spiked in with WSN infected MDCK")
countreg::rootogram(model_zinb, max = 100, main = "HC4 spiked in with WSN infected MDCK")
ks.test(pd$predict_infection_zinb[pd$species == "human"], pd$flu_gene_exp[pd$species == "human"])
ks.test(pd$predict_infection_zinb, pd$flu_gene_exp)
cor(pd$predict_infection_zinb[pd$species == "human"], pd$flu_gene_exp[pd$species == "human"], method = "spearman")
?cor

plot(fitted(model_zinb), residuals(model_zinb), pch = 20, main = "HC4 model")
cor(pd$predict_infection_zinb, pd$flu_gene_exp, method = "spearman")

ggplot(pd[pd$species == "human",],aes(flu_gene_exp + 1,predict_infection_zinb + 1)) +
  geom_point() +
  coord_trans(x= "log10", y = "log10")

ggplot(pd,aes(flu_gene_exp + 1,predict_infection_zinb + 1, color=species)) +
  geom_point() +
  coord_trans(x= "log10", y = "log10")


chisq.test(pd$flu_gene_exp[pd$species == "human"], 1- pd$pvals[pd$species == "human"])
plot(ecdf(1- pd$pvals[pd$species == "human"]))

qqrplot(model_zinb,main = "HC4 Hurdle Negative Binomial", pch = 20)
ks.test(x = pd$flu_gene_exp[pd$species == "human"],
        p = 1- pd$pvals[pd$species == "human"])

ggplot(pd) +
  stat_ecdf(aes(flu_gene_exp + 1)) +
  stat_ecdf(aes(predict_infection_zinb+ 1), col = "red") +
  scale_x_log10()

ggplot(pd) +
  geom_histogram(data = pd, aes(predict_infection_zinb+ 1, alpha = 0.5),  fill = "red", bins = 100) +
  geom_histogram(data = pd, aes(flu_gene_exp + 1, alpha = 0.5), bins = 100) +
  scale_x_log10() +
  facet_wrap(~species)

pd$infected_zinb <- "no"
pd[which(pd$padj_infection_zinb <= 0.05),"infected_zinb"] <- "yes"

pd$dog_flu_gene_exp <- pd$flu_gene_exp + pd$dog_gene_exp
pd$dog_flu_gene_frac <- pd$dog_flu_gene_exp /pd$UMI_sum
ggplot()+
  geom_histogram(data = pd[pd$infected_zinb == "yes" & pd$species == "human",], aes(human_gene_frac, alpha = 0.5), bins = 100) +
#  geom_histogram(data = pd[pd$infected_zinb == "no" & pd$species == "human",], aes(human_gene_frac, alpha = 0.5), bins = 100, fill = "red") +
  scale_y_log10()

ggplot(pd, aes(dog_gene_frac, human_gene_frac, color = -log10(padj_infection_zinb)))+
  geom_point() 

ggplot(pd, aes(dog_gene_exp +1, human_gene_exp + 1, color = -log10(padj_infection_zinb)))+
  geom_point() +
  scale_y_log10()+
  scale_x_log10()

pd <- pData(raw_filt)
pd$potential_doublets <- "no"
pd[which(pd$human_gene_frac<=0.9 & pd$species == "human"), "potential_doublets"] <- "doublets"

ggplot(pd[pd$UMI_sum<=5000,], aes(dog_gene_exp, human_gene_exp , color = potential_doublets))+
  geom_point(size = 0.5) +
  facet_wrap(~ donor +  infection_pred)
ggsave("~/nl/ycao/flu/NovaseqRuns/analysis/defineInfection/spikeindonors_prediction.eps")

ggplot(pd, aes((predict_infection_zinb),(flu_gene_exp), color=species))+
  geom_point(aes(alpha=0.2))+scale_y_log10()+scale_x_log10()+facet_wrap(~species)

ggplot(pd,aes(flu_gene_exp, color=species)) +
  geom_histogram()+
  facet_wrap(~infected_zinb)

######## zero inflated NB on human with various amount of infected dog cells ############
# 1% infected dog cells 
getdog <- sample(which(pd$species == "dog" & pd$infected_zinb == "yes"), 0.01 * 2152)
getall <- c(getdog, which(x$species == "human"))
x1 <- x[getall,]
model_zinb_x1 = hurdle((flu_gene_exp)~log(UMI_sum)|log(UMI_sum),data=x1, dist="negbin")
summary(model_zinb_x1)
pd$predict_infection_zinb_0.01dog <- predict(model_zinb_x1, pd, type = "response")
phat = predprob(model_zinb_x1,newdata=pd)[,1]
summary(phat)
pd$pvals_zinb_x1 = 1-pzanegbin(q=pd$flu_gene_exp,size=1/model_zinb_x1$theta,munb=pd$predict_infection_zinb_0.01dog,pobs0=phat)
pd$padj_infection_zinb_x1=p.adjust(pd$pvals_zinb_x1,method = "fdr")
table(pd[which(pd$padj_infection_zinb_x1 <= 0.05),"species"])

countreg::rootogram(model_zinb_x1,style="standing")
countreg::rootogram(model_zinb_x1)

#### 2% infected dog cells
getdog <- sample(which(pd$species == "dog" & pd$infected_zinb == "yes"), 0.02 * 2152)
getall <- c(getdog, which(x$species == "human"))
x2 <- x[getall,]
model_zinb_x2 = hurdle((flu_gene_exp)~log(UMI_sum)|log(UMI_sum),data=x2, dist="negbin")
summary(model_zinb_x2)
pd$predict_infection_zinb_0.02dog <- predict(model_zinb_x2, pd, type = "response")
phat = predprob(model_zinb_x2,newdata=pd)[,1]
summary(phat)
pd$pvals_zinb_x2 = 1-pzanegbin(q=pd$flu_gene_exp,size=1/model_zinb_x2$theta,munb=pd$predict_infection_zinb_0.02dog,pobs0=phat)
pd$padj_infection_zinb_x2=p.adjust(pd$pvals_zinb_x2,method = "fdr")
table(pd[which(pd$padj_infection_zinb_x2 <= 0.05),"species"])

countreg::rootogram(model_zinb_x2,style="standing")
countreg::rootogram(model_zinb_x2)

# 4% infected dog cells 
getdog <- sample(which(pd$species == "dog" & pd$infected_zinb == "yes"), 0.04 * 2152)
getall <- c(getdog, which(x$species == "human"))
x4 <- x[getall,]
model_zinb_x4 = hurdle((flu_gene_exp)~log(UMI_sum)|log(UMI_sum),data=x4, dist="negbin")
summary(model_zinb_x4)
pd$predict_infection_zinb_0.04dog <- predict(model_zinb_x4, pd, type = "response")
phat = predprob(model_zinb_x4,newdata=pd)[,1]
summary(phat)
pd$pvals_zinb_x4 = 1-pzanegbin(q=pd$flu_gene_exp,size=1/model_zinb_x4$theta,munb=pd$predict_infection_zinb_0.04dog,pobs0=phat)
pd$padj_infection_zinb_x4=p.adjust(pd$pvals_zinb_x4,method = "fdr")
table(pd[which(pd$padj_infection_zinb_x4 <= 0.05),"species"])

countreg::rootogram(model_zinb_x4,style="standing")
countreg::rootogram(model_zinb_x4)

# 8% infected dog cells 
getdog <- sample(which(pd$species == "dog" & pd$infected_zinb == "yes"), 0.08 * 2152)
getall <- c(getdog, which(x$species == "human"))
x8 <- x[getall,]
model_zinb_x8 = hurdle((flu_gene_exp)~log(UMI_sum)|log(UMI_sum),data=x8, dist="negbin")
summary(model_zinb_x8)
pd$predict_infection_zinb_0.08dog <- predict(model_zinb_x8, pd, type = "response")
phat = predprob(model_zinb_x8,newdata=pd)[,1]
summary(phat)
pd$pvals_zinb_x8 = 1-pzanegbin(q=pd$flu_gene_exp,size=1/model_zinb_x8$theta,munb=pd$predict_infection_zinb_0.08dog,pobs0=phat)
pd$padj_infection_zinb_x8=p.adjust(pd$pvals_zinb_x8,method = "fdr")
table(pd[which(pd$padj_infection_zinb_x8 <= 0.05),"species"])

countreg::rootogram(model_zinb_x8,style="standing")
countreg::rootogram(model_zinb_x8)

par(mfrow=c(2,3))
countreg::rootogram(model_zinb,style="standing")
countreg::rootogram(model_zinb_x1,style="standing")
countreg::rootogram(model_zinb_x2,style="standing")
countreg::rootogram(model_zinb_x4,style="standing")
countreg::rootogram(model_zinb_x8,style="standing")


countreg::rootogram(model_zinb)
countreg::rootogram(model_zinb_x1)
countreg::rootogram(model_zinb_x2)
countreg::rootogram(model_zinb_x4)
countreg::rootogram(model_zinb_x8)


table(pd[which(pd$padj_infection_zinb <= 0.05),"species"])
table(pd[which(pd$padj_infection_zinb_x1 <= 0.05),"species"])
table(pd[which(pd$padj_infection_zinb_x2 <= 0.05),"species"])
table(pd[which(pd$padj_infection_zinb_x4 <= 0.05),"species"])
table(pd[which(pd$padj_infection_zinb_x8 <= 0.05),"species"])


################ Donors ##########
data <- readRDS("~/nl/ycao/flu/NovaseqRuns/analysis/defineInfection/all_donors_pd.RData")

donors <- c("IAV1","IAV2","IAV3","IAV5","IAV6","IAV7", 
            "IBV1","IBV2","IBV3","IBV4","IBV5","IBV6")

for (i in 1:12){
  d <- donors[i] 
  model_zinb = hurdle((flu_exp)~log(UMI_sum)|log(UMI_sum),
                      data = data[data$Donor == d,], 
                      dist="negbin")
  summary(model_zinb)
  filen <- paste("~/nl/ycao/flu/NovaseqRuns/analysis/defineInfection/",d,"_zinb_model.RData", sep = "")
  print(filen)
 # saveRDS(model_zinb, file = filen)
  data$predict_infection_zinb[data$Donor == d] <- predict(model_zinb, 
                                                           data[data$Donor == d,],
                                                           type = "response")
  data$pzero[data$Donor == d] <- round(predprob(model_zinb,newdata=data[data$Donor == d,])[,1])
  phat = predprob(model_zinb, newdata = data[data$Donor == d,])[,1]
  summary(phat)
  data$pvals[data$Donor == d] = 1-pzanegbin(q = data$flu_exp[data$Donor == d],
                                             size = 1/model_zinb$theta,
                                             munb = data$predict_infection_zinb[data$Donor == d],
                                             pobs0 = phat)
  data$padj_infection_zinb[data$Donor == d] = p.adjust(data$pvals[data$Donor == d], method = "fdr")
}

data$pred_infection <- "no"
data$pred_infection[data$padj_infection_zinb <= 0.05] <- "yes"
table(data$Celltype_2, data$orig)
table(data$Celltype_2, data$pred_infection)
table(data$pred_infection, data$Donor)

model_IAV1 <- readRDS("~/nl/ycao/flu/NovaseqRuns/analysis/defineInfection/IAV1_zinb_model.RData")
model_IAV2 <- readRDS("~/nl/ycao/flu/NovaseqRuns/analysis/defineInfection/IAV2_zinb_model.RData")
model_IAV3 <- readRDS("~/nl/ycao/flu/NovaseqRuns/analysis/defineInfection/IAV3_zinb_model.RData")
model_IAV5 <- readRDS("~/nl/ycao/flu/NovaseqRuns/analysis/defineInfection/IAV5_zinb_model.RData")
model_IAV6 <- readRDS("~/nl/ycao/flu/NovaseqRuns/analysis/defineInfection/IAV6_zinb_model.RData")
model_IAV7 <- readRDS("~/nl/ycao/flu/NovaseqRuns/analysis/defineInfection/IAV7_zinb_model.RData")
model_IBV1 <- readRDS("~/nl/ycao/flu/NovaseqRuns/analysis/defineInfection/IBV1_zinb_model.RData")
model_IBV2 <- readRDS("~/nl/ycao/flu/NovaseqRuns/analysis/defineInfection/IBV2_zinb_model.RData")
model_IBV3 <- readRDS("~/nl/ycao/flu/NovaseqRuns/analysis/defineInfection/IBV3_zinb_model.RData")
model_IBV4 <- readRDS("~/nl/ycao/flu/NovaseqRuns/analysis/defineInfection/IBV4_zinb_model.RData")
model_IBV5 <- readRDS("~/nl/ycao/flu/NovaseqRuns/analysis/defineInfection/IBV5_zinb_model.RData")
model_IBV6 <- readRDS("~/nl/ycao/flu/NovaseqRuns/analysis/defineInfection/IBV6_zinb_model.RData")

summary(model_IAV1)
countreg::rootogram(model_IAV1,style="standing", max = 5000)
chisq.test(x = data[data$Donor == "IAV1","flu_exp"], p = 1- data[data$Donor == "IAV1","pvals"], rescale.p = T)
summary(model_IAV2)
countreg::rootogram(model_IAV2,style="standing", max = 3000)
summary(model_IAV3) # bad
summary(model_IAV5)
countreg::rootogram(model_IAV5,style="standing", max = 2000)
summary(model_IAV6) # bad
summary(model_IAV7)
countreg::rootogram(model_IAV7,style="standing", max = 3000)

summary(model_IBV1)
countreg::rootogram(model_IBV1,style="standing", max = 3000)

summary(model_IBV2)
countreg::rootogram(model_IBV2,style="standing", max = 10000)
summary(model_IBV3) # bad
summary(model_IBV4) # bad
summary(model_IBV5) # bad
summary(model_IBV6) # bad

donor_sub <- c("IAV1","IAV2","IAV5","IAV7","IBV1","IBV2")
before <- as.data.frame(prop.table(table(data$Celltype_2[data$Donor %in% donor_sub], data$orig[data$Donor %in% donor_sub]), 1))
after_IAV <- as.data.frame(prop.table(table(data$Celltype_2[data$Donor %in% donor_sub[1:4]], data$pred_infection[data$Donor %in% donor_sub[1:4]]), 1))
saveRDS(data, file = "~/nl/ycao/flu/NovaseqRuns/analysis/defineInfection/all_donors_pd.RData")

after_IAV$Var1 <- factor(after_IAV$Var1, levels = rev(c("CEP","GOB","MHCIIEPI","Squamous","NEU","EOS","MAST","MAC","pDC","cDC","CD4T","CD8T","proliferating_CD8T",
                                                        "gdT","B","NK")))
col = c("#BEBADA","#FB8072")
ggplot(after_IAV, aes(Var1, Freq, fill = Var2))+
  geom_bar(stat='identity') +
  coord_flip()+
  theme_classic()+
  scale_fill_manual(values = col) +
  ggtitle("IAV")
ggsave("~/nl/ycao/flu/NovaseqRuns/analysis/defineInfection/after_infection_bar_IAV.eps", height = 4, width = 5)

after_IBV <- as.data.frame(prop.table(table(data$Celltype_2[data$Donor %in% donor_sub[5:6]], data$pred_infection[data$Donor %in% donor_sub[5:6]]), 1))
after_IBV$Var1 <- factor(after_IBV$Var1, levels = rev(c("CEP","GOB","MHCIIEPI","Squamous","NEU","EOS","MAST","MAC","pDC","cDC","CD4T","CD8T","proliferating_CD8T",
                                                        "gdT","B","NK")))
ggplot(after_IBV, aes(Var1, Freq, fill = Var2))+
  geom_bar(stat='identity') +
  scale_x_discrete(drop=FALSE) +
  coord_flip()+
  theme_classic()+
  scale_fill_manual(values = col) +
  ggtitle("IBV")
ggsave("~/nl/ycao/flu/NovaseqRuns/analysis/defineInfection/after_infection_bar_IBV.eps", height = 4, width = 5)


