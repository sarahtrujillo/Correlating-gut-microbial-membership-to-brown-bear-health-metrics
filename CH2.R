####### Supplemental File 4 Script for Statistical Analysis in R #######

#ST 10/2021

packageVersion(qiime2R)

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("microbiome")
library(microbiome)
library("BiocManager")
## data analysis
install.packages("remotes")
remotes::install_github("jbisanz/qiime2R")
library("remotes")
library(qiime2R) # import data
library(phyloseq) # also the basis of data object. Data analysis and visualization
library(vegan) # some utility tools
library(data.table) # alternative to data.frame
library(dplyr) # data handling
library(tidyverse)
library(DT) ## interactive tables
library(ggpubr) ## plotting 
library(ggplot2)
devtools::install_github("leffj/mctoolsr")
library(mctoolsr)
library(picante) ## faith's PD
library(see)
library(Rmisc)## graphing
library(SRS)
library(cowplot)
library(shiny)
library(glue)
setwd("/Users/sarahtrujillo/Documents/CH2/R/")

#### Import & create phyloseq dataframe with qiime2R and QIIME2 artifacts #####
## Following Tutorial: Integrating QIIME2 and R for data visualization and analysis using qiime2R by J. Bisanz
## you will need
# 1.) Metafile.tsv will need to have the second row removed and the # infront of SampleID removed for it to read okay
# 2.) taxonomy.qza
# 3.) cleantableblahblah-rm.qza
# 4.) filterrooted.qza

## import artifacts & metadata file
metadata1<-read_tsv("brownbearmeta.tsv")
metadata2<-read_tsv("metadata.tsv")
head(metadata2)
metadata<-full_join(metadata1, metadata2)
head(metadata)
SVs<-read_qza("clean-brownbear-table-unassigned_Unknown_Arch-rm.qza")
head(SVs)
taxonomy<-read_qza("brownbear-taxonomy_renamed.qza")
taxtable<-taxonomy$data %>% as_tibble() %>% separate(Taxon, sep=";", c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")) #convert the table into a tabular split version
tree<-read_qza("filter-rooted-brownbear-tree.qza")

metadata$`Body Fat`<- cut(metadata$`Body Fat (%)`,3)
metadata$NetBodyMass<-cut(metadata$NetBodyMasskgs,3)
metadata$`Fat Mass`<-cut(metadata$`Fat Mass (kg)`,3)
metadata$`Lean Mass`<-cut(metadata$`Lean Mass (kg)`,3)

## Create the phyloseq object
phy_obj<-phyloseq(
  otu_table(SVs$data, taxa_are_rows = TRUE), 
  phy_tree(tree$data), 
  tax_table(as.data.frame(taxtable) %>% dplyr::select(-Confidence) %>% column_to_rownames("Feature.ID") %>% as.matrix()), #moving the taxonomy to the way phyloseq wants it
  sample_data(metadata %>% as.data.frame() %>% column_to_rownames("SampleID")))


## view data table
datatable(tax_table(phy_obj))

##### Clean Taxonomy table #####
## Rename NAs to last known group
tax.clean <- data.frame(tax_table(phy_obj))
for (i in 1:7){ tax.clean[,i] <- as.character(tax.clean[,i])}
tax.clean[is.na(tax.clean)] <- ""

for (i in 1:nrow(tax.clean)){
  if (tax.clean[i,2] == ""){
    kingdom <- paste("Kingdom_", tax.clean[i,1], sep = "")
    tax.clean[i, 2:7] <- kingdom
  } else if (tax.clean[i,3] == ""){
    phylum <- paste("Phylum_", tax.clean[i,2], sep = "")
    tax.clean[i, 3:7] <- phylum
  } else if (tax.clean[i,4] == ""){
    class <- paste("Class_", tax.clean[i,3], sep = "")
    tax.clean[i, 4:7] <- class
  } else if (tax.clean[i,5] == ""){
    order <- paste("Order_", tax.clean[i,4], sep = "")
    tax.clean[i, 5:7] <- order
  } else if (tax.clean[i,6] == ""){
    family <- paste("Family_", tax.clean[i,5], sep = "")
    tax.clean[i, 6:7] <- family
  } else if (tax.clean[i,7] == ""){
    tax.clean$Species[i] <- paste("Genus",tax.clean$Genus[i], sep = "_")
  }
}
## import new taxonomy table
tax_table(phy_obj) <- as.matrix(tax.clean)

datatable(tax_table(phy_obj))

###### Rename uncultured
tax.clean2 <- data.frame(tax_table(phy_obj))

for (i in 1:7){ tax.clean2[,i] <- as.character(tax.clean2[,i])}
for (i in 1:nrow(tax.clean2)){
  if (tax.clean2[i,2] == "uncultured"){
    kingdom <- paste("Kingdom_", tax.clean2[i,1], sep = "")
    tax.clean2[i, 2:7] <- kingdom
  } else if (tax.clean2[i,3] == "uncultured"){
    phylum <- paste("Phylum_", tax.clean2[i,2], sep = "")
    tax.clean2[i, 3:7] <- phylum
  } else if (tax.clean2[i,4] == "uncultured"){
    class <- paste("Class_", tax.clean2[i,3], sep = "")
    tax.clean2[i, 4:7] <- class
  } else if (tax.clean2[i,5] == "uncultured"){
    order <- paste("Order_", tax.clean2[i,4], sep = "")
    tax.clean2[i, 5:7] <- order
  } else if (tax.clean2[i,6] == "uncultured"){
    family <- paste("Family_", tax.clean2[i,5], sep = "")
    tax.clean2[i, 6:7] <- family
  } else if (tax.clean2[i,7] == ""){
    tax.clean2$Species[i] <- paste("Genus",tax.clean.bear$Genus[i], sep = "_")
  }
}

## import new taxonomy table
tax_table(phy_obj) <- as.matrix(tax.clean2)

datatable(tax_table(phy_obj))

saveRDS(phy_obj, "~/physeq.rds")

phy_obj<- readRDS("physeq.rds")

#SRS

set.seed(9242)
otu<-as.data.frame(otu_table(phy_obj))
SRS.shiny.app(otu)

#normalize
new_otu<-as.matrix(SRS(otu, Cmin=4087, set_seed = T, seed=9242), rownames=T)
# lost 4 samples
rownames(new_otu)<-rownames(otu)
srs_obj<-phyloseq(otu_table(new_otu, taxa_are_rows = T),
                  phy_tree(phy_obj), 
                  tax_table(phy_obj),
                  sample_data(phy_obj))

summary(sample_sums(srs_obj))
saveRDS(srs_obj, "~/physeq_srs.rds")
srs_obj<- readRDS("physeq_srs.rds")
                                

####  Alpha Diversity ####
##  pull metadata from physeq object
sam.meta <- meta(srs_obj)
sam.meta

## Add the rownames as a new colum for easy integration later.
sam.meta$sam_name <- rownames(sam.meta)

##### Non-phylogenetic diversities: Shannon ####
## calculated with microbiome package
div_shan<- microbiome::alpha(srs_obj, index = "shannon")

## Add the rownames to diversity table
div_shan$sam_name <- rownames(div_shan)
div_shan$sam_name<-factor(div_shan$sam_name)
write.csv(div_shan,'shannon_diversity.csv')
div_shan<- read.csv('shannon_diversity.csv')

##### Non-phylogenetic diversities: Simpson ####
## calculated with microbiome package
div_sim<- microbiome::alpha(srs_obj, index="diversity_inverse_simpson") 

## Add the rownames to diversity table
div_sim$sam_name <- rownames(div_sim)
as.factor("sam_name")

##### Phylogenetic diversity: Faith's PD #####
#Phylogenetic diversity is calculated using the picante package.

## pull ASV table
phyb.rar.asvtab <- as.data.frame(srs_obj@otu_table)

## pull tree
phyb.rar.tree <- srs_obj@phy_tree

##check if the tree is rooted or not 

srs_obj@phy_tree
###rooted so we are good to go

## Getting the data ready
div_pd <- pd(t(phyb.rar.asvtab), phyb.rar.tree,include.root=T) 
# t(ou_table) transposes the table for use in picante and the
#tree file comes from the first code  we used to read tree
#file (see making a phyloseq object section)

## Add the rownames to diversity table
div_pd$sam_name <- rownames(div_pd)

##merge all of the alphas into one file
merged_table2<-merge(div_pd,div_shan, by = "sam_name", all=T)
merged_table3<-merge(merged_table2,sam.meta, by = "sam_name", all=T)
alpha_table <- merge(merged_table3,div_sim, by = "sam_name", all=T)

datatable(alpha_table)

#Visualization
library(scico)

alpha_table$Body.Fat2<-cut(alpha_table$Body.Fat....,3, labels = c("Below Median", "Median", "Above Median"))
alpha_table$NetBodyMass2<- cut(alpha_table$NetBodyMasskgs, 3, labels = c("Below Median", "Median", "Above Median"))
alpha_table$Lean.Mass2<- cut(alpha_table$Lean.Mass..kg., 3,labels = c("Below Median", "Median", "Above Median"))
alpha_table$Fat.Mass2<- cut(alpha_table$Fat.Mass..kg., 3, labels = c("Below Median", "Median", "Above Median"))

alpha_table$Park<-factor(alpha_table$Park, levels=c( "KATM", "LACL", "GAAR"))
alpha_table$Body.Fat2<-factor(alpha_table$Body.Fat2, levels=c( "Below Median", "Median", "Above Median"))
alpha_table$NetBodyMass2<-factor(alpha_table$NetBodyMass2, levels=c( "Below Median", "Median", "Above Median"))
alpha_table$Lean.Mass2<-factor(alpha_table$Lean.Mass2, levels=c( "Below Median", "Median", "Above Median"))
alpha_table$Fat.Mass2<-factor(alpha_table$Fat.Mass2, levels=c( "Below Median", "Median", "Above Median"))

alpha_table<-alpha_table[complete.cases(alpha_table$Body.Fat2),]
alpha_table<-alpha_table[complete.cases(alpha_table$NetBodyMass2),]
alpha_table<-alpha_table[complete.cases(alpha_table$Lean.Mass2),]
alpha_table<-alpha_table[complete.cases(alpha_table$Fat.Mass2),]

bodyfat_pd<-ggplot(data = alpha_table, aes(x=Body.Fat2, y=PD))+
  geom_boxplot()+
  ylab("Faith's PD")+
  xlab("% Body Fat")+
  theme_bw()+
  theme(axis.text.x = element_text(size=8, color="black"),axis.text.y = element_text(size=8, color="black"),axis.ticks = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none")+
  scale_fill_manual(values=c("#FDB4B2", "#A58B2D","#185364"))
bodyfat_pd

netmass_pd<-ggplot(data = alpha_table, aes(x=NetBodyMass2, y=PD))+
  geom_boxplot()+
  ylab("Faith's PD")+
  xlab("Net Mass")+
  theme_bw()+
  theme(axis.text.x = element_text(size=8, color="black"),axis.text.y = element_text(size=8, color="black"),axis.ticks = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none")
netmass_pd

leanmass_pd<-ggplot(data = alpha_table, aes(x=Lean.Mass2, y=PD))+
  geom_boxplot()+
  ylab("Faith's PD")+
  xlab("Lean Mass")+
  theme_bw()+
  theme(axis.text.x = element_text(size=8, color="black"),axis.text.y = element_text(size=8, color="black"),axis.ticks = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none")
leanmass_pd

fatmass_PD<-ggplot(data = alpha_table, aes(x=Fat.Mass2, y=PD))+
  geom_boxplot()+
  ylab("Faith's PD")+
  xlab("Fat Mass")+
  theme_bw()+
  theme(axis.text.x = element_text(size=8, color="black"),axis.text.y = element_text(size=8, color="black"),legend.position = "none",panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(data = alpha_table, aes(x=Body.Fat2, y=diversity_shannon))+
  geom_boxplot()+
  ylab("")+
  xlab(" ")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none")

ggplot(data = alpha_table, aes(x=Body.Fat2, y=diversity_inverse_simpson))+
  geom_boxplot()+
  ylab("Inverse Simpson")+
  xlab(" ")+
  ylim(0,100)+
  theme_bw()+
  theme(legend.position = "none",panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 

ggplot(data = alpha_table, aes(x=NetBodyMass2, y=diversity_shannon))+
  geom_boxplot()+
  ylab("")+
  xlab(" ")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none")

ggplot(data = alpha_table, aes(x=NetBodyMass2, y=diversity_inverse_simpson))+
  geom_boxplot()+
  ylab("Inverse Simpson")+
  xlab(" ")+
  ylim(0,100)+
  theme_bw()+
  theme(legend.position = "none",panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 

ggplot(data = alpha_table, aes(x=Lean.Mass2, y=diversity_shannon))+
  geom_boxplot()+
  ylab("")+
  xlab(" ")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none")

ggplot(data = alpha_table, aes(x=Lean.Mass2, y=diversity_inverse_simpson))+
  geom_boxplot()+
  ylab("Inverse Simpson")+
  xlab(" ")+
  ylim(0,100)+
  theme_bw()+
  theme(legend.position = "none",panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 

ggplot(data = alpha_table, aes(x=Fat.Mass2, y=diversity_shannon))+
  geom_boxplot()+
  ylab("")+
  xlab(" ")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none")

ggplot(data = alpha_table, aes(x=Fat.Mass2, y=diversity_inverse_simpson))+
  geom_boxplot()+
  ylab("Inverse Simpson")+
  xlab(" ")+
  ylim(0,100)+
  theme_bw()+
  theme(legend.position = "none",panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 

summarySE(alpha_table, measurevar = "PD", groupvars =c("Body.Fat2"))
summarySE(alpha_table, measurevar = "PD", groupvars =c("NetBodyMass2"))
summarySE(alpha_table, measurevar = "PD", groupvars =c("Lean.Mass2"))
summarySE(alpha_table, measurevar = "PD", groupvars =c("Fat.Mass2"))
summarySE(alpha_table, measurevar = "diversity_inverse_simpson", groupvars=c("Body.Fat2"))
summarySE(alpha_table, measurevar = "diversity_inverse_simpson",groupvars =c("NetBodyMass2"))
summarySE(alpha_table, measurevar = "diversity_inverse_simpson", groupvars=c("Lean.Mass2"))
summarySE(alpha_table, measurevar = "diversity_inverse_simpson", groupvars=c("Fat.Mass2"))
summarySE(alpha_table, measurevar = "diversity_shannon", groupvars =c("Body.Fat2"))
summarySE(alpha_table, measurevar = "diversity_shannon", groupvars =c("NetBodyMass2"))
summarySE(alpha_table, measurevar = "diversity_shannon", groupvars =c("Lean.Mass2"))
summarySE(alpha_table, measurevar = "diversity_shannon", groupvars =c("Fat.Mass2"))

#####alpha diversity significant differences####

library(car)
library(RVAideMemoire)

byf.hist(data=alpha_table, PD~Body.Fat2, density = TRUE, sep = FALSE) #plot for every comparison

byf.shapiro(data=alpha_table, log(PD)~Body.Fat2)#pass
leveneTest(data=alpha_table, log(PD)~Body.Fat2)#pass
byf.shapiro(data=alpha_table, log(PD)~NetBodyMass2)#not enough above median samples
byf.shapiro(data=alpha_table, log(PD)~Lean.Mass2)#pass
leveneTest(data=alpha_table, log(PD)~Lean.Mass2)#pass
byf.shapiro(data=alpha_table, log(PD)~Fat.Mass2)#pass
leveneTest(data=alpha_table, log(PD)~Fat.Mass2)#pass

byf.shapiro(data=alpha_table, log(diversity_shannon)~Body.Fat2)#pass
leveneTest(data=alpha_table, log(diversity_shannon)~Body.Fat2)#pass
byf.shapiro(data=alpha_table, log(diversity_shannon)~NetBodyMass2)#not enough above median samples
byf.shapiro(data=alpha_table, log(diversity_shannon)~Lean.Mass2)#pass
leveneTest(data=alpha_table, log(diversity_shannon)~Lean.Mass2)#pass
byf.shapiro(data=alpha_table, log(diversity_shannon)~Fat.Mass2)#pass
leveneTest(data=alpha_table, log(diversity_shannon)~Fat.Mass2)#pass

byf.shapiro(data=alpha_table, log(diversity_inverse_simpson)~Body.Fat2)#pass
leveneTest(data=alpha_table, log(diversity_inverse_simpson)~Body.Fat2)#pass
byf.shapiro(data=alpha_table, log(diversity_inverse_simpson)~NetBodyMass2)#not enough above median samples
byf.shapiro(data=alpha_table, log(diversity_inverse_simpson)~Lean.Mass2)#below median not normal
leveneTest(data=alpha_table, log(diversity_inverse_simpson)~Lean.Mass2)#pass
byf.shapiro(data=alpha_table, log(diversity_inverse_simpson)~Fat.Mass2)#below median not normal
leveneTest(data=alpha_table, log(diversity_inverse_simpson)~Fat.Mass2)#pass

mod1<-aov(data=alpha_table, PD~Body.Fat2)
summary(mod1)
plot(mod1)
mod2<-aov(data=alpha_table, diversity_shannon~Body.Fat2)
summary(mod2)
with(alpha_table, par(mfrow=c(2,2)))
plot(mod2)
mod3<-aov(data=alpha_table, diversity_inverse_simpson~Body.Fat2)
summary(mod3)
with(alpha_table, par(mfrow=c(2,2)))
plot(mod3)

mod1b<-aov(data=alpha_table, PD~NetBodyMass2)
summary(mod1b)
plot(mod1b)
mod2b<-aov(data=alpha_table, diversity_shannon~NetBodyMass2)
summary(mod2b)
with(alpha_table, par(mfrow=c(2,2)))
plot(mod2b)
mod3b<-aov(data=alpha_table, diversity_inverse_simpson~NetBodyMass2)
summary(mod3b)
with(alpha_table, par(mfrow=c(2,2)))
plot(mod3b)

mod4<-aov(data=alpha_table, PD~Lean.Mass2)
summary(mod4)
plot(mod4)
mod5<-aov(data=alpha_table, diversity_shannon~Lean.Mass2)
summary(mod5)
with(alpha_table, par(mfrow=c(2,2)))
plot(mod5)
mod6<-aov(data=alpha_table, diversity_inverse_simpson~Lean.Mass2)
summary(mod6)
with(alpha_table, par(mfrow=c(2,2)))
plot(mod6)

mod7<-aov(data=alpha_table, PD~Fat.Mass2)
summary(mod7)
plot(mod7)
mod8<-aov(data=alpha_table, diversity_shannon~Fat.Mass2)
summary(mod8)
with(alpha_table, par(mfrow=c(2,2)))
plot(mod8)
mod9<-aov(data=alpha_table, diversity_inverse_simpson~Fat.Mass2)
summary(mod9)
with(alpha_table, par(mfrow=c(2,2)))
plot(mod9)

#####Alpha diversity correlation####

ggplot(alpha_table, aes(x=PD, y=Body.Fat....))+
  geom_point() #run for every cor.test

cor.test(alpha_table$PD, alpha_table$Body.Fat...., method="spearman", use="complete.obs", exact=FALSE)
cor.test(alpha_table$PD, alpha_table$NetBodyMasskgs, method="spearman", use="complete.obs", exact=FALSE)
cor.test(alpha_table$PD, alpha_table$Lean.Mass..kg., method="spearman", use="complete.obs", exact=FALSE)
cor.test(alpha_table$PD, alpha_table$Fat.Mass..kg., method="spearman", use="complete.obs", exact=FALSE)

cor.test(alpha_table$diversity_shannon, alpha_table$Body.Fat...., method="spearman", use="complete.obs", exact=FALSE)
cor.test(alpha_table$diversity_shannon, alpha_table$NetBodyMasskgs, method="spearman", use="complete.obs", exact=FALSE)
cor.test(alpha_table$diversity_shannon, alpha_table$Lean.Mass..kg., method="spearman", use="complete.obs", exact=FALSE)
cor.test(alpha_table$diversity_shannon, alpha_table$Fat.Mass..kg., method="spearman", use="complete.obs", exact=FALSE)

cor.test(alpha_table$diversity_inverse_simpson, alpha_table$Body.Fat...., method="spearman", use="complete.obs", exact=FALSE)
cor.test(alpha_table$diversity_inverse_simpson, alpha_table$NetBodyMasskgs, method="spearman", use="complete.obs", exact=FALSE)
cor.test(alpha_table$diversity_inverse_simpson, alpha_table$Lean.Mass..kg., method="spearman", use="complete.obs", exact=FALSE)
cor.test(alpha_table$diversity_inverse_simpson, alpha_table$Fat.Mass..kg., method="spearman", use="complete.obs", exact=FALSE)

#### Community composition  ######
## filter
# Remove taxa not seen more than 5 times in at least 20% of the samples
## relative abundance
pseq.rel <- microbiome::transform(srs_obj, "compositional")

## merge to phylum rank
phlyum <- tax_glom(pseq.rel, taxrank = "Phylum")
ntaxa(phlyum)
#44
## melt
phylum_melt<- psmelt(phlyum)

unique(phylum_melt$Phylum)
#44
## get summary statistics phyla 

summarySE(phylum_melt, measurevar = "Abundance", groupvars= "Phylum")

p_abund<-summarySE(phylum_melt, measurevar = "Abundance", groupvars =c("Phylum","Body.Fat"))
p_abund%>%
  group_by(Phylum,Body.Fat) %>%
  dplyr::summarise(sum(Abundance))%>%
  View()
p_abund$Abundance[p_abund$Abundance==0] <- NA
p_abund<-p_abund[complete.cases(p_abund$Abundance),]
p_abund<- p_abund %>% 
  mutate_if(is.numeric, round, digits = 5)
p_abund<-p_abund[complete.cases(p_abund$Body.Fat),]

p_abund_netmass<-summarySE(phylum_melt, measurevar = "Abundance", groupvars =c("Phylum", "NetBodyMass"))
p_abund_netmass%>%
  group_by(Phylum,NetBodyMass) %>%
  dplyr::summarise(sum(Abundance))%>%
  View()
p_abund_netmass$Abundance[p_abund_netmass$Abundance==0] <- NA
p_abund_netmass<-p_abund_netmass[complete.cases(p_abund_netmass$Abundance),]
p_abund_netmass<- p_abund_netmass %>% 
  mutate_if(is.numeric, round, digits = 5)
p_abund_netmass<-p_abund_netmass[complete.cases(p_abund_netmass$NetBodyMass),]

p_abund_lean<-summarySE(phylum_melt, measurevar = "Abundance", groupvars =c("Phylum", "Lean.Mass"))
summarySE(phylum_melt, measurevar = "Abundance", groupvars =c("Firmicutes"))
p_abund_lean$Abundance[p_abund_lean$Abundance==0] <- NA
p_abund_lean<-p_abund_lean[complete.cases(p_abund_lean$Abundance),]
p_abund_lean<- p_abund_lean %>% 
  mutate_if(is.numeric, round, digits = 5)
p_abund_lean<-p_abund_lean[complete.cases(p_abund_lean$Lean.Mass),]

p_abund_fat<-summarySE(phylum_melt, measurevar = "Abundance", groupvars =c("Phylum", "Fat.Mass"))
p_abund_fat%>%
  group_by(Phylum,Fat.Mass) %>%
  dplyr::summarise(sum(Abundance))%>%
  View()
p_abund_fat$Abundance[p_abund_fat$Abundance==0] <- NA
p_abund_fat<-p_abund_fat[complete.cases(p_abund_fat$Abundance),]
p_abund_fat<- p_abund_fat %>% 
  mutate_if(is.numeric, round, digits = 5)
p_abund_fat<-p_abund_fat[complete.cases(p_abund_fat$Fat.Mass),]

## genus 
## merge to Genus rank
rel_genus <- tax_glom(pseq.rel, taxrank = "Genus")
genus_melt<- psmelt(rel_genus)
ntaxa(rel_genus)

genus<-summarySE(genus_melt, measurevar = "Abundance", groupvars = c("Genus", "Phylum"))
genus$Abundance[genus$Abundance==0] <- NA
genus<-genus[complete.cases(genus$Abundance),]
genus<- genus %>% 
  mutate_if(is.numeric, round, digits = 5)
genus$Genus[genus$Abundance <= 0.01] <- "Minor"

g_abund<-summarySE(genus_melt, measurevar = "Abundance", groupvars =c("Genus", "Body.Fat"))

##remove 0 abundance
g_abund$Abundance[g_abund$Abundance==0] <- NA
g_abund<-g_abund[complete.cases(g_abund$Abundance),]
g_abund<- g_abund %>% 
  mutate_if(is.numeric, round, digits = 5)
g_abund<-g_abund[complete.cases(g_abund$Body.Fat),]
g_abund%>%
  group_by(Genus,Body.Fat) %>%
  dplyr::summarise(sum(Abundance))%>%
  View()

g_abund_netmass<-summarySE(genus_melt, measurevar = "Abundance", groupvars =c("Genus", "NetBodyMass"))
g_abund_netmass$Abundance[g_abund_netmass$Abundance==0] <- NA
g_abund_netmass<-g_abund_netmass[complete.cases(g_abund_netmass$Abundance),]
g_abund_netmass<- g_abund_netmass %>% 
  mutate_if(is.numeric, round, digits = 5)
g_abund_netmass<-g_abund_netmass[complete.cases(g_abund_netmass$NetBodyMass),]
g_abund_netmass%>%
  group_by(Genus,NetBodyMass) %>%
  dplyr::summarise(sum(Abundance))%>%
  View()

g_abund_lean<-summarySE(genus_melt, measurevar = "Abundance", groupvars =c("Genus", "Lean.Mass"))
g_abund_lean$Abundance[g_abund_lean$Abundance==0] <- NA
g_abund_lean<-g_abund_lean[complete.cases(g_abund_lean$Abundance),]
g_abund_lean<- g_abund_lean %>% 
  mutate_if(is.numeric, round, digits = 5)
g_abund_lean<-g_abund_lean[complete.cases(g_abund_lean$Lean.Mass),]
g_abund_lean%>%
  group_by(Genus,Lean.Mass) %>%
  dplyr::summarise(sum(Abundance))%>%
  View()

g_abund_fat<-summarySE(genus_melt, measurevar = "Abundance", groupvars =c("Genus", "Fat.Mass"))
g_abund_fat$Abundance[g_abund_fat$Abundance==0] <- NA
g_abund_fat<-g_abund_fat[complete.cases(g_abund_fat$Abundance),]
g_abund_fat<- g_abund_fat %>% 
  mutate_if(is.numeric, round, digits = 5)
g_abund_fat<-g_abund_fat[complete.cases(g_abund_fat$Fat.Mass),]
g_abund_fat%>%
  group_by(Genus,Fat.Mass) %>%
  dplyr::summarise(sum(Abundance))%>%
  View()

write.csv(g_abund,'gabund.csv')

#####Major phyla significant differences ####

supplement1a<-summarySE(P_kruskal, measurevar = "Abundance", groupvars = c("Phylum", "Body.Fat"))
supplement1b<-summarySE(P_kruskal, measurevar = "Abundance", groupvars = c("Phylum", "NetBodyMass"))
supplement1c<-summarySE(P_kruskal, measurevar = "Abundance", groupvars = c("Phylum", "Lean.Mass"))
supplement1d<-summarySE(P_kruskal, measurevar = "Abundance", groupvars = c("Phylum", "Fat.Mass"))

P_kruskal<-phylum_melt %>%
  select("Abundance", "Phylum", "Body.Fat", "NetBodyMass", "Lean.Mass", "Fat.Mass")
P_kruskal$Abundance[P_kruskal$Abundance==0] <- NA
P_kruskal<-P_kruskal[complete.cases(P_kruskal$Abundance),]
P_kruskal<- P_kruskal%>% 
  mutate_if(is.numeric, round, digits = 5)

actino_kw<- P_kruskal%>%
  filter(Phylum == "Actinobacteria")
bacteroidetes_kw<- P_kruskal%>%
  filter(Phylum == "Bacteroidetes")
epsilon_kw<- P_kruskal%>%
  filter(Phylum == "Epsilonbacteraeota")
firmicutes_kw<- P_kruskal%>%
  filter(Phylum == "Firmicutes")
proteo_kw<- P_kruskal%>%
  filter(Phylum == "Proteobacteria")
tenericutes_kw<- P_kruskal%>%
  filter(Phylum == "Tenericutes")

byf.hist(data=actino_kw, Abundance~Body.Fat, density = TRUE, sep = FALSE) #plot for every comparison

byf.shapiro(data=actino_kw, log(Abundance)~Body.Fat)#pass
leveneTest(data=actino_kw, log(Abundance)~Body.Fat)#pass
byf.shapiro(data=actino_kw, log(Abundance)~NetBodyMass)#not enough above median samples
byf.shapiro(data=actino_kw, log(Abundance)~Lean.Mass)#pass
leveneTest(data=actino_kw, log(Abundance)~Lean.Mass)#pass
byf.shapiro(data=actino_kw, log(Abundance)~Fat.Mass)#pass
leveneTest(data=actino_kw, log(Abundance)~Fat.Mass)#pass

byf.shapiro(data=bacteroidetes_kw, log(Abundance)~Body.Fat)#pass
leveneTest(data=bacteroidetes_kw, log(Abundance)~Body.Fat)#pass
byf.shapiro(data=bacteroidetes_kw, log(Abundance)~NetBodyMass)#pass
leveneTest(data=bacteroidetes_kw, log(Abundance)~NetBodyMass)#not enough above median samples
leveneTest(data=bacteroidetes_kw, log(Abundance)~Lean.Mass)#pass
byf.shapiro(data=bacteroidetes_kw, log(Abundance)~Fat.Mass)#pass
leveneTest(data=bacteroidetes_kw, log(Abundance)~Fat.Mass)#pass

byf.shapiro(data=epsilon_kw, log(Abundance)~Body.Fat)#pass
leveneTest(data=epsilon_kw, log(Abundance)~Body.Fat)#pass
byf.shapiro(data=epsilon_kw, log(Abundance)~NetBodyMass)#pass
leveneTest(data=epsilon_kw, log(Abundance)~NetBodyMass)#not enough above median samples
leveneTest(data=epsilon_kw, log(Abundance)~Lean.Mass)#pass
byf.shapiro(data=epsilon_kw, log(Abundance)~Fat.Mass)#pass
leveneTest(data=epsilon_kw, log(Abundance)~Fat.Mass) #pass

byf.shapiro(data=firmicutes_kw, log(Abundance)~Body.Fat)#fail
leveneTest(data=firmicutes_kw, log(Abundance)~Body.Fat)#pass
byf.shapiro(data=firmicutes_kw, log(Abundance)~NetBodyMass)#not enough above median samples
byf.shapiro(data=firmicutes_kw, log(Abundance)~Lean.Mass)#fail
leveneTest(data=firmicutes_kw, log(Abundance)~Lean.Mass)#fail
byf.shapiro(data=firmicutes_kw, log(Abundance)~Fat.Mass)#fail
leveneTest(data=firmicutes_kw, log(Abundance)~Fat.Mass)#pass

byf.shapiro(data=proteo_kw, log(Abundance)~Body.Fat)#fail
leveneTest(data=proteo_kw, log(Abundance)~Body.Fat)#pass
byf.shapiro(data=proteo_kw, log(Abundance)~NetBodyMass)#not enough above median samples
byf.shapiro(data=proteo_kw, Abundance~Lean.Mass)#fail
leveneTest(data=proteo_kw, log(Abundance)~Lean.Mass)#pass
byf.shapiro(data=proteo_kw, Abundance~Fat.Mass)#fail
leveneTest(data=proteo_kw, log(Abundance)~Fat.Mass)#pass

byf.shapiro(data=tenericutes_kw, log(Abundance)~Body.Fat)#not enough above median samples
byf.shapiro(data=tenericutes_kw, log(Abundance)~NetBodyMass)#not enough above median samples
byf.shapiro(data=tenericutes_kw, log(Abundance)~Lean.Mass)#pass
leveneTest(data=tenericutes_kw, log(Abundance)~Lean.Mass)#pass
byf.shapiro(data=tenericutes_kw, log(Abundance)~Fat.Mass)#not enough above median samples

mod10<-aov(data=actino_kw, Abundance~Body.Fat)
summary(mod10)
plot(mod10)
mod11<-aov(data=actino_kw, Abundance~NetBodyMass)
summary(mod11)
with(actino_kw, par(mfrow=c(2,2)))
plot(mod11)
mod12<-aov(data=actino_kw, Abundance~Lean.Mass)
summary(mod12)
with(actino_kw, par(mfrow=c(2,2)))
plot(mod12)
mod13<-aov(data=actino_kw, Abundance~Fat.Mass)
summary(mod13)
with(actino_kw, par(mfrow=c(2,2)))
plot(mod13)

mod14<-aov(data=bacteroidetes_kw, Abundance~Body.Fat)
summary(mod14)
with(bacteroidetes_kw, par(mfrow=c(2,2)))
plot(mod14)
mod15<-aov(data=bacteroidetes_kw, Abundance~NetBodyMass)
summary(mod15)
with(bacteroidetes_kw, par(mfrow=c(2,2)))
plot(mod15)
mod16<-aov(data=bacteroidetes_kw, Abundance~Lean.Mass)
summary(mod16)
with(bacteroidetes_kw, par(mfrow=c(2,2)))
plot(mod16)
mod17<-aov(data=bacteroidetes_kw, Abundance~Fat.Mass)
summary(mod17)
with(bacteroidetes_kw, par(mfrow=c(2,2)))
plot(mod17)

mod18<-aov(data=epsilon_kw, Abundance~Body.Fat)
summary(mod18)
with(epsilon_kw, par(mfrow=c(2,2)))
plot(mod18)
mod19<-aov(data=epsilon_kw, Abundance~NetBodyMass)
summary(mod19)
with(epsilon_kw, par(mfrow=c(2,2)))
plot(mod19)
mod20<-aov(data=epsilon_kw, Abundance~Lean.Mass)
summary(mod20)
with(epsilon_kw, par(mfrow=c(2,2)))
plot(mod20)
mod21<-aov(data=epsilon_kw, Abundance~Fat.Mass)
summary(mod21)
with(epsilon_kw, par(mfrow=c(2,2)))
plot(mod21)

mod22<-aov(data=tenericutes_kw, Abundance~Body.Fat)
summary(mod22)
with(tenericutes_kw, par(mfrow=c(2,2)))
plot(mod22)
mod23<-aov(data=tenericutes_kw, Abundance~NetBodyMass)
summary(mod23)
with(tenericutes_kw, par(mfrow=c(2,2)))
plot(mod23)
mod24<-aov(data=tenericutes_kw, Abundance~Lean.Mass)
summary(mod24)
with(tenericutes_kw, par(mfrow=c(2,2)))
plot(mod24)
mod25<-aov(data=tenericutes_kw, Abundance~Fat.Mass)
summary(mod25)
with(tenericutes_kw, par(mfrow=c(2,2)))
plot(mod25)

kruskal.test(data=firmicutes_kw, Abundance~Body.Fat)
#Kruskal-Wallis chi-squared =0.489, df = 2, p-value = 0.7831
kruskal.test(data=firmicutes_kw, Abundance~NetBodyMass)
#Kruskal-Wallis chi-squared = 1.6216, df = 2, p-value = 0.4445
kruskal.test(data=firmicutes_kw, Abundance~Lean.Mass)
#Kruskal-Wallis chi-squared = 3.4136, df = 2, p-value = 0.1814
kruskal.test(data=firmicutes_kw, Abundance~Fat.Mass)
#Kruskal-Wallis chi-squared =1.2877, df = 2, p-value = 0.5253
kruskal.test(data=proteo_kw, Abundance~Body.Fat)
#Kruskal-Wallis chi-squared =0.064987, df = 2, p-value = 0.968
kruskal.test(data=proteo_kw, Abundance~NetBodyMass)
#Kruskal-Wallis chi-squared =0.81357, df = 2, p-value = 0.6658
kruskal.test(data=proteo_kw, Abundance~Lean.Mass)
#Kruskal-Wallis chi-squared = 2.3159, df = 2, p-value = 0.3141
kruskal.test(data=proteo_kw, Abundance~Fat.Mass)
#Kruskal-Wallis chi-squared =0.18392, df = 2, p-value = 0.9121

#####Major phyla correlation#####

P_correlation<-phylum_melt %>%
  select( "Abundance","Phylum", "Body.Fat....", "NetBodyMasskgs", "Lean.Mass..kg.", "Fat.Mass..kg.")
P_correlation$Abundance[P_correlation$Abundance==0] <- NA
P_correlation<-P_correlation[complete.cases(P_correlation$Abundance),]
P_correlation<- P_correlation%>% 
  mutate_if(is.numeric, round, digits = 5)

actino<- P_correlation%>%
  filter(Phylum == "Actinobacteria")
bacteroidetes<- P_correlation%>%
  filter(Phylum == "Bacteroidetes")
epsilon<- P_correlation%>%
  filter(Phylum == "Epsilonbacteraeota")
firmicutes<- P_correlation%>%
  filter(Phylum == "Firmicutes")
proteo<- P_correlation%>%
  filter(Phylum == "Proteobacteria")
tenericutes<- P_correlation%>%
  filter(Phylum == "Tenericutes")

ggplot(actino, aes(x=Abundance, y=Body.Fat....))+
  geom_point() #run for every cor.test

cor.test(actino$Abundance, actino$Body.Fat...., method="spearman", use="complete.obs", exact=FALSE)
cor.test(actino$Abundance, actino$NetBodyMasskgs, method="spearman", use="complete.obs", exact=FALSE)
cor.test(actino$Abundance, actino$Lean.Mass..kg., method="spearman", use="complete.obs", exact=FALSE)
cor.test(actino$Abundance, actino$Fat.Mass..kg., method="spearman", use="complete.obs", exact=FALSE)

cor.test(bacteroidetes$Abundance, bacteroidetes$Body.Fat...., method="spearman", use="complete.obs", exact=FALSE)
cor.test(bacteroidetes$Abundance, bacteroidetes$NetBodyMasskgs, method="spearman", use="complete.obs", exact=FALSE)
cor.test(bacteroidetes$Abundance, bacteroidetes$Lean.Mass..kg., method="spearman", use="complete.obs", exact=FALSE)
cor.test(bacteroidetes$Abundance, bacteroidetes$Fat.Mass..kg., method="spearman", use="complete.obs", exact=FALSE)

cor.test(epsilon$Abundance, epsilon$Body.Fat...., method="spearman", use="complete.obs", exact=FALSE)
cor.test(epsilon$Abundance, epsilon$NetBodyMasskgs, method="spearman", use="complete.obs", exact=FALSE)
cor.test(epsilon$Abundance, epsilon$Lean.Mass..kg., method="spearman", use="complete.obs", exact=FALSE)
cor.test(epsilon$Abundance, epsilon$Fat.Mass..kg., method="spearman", use="complete.obs", exact=FALSE)
                            
cor.test(firmicutes$Abundance, firmicutes$Body.Fat...., method="spearman", use="complete.obs", exact=FALSE)
cor.test(firmicutes$Abundance, firmicutes$NetBodyMasskgs, method="spearman", use="complete.obs", exact=FALSE)
cor.test(firmicutes$Abundance, firmicutes$Lean.Mass..kg., method="spearman", use="complete.obs", exact=FALSE)
cor.test(firmicutes$Abundance, firmicutes$Fat.Mass..kg., method="spearman", use="complete.obs", exact=FALSE)

cor.test(proteo$Abundance, proteo$Body.Fat...., method="spearman", use="complete.obs", exact=FALSE)
cor.test(proteo$Abundance, proteo$NetBodyMasskgs, method="spearman", use="complete.obs", exact=FALSE)
cor.test(proteo$Abundance, proteo$Lean.Mass..kg., method="spearman", use="complete.obs", exact=FALSE)
cor.test(proteo$Abundance, proteo$Fat.Mass..kg., method="spearman", use="complete.obs", exact=FALSE)

cor.test(tenericutes$Abundance, tenericutes$Body.Fat...., method="spearman", use="complete.obs", exact=FALSE)
cor.test(tenericutes$Abundance, tenericutes$NetBodyMasskgs, method="spearman", use="complete.obs", exact=FALSE)
cor.test(tenericutes$Abundance, tenericutes$Lean.Mass..kg., method="spearman", use="complete.obs", exact=FALSE)
cor.test(tenericutes$Abundance, tenericutes$Fat.Mass..kg., method="spearman", use="complete.obs", exact=FALSE)

#### Community composition  Visualization ####

#rename phyla with < 1% abundance
p_abund$Phylum[p_abund$Abundance <= 0.01] <- "Minor"
unique(p_abund$Phylum)

write_csv(p_abund, "p_abund.csv")
## phyla order
p_abund$Phylum <- factor(p_abund$Phylum, levels = c( "Minor"     ,           "Bacteroidetes"  ,   "Actinobacteria"  ,  "Epsilonbacteraeota",
                                                        "Tenericutes"      ,   "Proteobacteria"  , "Firmicutes"      ))
spatial_plot <- ggplot(data=p_abund, aes(x=Body.Fat, y=Abundance, fill=Phylum, width=.5))

install.packages("ggthemes")
library(ggthemes)

#%body fat
p1<-spatial_plot + geom_bar(aes(),stat="identity", position="stack", width =.9) +
  scale_color_scico_d(palette = "batlow",
                      aesthetics = "fill")+
  theme_classic()+scale_y_continuous(position = "left", expand = c(0,0))+
  theme(legend.position="bottom", legend.key.size = unit(0.2, "cm"),
        legend.text = element_text(color="black", size=8),legend.title.align=0,
        legend.key.width=unit(0.2,'cm'),legend.spacing.x = unit(.2, 'cm'),
        legend.spacing.y=unit(.2, 'cm'),legend.title = element_blank(),
        axis.title=element_text(size=10),
        axis.text.x =element_text(color="black", size = 8),axis.ticks.x =element_blank(),
        axis.text.y =element_text(color="black", size = 8 ),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  xlab("% Body Fat")+ ylab("Phyla relative abundance")+ 
  scale_x_discrete(labels = c('Below\n','Median\n','Above\n'))
p1 

#Net Mass

p_abund_netmass$Phylum[p_abund_netmass$Abundance <= 0.01] <- "Minor"
p_abund_netmass$Phylum <- factor(p_abund_netmass$Phylum, levels = c( "Minor",  "Spirochaetes",     "Chloroflexi" ,    "Cyanobacteria"   ,    "Acidobacteria"  ,        "Bacteroidetes",     
 "Chlamydiae"  ,    "Planctomycetes",  "Dependentiae"   ,   "Actinobacteria" ,
"Epsilonbacteraeota"      ,     "Verrucomicrobia",    "Fusobacteria","Tenericutes", "Proteobacteria", "Firmicutes" ))

spatial_plot2 <- ggplot(data=p_abund_netmass, aes(x=NetBodyMass, y=Abundance, fill=Phylum, width=.5))


p2<-spatial_plot2 + geom_bar(aes(),stat="identity", position="stack", width =.9) +
  scale_color_scico_d(palette = "batlow",
                      aesthetics = "fill")+
  theme_classic()+
  theme(legend.position="bottom",legend.key.size = unit(0.2, "cm"),
        legend.text = element_text(color="black", size=8),legend.title.align=0,
        legend.key.width=unit(0.2,'cm'),legend.spacing.x = unit(.2, 'cm'),
        legend.spacing.y=unit(.2, 'cm'),legend.title = element_blank(),
        panel.background = element_blank(), axis.text.x =element_text(color="black", size=8 ),axis.ticks=element_blank(),axis.text.y = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  xlab("Net Body Mass (kg)")+ ylab("")+ scale_y_continuous(position = "left", expand = c(0,0)) +
  scale_x_discrete(labels = c('Below\n','Median\n','Above\n'))
p2


#Fat Mass

p_abund_fat$Phylum[p_abund_fat$Abundance <= 0.01] <- "Minor"
p_abund_fat$Phylum <- factor(p_abund_fat$Phylum, levels = c( "Minor"     ,           "Bacteroidetes"  ,   "Actinobacteria"  ,  "Epsilonbacteraeota",   "Fusobacteria" ,
                                                     "Tenericutes"      ,   "Proteobacteria"  , "Firmicutes"      ))

spatial_plot3 <- ggplot(data=p_abund_fat, aes(x=Fat.Mass, y=Abundance, fill=Phylum, width=.5))

p3<-spatial_plot3 + geom_bar(aes(),stat="identity", position="stack", width =.9) +
  scale_color_scico_d(palette = "batlow",
                      aesthetics = "fill")+
  theme_classic()+
  theme(legend.position="bottom", legend.key.size = unit(0.2, "cm"),
        legend.text = element_text(color="black", size=8),legend.title.align=0,
        legend.key.width=unit(0.2,'cm'),legend.spacing.x = unit(.2, 'cm'),
        legend.spacing.y=unit(.2, 'cm'),legend.title = element_blank(),
        axis.ticks.x =element_blank(), axis.text.x =element_text(color="black", size=8 ),
        axis.ticks.y = element_blank(), axis.text.y = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  xlab("Fat Mass (kg)")+ ylab("")+ guides(fill=guide_legend(nrow = 2)) +scale_y_continuous(position = "left", expand = c(0,0))+
  scale_x_discrete(labels = c('Below\n','Median\n','Above\n'))
p3

#Lean Mass
p_abund_lean$Phylum[p_abund_lean$Abundance <= 0.01] <- "Minor"
p_abund_lean$Phylum <- factor(p_abund_lean$Phylum, levels = c( "Minor"     ,   "Fusobacteria" ,       "Bacteroidetes"  ,   "Actinobacteria"  ,  "Epsilonbacteraeota",   
                                                             "Tenericutes"      ,   "Proteobacteria"  , "Firmicutes"      ))

spatial_plot4 <- ggplot(data=p_abund_lean, aes(x=Lean.Mass, y=Abundance, fill=Phylum, width=.5))

p4<-spatial_plot4 + geom_bar(aes(),stat="identity", position="stack", width =.9) +
  scale_color_scico_d(palette = "batlow",
                      aesthetics = "fill")+
  theme_classic()+
  theme(legend.position="bottom", legend.key.size = unit(0.2, "cm"),
        legend.text = element_text(color="black", size=8),legend.title.align=0,
        legend.key.width=unit(0.2,'cm'),legend.spacing.x = unit(.2, 'cm'),
        legend.spacing.y=unit(.2, 'cm'),legend.title = element_blank(),
        axis.ticks.x =element_blank(), axis.text.x =element_text(color="black", size=8 ),
        axis.text.y = element_text(color = "black"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  xlab("Lean Mass (kg)")+ ylab("Phyla relative abundance")+ guides(fill=guide_legend(nrow = 2)) +scale_y_continuous(position = "left", expand = c(0,0))+
  scale_x_discrete(labels = c('Below\n','Median\n','Above\n'))
p4

ggdraw() +
  draw_plot(p1, x = 0, y = .5, width = .5, height = .5) +
  draw_plot(p2, x = .5, y = .5, width = .5, height = .5) +
  draw_plot(p3, x = .5, y = 0, width = .5, height = .5) +
  draw_plot(p4, x = 0, y = 0, width = .5, height = 0.5) +
  draw_plot_label(label = c("A", "B", "C", "D"), size = 8,
                  x = c(0, 0.5, 0, 0.5), y = c(1, 1, 0.5, 0.5))

## genus

g_abund$Genus[g_abund$Abundance <= 0.01] <- "Minor"
write_csv(g_abund, "g_abund.csv")
unique(g_abund$Genus)
### put in order you want 
g_abund$Genus <- factor(g_abund$Genus, levels = c("Minor", "Actinobacillus",  "Bacteroides",                            
"Cellulosilyticum", "Clostridium sensu stricto 1" ,"Edwardsiella",       "Escherichia-Shigella" ,"Family_Enterobacteriaceae"  ,  "Family_Peptostreptococcaceae",                   
"Helicobacter" , "Mycoplasma" , "Order_Lactobacillales",  "Paenalcaligenes" , "Pseudomonas", "Romboutsia", "Streptococcus",  "Terrisporobacter" ,
"Turicibacter", "Ureaplasma" ,"Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium", "Lactobacillus",)) 
unique(g_abund_netmass$Genus)
g_abund_netmass$Genus <- factor(g_abund_netmass$Genus, levels = c( "Minor",  "Actinobacillus","Bacteroides", "Cellulosilyticum","Clostridium sensu stricto 1",  
"Edwardsiella", "Escherichia-Shigella", "Family_Enterobacteriaceae", "Family_Peptostreptococcaceae", "Helicobacter",  "Mycoplasma",
"Order_Lactobacillales", "Paenalcaligenes", "Pseudomonas", "Romboutsia", "Streptococcus", "Terrisporobacter", "Turicibacter", "Ureaplasma", "Family_Microbacteriaceae",  ))

#% body fat
spatial_plot5 <- ggplot(data=g_abund, aes(x=Body.Fat, y=Abundance, fill=Genus, width=.5))
g1<-spatial_plot5 + geom_bar(aes(),stat="identity", position="stack", width =.9) +
  theme_classic()+
  scale_color_scico_d(palette = "batlow",
                      aesthetics = "fill")+
  theme(legend.position="none",legend.key.size = unit(0.2, "cm"),
        legend.text = element_text(size=8),
        legend.key.width=unit(0.2,'cm'),legend.spacing.x = unit(0.2, 'cm'),legend.title.align=0,
        legend.spacing.y=unit(0, 'cm'),legend.title = element_blank(),
        axis.text =element_text(color="black", size=8),axis.ticks.x =element_blank(),  axis.title=element_text(color="black", size = 10 ),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  guides(fill=guide_legend(ncol=1))+
  xlab("% Body Fat")+ ylab("Genera relative abundance")+
  scale_y_continuous(position = "left", expand = c(0,0), limits=c(0, 1)) +
  scale_x_discrete(labels = c('Below\n','Median\n','Above\n'))
g1

 #Net Mass
 g_abund_netmass$Genus[g_abund_netmass$Abundance <= 0.01] <- "Minor"
 
 spatial_plot6 <- ggplot(data=g_abund_netmass, aes(x=NetBodyMass, y=Abundance, fill=Genus, width=.5))
 g2<-spatial_plot6 + geom_bar(aes(),stat="identity", position="stack", width =.9) +
   theme_classic()+scale_y_continuous(position = "left", expand = c(0,0), limits=c(0, 1))+
   scale_color_scico_d(palette = "batlow",
                       aesthetics = "fill")+
   theme(legend.position="none",legend.key.size = unit(0.2, "cm"),
         legend.text = element_text(size=8),
         legend.key.width=unit(0.2,'cm'),legend.spacing.x = unit(.2, 'cm'),legend.title.align=0,
         legend.spacing.y=unit(0, 'cm'),legend.title = element_blank(), axis.text.x =element_text(color="black", size = 8),
         axis.ticks=element_blank(), axis.text.y =element_blank(),
         panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
   guides(fill=guide_legend(ncol=1))+
   xlab("Net Body Mass (kg)")+ ylab("")+
 scale_x_discrete(labels = c('Bellow\n','Median\n','Above\n'))
 g2

 #Lean mass 
 g_abund_lean$Genus[g_abund_lean$Abundance <= 0.01] <- "Minor"
 unique(g_abund_lean$Genus)
 g_abund_lean$Genus <- factor(g_abund_lean$Genus, levels = c( "Minor", "Actinobacillus", "Bacteroides", "Cellulosilyticum", "Clostridium sensu stricto 1", "Edwardsiella",
"Escherichia-Shigella", "Family_Enterobacteriaceae", "Family_Peptostreptococcaceae", "Helicobacter", "Mycoplasma", 
"Order_Lactobacillales", "Paenalcaligenes", "Romboutsia", "Streptococcus", "Terrisporobacter", "Turicibacter", "Ureaplasma", "Family_Microbacteriaceae","Fusobacterium")) 
 
 spatial_plot7 <- ggplot(data=g_abund_lean, aes(x=Lean.Mass, y=Abundance, fill=Genus, width=.5))
 g3<-spatial_plot7 + geom_bar(aes(),stat="identity", position="stack", width =.9) +
   theme_classic()+scale_y_continuous(position = "left", expand = c(0,0), limits=c(0, 1))+
   scale_color_scico_d(palette = "batlow",
                       aesthetics = "fill")+
   theme(legend.position="none",legend.key.size = unit(0.2, "cm"),
         legend.text = element_text(size=8),
         legend.key.width=unit(0.2,'cm'),legend.spacing.x = unit(.2, 'cm'),legend.title.align=0,
         legend.spacing.y=unit(0, 'cm'),legend.title = element_blank(),  axis.text =element_text(color="black", size = 8),
         axis.text.y =element_blank(),axis.ticks.x =element_blank(),axis.ticks.y =element_blank(),
         panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
   guides(fill=guide_legend(ncol=1))+
   xlab("Lean Mass (kg)")+ ylab("")+
   scale_x_discrete(labels = c('Below\n','Median\n','Above\n'))
 g3

 #Fat Mass
 g_abund_fat$Genus[g_abund_fat$Abundance <= 0.01] <- "Minor"
 unique(g_abund_fat$Genus)
 g_abund_fat$Genus <- factor(g_abund_fat$Genus, levels = c("Minor", "Actinobacillus",  "Cellulosilyticum", "Clostridium sensu stricto 1" ,                     
"Edwardsiella", "Escherichia-Shigella","Family_Enterobacteriaceae",  "Family_Peptostreptococcaceae", "Helicobacter", "Mycoplasma",                                       
 "Order_Lactobacillales", "Romboutsia","Streptococcus", "Terrisporobacter", "Turicibacter", "Ureaplasma", 
"Ursidibacter", "Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium" ,"Bibersteinia","Family_Pasteurellaceae","Lactobacillus", "Order_Bacteroidales","Pseudomonas")) 
 
 
 spatial_plot8 <- ggplot(data=g_abund_fat, aes(x=Fat.Mass, y=Abundance, fill=Genus, width=.5))
 g4<-spatial_plot8 + geom_bar(aes(),stat="identity", position="stack", width =.9) +
   theme_classic()+scale_y_continuous(position = "left", expand = c(0,0), limits=c(0, 1))+
   scale_color_scico_d(palette = "batlow",
                       aesthetics = "fill")+
   theme(legend.position="none",legend.key.size = unit(0.2, "cm"),
         legend.text = element_text(size=8),
         legend.key.width=unit(0.2,'cm'),legend.spacing.x = unit(.2, 'cm'),legend.title.align=0,
         legend.spacing.y=unit(0, 'cm'),legend.title = element_blank(),  axis.text =element_text(color="black", size = 8),
         axis.text.y =element_blank(),axis.ticks.x =element_blank(),axis.ticks.y =element_blank(),
         panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
   guides(fill=guide_legend(ncol=1))+
   xlab("Fat Mass (kg)")+ ylab("")+
   scale_x_discrete(labels = c('Below\n','Median\n','Above\n'))
 g4
 
 write_csv(genus_melt, "genus.csv")
 
 ggdraw() +
   draw_plot(g1, x = 0, y = .5, width = .5, height = .5) +
   draw_plot(g2, x = .5, y = .5, width = .5, height = .5) +
   draw_plot(g4, x = .5, y = 0, width = .5, height = .5) +
   draw_plot(g3, x = 0, y = 0, width = .5, height = 0.5) +
   draw_plot_label(label = c("A", "B", "C", "D"), size = 8,
                   x = c(0, 0.5, 0, 0.5), y = c(1, 1, 0.5, 0.5))
