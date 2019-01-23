#source("http://bioconductor.org/biocLite.R")
#biocLite("illuminaHumanv3.db")
#biocLite("hgu133plus2.db")
#install.packages("gProfileR")
#biocLite("org.Hs.eg.db")
library(org.Hs.eg.db)
library(hgu133plus2.db)
library(annotate)
library(limma)
library(gProfileR)
############################################################
# Data by the user                       ###################
############################################################
exp_mat_path<-c("/media/rmejia/ADATA-HD710/boba-bk-postsismo/rmejia/Documents/Otros_Proyectos_academicos/NOTCH_Github_pelon/Data_Castillo/joinNormDisc_colapsed_d4828.txt")
Labels_path<-c("/media/rmejia/ADATA-HD710/boba-bk-postsismo/rmejia/Documents/Otros_Proyectos_academicos/NOTCH_Github_pelon/Data/Recep_Labels_Cases_And_Controls.txt")
results_path<-c("/media/rmejia/ADATA-HD710/boba-bk-postsismo/rmejia/Documents/Otros_Proyectos_academicos/NOTCH_Github_pelon/Results/")
# Don forget rename the colnames(design)

###############################################################
### Loading required libraries                  ###############
###############################################################
source("http://bioconductor.org/biocLite.R")
if (!require("org.Hs.eg.db")) {
  biocLite("org.Hs.eg.db", ask =FALSE)
  library(org.Hs.eg.db)
}
if (!require("hgu133plus2.db")){
  biocLite("hgu133plus2.db", ask =FALSE)
  library(hgu133plus2.db)
}
if (!require("annotate")){
  biocLite("annotate", ask =FALSE)
  library(annotate)
}
if (!require("limma")){
  biocLite("limma", ask =FALSE)
  library(limma)
}
#########################################################
####      Let's start                               #####
#########################################################

# Load the data frame (expression matrix)
#Exp_Mat <- read.table(exp_mat_path ,header=TRUE , sep="\t")

#Load the annotation data to extract subtypes
Labels<-read.table(Labels_path)


# convert into factors
samples<-Labels$Labels
# check factors have been assigned
samples
# set up the experimental design
design <- model.matrix(~0 + samples)
head(design)
colnames(design) <- c("Control","EnHnPn","EnHnPp","EnHpPn","EpHnPn",
                      "EpHnPp","EpHpPn","EpHpPp")

# fit the linear model to the your expression set "eset"
eset<-Exp_Mat
fit <- lmFit(eset, design)

# set up a contrast matrix to compare tissues v cell line
contrast.matrix <- makeContrasts(EnHnPn_Control = EnHnPn - Control,
                                 EnHnPp_Control = EnHnPp - Control,
                                 EnHpPn_Control = EnHpPn - Control,
                                 EpHnPn_Control = EpHnPn - Control,
                                 EpHnPp__Control = EpHnPp -  Control,  
                                 levels=design)
# check the contrast matrix
contrast.matrix
# Now the contrast matrix is combined with the per-probeset linear model fit.
huvec_fits <- contrasts.fit(fit, contrast.matrix)
huvec_ebFit <- eBayes(huvec_fits)

#####################################
##DGEs with  Affy ids no anotation###
#####################################
results_only_affys<-list(length(colnames(contrast.matrix)))
for(i in 1:length(colnames(contrast.matrix))){
  DGEs<-topTable(huvec_ebFit, coef=colnames(contrast.matrix)[i], number=10000, p.value = 10^-3, lfc =1)
  # Saving the matrices of DGE by subtype
  write.table(DGEs, paste(results_path,"DGE_MET_NormDisc_",colnames(contrast.matrix)[i],".txt",sep=""), sep="\t", quote=FALSE)
  # Saving with no _at
  write.table(DGEs[!grepl("ILMN_",rownames(DGEs)),], paste(results_path,"DGE_MET_NormDisc_",colnames(contrast.matrix)[i],"_ID_no_ILMN_.txt",sep=""), sep="\t", quote=FALSE)
  results_only_affys[[i]]<-DGEs
}
save(results_only_affys,file=paste(results_path,c("List_of_DGE_MET_NormDisc_with_at.RData"),sep=""))

rm(results_only_affys)

results_only_affys<-list(length(colnames(contrast.matrix)))
for(i in 1:length(colnames(contrast.matrix))){
  DGEs<-topTable(huvec_ebFit, coef=colnames(contrast.matrix)[i], number=10000, p.value = 10^-3, lfc =1)
  # filtering the ILMN_
  results_only_affys[[i]]<-DGEs[!grepl("ILMN_",rownames(DGEs)),]
}

save(results_only_affys,file=paste(results_path,c("List_of_DGE_MET_NormDisc_NO_ILMN_.RData"),sep=""))
