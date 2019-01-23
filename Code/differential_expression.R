setwd("/media/rmejia/ADATA-HD710/boba-bk-postsismo/rmejia/Documents/Otros_Proyectos_academicos/NOTCH_Github_pelon/NOTCH-in-Breast_Cancer")
# data given by the user
exp_mat_path<-c("/media/rmejia/ADATA-HD710/boba-bk-postsismo/rmejia/Documents/Otros_Proyectos_academicos/NOTCH_Github_pelon/Data_Castillo/joinNormDisc_colapsed_d4828.txt")
Exp_Mat <- read.table(exp_mat_path ,header=TRUE , sep="\t")

clin_table <-read.table("../Data/brca_metabric_clinical_data.tsv",header=TRUE,
                        sep = "\t", row.names = 1, stringsAsFactors = FALSE)


Num_rstatus<-data.frame(clin_table[,c(6,8,10)])
Num_rstatus[,"ER.Status"]<- gsub("Positive","4",Num_rstatus[,"ER.Status"])
Num_rstatus[,"ER.Status"] <- gsub("Negative","0",Num_rstatus[,"ER.Status"])
Num_rstatus[,"HER2.Status"]<- gsub("Positive","2",Num_rstatus[,"HER2.Status"])
Num_rstatus[,"HER2.Status"] <- gsub("Negative","0",Num_rstatus[,"HER2.Status"])
Num_rstatus[,"PR.Status"]<- gsub("Positive","1",Num_rstatus[,"PR.Status"])
Num_rstatus[,"PR.Status"] <- gsub("Negative","0",Num_rstatus[,"PR.Status"])
Num_rstatus <- Num_rstatus[!apply(Num_rstatus, 1, function(x) any(x=="")),] 
str(Num_rstatus)
Num_rstatus<- as.matrix(Num_rstatus)
mode(Num_rstatus)<-"numeric"
Sum <- Num_rstatus[,"ER.Status"] + Num_rstatus[,"HER2.Status"] + Num_rstatus[,"PR.Status"]
Num_rstatus<- cbind(Num_rstatus,Sum)
head(Num_rstatus)
Recep_Labels <- as.data.frame(Num_rstatus[,"Sum"])
colnames(Recep_Labels)<- "Labels"
Recep_Labels[,"Labels"] <- gsub(7,"EpHpPp",Recep_Labels[,"Labels"])
Recep_Labels[,"Labels"] <- gsub(6,"EpHpPn",Recep_Labels[,"Labels"])
Recep_Labels[,"Labels"] <- gsub(5,"EpHnPp",Recep_Labels[,"Labels"])
Recep_Labels[,"Labels"] <- gsub(4,"EpHnPn",Recep_Labels[,"Labels"])
Recep_Labels[,"Labels"] <- gsub(3,"EnHpPp",Recep_Labels[,"Labels"])
Recep_Labels[,"Labels"] <- gsub(2,"EnHpPn",Recep_Labels[,"Labels"])
Recep_Labels[,"Labels"] <- gsub(1,"EnHnPp",Recep_Labels[,"Labels"])
Recep_Labels[,"Labels"] <- gsub(0,"EnHnPn",Recep_Labels[,"Labels"])

head(Recep_Labels,10)

#write.table(Recep_Labels,file="/media/rmejia/ADATA-HD710/boba-bk-postsismo/rmejia/Documents/Otros_Proyectos_academicos/NOTCH_Github_pelon/Data/Recep_Labels.txt")
Dot_Labels<-read.table(file="/media/rmejia/ADATA-HD710/boba-bk-postsismo/rmejia/Documents/Otros_Proyectos_academicos/NOTCH_Github_pelon/Data/Recep_Labels.txt", stringsAsFactors = FALSE)
str(Dot_Labels)
head(Dot_Labels)
positions <- which(rownames(Dot_Labels) %in% colnames(Exp_Mat))
Recep <- data.frame(row.names = rownames(Dot_Labels)[positions],Dot_Labels[positions,])
?data.frame
colnames(Recep)<- "Labels"
head(Recep)
write.table(Recep,file = "/media/rmejia/ADATA-HD710/boba-bk-postsismo/rmejia/Documents/Otros_Proyectos_academicos/NOTCH_Github_pelon/Data/Recep_Labels_cases.txt" )
Labels_path<-c("/media/rmejia/ADATA-HD710/boba-bk-postsismo/rmejia/Documents/Otros_Proyectos_academicos/NOTCH_Github_pelon/Data/Labels_Ctrl_and_NL_Recal_separated.txt")
Labels<-read.table(Labels_path, stringsAsFactors = FALSE)
head(Labels)
Labels[143:145,]


EpHpPp <- rownames(Num_rstatus)[Num_rstatus[,"Sum"] == 7]
EpHpPn <- rownames(Num_rstatus)[Num_rstatus[,"Sum"] == 6]
EpHnPp <- rownames(Num_rstatus)[Num_rstatus[,"Sum"] == 5]
EpHnPn <- rownames(Num_rstatus)[Num_rstatus[,"Sum"] == 4]
EnHpPp <- rownames(Num_rstatus)[Num_rstatus[,"Sum"] == 3]
EnHpPn <- rownames(Num_rstatus)[Num_rstatus[,"Sum"] == 2]
EnHnPp <- rownames(Num_rstatus)[Num_rstatus[,"Sum"] == 1]
EnHnPn <- rownames(Num_rstatus)[Num_rstatus[,"Sum"] == 0]

# differential expression analysis

# differential expression results
write.table()
