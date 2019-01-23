setwd("/media/rmejia/ADATA-HD710/boba-bk-postsismo/rmejia/Documents/Otros_Proyectos_academicos/NOTCH_Github_pelon/NOTCH-in-Breast_Cancer")
# data given by the user

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

EpHpPp <- rownames(Num_rstatus)[Num_rstatus[,"Sum"] == 7]
EpHpPn <- rownames(Num_rstatus)[Num_rstatus[,"Sum"] == 6]
EpHnPp <- rownames(Num_rstatus)[Num_rstatus[,"Sum"] == 5]
EpHnPn <- rownames(Num_rstatus)[Num_rstatus[,"Sum"] == 4]
EnHpPp <- rownames(Num_rstatus)[Num_rstatus[,"Sum"] == 3]
EnHpPn <- rownames(Num_rstatus)[Num_rstatus[,"Sum"] == 2]
EnHnPp <- rownames(Num_rstatus)[Num_rstatus[,"Sum"] == 1]
EnHnPn <- rownames(Num_rstatus)[Num_rstatus[,"Sum"] == 0]


head(Num_rstatus, 30)

Num_rstatus[,"ER.Status"] <- gsub("","NA",Num_rstatus[,"ER.Status"])
Num_rstatus[,"HER2.Status"]<- gsub("","NA",Num_rstatus[,"HER2.Status"])
Num_rstatus[,"PR.Status"] <- gsub("",NA,Num_rstatus[,"PR.Status"])

head(Num_rstatus,30)
is.na(Num_rstatus)
<- Num_rstatus
Num_rstatus <- as.numeric(Num_rstatus)
is.na(Num_rstatus)
Num_rstatus[,"HER2.Status"]==""
Num_rstatus == "" 
is.na(clin_table,)
head(Num_rstatus, 20)
head(clin_table[,c(6,8,10)],30)
?gsub
str(Num_rstatus)
clin_table[,]
class(clin_table)
str(clin_table)
head(clin_table)
matrix <-read.table()
annotation <- read.table()

# differential expression analysis

# differential expression results
write.table()
