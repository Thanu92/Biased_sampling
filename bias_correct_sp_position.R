#Sep 19, 2023

#Sister species table
getwd()
library(foreach)
#install.packages("phytools")
library(phytools)
library(maps)
library(phangorn)
library(phylotools)
#Read fasta files of placed co1 on the bb tree
bias_co1_1=read.fasta(file = "CO1_df1.fasta")
bias_co1_2=read.fasta(file = "CO1_df2.fasta")
bias_co1_3=read.fasta(file = "CO1_df3.fasta")
bias_co1_4=read.fasta(file = "CO1_df4.fasta")
bias_co1_5=read.fasta(file = "CO1_df5.fasta")
bias_co1_6=read.fasta(file = "CO1_df6.fasta")
bias_co1_7=read.fasta(file = "CO1_df7.fasta")
bias_co1_8=read.fasta(file = "CO1_df8.fasta")
bias_co1_9=read.fasta(file = "CO1_df9.fasta")
bias_co1_10=read.fasta(file = "CO1_df10.fasta")

#read the refernce tree
ReferenceTree <- read.tree("RAxML_bestTree.FR100_new")

#To compare the 10 data-driven samples with random and stratified samples
# parseTreeComp <- read.tree("RAxML_bestTree.FR100_new")
fullyRandom20_1 <- read.tree("20_1_new")
fullyRandom40_1 <- read.tree("40_1_new")
fullyRandom60_1 <- read.tree("60_1_new")
fullyRandom80_1 <- read.tree("80_1_new")

fullyStrati20_1 <- read.tree("S20_1_new")
fullyStrati40_1 <- read.tree("S40_1_new")
fullyStrati60_1 <- read.tree("S60_1_new")
fullyStrati80_1 <- read.tree("S80_1_new")

#Read the 10 bias trees generated using bias samples
B1_tree<-read.tree("B1_new")
B2_tree<-read.tree("B2_new")
B3_tree<-read.tree("B3_new")
B4_tree<-read.tree("B4_new")
B5_tree<-read.tree("B5_new")
B6_tree<-read.tree("B6_new")
B7_tree<-read.tree("B7_new")
B8_tree<-read.tree("B8_new")
B9_tree<-read.tree("B9_new")
B10_tree<-read.tree("B10_new")
#------
bias_co1_1species <- bias_co1_1$seq.name
bias_co1_2species <- bias_co1_2$seq.name
bias_co1_3species <- bias_co1_3$seq.name
bias_co1_4species <- bias_co1_4$seq.name
bias_co1_5species <- bias_co1_5$seq.name
bias_co1_6species <- bias_co1_6$seq.name
bias_co1_7species <- bias_co1_7$seq.name
bias_co1_8species <- bias_co1_8$seq.name
bias_co1_9species <- bias_co1_9$seq.name
bias_co1_10species <- bias_co1_10$seq.name

RefTree_biasco1_1_List <- foreach(i=1:length(bias_co1_1species)) %do% getSisters(ReferenceTree,bias_co1_1species[i],mode="label")
biasco1_1_List <- foreach(i=1:length(bias_co1_1species)) %do% getSisters(B1_tree,bias_co1_1species[i],mode="label")

# RefTree_biasco1_1_List  <- foreach(i=1:length(bias_co1_1species)) %do% sister(ReferenceTree,bias_co1_1species[i],type="terminal",label=T)
# biasco1_1_List <- foreach(i=1:length(bias_co1_1species)) %do% sister(B1_tree,bias_co1_1species[i],type="terminal",label=T)
#get the number of species have similar sisters by comparing with reference tree
RefTree_biasco1_1_List [(RefTree_biasco1_1_List  %in% biasco1_1_List )] #158
#get the number of species doesn't have similar sisters by comparing with reference tree
RefTree_biasco1_1_List [!(RefTree_biasco1_1_List  %in% biasco1_1_List )] #3733

RefTree_biasco1_2_List <- foreach(i=1:length(bias_co1_2species)) %do% getSisters(ReferenceTree,bias_co1_2species[i],mode="label")
biasco1_2_List <- foreach(i=1:length(bias_co1_2species)) %do% getSisters(B2_tree,bias_co1_2species[i],mode="label")

#get the number of species have similar sisters by comparing with reference tree
RefTree_biasco1_2_List [(RefTree_biasco1_2_List  %in% biasco1_2_List )] #243
#get the number of species doesn't have similar sisters by comparing with reference tree
RefTree_biasco1_2_List [!(RefTree_biasco1_2_List  %in% biasco1_2_List )] #3349

RefTree_biasco1_3_List <- foreach(i=1:length(bias_co1_3species)) %do% getSisters(ReferenceTree,bias_co1_3species[i],mode="label")
biasco1_3_List <- foreach(i=1:length(bias_co1_3species)) %do% getSisters(B3_tree,bias_co1_3species[i],mode="label")

#get the number of species have similar sisters by comparing with reference tree
RefTree_biasco1_3_List [(RefTree_biasco1_3_List  %in% biasco1_3_List )] #371
#get the number of species doesn't have similar sisters by comparing with reference tree
RefTree_biasco1_3_List [!(RefTree_biasco1_3_List  %in% biasco1_3_List )] #2888

RefTree_biasco1_4_List <- foreach(i=1:length(bias_co1_4species)) %do% getSisters(ReferenceTree,bias_co1_4species[i],mode="label")
biasco1_4_List <- foreach(i=1:length(bias_co1_4species)) %do% getSisters(B4_tree,bias_co1_4species[i],mode="label")

#get the number of species have similar sisters by comparing with reference tree
RefTree_biasco1_4_List [(RefTree_biasco1_4_List  %in% biasco1_4_List )] #330
#get the number of species doesn't have similar sisters by comparing with reference tree
RefTree_biasco1_4_List [!(RefTree_biasco1_4_List  %in% biasco1_4_List )] #3221

RefTree_biasco1_5_List <- foreach(i=1:length(bias_co1_5species)) %do% getSisters(ReferenceTree,bias_co1_5species[i],mode="label")
biasco1_5_List <- foreach(i=1:length(bias_co1_5species)) %do% getSisters(B5_tree,bias_co1_5species[i],mode="label")

#get the number of species have similar sisters by comparing with reference tree
RefTree_biasco1_5_List [(RefTree_biasco1_5_List  %in% biasco1_5_List )] #233
#get the number of species doesn't have similar sisters by comparing with reference tree
RefTree_biasco1_5_List [!(RefTree_biasco1_5_List  %in% biasco1_5_List )] #3422

RefTree_biasco1_6_List <- foreach(i=1:length(bias_co1_6species)) %do% getSisters(ReferenceTree,bias_co1_6species[i],mode="label")
biasco1_6_List <- foreach(i=1:length(bias_co1_6species)) %do% getSisters(B6_tree,bias_co1_6species[i],mode="label")

#get the number of species have similar sisters by comparing with reference tree
RefTree_biasco1_6_List [(RefTree_biasco1_6_List  %in% biasco1_6_List )] #273
#get the number of species doesn't have similar sisters by comparing with reference tree
RefTree_biasco1_6_List [!(RefTree_biasco1_6_List  %in% biasco1_6_List )] #3601

RefTree_biasco1_7_List <- foreach(i=1:length(bias_co1_7species)) %do% getSisters(ReferenceTree,bias_co1_7species[i],mode="label")
biasco1_7_List <- foreach(i=1:length(bias_co1_7species)) %do% getSisters(B7_tree,bias_co1_7species[i],mode="label")

#get the number of species have similar sisters by comparing with reference tree
RefTree_biasco1_7_List [(RefTree_biasco1_7_List  %in% biasco1_7_List )] #297
#get the number of species doesn't have similar sisters by comparing with reference tree
RefTree_biasco1_7_List [!(RefTree_biasco1_7_List  %in% biasco1_7_List )] #3600

RefTree_biasco1_8_List <- foreach(i=1:length(bias_co1_8species)) %do% getSisters(ReferenceTree,bias_co1_8species[i],mode="label")
biasco1_8_List <- foreach(i=1:length(bias_co1_8species)) %do% getSisters(B8_tree,bias_co1_8species[i],mode="label")

#get the number of species have similar sisters by comparing with reference tree
RefTree_biasco1_8_List [(RefTree_biasco1_8_List  %in% biasco1_8_List )] #345
#get the number of species doesn't have similar sisters by comparing with reference tree
RefTree_biasco1_8_List [!(RefTree_biasco1_8_List  %in% biasco1_8_List )] #3250

RefTree_biasco1_9_List <- foreach(i=1:length(bias_co1_9species)) %do% getSisters(ReferenceTree,bias_co1_9species[i],mode="label")
biasco1_9_List <- foreach(i=1:length(bias_co1_9species)) %do% getSisters(B9_tree,bias_co1_9species[i],mode="label")

#get the number of species have similar sisters by comparing with reference tree
RefTree_biasco1_9_List [(RefTree_biasco1_9_List  %in% biasco1_9_List )] #389
#get the number of species doesn't have similar sisters by comparing with reference tree
RefTree_biasco1_9_List [!(RefTree_biasco1_9_List  %in% biasco1_9_List )] #3208

RefTree_biasco1_10_List <- foreach(i=1:length(bias_co1_10species)) %do% getSisters(ReferenceTree,bias_co1_10species[i],mode="label")
biasco1_10_List <- foreach(i=1:length(bias_co1_10species)) %do% getSisters(B10_tree,bias_co1_10species[i],mode="label")

#get the number of species have similar sisters by comparing with reference tree
RefTree_biasco1_10_List [(RefTree_biasco1_10_List  %in% biasco1_10_List )] #267
#get the number of species doesn't have similar sisters by comparing with reference tree
RefTree_biasco1_10_List [!(RefTree_biasco1_10_List  %in% biasco1_10_List )] #3505

#To compare the 10 data-driven samples with random and stratified samples
#For fully random 20_1
#read the refernce tree
ReferenceTree <- read.tree("20_1_new")

RefTree_biasco1_1_List <- foreach(i=1:length(bias_co1_1species)) %do% getSisters(ReferenceTree,bias_co1_1species[i],mode="label")
biasco1_1_List <- foreach(i=1:length(bias_co1_1species)) %do% getSisters(B1_tree,bias_co1_1species[i],mode="label")

#get the number of species have similar sisters by comparing with reference tree
RefTree_biasco1_1_List [(RefTree_biasco1_1_List  %in% biasco1_1_List )] #61
#get the number of species doesn't have similar sisters by comparing with reference tree
RefTree_biasco1_1_List [!(RefTree_biasco1_1_List  %in% biasco1_1_List )] #3830

RefTree_biasco1_2_List <- foreach(i=1:length(bias_co1_2species)) %do% getSisters(ReferenceTree,bias_co1_2species[i],mode="label")
biasco1_2_List <- foreach(i=1:length(bias_co1_2species)) %do% getSisters(B2_tree,bias_co1_2species[i],mode="label")

#get the number of species have similar sisters by comparing with reference tree
RefTree_biasco1_2_List [(RefTree_biasco1_2_List  %in% biasco1_2_List )] #79
#get the number of species doesn't have similar sisters by comparing with reference tree
RefTree_biasco1_2_List [!(RefTree_biasco1_2_List  %in% biasco1_2_List )] #3513

RefTree_biasco1_3_List <- foreach(i=1:length(bias_co1_3species)) %do% getSisters(ReferenceTree,bias_co1_3species[i],mode="label")
biasco1_3_List <- foreach(i=1:length(bias_co1_3species)) %do% getSisters(B3_tree,bias_co1_3species[i],mode="label")

#get the number of species have similar sisters by comparing with reference tree
RefTree_biasco1_3_List [(RefTree_biasco1_3_List  %in% biasco1_3_List )] #127
#get the number of species doesn't have similar sisters by comparing with reference tree
RefTree_biasco1_3_List [!(RefTree_biasco1_3_List  %in% biasco1_3_List )] #3132

RefTree_biasco1_4_List <- foreach(i=1:length(bias_co1_4species)) %do% getSisters(ReferenceTree,bias_co1_4species[i],mode="label")
biasco1_4_List <- foreach(i=1:length(bias_co1_4species)) %do% getSisters(B4_tree,bias_co1_4species[i],mode="label")

#get the number of species have similar sisters by comparing with reference tree
RefTree_biasco1_4_List [(RefTree_biasco1_4_List  %in% biasco1_4_List )] #119
#get the number of species doesn't have similar sisters by comparing with reference tree
RefTree_biasco1_4_List [!(RefTree_biasco1_4_List  %in% biasco1_4_List )] #3432

RefTree_biasco1_5_List <- foreach(i=1:length(bias_co1_5species)) %do% getSisters(ReferenceTree,bias_co1_5species[i],mode="label")
biasco1_5_List <- foreach(i=1:length(bias_co1_5species)) %do% getSisters(B5_tree,bias_co1_5species[i],mode="label")

#get the number of species have similar sisters by comparing with reference tree
RefTree_biasco1_5_List [(RefTree_biasco1_5_List  %in% biasco1_5_List )] #103
#get the number of species doesn't have similar sisters by comparing with reference tree
RefTree_biasco1_5_List [!(RefTree_biasco1_5_List  %in% biasco1_5_List )] #3552

RefTree_biasco1_6_List <- foreach(i=1:length(bias_co1_6species)) %do% getSisters(ReferenceTree,bias_co1_6species[i],mode="label")
biasco1_6_List <- foreach(i=1:length(bias_co1_6species)) %do% getSisters(B6_tree,bias_co1_6species[i],mode="label")

#get the number of species have similar sisters by comparing with reference tree
RefTree_biasco1_6_List [(RefTree_biasco1_6_List  %in% biasco1_6_List )] #100
#get the number of species doesn't have similar sisters by comparing with reference tree
RefTree_biasco1_6_List [!(RefTree_biasco1_6_List  %in% biasco1_6_List )] #3774

RefTree_biasco1_7_List <- foreach(i=1:length(bias_co1_7species)) %do% getSisters(ReferenceTree,bias_co1_7species[i],mode="label")
biasco1_7_List <- foreach(i=1:length(bias_co1_7species)) %do% getSisters(B7_tree,bias_co1_7species[i],mode="label")

#get the number of species have similar sisters by comparing with reference tree
RefTree_biasco1_7_List [(RefTree_biasco1_7_List  %in% biasco1_7_List )] #140
#get the number of species doesn't have similar sisters by comparing with reference tree
RefTree_biasco1_7_List [!(RefTree_biasco1_7_List  %in% biasco1_7_List )] #3757

RefTree_biasco1_8_List <- foreach(i=1:length(bias_co1_8species)) %do% getSisters(ReferenceTree,bias_co1_8species[i],mode="label")
biasco1_8_List <- foreach(i=1:length(bias_co1_8species)) %do% getSisters(B8_tree,bias_co1_8species[i],mode="label")

#get the number of species have similar sisters by comparing with reference tree
RefTree_biasco1_8_List [(RefTree_biasco1_8_List  %in% biasco1_8_List )] #160
#get the number of species doesn't have similar sisters by comparing with reference tree
RefTree_biasco1_8_List [!(RefTree_biasco1_8_List  %in% biasco1_8_List )] #3435

RefTree_biasco1_9_List <- foreach(i=1:length(bias_co1_9species)) %do% getSisters(ReferenceTree,bias_co1_9species[i],mode="label")
biasco1_9_List <- foreach(i=1:length(bias_co1_9species)) %do% getSisters(B9_tree,bias_co1_9species[i],mode="label")

#get the number of species have similar sisters by comparing with reference tree
RefTree_biasco1_9_List [(RefTree_biasco1_9_List  %in% biasco1_9_List )] #169
#get the number of species doesn't have similar sisters by comparing with reference tree
RefTree_biasco1_9_List [!(RefTree_biasco1_9_List  %in% biasco1_9_List )] #3428

RefTree_biasco1_10_List <- foreach(i=1:length(bias_co1_10species)) %do% getSisters(ReferenceTree,bias_co1_10species[i],mode="label")
biasco1_10_List <- foreach(i=1:length(bias_co1_10species)) %do% getSisters(B10_tree,bias_co1_10species[i],mode="label")

#get the number of species have similar sisters by comparing with reference tree
RefTree_biasco1_10_List [(RefTree_biasco1_10_List  %in% biasco1_10_List )] #96
#get the number of species doesn't have similar sisters by comparing with reference tree
RefTree_biasco1_10_List [!(RefTree_biasco1_10_List  %in% biasco1_10_List )] #3676

#To compare the 10 data-driven samples with random and stratified samples
#For fully random 20_1
#read the refernce tree
ReferenceTree <- read.tree("40_1_new")

RefTree_biasco1_1_List <- foreach(i=1:length(bias_co1_1species)) %do% getSisters(ReferenceTree,bias_co1_1species[i],mode="label")
biasco1_1_List <- foreach(i=1:length(bias_co1_1species)) %do% getSisters(B1_tree,bias_co1_1species[i],mode="label")

#get the number of species have similar sisters by comparing with reference tree
RefTree_biasco1_1_List [(RefTree_biasco1_1_List  %in% biasco1_1_List )] #104
#get the number of species doesn't have similar sisters by comparing with reference tree
RefTree_biasco1_1_List [!(RefTree_biasco1_1_List  %in% biasco1_1_List )] #3787

RefTree_biasco1_2_List <- foreach(i=1:length(bias_co1_2species)) %do% getSisters(ReferenceTree,bias_co1_2species[i],mode="label")
biasco1_2_List <- foreach(i=1:length(bias_co1_2species)) %do% getSisters(B2_tree,bias_co1_2species[i],mode="label")

#get the number of species have similar sisters by comparing with reference tree
RefTree_biasco1_2_List [(RefTree_biasco1_2_List  %in% biasco1_2_List )] #186
#get the number of species doesn't have similar sisters by comparing with reference tree
RefTree_biasco1_2_List [!(RefTree_biasco1_2_List  %in% biasco1_2_List )] #3406

RefTree_biasco1_3_List <- foreach(i=1:length(bias_co1_3species)) %do% getSisters(ReferenceTree,bias_co1_3species[i],mode="label")
biasco1_3_List <- foreach(i=1:length(bias_co1_3species)) %do% getSisters(B3_tree,bias_co1_3species[i],mode="label")

#get the number of species have similar sisters by comparing with reference tree
RefTree_biasco1_3_List [(RefTree_biasco1_3_List  %in% biasco1_3_List )] #236
#get the number of species doesn't have similar sisters by comparing with reference tree
RefTree_biasco1_3_List [!(RefTree_biasco1_3_List  %in% biasco1_3_List )] #3023

RefTree_biasco1_4_List <- foreach(i=1:length(bias_co1_4species)) %do% getSisters(ReferenceTree,bias_co1_4species[i],mode="label")
biasco1_4_List <- foreach(i=1:length(bias_co1_4species)) %do% getSisters(B4_tree,bias_co1_4species[i],mode="label")

#get the number of species have similar sisters by comparing with reference tree
RefTree_biasco1_4_List [(RefTree_biasco1_4_List  %in% biasco1_4_List )] #210
#get the number of species doesn't have similar sisters by comparing with reference tree
RefTree_biasco1_4_List [!(RefTree_biasco1_4_List  %in% biasco1_4_List )] #3341

RefTree_biasco1_5_List <- foreach(i=1:length(bias_co1_5species)) %do% getSisters(ReferenceTree,bias_co1_5species[i],mode="label")
biasco1_5_List <- foreach(i=1:length(bias_co1_5species)) %do% getSisters(B5_tree,bias_co1_5species[i],mode="label")

#get the number of species have similar sisters by comparing with reference tree
RefTree_biasco1_5_List [(RefTree_biasco1_5_List  %in% biasco1_5_List )] #174
#get the number of species doesn't have similar sisters by comparing with reference tree
RefTree_biasco1_5_List [!(RefTree_biasco1_5_List  %in% biasco1_5_List )] #3481

RefTree_biasco1_6_List <- foreach(i=1:length(bias_co1_6species)) %do% getSisters(ReferenceTree,bias_co1_6species[i],mode="label")
biasco1_6_List <- foreach(i=1:length(bias_co1_6species)) %do% getSisters(B6_tree,bias_co1_6species[i],mode="label")

#get the number of species have similar sisters by comparing with reference tree
RefTree_biasco1_6_List [(RefTree_biasco1_6_List  %in% biasco1_6_List )] #218
#get the number of species doesn't have similar sisters by comparing with reference tree
RefTree_biasco1_6_List [!(RefTree_biasco1_6_List  %in% biasco1_6_List )] #3656

RefTree_biasco1_7_List <- foreach(i=1:length(bias_co1_7species)) %do% getSisters(ReferenceTree,bias_co1_7species[i],mode="label")
biasco1_7_List <- foreach(i=1:length(bias_co1_7species)) %do% getSisters(B7_tree,bias_co1_7species[i],mode="label")

#get the number of species have similar sisters by comparing with reference tree
RefTree_biasco1_7_List [(RefTree_biasco1_7_List  %in% biasco1_7_List )] #198
#get the number of species doesn't have similar sisters by comparing with reference tree
RefTree_biasco1_7_List [!(RefTree_biasco1_7_List  %in% biasco1_7_List )] #3699

RefTree_biasco1_8_List <- foreach(i=1:length(bias_co1_8species)) %do% getSisters(ReferenceTree,bias_co1_8species[i],mode="label")
biasco1_8_List <- foreach(i=1:length(bias_co1_8species)) %do% getSisters(B8_tree,bias_co1_8species[i],mode="label")

#get the number of species have similar sisters by comparing with reference tree
RefTree_biasco1_8_List [(RefTree_biasco1_8_List  %in% biasco1_8_List )] #273
#get the number of species doesn't have similar sisters by comparing with reference tree
RefTree_biasco1_8_List [!(RefTree_biasco1_8_List  %in% biasco1_8_List )] #3322

RefTree_biasco1_9_List <- foreach(i=1:length(bias_co1_9species)) %do% getSisters(ReferenceTree,bias_co1_9species[i],mode="label")
biasco1_9_List <- foreach(i=1:length(bias_co1_9species)) %do% getSisters(B9_tree,bias_co1_9species[i],mode="label")

#get the number of species have similar sisters by comparing with reference tree
RefTree_biasco1_9_List [(RefTree_biasco1_9_List  %in% biasco1_9_List )] #266
#get the number of species doesn't have similar sisters by comparing with reference tree
RefTree_biasco1_9_List [!(RefTree_biasco1_9_List  %in% biasco1_9_List )] #3331

RefTree_biasco1_10_List <- foreach(i=1:length(bias_co1_10species)) %do% getSisters(ReferenceTree,bias_co1_10species[i],mode="label")
biasco1_10_List <- foreach(i=1:length(bias_co1_10species)) %do% getSisters(B10_tree,bias_co1_10species[i],mode="label")

#get the number of species have similar sisters by comparing with reference tree
RefTree_biasco1_10_List [(RefTree_biasco1_10_List  %in% biasco1_10_List )] #191
#get the number of species doesn't have similar sisters by comparing with reference tree
RefTree_biasco1_10_List [!(RefTree_biasco1_10_List  %in% biasco1_10_List )] #3581

#To compare the 10 data-driven samples with random and stratified samples
#For fully random 20_1
#read the refernce tree
ReferenceTree <- read.tree("60_1_new")

RefTree_biasco1_1_List <- foreach(i=1:length(bias_co1_1species)) %do% getSisters(ReferenceTree,bias_co1_1species[i],mode="label")
biasco1_1_List <- foreach(i=1:length(bias_co1_1species)) %do% getSisters(B1_tree,bias_co1_1species[i],mode="label")

#get the number of species have similar sisters by comparing with reference tree
RefTree_biasco1_1_List [(RefTree_biasco1_1_List  %in% biasco1_1_List )] #150
#get the number of species doesn't have similar sisters by comparing with reference tree
RefTree_biasco1_1_List [!(RefTree_biasco1_1_List  %in% biasco1_1_List )] #3741

RefTree_biasco1_2_List <- foreach(i=1:length(bias_co1_2species)) %do% getSisters(ReferenceTree,bias_co1_2species[i],mode="label")
biasco1_2_List <- foreach(i=1:length(bias_co1_2species)) %do% getSisters(B2_tree,bias_co1_2species[i],mode="label")

#get the number of species have similar sisters by comparing with reference tree
RefTree_biasco1_2_List [(RefTree_biasco1_2_List  %in% biasco1_2_List )] #223
#get the number of species doesn't have similar sisters by comparing with reference tree
RefTree_biasco1_2_List [!(RefTree_biasco1_2_List  %in% biasco1_2_List )] #3369

RefTree_biasco1_3_List <- foreach(i=1:length(bias_co1_3species)) %do% getSisters(ReferenceTree,bias_co1_3species[i],mode="label")
biasco1_3_List <- foreach(i=1:length(bias_co1_3species)) %do% getSisters(B3_tree,bias_co1_3species[i],mode="label")

#get the number of species have similar sisters by comparing with reference tree
RefTree_biasco1_3_List [(RefTree_biasco1_3_List  %in% biasco1_3_List )] #318
#get the number of species doesn't have similar sisters by comparing with reference tree
RefTree_biasco1_3_List [!(RefTree_biasco1_3_List  %in% biasco1_3_List )] #2941

RefTree_biasco1_4_List <- foreach(i=1:length(bias_co1_4species)) %do% getSisters(ReferenceTree,bias_co1_4species[i],mode="label")
biasco1_4_List <- foreach(i=1:length(bias_co1_4species)) %do% getSisters(B4_tree,bias_co1_4species[i],mode="label")

#get the number of species have similar sisters by comparing with reference tree
RefTree_biasco1_4_List [(RefTree_biasco1_4_List  %in% biasco1_4_List )] #280
#get the number of species doesn't have similar sisters by comparing with reference tree
RefTree_biasco1_4_List [!(RefTree_biasco1_4_List  %in% biasco1_4_List )] #3271

RefTree_biasco1_5_List <- foreach(i=1:length(bias_co1_5species)) %do% getSisters(ReferenceTree,bias_co1_5species[i],mode="label")
biasco1_5_List <- foreach(i=1:length(bias_co1_5species)) %do% getSisters(B5_tree,bias_co1_5species[i],mode="label")

#get the number of species have similar sisters by comparing with reference tree
RefTree_biasco1_5_List [(RefTree_biasco1_5_List  %in% biasco1_5_List )] #230
#get the number of species doesn't have similar sisters by comparing with reference tree
RefTree_biasco1_5_List [!(RefTree_biasco1_5_List  %in% biasco1_5_List )] #3425

RefTree_biasco1_6_List <- foreach(i=1:length(bias_co1_6species)) %do% getSisters(ReferenceTree,bias_co1_6species[i],mode="label")
biasco1_6_List <- foreach(i=1:length(bias_co1_6species)) %do% getSisters(B6_tree,bias_co1_6species[i],mode="label")

#get the number of species have similar sisters by comparing with reference tree
RefTree_biasco1_6_List [(RefTree_biasco1_6_List  %in% biasco1_6_List )] #255
#get the number of species doesn't have similar sisters by comparing with reference tree
RefTree_biasco1_6_List [!(RefTree_biasco1_6_List  %in% biasco1_6_List )] #3619

RefTree_biasco1_7_List <- foreach(i=1:length(bias_co1_7species)) %do% getSisters(ReferenceTree,bias_co1_7species[i],mode="label")
biasco1_7_List <- foreach(i=1:length(bias_co1_7species)) %do% getSisters(B7_tree,bias_co1_7species[i],mode="label")

#get the number of species have similar sisters by comparing with reference tree
RefTree_biasco1_7_List [(RefTree_biasco1_7_List  %in% biasco1_7_List )] #262
#get the number of species doesn't have similar sisters by comparing with reference tree
RefTree_biasco1_7_List [!(RefTree_biasco1_7_List  %in% biasco1_7_List )] #3635

RefTree_biasco1_8_List <- foreach(i=1:length(bias_co1_8species)) %do% getSisters(ReferenceTree,bias_co1_8species[i],mode="label")
biasco1_8_List <- foreach(i=1:length(bias_co1_8species)) %do% getSisters(B8_tree,bias_co1_8species[i],mode="label")

#get the number of species have similar sisters by comparing with reference tree
RefTree_biasco1_8_List [(RefTree_biasco1_8_List  %in% biasco1_8_List )] #316
#get the number of species doesn't have similar sisters by comparing with reference tree
RefTree_biasco1_8_List [!(RefTree_biasco1_8_List  %in% biasco1_8_List )] #3279

RefTree_biasco1_9_List <- foreach(i=1:length(bias_co1_9species)) %do% getSisters(ReferenceTree,bias_co1_9species[i],mode="label")
biasco1_9_List <- foreach(i=1:length(bias_co1_9species)) %do% getSisters(B9_tree,bias_co1_9species[i],mode="label")

#get the number of species have similar sisters by comparing with reference tree
RefTree_biasco1_9_List [(RefTree_biasco1_9_List  %in% biasco1_9_List )] #340
#get the number of species doesn't have similar sisters by comparing with reference tree
RefTree_biasco1_9_List [!(RefTree_biasco1_9_List  %in% biasco1_9_List )] #3257

RefTree_biasco1_10_List <- foreach(i=1:length(bias_co1_10species)) %do% getSisters(ReferenceTree,bias_co1_10species[i],mode="label")
biasco1_10_List <- foreach(i=1:length(bias_co1_10species)) %do% getSisters(B10_tree,bias_co1_10species[i],mode="label")

#get the number of species have similar sisters by comparing with reference tree
RefTree_biasco1_10_List [(RefTree_biasco1_10_List  %in% biasco1_10_List )] #232
#get the number of species doesn't have similar sisters by comparing with reference tree
RefTree_biasco1_10_List [!(RefTree_biasco1_10_List  %in% biasco1_10_List )] #3540

#To compare the 10 data-driven samples with random and stratified samples
#For fully random 20_1
#read the refernce tree
ReferenceTree <- read.tree("80_1_new")

RefTree_biasco1_1_List <- foreach(i=1:length(bias_co1_1species)) %do% getSisters(ReferenceTree,bias_co1_1species[i],mode="label")
biasco1_1_List <- foreach(i=1:length(bias_co1_1species)) %do% getSisters(B1_tree,bias_co1_1species[i],mode="label")

#get the number of species have similar sisters by comparing with reference tree
RefTree_biasco1_1_List [(RefTree_biasco1_1_List  %in% biasco1_1_List )] #163
#get the number of species doesn't have similar sisters by comparing with reference tree
RefTree_biasco1_1_List [!(RefTree_biasco1_1_List  %in% biasco1_1_List )] #3728

RefTree_biasco1_2_List <- foreach(i=1:length(bias_co1_2species)) %do% getSisters(ReferenceTree,bias_co1_2species[i],mode="label")
biasco1_2_List <- foreach(i=1:length(bias_co1_2species)) %do% getSisters(B2_tree,bias_co1_2species[i],mode="label")

#get the number of species have similar sisters by comparing with reference tree
RefTree_biasco1_2_List [(RefTree_biasco1_2_List  %in% biasco1_2_List )] #230
#get the number of species doesn't have similar sisters by comparing with reference tree
RefTree_biasco1_2_List [!(RefTree_biasco1_2_List  %in% biasco1_2_List )] #3362

RefTree_biasco1_3_List <- foreach(i=1:length(bias_co1_3species)) %do% getSisters(ReferenceTree,bias_co1_3species[i],mode="label")
biasco1_3_List <- foreach(i=1:length(bias_co1_3species)) %do% getSisters(B3_tree,bias_co1_3species[i],mode="label")

#get the number of species have similar sisters by comparing with reference tree
RefTree_biasco1_3_List [(RefTree_biasco1_3_List  %in% biasco1_3_List )] #344
#get the number of species doesn't have similar sisters by comparing with reference tree
RefTree_biasco1_3_List [!(RefTree_biasco1_3_List  %in% biasco1_3_List )] #2915

RefTree_biasco1_4_List <- foreach(i=1:length(bias_co1_4species)) %do% getSisters(ReferenceTree,bias_co1_4species[i],mode="label")
biasco1_4_List <- foreach(i=1:length(bias_co1_4species)) %do% getSisters(B4_tree,bias_co1_4species[i],mode="label")

#get the number of species have similar sisters by comparing with reference tree
RefTree_biasco1_4_List [(RefTree_biasco1_4_List  %in% biasco1_4_List )] #315
#get the number of species doesn't have similar sisters by comparing with reference tree
RefTree_biasco1_4_List [!(RefTree_biasco1_4_List  %in% biasco1_4_List )] #3236

RefTree_biasco1_5_List <- foreach(i=1:length(bias_co1_5species)) %do% getSisters(ReferenceTree,bias_co1_5species[i],mode="label")
biasco1_5_List <- foreach(i=1:length(bias_co1_5species)) %do% getSisters(B5_tree,bias_co1_5species[i],mode="label")

#get the number of species have similar sisters by comparing with reference tree
RefTree_biasco1_5_List [(RefTree_biasco1_5_List  %in% biasco1_5_List )] #222
#get the number of species doesn't have similar sisters by comparing with reference tree
RefTree_biasco1_5_List [!(RefTree_biasco1_5_List  %in% biasco1_5_List )] #3433

RefTree_biasco1_6_List <- foreach(i=1:length(bias_co1_6species)) %do% getSisters(ReferenceTree,bias_co1_6species[i],mode="label")
biasco1_6_List <- foreach(i=1:length(bias_co1_6species)) %do% getSisters(B6_tree,bias_co1_6species[i],mode="label")

#get the number of species have similar sisters by comparing with reference tree
RefTree_biasco1_6_List [(RefTree_biasco1_6_List  %in% biasco1_6_List )] #247
#get the number of species doesn't have similar sisters by comparing with reference tree
RefTree_biasco1_6_List [!(RefTree_biasco1_6_List  %in% biasco1_6_List )] #3627

RefTree_biasco1_7_List <- foreach(i=1:length(bias_co1_7species)) %do% getSisters(ReferenceTree,bias_co1_7species[i],mode="label")
biasco1_7_List <- foreach(i=1:length(bias_co1_7species)) %do% getSisters(B7_tree,bias_co1_7species[i],mode="label")

#get the number of species have similar sisters by comparing with reference tree
RefTree_biasco1_7_List [(RefTree_biasco1_7_List  %in% biasco1_7_List )] #288
#get the number of species doesn't have similar sisters by comparing with reference tree
RefTree_biasco1_7_List [!(RefTree_biasco1_7_List  %in% biasco1_7_List )] #3609

RefTree_biasco1_8_List <- foreach(i=1:length(bias_co1_8species)) %do% getSisters(ReferenceTree,bias_co1_8species[i],mode="label")
biasco1_8_List <- foreach(i=1:length(bias_co1_8species)) %do% getSisters(B8_tree,bias_co1_8species[i],mode="label")

#get the number of species have similar sisters by comparing with reference tree
RefTree_biasco1_8_List [(RefTree_biasco1_8_List  %in% biasco1_8_List )] #331
#get the number of species doesn't have similar sisters by comparing with reference tree
RefTree_biasco1_8_List [!(RefTree_biasco1_8_List  %in% biasco1_8_List )] #3264

RefTree_biasco1_9_List <- foreach(i=1:length(bias_co1_9species)) %do% getSisters(ReferenceTree,bias_co1_9species[i],mode="label")
biasco1_9_List <- foreach(i=1:length(bias_co1_9species)) %do% getSisters(B9_tree,bias_co1_9species[i],mode="label")

#get the number of species have similar sisters by comparing with reference tree
RefTree_biasco1_9_List [(RefTree_biasco1_9_List  %in% biasco1_9_List )] #366
#get the number of species doesn't have similar sisters by comparing with reference tree
RefTree_biasco1_9_List [!(RefTree_biasco1_9_List  %in% biasco1_9_List )] #3231

RefTree_biasco1_10_List <- foreach(i=1:length(bias_co1_10species)) %do% getSisters(ReferenceTree,bias_co1_10species[i],mode="label")
biasco1_10_List <- foreach(i=1:length(bias_co1_10species)) %do% getSisters(B10_tree,bias_co1_10species[i],mode="label")

#get the number of species have similar sisters by comparing with reference tree
RefTree_biasco1_10_List [(RefTree_biasco1_10_List  %in% biasco1_10_List )] #269
#get the number of species doesn't have similar sisters by comparing with reference tree
RefTree_biasco1_10_List [!(RefTree_biasco1_10_List  %in% biasco1_10_List )] #3503

#To compare the 10 data-driven samples with random and stratified samples
#For fully random 20_1
#read the refernce tree
ReferenceTree <- read.tree("S20_1_new")

RefTree_biasco1_1_List <- foreach(i=1:length(bias_co1_1species)) %do% getSisters(ReferenceTree,bias_co1_1species[i],mode="label")
biasco1_1_List <- foreach(i=1:length(bias_co1_1species)) %do% getSisters(B1_tree,bias_co1_1species[i],mode="label")

#get the number of species have similar sisters by comparing with reference tree
RefTree_biasco1_1_List [(RefTree_biasco1_1_List  %in% biasco1_1_List )] #60
#get the number of species doesn't have similar sisters by comparing with reference tree
RefTree_biasco1_1_List [!(RefTree_biasco1_1_List  %in% biasco1_1_List )] #3831

RefTree_biasco1_2_List <- foreach(i=1:length(bias_co1_2species)) %do% getSisters(ReferenceTree,bias_co1_2species[i],mode="label")
biasco1_2_List <- foreach(i=1:length(bias_co1_2species)) %do% getSisters(B2_tree,bias_co1_2species[i],mode="label")

#get the number of species have similar sisters by comparing with reference tree
RefTree_biasco1_2_List [(RefTree_biasco1_2_List  %in% biasco1_2_List )] #90
#get the number of species doesn't have similar sisters by comparing with reference tree
RefTree_biasco1_2_List [!(RefTree_biasco1_2_List  %in% biasco1_2_List )] #3502

RefTree_biasco1_3_List <- foreach(i=1:length(bias_co1_3species)) %do% getSisters(ReferenceTree,bias_co1_3species[i],mode="label")
biasco1_3_List <- foreach(i=1:length(bias_co1_3species)) %do% getSisters(B3_tree,bias_co1_3species[i],mode="label")

#get the number of species have similar sisters by comparing with reference tree
RefTree_biasco1_3_List [(RefTree_biasco1_3_List  %in% biasco1_3_List )] #146
#get the number of species doesn't have similar sisters by comparing with reference tree
RefTree_biasco1_3_List [!(RefTree_biasco1_3_List  %in% biasco1_3_List )] #3113

RefTree_biasco1_4_List <- foreach(i=1:length(bias_co1_4species)) %do% getSisters(ReferenceTree,bias_co1_4species[i],mode="label")
biasco1_4_List <- foreach(i=1:length(bias_co1_4species)) %do% getSisters(B4_tree,bias_co1_4species[i],mode="label")

#get the number of species have similar sisters by comparing with reference tree
RefTree_biasco1_4_List [(RefTree_biasco1_4_List  %in% biasco1_4_List )] #149
#get the number of species doesn't have similar sisters by comparing with reference tree
RefTree_biasco1_4_List [!(RefTree_biasco1_4_List  %in% biasco1_4_List )] #3402

RefTree_biasco1_5_List <- foreach(i=1:length(bias_co1_5species)) %do% getSisters(ReferenceTree,bias_co1_5species[i],mode="label")
biasco1_5_List <- foreach(i=1:length(bias_co1_5species)) %do% getSisters(B5_tree,bias_co1_5species[i],mode="label")

#get the number of species have similar sisters by comparing with reference tree
RefTree_biasco1_5_List [(RefTree_biasco1_5_List  %in% biasco1_5_List )] #108
#get the number of species doesn't have similar sisters by comparing with reference tree
RefTree_biasco1_5_List [!(RefTree_biasco1_5_List  %in% biasco1_5_List )] #3547

RefTree_biasco1_6_List <- foreach(i=1:length(bias_co1_6species)) %do% getSisters(ReferenceTree,bias_co1_6species[i],mode="label")
biasco1_6_List <- foreach(i=1:length(bias_co1_6species)) %do% getSisters(B6_tree,bias_co1_6species[i],mode="label")

#get the number of species have similar sisters by comparing with reference tree
RefTree_biasco1_6_List [(RefTree_biasco1_6_List  %in% biasco1_6_List )] #145
#get the number of species doesn't have similar sisters by comparing with reference tree
RefTree_biasco1_6_List [!(RefTree_biasco1_6_List  %in% biasco1_6_List )] #3729

RefTree_biasco1_7_List <- foreach(i=1:length(bias_co1_7species)) %do% getSisters(ReferenceTree,bias_co1_7species[i],mode="label")
biasco1_7_List <- foreach(i=1:length(bias_co1_7species)) %do% getSisters(B7_tree,bias_co1_7species[i],mode="label")

#get the number of species have similar sisters by comparing with reference tree
RefTree_biasco1_7_List [(RefTree_biasco1_7_List  %in% biasco1_7_List )] #132
#get the number of species doesn't have similar sisters by comparing with reference tree
RefTree_biasco1_7_List [!(RefTree_biasco1_7_List  %in% biasco1_7_List )] #3765

RefTree_biasco1_8_List <- foreach(i=1:length(bias_co1_8species)) %do% getSisters(ReferenceTree,bias_co1_8species[i],mode="label")
biasco1_8_List <- foreach(i=1:length(bias_co1_8species)) %do% getSisters(B8_tree,bias_co1_8species[i],mode="label")

#get the number of species have similar sisters by comparing with reference tree
RefTree_biasco1_8_List [(RefTree_biasco1_8_List  %in% biasco1_8_List )] #148
#get the number of species doesn't have similar sisters by comparing with reference tree
RefTree_biasco1_8_List [!(RefTree_biasco1_8_List  %in% biasco1_8_List )] #3447

RefTree_biasco1_9_List <- foreach(i=1:length(bias_co1_9species)) %do% getSisters(ReferenceTree,bias_co1_9species[i],mode="label")
biasco1_9_List <- foreach(i=1:length(bias_co1_9species)) %do% getSisters(B9_tree,bias_co1_9species[i],mode="label")

#get the number of species have similar sisters by comparing with reference tree
RefTree_biasco1_9_List [(RefTree_biasco1_9_List  %in% biasco1_9_List )] #188
#get the number of species doesn't have similar sisters by comparing with reference tree
RefTree_biasco1_9_List [!(RefTree_biasco1_9_List  %in% biasco1_9_List )] #3409

RefTree_biasco1_10_List <- foreach(i=1:length(bias_co1_10species)) %do% getSisters(ReferenceTree,bias_co1_10species[i],mode="label")
biasco1_10_List <- foreach(i=1:length(bias_co1_10species)) %do% getSisters(B10_tree,bias_co1_10species[i],mode="label")

#get the number of species have similar sisters by comparing with reference tree
RefTree_biasco1_10_List [(RefTree_biasco1_10_List  %in% biasco1_10_List )] #148
#get the number of species doesn't have similar sisters by comparing with reference tree
RefTree_biasco1_10_List [!(RefTree_biasco1_10_List  %in% biasco1_10_List )] #3624

#To compare the 10 data-driven samples with random and stratified samples
#For fully random 20_1
#read the refernce tree
ReferenceTree <- read.tree("S40_1_new")

RefTree_biasco1_1_List <- foreach(i=1:length(bias_co1_1species)) %do% getSisters(ReferenceTree,bias_co1_1species[i],mode="label")
biasco1_1_List <- foreach(i=1:length(bias_co1_1species)) %do% getSisters(B1_tree,bias_co1_1species[i],mode="label")

#get the number of species have similar sisters by comparing with reference tree
RefTree_biasco1_1_List [(RefTree_biasco1_1_List  %in% biasco1_1_List )] #113
#get the number of species doesn't have similar sisters by comparing with reference tree
RefTree_biasco1_1_List [!(RefTree_biasco1_1_List  %in% biasco1_1_List )] #3778

RefTree_biasco1_2_List <- foreach(i=1:length(bias_co1_2species)) %do% getSisters(ReferenceTree,bias_co1_2species[i],mode="label")
biasco1_2_List <- foreach(i=1:length(bias_co1_2species)) %do% getSisters(B2_tree,bias_co1_2species[i],mode="label")

#get the number of species have similar sisters by comparing with reference tree
RefTree_biasco1_2_List [(RefTree_biasco1_2_List  %in% biasco1_2_List )] #169
#get the number of species doesn't have similar sisters by comparing with reference tree
RefTree_biasco1_2_List [!(RefTree_biasco1_2_List  %in% biasco1_2_List )] #3423

RefTree_biasco1_3_List <- foreach(i=1:length(bias_co1_3species)) %do% getSisters(ReferenceTree,bias_co1_3species[i],mode="label")
biasco1_3_List <- foreach(i=1:length(bias_co1_3species)) %do% getSisters(B3_tree,bias_co1_3species[i],mode="label")

#get the number of species have similar sisters by comparing with reference tree
RefTree_biasco1_3_List [(RefTree_biasco1_3_List  %in% biasco1_3_List )] #229
#get the number of species doesn't have similar sisters by comparing with reference tree
RefTree_biasco1_3_List [!(RefTree_biasco1_3_List  %in% biasco1_3_List )] #3030

RefTree_biasco1_4_List <- foreach(i=1:length(bias_co1_4species)) %do% getSisters(ReferenceTree,bias_co1_4species[i],mode="label")
biasco1_4_List <- foreach(i=1:length(bias_co1_4species)) %do% getSisters(B4_tree,bias_co1_4species[i],mode="label")

#get the number of species have similar sisters by comparing with reference tree
RefTree_biasco1_4_List [(RefTree_biasco1_4_List  %in% biasco1_4_List )] #231
#get the number of species doesn't have similar sisters by comparing with reference tree
RefTree_biasco1_4_List [!(RefTree_biasco1_4_List  %in% biasco1_4_List )] #3320

RefTree_biasco1_5_List <- foreach(i=1:length(bias_co1_5species)) %do% getSisters(ReferenceTree,bias_co1_5species[i],mode="label")
biasco1_5_List <- foreach(i=1:length(bias_co1_5species)) %do% getSisters(B5_tree,bias_co1_5species[i],mode="label")

#get the number of species have similar sisters by comparing with reference tree
RefTree_biasco1_5_List [(RefTree_biasco1_5_List  %in% biasco1_5_List )] #178
#get the number of species doesn't have similar sisters by comparing with reference tree
RefTree_biasco1_5_List [!(RefTree_biasco1_5_List  %in% biasco1_5_List )] #3477

RefTree_biasco1_6_List <- foreach(i=1:length(bias_co1_6species)) %do% getSisters(ReferenceTree,bias_co1_6species[i],mode="label")
biasco1_6_List <- foreach(i=1:length(bias_co1_6species)) %do% getSisters(B6_tree,bias_co1_6species[i],mode="label")

#get the number of species have similar sisters by comparing with reference tree
RefTree_biasco1_6_List [(RefTree_biasco1_6_List  %in% biasco1_6_List )] #192
#get the number of species doesn't have similar sisters by comparing with reference tree
RefTree_biasco1_6_List [!(RefTree_biasco1_6_List  %in% biasco1_6_List )] #3682

RefTree_biasco1_7_List <- foreach(i=1:length(bias_co1_7species)) %do% getSisters(ReferenceTree,bias_co1_7species[i],mode="label")
biasco1_7_List <- foreach(i=1:length(bias_co1_7species)) %do% getSisters(B7_tree,bias_co1_7species[i],mode="label")

#get the number of species have similar sisters by comparing with reference tree
RefTree_biasco1_7_List [(RefTree_biasco1_7_List  %in% biasco1_7_List )] #234
#get the number of species doesn't have similar sisters by comparing with reference tree
RefTree_biasco1_7_List [!(RefTree_biasco1_7_List  %in% biasco1_7_List )] #3663

RefTree_biasco1_8_List <- foreach(i=1:length(bias_co1_8species)) %do% getSisters(ReferenceTree,bias_co1_8species[i],mode="label")
biasco1_8_List <- foreach(i=1:length(bias_co1_8species)) %do% getSisters(B8_tree,bias_co1_8species[i],mode="label")

#get the number of species have similar sisters by comparing with reference tree
RefTree_biasco1_8_List [(RefTree_biasco1_8_List  %in% biasco1_8_List )] #238
#get the number of species doesn't have similar sisters by comparing with reference tree
RefTree_biasco1_8_List [!(RefTree_biasco1_8_List  %in% biasco1_8_List )] #3357

RefTree_biasco1_9_List <- foreach(i=1:length(bias_co1_9species)) %do% getSisters(ReferenceTree,bias_co1_9species[i],mode="label")
biasco1_9_List <- foreach(i=1:length(bias_co1_9species)) %do% getSisters(B9_tree,bias_co1_9species[i],mode="label")

#get the number of species have similar sisters by comparing with reference tree
RefTree_biasco1_9_List [(RefTree_biasco1_9_List  %in% biasco1_9_List )] #300
#get the number of species doesn't have similar sisters by comparing with reference tree
RefTree_biasco1_9_List [!(RefTree_biasco1_9_List  %in% biasco1_9_List )] #3297

RefTree_biasco1_10_List <- foreach(i=1:length(bias_co1_10species)) %do% getSisters(ReferenceTree,bias_co1_10species[i],mode="label")
biasco1_10_List <- foreach(i=1:length(bias_co1_10species)) %do% getSisters(B10_tree,bias_co1_10species[i],mode="label")

#get the number of species have similar sisters by comparing with reference tree
RefTree_biasco1_10_List [(RefTree_biasco1_10_List  %in% biasco1_10_List )] #210
#get the number of species doesn't have similar sisters by comparing with reference tree
RefTree_biasco1_10_List [!(RefTree_biasco1_10_List  %in% biasco1_10_List )] #3562

#To compare the 10 data-driven samples with random and stratified samples
#For fully random 20_1
#read the refernce tree
ReferenceTree <- read.tree("S60_1_new")

RefTree_biasco1_1_List <- foreach(i=1:length(bias_co1_1species)) %do% getSisters(ReferenceTree,bias_co1_1species[i],mode="label")
biasco1_1_List <- foreach(i=1:length(bias_co1_1species)) %do% getSisters(B1_tree,bias_co1_1species[i],mode="label")

#get the number of species have similar sisters by comparing with reference tree
RefTree_biasco1_1_List [(RefTree_biasco1_1_List  %in% biasco1_1_List )] #142
#get the number of species doesn't have similar sisters by comparing with reference tree
RefTree_biasco1_1_List [!(RefTree_biasco1_1_List  %in% biasco1_1_List )] #3749

RefTree_biasco1_2_List <- foreach(i=1:length(bias_co1_2species)) %do% getSisters(ReferenceTree,bias_co1_2species[i],mode="label")
biasco1_2_List <- foreach(i=1:length(bias_co1_2species)) %do% getSisters(B2_tree,bias_co1_2species[i],mode="label")

#get the number of species have similar sisters by comparing with reference tree
RefTree_biasco1_2_List [(RefTree_biasco1_2_List  %in% biasco1_2_List )] #210
#get the number of species doesn't have similar sisters by comparing with reference tree
RefTree_biasco1_2_List [!(RefTree_biasco1_2_List  %in% biasco1_2_List )] #3382

RefTree_biasco1_3_List <- foreach(i=1:length(bias_co1_3species)) %do% getSisters(ReferenceTree,bias_co1_3species[i],mode="label")
biasco1_3_List <- foreach(i=1:length(bias_co1_3species)) %do% getSisters(B3_tree,bias_co1_3species[i],mode="label")

#get the number of species have similar sisters by comparing with reference tree
RefTree_biasco1_3_List [(RefTree_biasco1_3_List  %in% biasco1_3_List )] #301
#get the number of species doesn't have similar sisters by comparing with reference tree
RefTree_biasco1_3_List [!(RefTree_biasco1_3_List  %in% biasco1_3_List )] #2958

RefTree_biasco1_4_List <- foreach(i=1:length(bias_co1_4species)) %do% getSisters(ReferenceTree,bias_co1_4species[i],mode="label")
biasco1_4_List <- foreach(i=1:length(bias_co1_4species)) %do% getSisters(B4_tree,bias_co1_4species[i],mode="label")

#get the number of species have similar sisters by comparing with reference tree
RefTree_biasco1_4_List [(RefTree_biasco1_4_List  %in% biasco1_4_List )] #292
#get the number of species doesn't have similar sisters by comparing with reference tree
RefTree_biasco1_4_List [!(RefTree_biasco1_4_List  %in% biasco1_4_List )] #3259

RefTree_biasco1_5_List <- foreach(i=1:length(bias_co1_5species)) %do% getSisters(ReferenceTree,bias_co1_5species[i],mode="label")
biasco1_5_List <- foreach(i=1:length(bias_co1_5species)) %do% getSisters(B5_tree,bias_co1_5species[i],mode="label")

#get the number of species have similar sisters by comparing with reference tree
RefTree_biasco1_5_List [(RefTree_biasco1_5_List  %in% biasco1_5_List )] #199
#get the number of species doesn't have similar sisters by comparing with reference tree
RefTree_biasco1_5_List [!(RefTree_biasco1_5_List  %in% biasco1_5_List )] #3456

RefTree_biasco1_6_List <- foreach(i=1:length(bias_co1_6species)) %do% getSisters(ReferenceTree,bias_co1_6species[i],mode="label")
biasco1_6_List <- foreach(i=1:length(bias_co1_6species)) %do% getSisters(B6_tree,bias_co1_6species[i],mode="label")

#get the number of species have similar sisters by comparing with reference tree
RefTree_biasco1_6_List [(RefTree_biasco1_6_List  %in% biasco1_6_List )] #260
#get the number of species doesn't have similar sisters by comparing with reference tree
RefTree_biasco1_6_List [!(RefTree_biasco1_6_List  %in% biasco1_6_List )] #3614

RefTree_biasco1_7_List <- foreach(i=1:length(bias_co1_7species)) %do% getSisters(ReferenceTree,bias_co1_7species[i],mode="label")
biasco1_7_List <- foreach(i=1:length(bias_co1_7species)) %do% getSisters(B7_tree,bias_co1_7species[i],mode="label")

#get the number of species have similar sisters by comparing with reference tree
RefTree_biasco1_7_List [(RefTree_biasco1_7_List  %in% biasco1_7_List )] #257
#get the number of species doesn't have similar sisters by comparing with reference tree
RefTree_biasco1_7_List [!(RefTree_biasco1_7_List  %in% biasco1_7_List )] #3640

RefTree_biasco1_8_List <- foreach(i=1:length(bias_co1_8species)) %do% getSisters(ReferenceTree,bias_co1_8species[i],mode="label")
biasco1_8_List <- foreach(i=1:length(bias_co1_8species)) %do% getSisters(B8_tree,bias_co1_8species[i],mode="label")

#get the number of species have similar sisters by comparing with reference tree
RefTree_biasco1_8_List [(RefTree_biasco1_8_List  %in% biasco1_8_List )] #301
#get the number of species doesn't have similar sisters by comparing with reference tree
RefTree_biasco1_8_List [!(RefTree_biasco1_8_List  %in% biasco1_8_List )] #3294

RefTree_biasco1_9_List <- foreach(i=1:length(bias_co1_9species)) %do% getSisters(ReferenceTree,bias_co1_9species[i],mode="label")
biasco1_9_List <- foreach(i=1:length(bias_co1_9species)) %do% getSisters(B9_tree,bias_co1_9species[i],mode="label")

#get the number of species have similar sisters by comparing with reference tree
RefTree_biasco1_9_List [(RefTree_biasco1_9_List  %in% biasco1_9_List )] #342
#get the number of species doesn't have similar sisters by comparing with reference tree
RefTree_biasco1_9_List [!(RefTree_biasco1_9_List  %in% biasco1_9_List )] #3255

RefTree_biasco1_10_List <- foreach(i=1:length(bias_co1_10species)) %do% getSisters(ReferenceTree,bias_co1_10species[i],mode="label")
biasco1_10_List <- foreach(i=1:length(bias_co1_10species)) %do% getSisters(B10_tree,bias_co1_10species[i],mode="label")

#get the number of species have similar sisters by comparing with reference tree
RefTree_biasco1_10_List [(RefTree_biasco1_10_List  %in% biasco1_10_List )] #243
#get the number of species doesn't have similar sisters by comparing with reference tree
RefTree_biasco1_10_List [!(RefTree_biasco1_10_List  %in% biasco1_10_List )] #3529


#To compare the 10 data-driven samples with random and stratified samples
#For fully random 20_1
#read the refernce tree
ReferenceTree <- read.tree("S80_1_new")

RefTree_biasco1_1_List <- foreach(i=1:length(bias_co1_1species)) %do% getSisters(ReferenceTree,bias_co1_1species[i],mode="label")
biasco1_1_List <- foreach(i=1:length(bias_co1_1species)) %do% getSisters(B1_tree,bias_co1_1species[i],mode="label")

#get the number of species have similar sisters by comparing with reference tree
RefTree_biasco1_1_List [(RefTree_biasco1_1_List  %in% biasco1_1_List )] #136
#get the number of species doesn't have similar sisters by comparing with reference tree
RefTree_biasco1_1_List [!(RefTree_biasco1_1_List  %in% biasco1_1_List )] #2755

RefTree_biasco1_2_List <- foreach(i=1:length(bias_co1_2species)) %do% getSisters(ReferenceTree,bias_co1_2species[i],mode="label")
biasco1_2_List <- foreach(i=1:length(bias_co1_2species)) %do% getSisters(B2_tree,bias_co1_2species[i],mode="label")

#get the number of species have similar sisters by comparing with reference tree
RefTree_biasco1_2_List [(RefTree_biasco1_2_List  %in% biasco1_2_List )] #240
#get the number of species doesn't have similar sisters by comparing with reference tree
RefTree_biasco1_2_List [!(RefTree_biasco1_2_List  %in% biasco1_2_List )] #3352

RefTree_biasco1_3_List <- foreach(i=1:length(bias_co1_3species)) %do% getSisters(ReferenceTree,bias_co1_3species[i],mode="label")
biasco1_3_List <- foreach(i=1:length(bias_co1_3species)) %do% getSisters(B3_tree,bias_co1_3species[i],mode="label")

#get the number of species have similar sisters by comparing with reference tree
RefTree_biasco1_3_List [(RefTree_biasco1_3_List  %in% biasco1_3_List )] #353
#get the number of species doesn't have similar sisters by comparing with reference tree
RefTree_biasco1_3_List [!(RefTree_biasco1_3_List  %in% biasco1_3_List )] #2906

RefTree_biasco1_4_List <- foreach(i=1:length(bias_co1_4species)) %do% getSisters(ReferenceTree,bias_co1_4species[i],mode="label")
biasco1_4_List <- foreach(i=1:length(bias_co1_4species)) %do% getSisters(B4_tree,bias_co1_4species[i],mode="label")

#get the number of species have similar sisters by comparing with reference tree
RefTree_biasco1_4_List [(RefTree_biasco1_4_List  %in% biasco1_4_List )] #312
#get the number of species doesn't have similar sisters by comparing with reference tree
RefTree_biasco1_4_List [!(RefTree_biasco1_4_List  %in% biasco1_4_List )] #3239

RefTree_biasco1_5_List <- foreach(i=1:length(bias_co1_5species)) %do% getSisters(ReferenceTree,bias_co1_5species[i],mode="label")
biasco1_5_List <- foreach(i=1:length(bias_co1_5species)) %do% getSisters(B5_tree,bias_co1_5species[i],mode="label")

#get the number of species have similar sisters by comparing with reference tree
RefTree_biasco1_5_List [(RefTree_biasco1_5_List  %in% biasco1_5_List )] #225
#get the number of species doesn't have similar sisters by comparing with reference tree
RefTree_biasco1_5_List [!(RefTree_biasco1_5_List  %in% biasco1_5_List )] #3430

RefTree_biasco1_6_List <- foreach(i=1:length(bias_co1_6species)) %do% getSisters(ReferenceTree,bias_co1_6species[i],mode="label")
biasco1_6_List <- foreach(i=1:length(bias_co1_6species)) %do% getSisters(B6_tree,bias_co1_6species[i],mode="label")

#get the number of species have similar sisters by comparing with reference tree
RefTree_biasco1_6_List [(RefTree_biasco1_6_List  %in% biasco1_6_List )] #265
#get the number of species doesn't have similar sisters by comparing with reference tree
RefTree_biasco1_6_List [!(RefTree_biasco1_6_List  %in% biasco1_6_List )] #3609

RefTree_biasco1_7_List <- foreach(i=1:length(bias_co1_7species)) %do% getSisters(ReferenceTree,bias_co1_7species[i],mode="label")
biasco1_7_List <- foreach(i=1:length(bias_co1_7species)) %do% getSisters(B7_tree,bias_co1_7species[i],mode="label")

#get the number of species have similar sisters by comparing with reference tree
RefTree_biasco1_7_List [(RefTree_biasco1_7_List  %in% biasco1_7_List )] #302
#get the number of species doesn't have similar sisters by comparing with reference tree
RefTree_biasco1_7_List [!(RefTree_biasco1_7_List  %in% biasco1_7_List )] #3595

RefTree_biasco1_8_List <- foreach(i=1:length(bias_co1_8species)) %do% getSisters(ReferenceTree,bias_co1_8species[i],mode="label")
biasco1_8_List <- foreach(i=1:length(bias_co1_8species)) %do% getSisters(B8_tree,bias_co1_8species[i],mode="label")

#get the number of species have similar sisters by comparing with reference tree
RefTree_biasco1_8_List [(RefTree_biasco1_8_List  %in% biasco1_8_List )] #345
#get the number of species doesn't have similar sisters by comparing with reference tree
RefTree_biasco1_8_List [!(RefTree_biasco1_8_List  %in% biasco1_8_List )] #3250

RefTree_biasco1_9_List <- foreach(i=1:length(bias_co1_9species)) %do% getSisters(ReferenceTree,bias_co1_9species[i],mode="label")
biasco1_9_List <- foreach(i=1:length(bias_co1_9species)) %do% getSisters(B9_tree,bias_co1_9species[i],mode="label")

#get the number of species have similar sisters by comparing with reference tree
RefTree_biasco1_9_List [(RefTree_biasco1_9_List  %in% biasco1_9_List )] #372
#get the number of species doesn't have similar sisters by comparing with reference tree
RefTree_biasco1_9_List [!(RefTree_biasco1_9_List  %in% biasco1_9_List )] #3225

RefTree_biasco1_10_List <- foreach(i=1:length(bias_co1_10species)) %do% getSisters(ReferenceTree,bias_co1_10species[i],mode="label")
biasco1_10_List <- foreach(i=1:length(bias_co1_10species)) %do% getSisters(B10_tree,bias_co1_10species[i],mode="label")

#get the number of species have similar sisters by comparing with reference tree
RefTree_biasco1_10_List [(RefTree_biasco1_10_List  %in% biasco1_10_List )] #274
#get the number of species doesn't have similar sisters by comparing with reference tree
RefTree_biasco1_10_List [!(RefTree_biasco1_10_List  %in% biasco1_10_List )] #3498


#checking
a <- list(a = 1:8,b= 2:3)
b <- list(a = 2:9,c=2:3)
a[(a %in% b )] #158
intersect(a,b)
#speciesList  [(speciesList   %in% RefTree_speciesList)]

class(Species_df)
Species_df
data <- c("Abalistes_stellaris","Abbottina_rivularis")
getSisters(tree,data,mode="label")


## starting here
dd<-lapply(1:tree$Nnode+Ntip(tree),function(n,t)
  Descendants(t,n)[[1]],t=tree)
nodes<-c(1:tree$Nnode+Ntip(tree))[which(sapply(dd,
                                               length)==3)]
sisters<-t(sapply(nodes,function(n,t) 
  t$tip.label[Descendants(t,n)[[1]]],t=tree))
rownames(sisters)<-nodes
sisters
#-------------
dd<-lapply(1:tree$Nnode+Ntip(tree),function(n,t)
  Descendants(t,n)[[1]],t=tree)
nodes<-c(1:tree$Nnode+Ntip(tree))[which(sapply(dd,length)==2)]
sisters<-t(sapply(nodes,function(n,t)
  t$tip.label[Descendants(t,n)[[1]]],t=tree))
#-------------

#-------------

