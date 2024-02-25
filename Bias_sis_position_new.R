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

RefTree_biasco1_1_List  <- foreach(i=1:length(bias_co1_1species)) %do% sister(ReferenceTree,bias_co1_1species[i],type="terminal",label=T)
biasco1_1_List <- foreach(i=1:length(bias_co1_1species)) %do% sister(B1_tree,bias_co1_1species[i],type="terminal",label=T)
#get the number of species have similar sisters by comparing with reference tree
RefTree_biasco1_1_List [(RefTree_biasco1_1_List  %in% biasco1_1_List )] #138
#get the number of species doesn't have similar sisters by comparing with reference tree
RefTree_biasco1_1_List [!(RefTree_biasco1_1_List  %in% biasco1_1_List )] #3753

RefTree_biasco1_2_List  <- foreach(i=1:length(bias_co1_2species)) %do% sister(ReferenceTree,bias_co1_2species[i],type="terminal",label=T)
biasco1_2_List <- foreach(i=1:length(bias_co1_2species)) %do% sister(B2_tree,bias_co1_2species[i],type="terminal",label=T)

#get the number of species have similar sisters by comparing with reference tree
RefTree_biasco1_2_List [(RefTree_biasco1_2_List  %in% biasco1_2_List )] #222
#get the number of species doesn't have similar sisters by comparing with reference tree
RefTree_biasco1_2_List [!(RefTree_biasco1_2_List  %in% biasco1_2_List )] #3370

RefTree_biasco1_3_List  <- foreach(i=1:length(bias_co1_3species)) %do% sister(ReferenceTree,bias_co1_3species[i],type="terminal",label=T)
biasco1_3_List <- foreach(i=1:length(bias_co1_3species)) %do% sister(B3_tree,bias_co1_3species[i],type="terminal",label=T)

#get the number of species have similar sisters by comparing with reference tree
RefTree_biasco1_3_List [(RefTree_biasco1_3_List  %in% biasco1_3_List )] #321
#get the number of species doesn't have similar sisters by comparing with reference tree
RefTree_biasco1_3_List [!(RefTree_biasco1_3_List  %in% biasco1_3_List )] #2938

RefTree_biasco1_4_List  <- foreach(i=1:length(bias_co1_4species)) %do% sister(ReferenceTree,bias_co1_4species[i],type="terminal",label=T)
biasco1_4_List <- foreach(i=1:length(bias_co1_4species)) %do% sister(B4_tree,bias_co1_4species[i],type="terminal",label=T)

#get the number of species have similar sisters by comparing with reference tree
RefTree_biasco1_4_List [(RefTree_biasco1_4_List  %in% biasco1_4_List )] #257
#get the number of species doesn't have similar sisters by comparing with reference tree
RefTree_biasco1_4_List [!(RefTree_biasco1_4_List  %in% biasco1_4_List )] #3294

RefTree_biasco1_5_List  <- foreach(i=1:length(bias_co1_5species)) %do% sister(ReferenceTree,bias_co1_5species[i],type="terminal",label=T)
biasco1_5_List <- foreach(i=1:length(bias_co1_5species)) %do% sister(B5_tree,bias_co1_5species[i],type="terminal",label=T)

#get the number of species have similar sisters by comparing with reference tree
RefTree_biasco1_5_List [(RefTree_biasco1_5_List  %in% biasco1_5_List )] #179
#get the number of species doesn't have similar sisters by comparing with reference tree
RefTree_biasco1_5_List [!(RefTree_biasco1_5_List  %in% biasco1_5_List )] #3476

RefTree_biasco1_6_List  <- foreach(i=1:length(bias_co1_6species)) %do% sister(ReferenceTree,bias_co1_6species[i],type="terminal",label=T)
biasco1_6_List <- foreach(i=1:length(bias_co1_6species)) %do% sister(B6_tree,bias_co1_6species[i],type="terminal",label=T)

#get the number of species have similar sisters by comparing with reference tree
RefTree_biasco1_6_List [(RefTree_biasco1_6_List  %in% biasco1_6_List )] #196
#get the number of species doesn't have similar sisters by comparing with reference tree
RefTree_biasco1_6_List [!(RefTree_biasco1_6_List  %in% biasco1_6_List )] #3678

RefTree_biasco1_7_List  <- foreach(i=1:length(bias_co1_7species)) %do% sister(ReferenceTree,bias_co1_7species[i],type="terminal",label=T)
biasco1_7_List <- foreach(i=1:length(bias_co1_7species)) %do% sister(B7_tree,bias_co1_7species[i],type="terminal",label=T)

#get the number of species have similar sisters by comparing with reference tree
RefTree_biasco1_7_List [(RefTree_biasco1_7_List  %in% biasco1_7_List )] #239
#get the number of species doesn't have similar sisters by comparing with reference tree
RefTree_biasco1_7_List [!(RefTree_biasco1_7_List  %in% biasco1_7_List )] #3658

RefTree_biasco1_8_List  <- foreach(i=1:length(bias_co1_8species)) %do% sister(ReferenceTree,bias_co1_8species[i],type="terminal",label=T)
biasco1_8_List <- foreach(i=1:length(bias_co1_8species)) %do% sister(B8_tree,bias_co1_8species[i],type="terminal",label=T)

#get the number of species have similar sisters by comparing with reference tree
RefTree_biasco1_8_List [(RefTree_biasco1_8_List  %in% biasco1_8_List )] #268
#get the number of species doesn't have similar sisters by comparing with reference tree
RefTree_biasco1_8_List [!(RefTree_biasco1_8_List  %in% biasco1_8_List )] #3327

RefTree_biasco1_9_List  <- foreach(i=1:length(bias_co1_9species)) %do% sister(ReferenceTree,bias_co1_9species[i],type="terminal",label=T)
biasco1_9_List <- foreach(i=1:length(bias_co1_9species)) %do% sister(B9_tree,bias_co1_9species[i],type="terminal",label=T)

#get the number of species have similar sisters by comparing with reference tree
RefTree_biasco1_9_List [(RefTree_biasco1_9_List  %in% biasco1_9_List )] #311
#get the number of species doesn't have similar sisters by comparing with reference tree
RefTree_biasco1_9_List [!(RefTree_biasco1_9_List  %in% biasco1_9_List )] #3286

RefTree_biasco1_10_List  <- foreach(i=1:length(bias_co1_10species)) %do% sister(ReferenceTree,bias_co1_10species[i],type="terminal",label=T)
biasco1_10_List <- foreach(i=1:length(bias_co1_10species)) %do% sister(B10_tree,bias_co1_10species[i],type="terminal",label=T)

#get the number of species have similar sisters by comparing with reference tree
RefTree_biasco1_10_List [(RefTree_biasco1_10_List  %in% biasco1_10_List )] #217
#get the number of species doesn't have similar sisters by comparing with reference tree
RefTree_biasco1_10_List [!(RefTree_biasco1_10_List  %in% biasco1_10_List )] #3555


