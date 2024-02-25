getwd()
#Function to combine two rows in a dataframe
combine_rows <- function(data, row1, row2) {
  data[row2, ] <- data[row1, ] + data[row2, ]
  data[-row1, ]
}

#Read the fish dataset in R
fish <- read.table("Fish_Realigned_multigene.phy")

#remove the 1st row with numbers
fish <- fish[-1,]


#LOAD PACKAGES----
# install.packages("tidyverse")
library(tidyverse)
library(dplyr)
#install.packages("rfishbase")                   
library(rfishbase)
library(stats)
library(grid)
library(futile.logger)
library(VennDiagram)
library(dplyr)
library(MASS)
library(ggplot2)
library(ape)
library(phylotools)
library(foreach)
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install(version = "3.16")
# BiocManager::install("Biostrings")
library(BiocManager)
library(Biostrings)
#STEP1----
#Filter out sequences without CO1 gene as we need only the multi-gene sequences including CO1 for downstream analysis 
#Rename columns
colnames(fish) <- c("species_name","seq")

#check the dataset
dim(fish)
length(unique(fish$species_name))
# cal_bif <- fish_multigene[902,2]
# class(cal_bif)
# as.list(cal_bif)
# 
# write_file(cal_bif,file= "Callogobius_bifasciatus.fa")
# 
# write.FASTA()

#Extract CO1 nucleotides (This step is to extract sequences with CO1 gene. So that we can remove sequences without CO1 gene)
fishCO1 <- substr(fish$seq,2710,3391)

#Check the class of fishCO1 data set
class(fishCO1)

#Convert to dataframe
fishCO1 <- data.frame(fishCO1)

#Check the converted data set
class(fishCO1)

#Checking the ACTG count of the CO1 sequences----
fishco1_new <- fishCO1
fishco1_new$fishCO1 <- DNAStringSet(fishco1_new$fishCO1)
fishCO1NucFreq <-as.data.frame(letterFrequency(fishco1_new$fishCO1,letters = c("A","T","C","G")))

#Bind the column "species name" to the fishCO1 dataframe (Co1 sequences to be used in EPA)
fishCO1_with_speciesname <- cbind(fish$species_name, fishCO1)

#Replace "-" with nothing. with this we can get rid of CO1 gene 
df_0 <- gsub('-','',fishCO1_with_speciesname$fishCO1)

#Convert to dataframe
df_0 <- as.data.frame(df_0)

#Bind the column "species name" of fish dataset with new dataframe
df_1 <- cbind(fish$species_name,df_0)

#check the dataset
dim(df_1)

#Remove rows with empty CO1 regions (For the phylogenetic placement we do need CO1 region. Therefore the empty CO1 regions are useless for further analysis)
df_2 <- subset(df_1, df_1$df_0 != "")

#check the dataset
dim(df_2)

#Rename column names
colnames(df_2) <- c("species_name","seq")

#Keep only rows with sequences
fish_multigene <- merge(fish,df_2, by= "species_name")

#check the dataset
dim(fish_multigene)
names(fish_multigene)

#New dataframe with first two columns (Get the dataframe which we need for the down stream analysis)
#This dataframe consists of species_name and multigene sequences with CO1 gene 
fish_multigene <- fish_multigene[, 1:2]

#check the dataset
dim(fish_multigene)
names(fish_multigene)
length(unique(fish_multigene$species_name))

#Rename columns
colnames(fish_multigene) <- c("species_name","seq")
class(fish_multigene)

#convert datafame to phylip format to use in RAxML
# this file has full multigene sequences (4520 species)
dat2phylip(fish_multigene,outfile= "fish_multigene.phy")
getwd()
#STEP2----
#Get the species count in each family

#As we need the family name we have to use PFCtaxonomy .csv in FishTreeOfLife database downloads
#Read the fish taxonomy spread sheet in R
fish_taxonomy <- read.csv("/home/thanu/Desktop/FishData/Stratified/PFCtaxonomy.csv")

#Extract family column and genus_species column for downstream analysis
df_fish_family_species <- fish_taxonomy[,c("family","genus.species")]

#Check data set
dim(df_fish_family_species) #31516 rows 2 columns
length(unique(df_fish_family_species$family))#490 families

#rename columns
colnames(df_fish_family_species) <- c("family","species_name")

#check the dataframe
head(df_fish_family_species) #In the down stream analysis we have to merge the species name column. So the format should be same

#Replace the space with "_" and covert the resultant list into a df
df_fish_family_species1 <- as.data.frame(gsub(' ','_',df_fish_family_species$species_name))

#check the dataframe
head(df_fish_family_species1)

#Bind the column "family" of fish_family_species data set with new dataframe df_fish_family_species1
newdf_fish_taxon <- cbind(df_fish_family_species1,df_fish_family_species$family)

dim(newdf_fish_taxon)
names(newdf_fish_taxon)
#rename columns
colnames(newdf_fish_taxon) <- c("species_name","family")

class(newdf_fish_taxon)
#merge fish taxon and sequences
#df_intersect <- intersect(newdf_fish_taxon,fish_multigene)
df_fish <- merge(newdf_fish_taxon,fish_multigene,by="species_name")
names(df_fish)
length(df_fish$species_name)
length(fish_multigene$species_name)
# #merge fix taxon with co1 sequences (important for EPA)
# df_fish_co1 <- merge(newdf_fish_taxon,df_2,by="species_name")
# names(df_fish_co1)# "species_name" "family" "seq"  
# dim(df_fish_co1)#6796 rows and 3 columns

#----
#Extract the range of COI gene to be used in EPA
fishsamplefull_co1 <- substr(fish_multigene$seq,2710,3391)
#Check the class
class(fishsamplefull_co1)
#Convert to dataframe
fishsamplefull_co1 <- data.frame(fishsamplefull_co1)
#Check the class
class(fishsamplefull_co1)
dim(fishsamplefull_co1)
#Combine columns of COI dataframe and fishsample100 dataframe to get the sepecies name
fishsamplefull_co1_with_speciesname <- cbind(fish_multigene$species_name,fishsamplefull_co1)
dim(fishsamplefull_co1_with_speciesname)
#Assign column names
colnames(fishsamplefull_co1_with_speciesname) <- c("species_name","seq")
#Checking the ACTG count of the CO1 sequences----
fishsamplefull_co1_with_speciesname_new <- fishsamplefull_co1_with_speciesname
fishsamplefull_co1_with_speciesname_new$seq <- DNAStringSet(fishsamplefull_co1_with_speciesname_new$seq)
fishsamplefull_co1_with_speciesnameNucFreq <-as.data.frame(letterFrequency(fishsamplefull_co1_with_speciesname_new$seq,letters = c("A","T","C","G")))

#Convert to phylip format
dat2phylip(fishsamplefull_co1_with_speciesname,outfile= "fishsamplefull_co1_with_speciesname.phy")

#----

#Checking the families with the number of species in each family and convert to a dataframe
df_tt <- data.frame(table(df_fish$family))
head(df_tt)

#Rename columns
colnames(df_tt) <- c("family","species_count")

#Check the umber of families in fish sample with co1 genes
length(unique(df_fish$family))#358 families

#STEP 3----
#Fishbase data set matches against the information obtained from FishTreeofLife data set 
library(rfishbase)
#Extract all of the species names that are available on FishBase
# fishbase_newversion <- load_taxa()#The data-object 'fishbase' is no longer included in R version 4.0.0
# AllFish <- fishbase
# dim(AllFish)
# dim(fishbase_newversion)
# getwd()
#write.csv(AllFish,file="fishbase.csv")
AllFish <- read.csv("fishbase.csv")
#Paste Genus and Species columns together to get the species name 
AllFish$FullName <- paste(AllFish$Genus, AllFish$Species) 
FishBaseSpecies <- AllFish$FullName # 33104 species names
#Match the species labels from FishTreeofLife dataset with the species names from FishBase

#First make FishBaseSpecies into a dataframe. So, it can be merged with the family column of AllFish dataset

dfFishBaseSpecies <- data.frame(FishBaseSpecies)
#names(dfFishBaseSpecies)
#Checking the columns of fishBase dataset
names(AllFish)
#names(f)
#combine the species name and family columns together in fishbase data
FB_family_species <- cbind(AllFish$Family,dfFishBaseSpecies)
names(FB_family_species)
dim(FB_family_species)

#Rename columns
colnames(FB_family_species) <- c("family","species_name")
class(FB_family_species)
length(unique(FB_family_species$species_name))#33104

#Checking the families with the number of species in each family
tt_FB <- table(FB_family_species$family)
class(tt_FB)
df_tt_FB <- data.frame(tt_FB <- table(FB_family_species$family))

#Rename columns
colnames(df_tt_FB) <- c("family","species_count")

#Check the number of families in fish sample
length(unique(df_tt_FB$family))#549

#Merge the two dataframes of rfishbase and fishTreeofLife dataset by "family" column
df_merged <- merge(df_tt, df_tt_FB, by="family")
dim(df_merged)
names(df_merged)

#Rename columns
colnames(df_merged) <- c("family","FishTreeofLife_species_count","FishBase_species_count")
length(df_merged$FishBase_species_count)#344 families

library(dplyr)
#Usually, fishTree data set should be a subset of Fishbase.
#So, species count in Fish base is greater than fish tree. I am checking whether there are counts greater in fishtree than fishbase

Species_count_greater <- filter(df_merged,FishTreeofLife_species_count>FishBase_species_count)
#There are 344 familes common in fishTree and FishBase
#By looking at this data set, I didn't find any species count greater in fishtree than fishbase. 

#Calculations
#Get the species count percentage per family by dividing the FishTreeofLife species count by FishBase species count
df_merged_1<- df_merged %>% mutate(present_species_count_percentage = (FishTreeofLife_species_count/ FishBase_species_count)*100)
dim(df_merged_1)#344
names(df_merged_1)

#df_merged_1
#missing data
df_merged_2<- df_merged_1 %>% mutate(missing_percentage=100-(FishTreeofLife_species_count/ FishBase_species_count)*100)
#6762/30435
dim(df_merged_2)
names(df_merged_2)

#Let's check a venn diagram,
# Chart
venn.diagram(
  x = list(df_fish$family,FB_family_species$family),
  fill = c("green", "red"),
  category.names = c("FishTree", "FishBase"),
  filename = 'venn_diagram',
  output=TRUE
)

#According to this 344 families are common in both datasets.
#14 familes are there only in fishTree and not in fishbase(Which cannot be acceptable according to our assumption: fishTree is a subset of FishBase).
#I checked these families deeply and found descripency of taxonomy between two datasets
#Hence, we can exclude these 14 for downstream analysis.
#205 familes can be seen in fishbase but not in fishTree (Which can be acceptable)
#So, I'm going to include these 205 families for downstream analysis.
#Now there are 549 families together for the missing data analysis 
#(present species percentage for these 205 families are 0% therefore, the missing species percentage is 100%)
#Number of families in the genetic dataset is 358 (in fish tree of life dataset).

#Union of the two dataset
unio <- union(df_fish$family,FB_family_species$family)
class(unio)

#get the families in fishbase but not in fishTree using the setdiff function.
setdif_FishBase <- data.frame(setdiff(unio,df_fish$family))
length(setdif_FishBase$setdiff.unio..df_fish.family.)
dim(setdif_FishBase)
names(setdif_FishBase)
#Add the column the missing species percentage which is 100% (100% missing from fish tree data set). 
#Now in the setdif_FishBase df has only the family colum. Those families are 100% missing from Fish tree dataset.
#So add a 100% missing percentage column to this dataframe
setdif_FishBase$missing_percentage <- 100

#549-344 = 205 (fish families in fishBase but not in fishTree)
#assign column names for the setdif_FishBase dataframe
colnames(setdif_FishBase) <- c("family","missing_percentage")
dim(setdif_FishBase)
names(setdif_FishBase)
head(setdif_FishBase)
#make a new dataframe by extracting family and missing percentage columns
df_merged_3 <- df_merged_2[,c("family","missing_percentage")]
dim(df_merged_3)#344

#Adding 205 families with 100% missing species percentage to 344 families, 
#the percentages are zeros here for present species
new <- rbind(df_merged_3,setdif_FishBase)
class(new)
#checking the column names of the new dataframe
names(new)

dim(new)
length(unique(new$family))
length(new$family)
#get the families in fishbase but not in fishTree using the setdiff function.
setdif_FishBase <- data.frame(setdiff(unio,df_fish$family))

#Add the column the present species percentage which is 0%
setdif_FishBase$present_species_count_percentage <- 0

#563-378 = 185 (fish families in fishTree but not in fish base)
#assign column names for the setdif_FishBase dataframe
colnames(setdif_FishBase) <- c("family","present_species_count_percentage")

#make a new dataframe by extracting family and missing percentage columns
#This dataframe does not contain 100% missing families
#In the downstream analysis combine this new dataframe with 100% missing families dataframe(setdif_FishBase)
df_merged_4 <- df_merged_2[,c("family","present_species_count_percentage")]
length(df_merged_4$family)
#Adding 185 families with 100% missing species percentage to 362 families, 
#the percentages are zeros here for present species
new_present <- rbind(df_merged_4,setdif_FishBase)
class(new_present)
#check the bottom of the dataset
tail(new_present)
dim(new_present)#549

attach(new)
#Histogram
hist(missing_percentage,main = "Histogram of percent missing species per family", xlab = "Missing species per family (%)")

#Statistics
sd(missing_percentage)#27.59942
mean(missing_percentage)#80.91226
median(missing_percentage)#93.75
var(missing_percentage)#761.7281
length(missing_percentage)#549

detach(new)
attach(new_present)
head(new_present)


#Histogram
hist(present_species_count_percentage ,main = "Histogram of percent present species per family", xlab = "Present species per family (%)")
abline(v = median(present_species_count_percentage),                     # Add line for median
       col = "red",
       lty=2,
       lwd = 2)
text(x = median(present_species_count_percentage) * 5,                 # Add text for median
     y = median(present_species_count_percentage) * 15,
     paste("Median =", median(present_species_count_percentage)),
     col = "red",
     cex = 1.2)

#Statistics
sd(present_species_count_percentage)#27.59942
mean(present_species_count_percentage)#19.08774
median(present_species_count_percentage)#6.25
var(present_species_count_percentage)#761.7281
length(present_species_count_percentage)#549

#round off the present percentage
new_present_round <- new_present %>% 
  mutate(round(new_present[,2]))
head(new_present_round)

colnames(new_present_round) <- c("family", "percentage", "presentspecies_percentage")
dim(new_present_round)#547
#sample(new_present_round$presentspecies_percentage)

set.seed(200)
random_sample_10 <- data.frame(replicate(10,sample(new_present_round$presentspecies_percentage)))
head(random_sample_10)

#now using the 10 sets of random numbers do the distribution
for (col in 2:ncol(random_sample_10)) {
  hist((random_sample_10[,col]),main = "Random percentage distribution for present species", xlab = "Random present percentages")
}
class(random_sample_10)

Biased_10samples <- cbind(new_present_round$family,random_sample_10)
head(Biased_10samples)

names(Biased_10samples)[1] <- "family"
head(Biased_10samples)
# #join two data frames together. So now sequences also in the same data frame
# joined_df <- inner_join(Biased_10samples,df_fish,by="family") 
# head(joined_df)
# 
# #Bind the dataframe with random numbers (sample_10_present) to the dataframe consists of family name and species count
# binded_df <- cbind(new_present,Biased_10samples)
# head(binded_df)

#get species count to the dataframe (Here,the families with zero species are left out)
#Here we extract species from the resultant merged df. So, families without species are removed.
#df_tt is the table of families with number of species count in each family
merged_df <- merge(df_tt,Biased_10samples,by="family")
head(df_tt)
class(merged_df)
dim(merged_df)#344 rows 12 cols

# #get species count to the dataframe (Here,the families with zero species are left out)
# merged_df <- merge(df_tt,joined_df,by="family")
#------------------------------------------------------------------
# 
# #Checking the best distribution which is closest to the real dataset
# set.seed(123)
# x <- rgamma(n=547, shape= 0.2,  rate = 1/25.50842)
# gamma <- fitdistr(x, "gamma")
# summary(gamma)
# hist(x,main= "Gamma Distribution")
# gamma$loglik
# AIC(gamma)
# 
# #To get the best shape value check the likelihood value
# set.seed(123)
# x <- rgamma(n=547, shape= 0.2,  rate = 1/25.50842)
# gamma <- fitdistr(x, "gamma")
# gamma$loglik
# 
# 
# set.seed(123)
# x4 <- rnegbin(547, mu = 25.50842, theta = 1)
# neg_bino <- fitdistr(x4, "Negative Binomial")
# neg_bino$loglik
# AIC(neg_bino)
# summary(neg_bino)
# hist(x4, main= "Negative binomial Distribution")
# 
# set.seed(123)
# x3 <- rweibull(547, shape = 1, scale = 25.50842)
# wbul <- fitdistr(x3, "weibull")
# wbul$loglik
# AIC(wbul)
# hist(x3,main= "Weibull Distribution")
# summary(wbul,main= "Distribution of Weibull")
# 
# x5 <- rexp(n=547, rate=1/25.50842)
# exp <- fitdistr(x5,"exponential")
# hist(x5,main= "Exponential distribution")
# exp$loglik
# AIC(exp)
# summary(exp)
# 
# #AIC comparison of models
# AIC(neg_bino,wbul,exp,gamma)
# 
# #I'm using gamma to simulate missingness. As it's the best model for that according to the AIC comparisson
# #random gamma distribution (rate=(1/mean))
# 
# set.seed(1002)
# #get 10 sets of 547 values 
# sample_10_gamma <- data.frame(replicate(10,rgamma(n=547, shape =0.2, rate = 1/25.50842)))
# head(sample_10_gamma)
# 
# #change values over 100 to 100
# sample_10_gamma[sample_10_gamma > 100] <- 100
# class(sample_10_gamma)
# # #In the stratified sampling we had 378 families. So here I'm gonna use 378 random values based on gamma model
# # sample_10_gamma_1 <- data.frame(sample(nrow(sample_10_gamma),378))
# # sample_10_gamma <- as.data.frame(sample_10_gamma)
# 
# #now using the 10 sets of random numbers do the distribution
# for (col in 2:ncol(sample_10_gamma)) {
#   hist((sample_10_gamma[,col]),main = "Gamma distribution for present species", xlab = "Gamma random numbers")
# }
# class(sample_10_gamma)
# #Bind the dataframe with random numbers (sample_10_present) to the dataframe consists of family name and species count
# binded_df <- cbind(new_present,sample_10_gamma)
# tail(binded_df)
# 
# #get species count to the dataframe (Here,the families with zero species are left out)
# merged_df <- merge(df_tt,binded_df,by="family")
# 
# class(merged_df)
# detach(new_present)
# attach(merged_df)
#-------------------------------------------------------------------------------
head(merged_df)
#Calculate the sample species to extract from the dataset
#merged_df1<- merged_df %>% mutate(sampled_species = (round(species_count*X1/100)))

#merged_df1<- merged_df %>% mutate(sampled_species = (ceiling(species_count*X1/100)))
#head(merged_df1)

merged_df1_present <- round(merged_df[,2]*merged_df[,3:12]/100)
head(merged_df1_present)
merged_df1_missing <- merged_df[,2]- merged_df1_present[,1:10]
head(merged_df1_missing)
tail(merged_df1_missing)

miss_wt_families <- cbind(merged_df$family,merged_df1_missing)
head(miss_wt_families)
#first we need present species to build trees
present_wt_families <- cbind(merged_df$family,merged_df1_present)
head(present_wt_families)

#Rename column names
colnames(present_wt_families) <- c("family","C1","C2","C3","C4","C5","C6","C7","C8","C9","C10")

#checking
head(present_wt_families)
class(present_wt_families)
# #get the first biased sample with family name
# present_wt_families_1 <- present_wt_families[,1:2]

#create a function to split the large data frame into small data frames

#STEP 4----
#Arrange 10 bias data sets
#Function to split the full dataframe into 10 (now all 10 bias sampling counts are in one dataframe, we need to split it)
split <- function(column){
  df_split <- present_wt_families[,c(1,column)]
  return(df_split)
}

# #A list of columns
# col_vector <- c(2:11)
# 
# #we can use foreach loop to iterate through each data frame (There were 40 data frames)
# split_col <- foreach(i=1:length(col_vector)) %do% present_wt_families[,c(1,col_vector[[i]])]
#   
# head(split_col)



#call the split function to generate 10 data frames splitting a large data frame
df_01 <- split(2)
df_02 <- split(3)
df_03 <- split(4)
df_04 <- split(5)
df_05 <- split(6)
df_06 <- split(7)
df_07 <- split(8)
df_08 <- split(9)
df_09 <- split(10)
df_10 <- split(11)

#The real bias species dataframe
df_bias <- merged_df[,c("family","species_count")]

#checking
head(df_10)
class(df_10)
class(df_bias)

#function to remove families with zero species
fam_no_zero <- function(df_new,C){
  new_df<- df_new[df_new$C >0,]
  return(new_df)
}


# call the fam_no_zero funtion to generate backbone trees that represent present species according to biased realistic samples
d01 <- fam_no_zero(df_01,"C1")
d02 <- fam_no_zero(df_02,"C2")
d03 <- fam_no_zero(df_03,"C3")
d04 <- fam_no_zero(df_04,"C4")
d05 <- fam_no_zero(df_05,"C5")
d06 <- fam_no_zero(df_06,"C6")
d07 <- fam_no_zero(df_07,"C7")
d08 <- fam_no_zero(df_08,"C8")
d09 <- fam_no_zero(df_09,"C9")
d10 <- fam_no_zero(df_10,"C10")
class(C2)
C3
colnames(df_bias) <- c("family", "C0")
#for the real bias species
dBias <- fam_no_zero(df_bias,"C0")

#pull function of the dplyr package to convert a column of a data frame into a vector
#C1 <- pull(present_wt_families_1,C1)
C1 <- pull(d01,C1)
C2 <- pull(d02,C2)
C3 <- pull(d03,C3)
C4 <- pull(d04,C4)
C5 <- pull(d05,C5)
C6 <- pull(d06,C6)
C7 <- pull(d07,C7)
C8 <- pull(d08,C8)
C9 <- pull(d09,C9)
C10 <- pull(d10,C10)
class(C1)
#Real bias species
C0 <- pull (dBias,C0)
# C1 <- as.numeric(C1)
# sample_sizes <- data.frame(
#   families = unique(d01$family),
#   n_to_sample = C1
# )
# 
# sample_size <- function (df,C){
#   sample_df <- data.frame(families = unique(df$family), n_to_sample = C)
#  return(sample_df)
# }
# 
# df_1s <- sample_size(d01,C1)

# #join two data frames together. So now sequences also in the same data frame
# #inner_join can use to drop the rows from both data where column values don't match (Here the column I used is family)
# join_d01 <- inner_join(d01,df_fish,by="family") 
# class(join_d01)
# tail(join_d01)
# names(join_d01)
# dim(join_d01)
# names(d01)
# dim(d01)
# join_d02 <- inner_join(d02,df_fish,by="family") 
# join_d03 <- inner_join(d03,df_fish,by="family") 
# join_d04 <- inner_join(d04,df_fish,by="family") 
# join_d05 <- inner_join(d05,df_fish,by="family") 
# join_d06 <- inner_join(d06,df_fish,by="family") 
# join_d07 <- inner_join(d07,df_fish,by="family") 
# join_d08 <- inner_join(d08,df_fish,by="family") 
# join_d09 <- inner_join(d09,df_fish,by="family") 
# join_d10 <- inner_join(d10,df_fish,by="family") 

#merge two data frames together. So now sequences also in the same data frame
#inner_join can use to drop the rows from both data where column values don't match (Here the column I used is family)
join_d01_new <- merge(d01,df_fish,by="family") 
join_d02_new <- merge(d02,df_fish,by="family") 
join_d03_new <- merge(d03,df_fish,by="family") 
join_d04_new <- merge(d04,df_fish,by="family") 
join_d05_new <- merge(d05,df_fish,by="family") 
join_d06_new <- merge(d06,df_fish,by="family") 
join_d07_new <- merge(d07,df_fish,by="family") 
join_d08_new <- merge(d08,df_fish,by="family") 
join_d09_new <- merge(d09,df_fish,by="family") 
join_d10_new <- merge(d10,df_fish,by="family") 


#----------------------------
#To get number of families in each bias sample
colnames(join_d01_new)
length(unique(join_d10_new$family))
#---------------------------
#Real bias species
join_dBias <- inner_join(dBias,df_fish,by="family") 
#join_dBias <- inner_join(dBias,df_fish_co1,by="family") 
# length(unique(join_dBias$family))
# join <- function(df){
#   df_seq <- inner_join(df,df_fish,by="family")
#   return(df_seq)
# }
# 
# group <- function(df_j,C){
#   df_join <- df_j %>%
#     group_by(family) %>%
#     sample_n(C)
#   return(df_join)
# 
# }
# 
# j1 <- join (d01)
# g <- group(j1,C1)

#get different samples based on family. This contains the sequences as well
df_g1 <- join_d01_new %>% 
  group_by(family) %>% 
  sample_n(C1[1])
class(df_g1)
df_g2 <- join_d02_new %>% 
  group_by(family) %>% 
  sample_n(C2[1])

df_g3 <- join_d03_new %>% 
  group_by(family) %>% 
  sample_n(C3[1])

df_g4 <- join_d04_new %>% 
  group_by(family) %>% 
  sample_n(C4[1])
df_g5 <- join_d05_new %>% 
  group_by(family) %>% 
  sample_n(C5[1])

df_g6 <- join_d06_new %>% 
  group_by(family) %>% 
  sample_n(C6[1])
df_g7 <- join_d07_new %>% 
  group_by(family) %>% 
  sample_n(C7[1])

df_g8 <- join_d08_new %>% 
  group_by(family) %>% 
  sample_n(C8[1])
df_g9 <- join_d09_new %>% 
  group_by(family) %>% 
  sample_n(C9[1])

df_g10 <- join_d10_new %>% 
  group_by(family) %>% 
  sample_n(C10[1])

df_g0 <- join_dBias %>% 
  group_by(family) %>% 
  sample_n(C0[1])

# df_g1 <- joined_df %>% 
#   group_by(family) %>% 
#   sample_n(X1, replace = TRUE) 
# 
# df_g2 <- joined_df %>% 
#   group_by(family) %>% 
#   sample_n(X2, replace = TRUE)
# 
# df_g3 <- joined_df %>% 
#   group_by(family) %>% 
#   sample_n(X3, replace = TRUE)
# 
# df_g4 <- joined_df %>% 
#   group_by(family) %>% 
#   sample_n(X4, replace = TRUE)
# df_g5 <- joined_df %>% 
#   group_by(family) %>% 
#   sample_n(X5, replace = TRUE)
# 
# df_g6 <- joined_df %>% 
#   group_by(family) %>% 
#   sample_n(X6, replace = TRUE)
# df_g7 <- joined_df %>% 
#   group_by(family) %>% 
#   sample_n(X7, replace = TRUE)
# 
# df_g8 <- joined_df %>% 
#   group_by(family) %>% 
#   sample_n(X8, replace = TRUE)
# df_g9 <- joined_df %>% 
#   group_by(family) %>% 
#   sample_n(X9, replace = TRUE)
# 
# df_g10 <- joined_df %>% 
#   group_by(family) %>% 
#   sample_n(X10, replace = TRUE)

class(df_g1)



#The data type should be dataf.frame to use dat2phylip (). otherwise the result give only one sequence in phylip output file. Soconvert to dataframe
df_g1 <- as.data.frame(df_g1)
class(df_g1)
df_g2 <- as.data.frame(df_g2)
class(df_g2)
df_g3 <- as.data.frame(df_g3)
df_g4 <- as.data.frame(df_g4)
df_g5 <- as.data.frame(df_g5)
df_g6 <- as.data.frame(df_g6)
df_g7 <- as.data.frame(df_g7)
df_g8 <- as.data.frame(df_g8)
df_g9 <- as.data.frame(df_g9)
df_g10 <- as.data.frame(df_g10)

#real bias species
df_g0 <- as.data.frame(df_g0)
#create a function to subset data frames 

# df_sub<- function(df){
#   df_subset <- df[,c(12:13)]
#   df_subset1 <- unique(df_subset)
#   return(df_subset1)
# }

df_sub<- function(df){
  df_subset <- df[,c(3:4)]
  df_subset1 <- unique(df_subset)
  return(df_subset1)
}
dim(df_g1)

#cc <- unique(df_g1[,c(12,13)])

df1 <- df_sub(df_g1)
class(df1)
df2 <- df_sub(df_g2)
class(df2)
df3 <- df_sub(df_g3)
df4 <- df_sub(df_g4)
df5 <- df_sub(df_g5)
df6 <- df_sub(df_g6)
df7 <- df_sub(df_g7)
df8 <- df_sub(df_g8)
df9 <- df_sub(df_g9)
df10 <- df_sub(df_g10)
colnames(df1)
df0 <- df_sub(df_g0)

# #Extract CO1 nucleotides (This step is to extract sequences with CO1 gene. So that we can remove sequences without CO1 gene)
# 
# df_substr <- function(dataf){
#   df_substring <- substr(dataf$seq,2292,2973)
#   return(df_substring)
# }
# 
# 
# df0_co1 <- df_substr(df0)
#df0_co1 <- substr(df0$seq,2292,2973)
# df1_co1 <- df_substr(df1)
# df2_co1 <- df_substr(df2)
# df3_co1 <- df_substr(df3)
# df4_co1 <- df_substr(df4)
# df5_co1 <- df_substr(df5)
# df6_co1 <- df_substr(df6)
# df7_co1 <- df_substr(df7)
# df8_co1 <- df_substr(df8)
# df9_co1 <- df_substr(df9)
# df10_co1 <- df_substr(df10)

#Check the class of df0_co1 dataset
#class(df0_co1)

#Convert to dataframe

#df0_co1 <- data.frame(df0_co1)
# df1_co1 <- data.frame(df1_co1)
# df2_co1 <- data.frame(df2_co1)
# df3_co1 <- data.frame(df3_co1)
# df4_co1 <- data.frame(df4_co1)
# df5_co1 <- data.frame(df5_co1)
# df6_co1 <- data.frame(df6_co1)
# df7_co1 <- data.frame(df7_co1)
# df8_co1 <- data.frame(df8_co1)
# df9_co1 <- data.frame(df9_co1)
# df10_co1 <- data.frame(df10_co1)

#Check the converted dataset
#class(df0_co1)

#Bind the column "species name" to the fishCO1 dataframe (Co1 sequences to be used in EPA)
#df0_co1_with_speciesname <- cbind(df0$species_name,df0_co1)

#Use the function to generate phylip files from the dataframe
dff <- function(xx){
  return(dat2phylip(xx,outfile = "bias.phy"))
}
#change from df1 to df10 to get 10 phylip files
dff(df0)

#fasta #we need this when using mafft to align co1 and biased sample fasta files
dataframe2fas(df10, file = "fishsample_df10.fasta")

#dff(df0_co1_with_speciesname)
getwd()

#STEP 5--------------------------------------------------------
#df1 (1st biased sample) and fishsamplefull_co1_with_speciesname (Co1 100% marker) antijoin to get species not in biased sample of co1 that exist only in co1 100% not in biased sample
CO1_df1 <- anti_join(fishsamplefull_co1_with_speciesname,df1,by="species_name")
CO1_df2 <- anti_join(fishsamplefull_co1_with_speciesname,df2,by="species_name")
CO1_df3 <- anti_join(fishsamplefull_co1_with_speciesname,df3,by="species_name")
CO1_df4 <- anti_join(fishsamplefull_co1_with_speciesname,df4,by="species_name")
CO1_df5 <- anti_join(fishsamplefull_co1_with_speciesname,df5,by="species_name")
CO1_df6 <- anti_join(fishsamplefull_co1_with_speciesname,df6,by="species_name")
CO1_df7 <- anti_join(fishsamplefull_co1_with_speciesname,df7,by="species_name")
CO1_df8 <- anti_join(fishsamplefull_co1_with_speciesname,df8,by="species_name")
CO1_df9 <- anti_join(fishsamplefull_co1_with_speciesname,df9,by="species_name")
CO1_df10 <- anti_join(fishsamplefull_co1_with_speciesname,df10,by="species_name")
CO1_df0 <- anti_join(fishsamplefull_co1_with_speciesname,df0,by="species_name")
names(df2)
names(fishsamplefull_co1_with_speciesname)
# colnames(newdfF2) <- c("species_name","seq")
# colnames(newdfF3) <- c("species_name","seq")
# colnames(newdfF4) <- c("species_name","seq")
# colnames(newdfF5) <- c("species_name","seq")
# colnames(newdfF6) <- c("species_name","seq")
# colnames(newdfF7) <- c("species_name","seq")
# colnames(newdfF8) <- c("species_name","seq")
# colnames(newdfF9) <- c("species_name","seq")
# colnames(newdfF10) <- c("species_name","seq")

# 
# colnames(fishsample60_1) <- c("species_name","seq")
# colnames(fishsample60_2) <- c("species_name","seq")
# colnames(fishsample60_3) <- c("species_name","seq")
# colnames(fishsample60_4) <- c("species_name","seq")
# colnames(fishsample60_5) <- c("species_name","seq")
# colnames(fishsample60_6) <- c("species_name","seq")
# colnames(fishsample60_7) <- c("species_name","seq")
# colnames(fishsample60_8) <- c("species_name","seq")
# colnames(fishsample60_9) <- c("species_name","seq")
# colnames(fishsample60_10) <- c("species_name","seq")

getwd()
#Genearte phylip files for aligned co1 dataset 
#Just change the dataset name (CO1_20_1...CO1_80_10) for the function dff (function dff leads to datarame to phylip format)
#Use the function to genetare phylip files from the dataframe
dff <- function(xx){
  
  return(dat2phylip(xx,outfile = "newphy.phy"))
  
}

dff(CO1_df0)
library(seqRFLP)
#Generate fasta files for aligned co1 dataset
dataframe2fas(CO1_df1, file = "CO1_df1.fasta")
dataframe2fas(CO1_df2, file = "CO1_df2.fasta")
dataframe2fas(CO1_df3, file = "CO1_df3.fasta")
dataframe2fas(CO1_df4, file = "CO1_df4.fasta")
dataframe2fas(CO1_df5, file = "CO1_df5.fasta")
dataframe2fas(CO1_df6, file = "CO1_df6.fasta")
dataframe2fas(CO1_df7, file = "CO1_df7.fasta")
dataframe2fas(CO1_df8, file = "CO1_df8.fasta")
dataframe2fas(CO1_df9, file = "CO1_df9.fasta")
dataframe2fas(CO1_df10, file = "CO1_df10.fasta")
dataframe2fas(CO1_df0, file = "CO1_df0.fasta")

#Convert fas to phylip format
#install.packages("seqmagick")
library(seqmagick)
library(seqinr)
library(phylotools)
#Here fish20_co1_80_1.fasta is the reference aligned sequence set which is an output of multigene 20% and co1 80% aligning by using muscle in Linux terminal
getwd()

# fish20_co1_80_1=read.fasta(file = "fish20_co1_80_1.fasta")
# class(fish20_co1_80_1)
# 
# fish20_co1_80_2=read.fasta(file = "fish20_co1_80_2.fasta")
# fish20_co1_80_3=read.fasta(file = "fish20_co1_80_3.fasta")
# fish20_co1_80_4=read.fasta(file = "fish20_co1_80_4.fasta")
# fish20_co1_80_5=read.fasta(file = "fish20_co1_80_5.fasta")
# fish20_co1_80_6=read.fasta(file = "fish20_co1_80_6.fasta")
# fish20_co1_80_7=read.fasta(file = "fish20_co1_80_7.fasta")
# fish20_co1_80_8=read.fasta(file = "fish20_co1_80_8.fasta")
# fish20_co1_80_9=read.fasta(file = "fish20_co1_80_9.fasta")
# fish20_co1_80_10=read.fasta(file = "fish20_co1_80_10.fasta")
# 
# fish40_co1_60_1=read.fasta(file = "fish40_co1_60_1.fasta")
# fish40_co1_60_2=read.fasta(file = "fish40_co1_60_2.fasta")
# fish40_co1_60_3=read.fasta(file = "fish40_co1_60_3.fasta")
# fish40_co1_60_4=read.fasta(file = "fish40_co1_60_4.fasta")
# fish40_co1_60_5=read.fasta(file = "fish40_co1_60_5.fasta")
# fish40_co1_60_6=read.fasta(file = "fish40_co1_60_6.fasta")
# fish40_co1_60_7=read.fasta(file = "fish40_co1_60_7.fasta")
# fish40_co1_60_8=read.fasta(file = "fish40_co1_60_8.fasta")
# fish40_co1_60_9=read.fasta(file = "fish40_co1_60_9.fasta")
# fish40_co1_60_10=read.fasta(file = "fish40_co1_60_10.fasta")
# 
# fish60_co1_40_1=read.fasta(file = "fish60_co1_40_1.fasta")
# fish60_co1_40_2=read.fasta(file = "fish60_co1_40_2.fasta")
# fish60_co1_40_3=read.fasta(file = "fish60_co1_40_3.fasta")
# fish60_co1_40_4=read.fasta(file = "fish60_co1_40_4.fasta")
# fish60_co1_40_5=read.fasta(file = "fish60_co1_40_5.fasta")
# fish60_co1_40_6=read.fasta(file = "fish60_co1_40_6.fasta")
# fish60_co1_40_7=read.fasta(file = "fish60_co1_40_7.fasta")
# fish60_co1_40_8=read.fasta(file = "fish60_co1_40_8.fasta")
# fish60_co1_40_9=read.fasta(file = "fish60_co1_40_9.fasta")
# fish60_co1_40_10=read.fasta(file = "fish60_co1_40_10.fasta")
# 
# fish80_co1_20_1=read.fasta(file = "fish80_co1_20.fasta")
# fish80_co1_20_2=read.fasta(file = "fish80_co1_20_2.fasta")
# fish80_co1_20_3=read.fasta(file = "fish80_co1_20_3.fasta")
# fish80_co1_20_4=read.fasta(file = "fish80_co1_20_4.fasta")
# fish80_co1_20_5=read.fasta(file = "fish80_co1_20_5.fasta")
# fish80_co1_20_6=read.fasta(file = "fish80_co1_20_6.fasta")
# fish80_co1_20_7=read.fasta(file = "fish80_co1_20_7.fasta")
# fish80_co1_20_8=read.fasta(file = "fish80_co1_20_8.fasta")
# fish80_co1_20_9=read.fasta(file = "fish80_co1_20_9.fasta")
# fish80_co1_20_10=read.fasta(file = "fish80_co1_20_10.fasta")
# 
# #Just change the dataset name (fish20_co1_80_1...fish80_co1_20_10) for the function dff (function dff leads to datarame to phylip format)
# dff(fish80_co1_20_10)
# 
# #STEP 6--------------------------------------------------------

# Read trees generated from RAxML into R 
# read all of the trees at once ensuring we have the correct path to our small_dataset file 
# So for me the path is /home/thanu/Desktop/FishData/missing_pattern/NewBiasPhylipfiles_Oct/ and the pattern is a regular expression to extract out the files we need
# Then I'm using foreach package to just iterate through each file
#install.packages("foreach")
library(ape)
library(foreach)
library(phangorn)
library(ggplot2)
#files <- list.files(path="/home/thanu/Desktop/FishData/missing_pattern/NewBiasPhylipfiles_Oct", pattern="BS[0-9]{1,2}_new", full.names=TRUE, recursive=FALSE)
files <- list.files(path="/home/thanu/Desktop/FishData/New_dataset/BIASED", pattern="B[0-9]{1,2}_new", full.names=TRUE, recursive=FALSE)

treeList <- foreach(i=1:length(files)) %do% read.tree(files[i])

# Now I have a list with each element being its own separate tree, so treeList[[1]] is my first tree 20_1_new, treeList[[2]] is 20_2_new etc.

# Also we can name each list element by substituting out our path from the names of each file
names(treeList) <- gsub("/home/thanu/Desktop/FishData/New_dataset/BIASED/", "", files)

# This will let us keep track of which list element corresponds to which tree file
names(treeList) 

#Standardize the branch lengths using compute.brlen (We did this because we couldn't get a good result in path distance. 
#So, we assumed that the reason could be the measures of path distance 
#(not only topology but also branch lengths take into account to measure the path length) )
#treeList<- foreach(i=1:length(treeList)) %do% compute.brlen(treeList[[i]],method = "Grafen", power = 1)
###----------------------------



#aa <- compute.brlen(treeList$PB1_new, method = "Grafen", power = 1)
#compute.brlen(treeList$PB1_new)

# # Then I keep the complete full tree ("RAxML_parsimonyTree.complete100") a separate variable since I'm using it for all the analyses below
# # I call it parseTreeComp
# parseTreeComp <- read.tree("RAxML_parsimonyTree.complete100")

# For Robinson-Foulds values, we can use foreach loop to iterate through each tree (There were 10 trees)

#I compare these 10 bias trees with the fully random complete full tree ("RAxML_bestTree.mgTree").  
# I call it fullyRandomparseTreeComp, a separate variable since I'm using it for all the analyses below
fullyRandomparseTreeComp <- read.tree("RAxML_bestTree.FR100_new")
plot(fullyRandomparseTreeComp)
fullyRandom20_1 <- read.tree("20_1_new")
fullyRandom40_1 <- read.tree("40_1_new")
fullyRandom60_1 <- read.tree("60_1_new")
fullyRandom80_1 <- read.tree("80_1_new")

fullyStrati20_1 <- read.tree("S20_1_new")
fullyStrati40_1 <- read.tree("S40_1_new")
fullyStrati60_1 <- read.tree("S60_1_new")
fullyStrati80_1 <- read.tree("S80_1_new")
#check whether unusual branches exsist in the tree
plot(fullyRandom20_1)
# fullyRandom20_1<- compute.brlen(fullyRandom20_1, method = "Grafen", power = 1)
# fullyRandom40_1<- compute.brlen(fullyRandom40_1, method = "Grafen", power = 1)
# fullyRandom60_1<- compute.brlen(fullyRandom60_1, method = "Grafen", power = 1)
# fullyRandom80_1<- compute.brlen(fullyRandom80_1, method = "Grafen", power = 1)


# plot(fullyRandom20_1)
# aa <- compute.brlen(fullyRandom20_1, method = "Grafen", power = 1)
# plot(aa)

RF_List_fullyRandomTree <- foreach(i=1:length(treeList)) %do% RF.dist(fullyRandomparseTreeComp,treeList[[i]],normalize = TRUE,check.labels = FALSE)
# RF_List_fullyRandomTree

#branch.sum <- tree@edge.length and then set tree@edge.length <- tree@edge.length/branch.sum
# treeList_new <- treeList
# 
# #call each tree in the tree list and assigned a name for each tree
# t1 <- treeList_new$RAxML_bestTree.b_df1
# t2 <- treeList_new$RAxML_bestTree.b_df2
# t3 <- treeList_new$RAxML_bestTree.b_df3
# t4 <- treeList_new$RAxML_bestTree.b_df4
# t5 <- treeList_new$RAxML_bestTree.b_df5
# t6 <- treeList_new$RAxML_bestTree.b_df6
# t7 <- treeList_new$RAxML_bestTree.b_df7
# t8 <- treeList_new$RAxML_bestTree.b_df8
# t9 <- treeList_new$RAxML_bestTree.b_df9
# t10 <- treeList_new$RAxML_bestTree.b_df10
# #standardize the tree (the sum of branch length to 1)
# t1$edge.length <- t1$edge.length/sum(t1$edge.length)
# sum(t1$edge.length)
# t2$edge.length <- t2$edge.length/sum(t2$edge.length)
# t3$edge.length <- t3$edge.length/sum(t3$edge.length)
# t4$edge.length <- t4$edge.length/sum(t4$edge.length)
# t5$edge.length <- t5$edge.length/sum(t5$edge.length)
# t6$edge.length <- t6$edge.length/sum(t6$edge.length)
# t7$edge.length <- t7$edge.length/sum(t7$edge.length)
# t8$edge.length <- t8$edge.length/sum(t8$edge.length)
# t9$edge.length <- t9$edge.length/sum(t9$edge.length)
# t10$edge.length <- t10$edge.length/sum(t10$edge.length)
# New_treeList <- list(t1,t2,t3,t4,t5,t6,t7,t8,t9,t10)
# #assign new names for the trees
# fullyRandom20_1_new <- fullyRandom20_1
# fullyRandom40_1_new <- fullyRandom40_1
# fullyRandom60_1_new <- fullyRandom60_1
# fullyRandom80_1_new <- fullyRandom80_1
# 
# fullyStrati20_1_new <- fullyStrati20_1
# fullyStrati40_1_new <- fullyStrati40_1
# fullyStrati60_1_new <- fullyStrati60_1
# fullyStrati80_1_new <- fullyStrati80_1
# #standardize the tree (the sum of branch length to 1)
# fullyRandom20_1_new$edge.length <- fullyRandom20_1$edge.length/sum(fullyRandom20_1$edge.length)
# fullyRandom40_1_new$edge.length <- fullyRandom40_1$edge.length/sum(fullyRandom40_1$edge.length)
# fullyRandom60_1_new$edge.length <- fullyRandom60_1$edge.length/sum(fullyRandom60_1$edge.length)
# fullyRandom80_1_new$edge.length <- fullyRandom80_1$edge.length/sum(fullyRandom80_1$edge.length)
# fullyStrati20_1_new$edge.length <- fullyStrati20_1$edge.length/sum(fullyStrati20_1$edge.length)
# fullyStrati40_1_new$edge.length <- fullyStrati40_1$edge.length/sum(fullyStrati40_1$edge.length)
# fullyStrati60_1_new$edge.length <- fullyStrati60_1$edge.length/sum(fullyStrati60_1$edge.length)
# fullyStrati80_1_new$edge.length <- fullyStrati80_1$edge.length/sum(fullyStrati80_1$edge.length)
# #Test whether the sum of the branch length is =1
# sum(fullyStrati80_1_new$edge.length)


# For Robinson-Foulds values, we can use foreach loop to iterate through each tree
RF_List_RandomTree20 <- foreach(i=1:length(treeList)) %do% RF.dist(fullyRandom20_1,treeList[[i]],normalize = TRUE)
RF_List_RandomTree40 <- foreach(i=1:length(treeList)) %do% RF.dist(fullyRandom40_1,treeList[[i]],normalize = TRUE)
RF_List_RandomTree60 <- foreach(i=1:length(treeList)) %do% RF.dist(fullyRandom60_1,treeList[[i]],normalize = TRUE)
RF_List_RandomTree80 <- foreach(i=1:length(treeList)) %do% RF.dist(fullyRandom80_1,treeList[[i]],normalize = TRUE)

# #playing around-----
# RF_List_RandomTree20 <- foreach(i=1:length(treeList)) %do% RF.dist(cF,treeList[[i]],normalize = TRUE,check.labels = FALSE)
# RF_List_RandomTree40 <- foreach(i=1:length(treeList)) %do% RF.dist(mF,treeList[[i]],normalize = TRUE,check.labels = FALSE)
# RF_List_RandomTree60 <- foreach(i=1:length(treeList)) %do% RF.dist(wF,treeList[[i]],normalize = TRUE,check.labels = FALSE)
# RF_List_RandomTree80 <- foreach(i=1:length(treeList)) %do% RF.dist(ggF,treeList[[i]],normalize = TRUE,check.labels = FALSE)
# RF_List_RandomTree20 <- RF.dist(aF,fullyRandomparseTreeComp,normalize = TRUE,check.labels = FALSE)
# RF_List_RandomTree40 <- RF.dist(lF,fullyRandomparseTreeComp,normalize = TRUE,check.labels = FALSE)
# RF_List_RandomTree60 <- RF.dist(vF,fullyRandomparseTreeComp,normalize = TRUE,check.labels = FALSE)
# RF_List_RandomTree80 <- RF.dist(ffF,fullyRandomparseTreeComp,normalize = TRUE,check.labels = FALSE)

#unlist the list for downstream analysis
Unlist_RF_List_fullyRandomTree <- unlist(RF_List_fullyRandomTree)

Unlist_RF_List_fullyRandomTree20 <- unlist(RF_List_RandomTree20)
Unlist_RF_List_fullyRandomTree40 <- unlist(RF_List_RandomTree40)
Unlist_RF_List_fullyRandomTree60 <- unlist(RF_List_RandomTree60)
Unlist_RF_List_fullyRandomTree80 <- unlist(RF_List_RandomTree80)

# Same thing for path distances 
P_List_fullyRandomTree <- foreach(i=1:length(treeList)) %do% path.dist(fullyRandomparseTreeComp,treeList[[i]],check.labels = FALSE)


p_List_RandomTree20 <- foreach(i=1:length(treeList)) %do% path.dist(fullyRandom20_1,treeList[[i]])
p_List_RandomTree40 <- foreach(i=1:length(treeList)) %do% path.dist(fullyRandom40_1,treeList[[i]])
p_List_RandomTree60 <- foreach(i=1:length(treeList)) %do% path.dist(fullyRandom60_1,treeList[[i]])
p_List_RandomTree80 <- foreach(i=1:length(treeList)) %do% path.dist(fullyRandom80_1,treeList[[i]])
#p_List_RandomTree20 <- foreach(i=1:length(treeList)) %do% SPRDist(fullyRandom20_1,treeList[[i]],check.labels = FALSE)


# #playing around2---
# 
# p_List_RandomTree20 <- path.dist(fullyRandom20_1,fullyRandomparseTreeComp,check.labels = FALSE)
# p_List_RandomTree40 <- path.dist(fullyRandom40_1,fullyRandomparseTreeComp,check.labels = FALSE)
# p_List_RandomTree60 <- path.dist(fullyRandom60_1,fullyRandomparseTreeComp,check.labels = FALSE)
# p_List_RandomTree80 <- path.dist(fullyRandom80_1,fullyRandomparseTreeComp,check.labels = FALSE)


#unlist the list for downstream analysis
Unlist_P_List_fullyRandomTree<- unlist(P_List_fullyRandomTree)

Unlist_P_List_fullyRandomTree20<- unlist(p_List_RandomTree20)
Unlist_P_List_fullyRandomTree40<- unlist(p_List_RandomTree40)
Unlist_P_List_fullyRandomTree60<- unlist(p_List_RandomTree60)
Unlist_P_List_fullyRandomTree80<- unlist(p_List_RandomTree80)
# Then I compare these 10 bias trees with the stratified complete full tree ("RAxML_bestTree.StratTree").  
# I call it fullyRandomparseTreeComp, a separate variable since I'm using it for all the analyses below
StratifiedparseTreeComp <- read.tree("RAxML_bestTree.FR100_new")




# For Robinson-Foulds values, we can use foreach loop to iterate through each tree (There were 40 trees)
RF_List_stratifiedTree <- foreach(i=1:length(treeList)) %do% RF.dist(StratifiedparseTreeComp,treeList[[i]],normalize = TRUE,check.labels = FALSE)
# RF_List_stratifiedTree 

RF_List_stratifiedTree20 <- foreach(i=1:length(treeList)) %do% RF.dist(fullyStrati20_1,treeList[[i]],normalize = TRUE)
RF_List_stratifiedTree40 <- foreach(i=1:length(treeList)) %do% RF.dist(fullyStrati40_1,treeList[[i]],normalize = TRUE)
RF_List_stratifiedTree60 <- foreach(i=1:length(treeList)) %do% RF.dist(fullyStrati60_1,treeList[[i]],normalize = TRUE)
RF_List_stratifiedTree80 <- foreach(i=1:length(treeList)) %do% RF.dist(fullyStrati80_1,treeList[[i]],normalize = TRUE)


#unlist the list for downstream analysis
Unlist_RF_List_stratifiedTree<- unlist(RF_List_stratifiedTree)

Unlist_RF_List_stratifiedTree20<- unlist(RF_List_stratifiedTree20)
Unlist_RF_List_stratifiedTree40<- unlist(RF_List_stratifiedTree40)
Unlist_RF_List_stratifiedTree60<- unlist(RF_List_stratifiedTree60)
Unlist_RF_List_stratifiedTree80<- unlist(RF_List_stratifiedTree80)


# Same thing for path distances 
P_List_stratifiedTree <- foreach(i=1:length(treeList)) %do% path.dist(StratifiedparseTreeComp,treeList[[i]],check.labels = FALSE)


P_List_stratifiedTree20 <- foreach(i=1:length(treeList)) %do% path.dist(fullyStrati20_1,treeList[[i]])
P_List_stratifiedTree40 <- foreach(i=1:length(treeList)) %do% path.dist(fullyStrati40_1,treeList[[i]])
P_List_stratifiedTree60 <- foreach(i=1:length(treeList)) %do% path.dist(fullyStrati60_1,treeList[[i]])
P_List_stratifiedTree80 <- foreach(i=1:length(treeList)) %do% path.dist(fullyStrati80_1,treeList[[i]])



#unlist the list for downstream analysis
Unlist_P_List_stratifiedTree<- unlist(P_List_stratifiedTree)

Unlist_P_List_stratifiedTree20<- unlist(P_List_stratifiedTree20)
Unlist_P_List_stratifiedTree40<- unlist(P_List_stratifiedTree40)
Unlist_P_List_stratifiedTree60<- unlist(P_List_stratifiedTree60)
Unlist_P_List_stratifiedTree80<- unlist(P_List_stratifiedTree80)

#I compare these 10 bias trees with the complete bias tree ("RAxML_bestTree.biasTree").  
# I call it fullyRandomparseTreeComp, a separate variable since I'm using it for all the analyses below
#biasTree <- read.tree("RAxML_bestTree.biasTree")

#Compare the 10 bias trees with  the complete original tree of FishTree of life dataset
biasTree  <- read.tree("Bias_new")


# #standardize the tree (the sum of branch length to 1)
# biasTree_new$edge.length <- biasTree$edge.length/sum(biasTree$edge.length)
# sum(biasTree_new$edge.length)
#biasTree<- compute.brlen(biasTree, method = "Grafen", power = 1)
# For Robinson-Foulds values, we can use foreach loop to iterate through each tree 
RF_List_biasTree <- foreach(i=1:length(treeList)) %do% RF.dist(biasTree,treeList[[i]],normalize = TRUE)
RF_List_biasTree

# Same thing for path distances 
P_List_biasTree <- foreach(i=1:length(treeList)) %do% path.dist(biasTree,treeList[[i]])
P_List_biasTree


#unlist the list for downstream analysis
Unlist_RF_List_biasTree<- unlist(RF_List_biasTree)

#unlist the list for downstream analysis
Unlist_P_List_biasTree<- unlist(P_List_biasTree)

#Newick to distance matrix. here I have used foreach loop to simplify the code 
#DM_List <- foreach(i=1:length(treeList)) %do% (cophenetic(treeList[[i]])/max(cophenetic(treeList[[i]]))) 
#All matrices must have the same number of objects to perform CADM. Here the matrices have different number of objects. So, no way of conducting CADM for bias trees
#Distance matrix of RAxML_bestTree.mgTree & RAxML_bestTree.StratTree


#create unstack data (dataframe) by combining Unlist_RF_List_fullyRandomTree and Unlist_RF_List_stratifiedTree to generate boxplot to see the difference between Random complete tree and dtratified complete tree based on biased samples
# RF <- data.frame(Unlist_RF_List_fullyRandomTree,Unlist_RF_List_stratifiedTree,Unlist_RF_List_biasTree)
RF1 <- data.frame(Unlist_RF_List_fullyRandomTree20,Unlist_RF_List_fullyRandomTree40,Unlist_RF_List_fullyRandomTree60,Unlist_RF_List_fullyRandomTree80,Unlist_RF_List_stratifiedTree20,Unlist_RF_List_stratifiedTree40,Unlist_RF_List_stratifiedTree60,Unlist_RF_List_stratifiedTree80,Unlist_RF_List_biasTree)
#Stack the dataframe with stack() to generate boxplot
stacked_RF <- stack(RF1)

#Generate the boxplot for stacked_RF

#ColorBrewer is a site that helps to choose different color palettes and there is a R package called RColorBrewer that helps to do this in R itself.
#install.packages("RColorBrewer")
library(RColorBrewer)
#brewer.pal.info
#display.brewer.all()

library(dplyr)
#Combine colors
#make my own color pallet  by grouping (Random samples in different intensities in blue, stratified samples in green, data-driven samples in violet)
#The following is just to get the color id
col1 <-  brewer.pal(11, "PiYG")
col2 <- brewer.pal(9, "Blues")
display.brewer.pal("all") 
# brewer.pal(11, "RdGy")
# display.brewer.pal(11, "RdGy")
#col_concat <- c(col1,col2,col3)
# display.brewer.all()
# colors_9 <- c("#6BAED6" ,"#4292C6", "#2171B5", "#084594","#B8E186", "#7FBC41", "#4D9221", "#276419","#e6550d")
#colors_9 <- c("#6BAED6" ,"#4292C6", "#2171B5", "#084594","#fc9272","#D6604D","#B2182B" , "#67001F","#01665E")
colors_9 <- c("#6BAED6" ,"#4292C6", "#2171B5", "#084594","#F1B6DA","#DE77AE","#C51B7D","#8E0152","#E6F5D0")

#testing
boxplot(values~ind, data = stacked_RF, outline=FALSE,
        xlab = "Sampling Method",
        ylab= "RF distance",
        ylim = range(stacked_RF$values),
        main= "RF distance difference among Random sampling, Stratified sampling and data-driven sampling",
        col= colors_9,
        names= c("20% Random","40% Random", "60% Random", "80% Random","20% Stratified","40% Stratified", "60% Stratified", "80% Stratified", "Data-driven sampling"))

# Points
stripchart(values ~ ind,
           data = stacked_RF,
           method = "jitter",
           pch = 19,
           col = "black",
           vertical = TRUE,
           add = TRUE)

#create unstack data (dataframe) by combining Unlist_P_List_fullyRandomTree and Unlist_P_List_stratifiedTree to generate boxplot to see the difference between Random complete tree and dtratified complete tree based on biased samples
# Path <- data.frame(Unlist_P_List_fullyRandomTree,Unlist_P_List_stratifiedTree)
Path1 <- data.frame(Unlist_P_List_fullyRandomTree20,Unlist_P_List_fullyRandomTree40,Unlist_P_List_fullyRandomTree60,Unlist_P_List_fullyRandomTree80,Unlist_P_List_stratifiedTree20,Unlist_P_List_stratifiedTree40,Unlist_P_List_stratifiedTree60,Unlist_P_List_stratifiedTree80,Unlist_P_List_biasTree)

#Stack the dataframe with stack() to generate boxplot
stacked_Path <- stack(Path1)

#Generate the boxplot for stacked_RF

# boxplot(values~ind, data = stacked_Path,
#         xlab = "Sampling Method",
#         ylab= "Path distance",
#         main= "Path distance difference between Random sampling tree and Stratified sampling tree based on biased samples",
#         col= c("light blue","purple"),
#         names= c("Complete Random Tree", "Complete Stratified Tree"))
boxplot(values~ind, data = stacked_Path,outline=FALSE,
        xlab = "Sampling Method",
        ylab= "Path distance",
        #ylim = rev(range(stacked_Path$values)),
        ylim = range(stacked_Path$values),
        main= "Path distance difference among Random sampling, Stratified sampling and data-driven sampling",
        col= colors_9,
        names= c("20% Random","40% Random", "60% Random", "80% Random","20% Stratified","40% Stratified", "60% Stratified", "80% Stratified","Data-driven sampling"))

# Points
stripchart(values ~ ind,
           data = stacked_Path,
           method = "jitter",
           pch = 19,
           col = "black",
           vertical = TRUE,
           add = TRUE)



#----------------
library(dplyr)
#Combine colors
#make my own color pallet  by grouping (Random samples in different intensities in blue, stratified samples in green, data-driven samples in violet)
#The following is just to get the color id
col1 <-  brewer.pal(11, "PiYG")
col2 <- brewer.pal(9, "Blues")

#col_concat <- c(col1,col2,col3)

colors_9 <- c("#6BAED6" ,"#4292C6", "#2171B5", "#084594","#F1B6DA","#DE77AE","#C51B7D","#8E0152","#B8E186")



#Newick to distance matrix. here I have used foreach loop to simplify the code 
DM_List <- foreach(i=1:length(treeList)) %do% (cophenetic(treeList[[i]])/max(cophenetic(treeList[[i]]))) 

#Distance matrix of 100% tree
DM_100 <- cophenetic(fullyRandomparseTreeComp)/max(cophenetic(fullyRandomparseTreeComp))  

#function to generate CADM results
CADM <- function(DM){
  #add two matrics together
  AB <- rbind(DM_100, DM)
  cadm_global <- CADM.global(AB, 2, 4520)$congruence_analysis[1,]
  return(cadm_global)
}

#Use CADM function to compare backbone trees with the complete tree (foreach loop to iterate through CADM list)
CADM_List <- foreach(i=1:length(treeList)) %do% CADM(DM_List[[i]])

#Using the metric(), generate CADM dataframe
CADM_metric <- metric(CADM_List)

# Same thing for path distances 
P_List_stratifiedTree <- foreach(i=1:length(treeList)) %do% path.dist(StratifiedparseTreeComp,treeList[[i]],check.labels = FALSE)


P_List_stratifiedTree20 <- foreach(i=1:length(treeList)) %do% path.dist(fullyStrati20_1,treeList[[i]])
P_List_stratifiedTree40 <- foreach(i=1:length(treeList)) %do% path.dist(fullyStrati40_1,treeList[[i]])
P_List_stratifiedTree60 <- foreach(i=1:length(treeList)) %do% path.dist(fullyStrati60_1,treeList[[i]])
P_List_stratifiedTree80 <- foreach(i=1:length(treeList)) %do% path.dist(fullyStrati80_1,treeList[[i]])



#unlist the list for downstream analysis
Unlist_P_List_stratifiedTree<- unlist(P_List_stratifiedTree)

Unlist_P_List_stratifiedTree20<- unlist(P_List_stratifiedTree20)
Unlist_P_List_stratifiedTree40<- unlist(P_List_stratifiedTree40)
Unlist_P_List_stratifiedTree60<- unlist(P_List_stratifiedTree60)
Unlist_P_List_stratifiedTree80<- unlist(P_List_stratifiedTree80)

















--------------------------------
  # Then reorder the lists such that the ordering goes 20_1, 40_1, 60_1, 80_1, 20_2, 40_2, 60_2, 80_2 etc. 
  order <- c(2, 12, 22, 32,# 20_1, 40_1, 60_1, 80_1
             3, 13, 23, 33, # 20_2, 40_2, 60_2, 80_2
             4, 14, 24, 34, 
             5, 15, 25, 35,
             6, 16, 26, 36,
             7, 17, 27, 37,
             8, 18, 28, 38,
             9, 19, 29, 39,
             10, 20, 30, 40,
             1, 11, 21, 31) # 20_10, 40_10, 60_10, 80_10

# When plotting we need the percentage of backbone tree as the x axis. For that we can use,
GRF <- c(20, 40, 60, 80)

# Here I have used a function to reorder the metric lists, group the values and assign values to variable names.
# From metric list (L_list) I'm going to group metric values (Robinson-Foulds values, CADM values, path distances etc.) (eg; 20_1, 40_1, 60_1, 80_1 as L_1). Then unlist it to get a vector and assign variable names for down stream analysis

metric <- function(L_List){
  #reordering
  L_List <- L_List[order]
  #group the values and assign values to variable names
  L_1 <- unlist(L_List[1:4])
  L_2 <- unlist(L_List[5:8])
  L_3 <- unlist(L_List[9:12])
  L_4 <- unlist(L_List[13:16])
  L_5 <- unlist(L_List[17:20])
  L_6 <- unlist(L_List[21:24])
  L_7 <- unlist(L_List[25:28])
  L_8 <- unlist(L_List[29:32])
  L_9 <- unlist(L_List[33:36])
  L_10 <- unlist(L_List[37:40])
  #generate a dataframe using backbone % and metric values
  results <- data.frame(GRF,L_1,L_2,L_3,L_4,L_5,L_6,L_7,L_8,L_9,L_10)
}

# Now we can use the above function to generate dataframes in different metrics
# First, Robinson-Foulds metric,
RF_metric<- metric(RF_List)

# Plot Robinson-Fould distances 
ggplot(RF_metric , aes(GRF)) +
  scale_y_reverse() + xlab("The level of backbone(%)")+ylab("RF Distance")+labs(title="The plot of RF distance vs. the level of backbone")+
  geom_line(aes(y=L_1),
            colour="red") +
  geom_line(aes(y=L_2),
            colour="green") +
  geom_line(aes(y=L_3),
            colour="blue") +
  geom_line(aes(y=L_4),
            colour="yellow") +
  geom_line(aes(y=L_5),
            colour="orange") +
  geom_line(aes(y=L_6),
            colour="black")+
  geom_line(aes(y=L_7),
            colour="purple") +
  geom_line(aes(y=L_8),
            colour="pink") +
  geom_line(aes(y=L_9),
            colour="brown") +
  geom_line(aes(y=L_10),
            colour="violet")

