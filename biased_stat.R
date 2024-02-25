#Statistical test results:

#This code is to check the results of phylogenetic metrics using statistical tests.

#Coded by M. A. Thanuja M. Fernando

#Last updated on Sep 06, 2023

#Acknowledgement: Idea of using "foreach package" in this  code by Matthew Orton 

#Load packages
#install.packages("foreach")
library(foreach)
library(phangorn)
library(reshape) #melt
library(ez)
library(ggpubr)
library(rstatix)#outliers
library(lme4)
library(rmcorr)
#install.packages("lmerTest")
library(lmerTest)


#Read trees generated from RAxML into R. Read all of the trees at once ensuring I have the correct path to our small_dataset file. Here the path is /home/thanu/Desktop/FishData/small_dataset and the pattern is a regular expression to extract out the files I need

files <- list.files(path="/home/thanu/Desktop/FishData/New_dataset/BIASED", pattern="B[0-9]{1,2}_new", full.names=TRUE, recursive=FALSE)
#foreach package to just iterate through each file
treeList <- foreach(i=1:length(files)) %do% read.tree(files[i])

#Now I have a list with each element being its own separate tree, so treeList[[1]] is my first tree 20_1_new, treeList[[2]] is 20_2_new etc. Now I name each list element by substituting out the path from the names of each file

names(treeList) <- gsub("/home/thanu/Desktop/FishData/New_dataset/BIASED/", "", files)

#This will let us keep track of which list element corresponds to which tree file

names(treeList)

#Then I keep the random, stratified and dat-driven samples as ("RAxML_bestTree.complete100") separate variables since I'm using it for all the analyses below. 

parseTreeComp <- read.tree("RAxML_bestTree.FR100_new")
fullyRandom20_1 <- read.tree("20_1_new")
fullyRandom40_1 <- read.tree("40_1_new")
fullyRandom60_1 <- read.tree("60_1_new")
fullyRandom80_1 <- read.tree("80_1_new")

fullyStrati20_1 <- read.tree("S20_1_new")
fullyStrati40_1 <- read.tree("S40_1_new")
fullyStrati60_1 <- read.tree("S60_1_new")
fullyStrati80_1 <- read.tree("S80_1_new")

#Compare the 10 bias trees with  the complete original tree of FishTree of life dataset
biasTree  <- read.tree("Bias_new")

#For Robinson-Foulds values, I use foreach loop to iterate through each tree (There were 40 trees)

RF_List <- foreach(i=1:length(treeList)) %do% RF.dist(parseTreeComp,treeList[[i]],normalize = TRUE)
RF_metric <- unlist(RF_List)
class(RF_metric)
RF_metric <- data.frame(RF_metric)
RF_metric <- stack(RF_metric)

# For Robinson-Foulds values, we can use foreach loop to iterate through each tree
RF_List_RandomTree20 <- foreach(i=1:length(treeList)) %do% RF.dist(fullyRandom20_1,treeList[[i]],normalize = TRUE)
RF_List_RandomTree40 <- foreach(i=1:length(treeList)) %do% RF.dist(fullyRandom40_1,treeList[[i]],normalize = TRUE)
RF_List_RandomTree60 <- foreach(i=1:length(treeList)) %do% RF.dist(fullyRandom60_1,treeList[[i]],normalize = TRUE)
RF_List_RandomTree80 <- foreach(i=1:length(treeList)) %do% RF.dist(fullyRandom80_1,treeList[[i]],normalize = TRUE)

Unlist_RF_List_fullyRandomTree20 <- unlist(RF_List_RandomTree20)
Unlist_RF_List_fullyRandomTree40 <- unlist(RF_List_RandomTree40)
Unlist_RF_List_fullyRandomTree60 <- unlist(RF_List_RandomTree60)
Unlist_RF_List_fullyRandomTree80 <- unlist(RF_List_RandomTree80)


# Same thing for path distances 

p_List_RandomTree20 <- foreach(i=1:length(treeList)) %do% path.dist(fullyRandom20_1,treeList[[i]])
p_List_RandomTree40 <- foreach(i=1:length(treeList)) %do% path.dist(fullyRandom40_1,treeList[[i]])
p_List_RandomTree60 <- foreach(i=1:length(treeList)) %do% path.dist(fullyRandom60_1,treeList[[i]])
p_List_RandomTree80 <- foreach(i=1:length(treeList)) %do% path.dist(fullyRandom80_1,treeList[[i]])

Unlist_P_List_fullyRandomTree20<- unlist(p_List_RandomTree20)
Unlist_P_List_fullyRandomTree40<- unlist(p_List_RandomTree40)
Unlist_P_List_fullyRandomTree60<- unlist(p_List_RandomTree60)
Unlist_P_List_fullyRandomTree80<- unlist(p_List_RandomTree80)

# RF_List_stratifiedTree 

RF_List_stratifiedTree20 <- foreach(i=1:length(treeList)) %do% RF.dist(fullyStrati20_1,treeList[[i]],normalize = TRUE)
RF_List_stratifiedTree40 <- foreach(i=1:length(treeList)) %do% RF.dist(fullyStrati40_1,treeList[[i]],normalize = TRUE)
RF_List_stratifiedTree60 <- foreach(i=1:length(treeList)) %do% RF.dist(fullyStrati60_1,treeList[[i]],normalize = TRUE)
RF_List_stratifiedTree80 <- foreach(i=1:length(treeList)) %do% RF.dist(fullyStrati80_1,treeList[[i]],normalize = TRUE)

Unlist_RF_List_stratifiedTree20<- unlist(RF_List_stratifiedTree20)
Unlist_RF_List_stratifiedTree40<- unlist(RF_List_stratifiedTree40)
Unlist_RF_List_stratifiedTree60<- unlist(RF_List_stratifiedTree60)
Unlist_RF_List_stratifiedTree80<- unlist(RF_List_stratifiedTree80)


P_List_stratifiedTree20 <- foreach(i=1:length(treeList)) %do% path.dist(fullyStrati20_1,treeList[[i]])
P_List_stratifiedTree40 <- foreach(i=1:length(treeList)) %do% path.dist(fullyStrati40_1,treeList[[i]])
P_List_stratifiedTree60 <- foreach(i=1:length(treeList)) %do% path.dist(fullyStrati60_1,treeList[[i]])
P_List_stratifiedTree80 <- foreach(i=1:length(treeList)) %do% path.dist(fullyStrati80_1,treeList[[i]])

Unlist_P_List_stratifiedTree20<- unlist(P_List_stratifiedTree20)
Unlist_P_List_stratifiedTree40<- unlist(P_List_stratifiedTree40)
Unlist_P_List_stratifiedTree60<- unlist(P_List_stratifiedTree60)
Unlist_P_List_stratifiedTree80<- unlist(P_List_stratifiedTree80)

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


#RF metric
#Now we can use the above function to generate dataframes in different metrics
#First, Robinson-Foulds metric,
#transpose the dataframe
trans_RF_metric_20R <- data.frame(Unlist_RF_List_fullyRandomTree20)
trans_RF_metric_40R <- data.frame(Unlist_RF_List_fullyRandomTree40)
trans_RF_metric_60R <- data.frame(Unlist_RF_List_fullyRandomTree60)
trans_RF_metric_80R <- data.frame(Unlist_RF_List_fullyRandomTree80)
trans_RF_metric_20S <- data.frame(Unlist_RF_List_stratifiedTree20)
trans_RF_metric_40S <- data.frame(Unlist_RF_List_stratifiedTree40)
trans_RF_metric_60S <- data.frame(Unlist_RF_List_stratifiedTree60)
trans_RF_metric_80S <- data.frame(Unlist_RF_List_stratifiedTree80)
trans_RF_metric_Bias <- data.frame(Unlist_RF_List_biasTree)
#delete the 1st row of the dataframe
# trans_RF_metric_20R = trans_RF_metric_20R[-1,]
# dim(trans_RF_metric_20R)
id <- c(1:10)
RF_metric_Random <- cbind(id,trans_RF_metric_20R,trans_RF_metric_40R,trans_RF_metric_60R,trans_RF_metric_80R)
RF_metric_Strat <- cbind(id,trans_RF_metric_20S,trans_RF_metric_40S,trans_RF_metric_60S,trans_RF_metric_80S)
#playing AROUND
RF_metric_RSB <- cbind(id,trans_RF_metric_20R,trans_RF_metric_40R,trans_RF_metric_60R,trans_RF_metric_80R,trans_RF_metric_20S,trans_RF_metric_40S,trans_RF_metric_60S,trans_RF_metric_80S,trans_RF_metric_Bias)
colnames(RF_metric_RSB) <-c("id","0.2R","0.4R","0.6R","0.8R","0.2S","0.4S","0.6S","0.8S","Biased") 
RF_metric_RSB
#Melt the dataframe into a form suitable for easy casting
longdata_RSB <- melt(RF_metric_RSB,
                 id="id",
                 measure=c("0.2R","0.4R","0.6R","0.8R","0.2S","0.4S","0.6S","0.8S","Biased"))
#assign column names
colnames(longdata_RSB)=c("id","sample_type","RF_distance")

# #delete the 1st row of the dataframe
# RF_metric_Random  = RF_metric_Random [-1,]
# #adding unique id
# RF_metric_Random$id <- cumsum(!duplicated(RF_metric_Random[1:2]))
colnames(RF_metric_Random) <- c("id","0.2","0.4","0.6","0.8")
colnames(RF_metric_Strat) <- c("id","0.2","0.4","0.6","0.8")
dim(RF_metric_Random)
class(RF_metric_Random)
#Melt the dataframe into a form suitable for easy casting
longdata <- melt(RF_metric_Random,
                 id="id",
                 measure=c("0.2","0.4","0.6","0.8"))
longdata_strat <- melt(RF_metric_Strat,
                 id="id",
                 measure=c("0.2","0.4","0.6","0.8"))
#assign column names
colnames(longdata)=c("id","sample_completeness","RF_distance")
colnames(longdata_strat)=c("id","sample_completeness","RF_distance")
longdata
longdata_strat
#ANOVA test. dv means dependant variable, The one-way repeated measures ANOVA can be used to determine whether the means of longdata RF distances are significantly different among four backbone levels(sample completeness 20%,40%,60%,80%)
ezANOVA(data = longdata, dv = RF_distance, wid = id, within = sample_completeness)
ezANOVA(data = longdata_strat, dv = RF_distance, wid = id, within = sample_completeness)
ezANOVA(data = longdata_RSB, dv = RF_distance, wid = id, within = sample_type)
library(ggpubr)
ggboxplot(longdata_RSB,x="sample_type",y="RF_distance", color = "sample_type",palette = c("#6BAED6" ,"#4292C6", "#2171B5", "#084594","#F1B6DA","#DE77AE","#C51B7D","#8E0152","#B8E186")
,main="Boxplot of RF distances based on sample types")
#boxplot(RF_distance~sample_type, data= longdata_RSB,main="Boxplot of RF distances based on sample types")
RF_anova <- aov(RF_distance~sample_type, data= longdata_RSB)
summary(RF_anova)
#Tukey HSD to perfoem pair-wise-comparison between the means of groups.
#The TukeyHD() takes the fitted ANOVA as an argument
TukeyHSD(RF_anova)
#if the 0 lie between lwr and upr end points of the confidence interval at 95% then that groups differences are significant with an adjusted p value of (adjusted p val is mention in the output)
#Check ANOVA assumptions
#The ANOVA test assumes that the data are normally dsitributed and the variance across groups are homogeneous
#We can check this with some diagnostic plots
#The residual vs. fits plotds acn be used to check the homogeneity of variences
par(mfrow=c(2,2))
plot(RF_anova)
#In this homogeneity of variances plot for RF distance, has no evident relationships between residuas and fitted values(the mean of each groups), which is good.
#So, we can assume the homogeneity of variances

#It is alos possible to use Bartlett's test or Levene's test (recommend as it is less sensitive to departures from normal distribution) to check the homogeneity of variances
library(car)
leveneTest(RF_distance~sample_type, data= longdata_RSB)
#If the p value is not lesss than the significance level of 0.05, then no evidence to suggest that the variance across the groups is statistically significant.
#Therefore, we can assume the homogeneity of variances in the different treatmet groups.

#One-way ANOVA test requiresequal variances for all groups.
#If the Levene test is not significant, then to save ANOVAtest, in a situation where the homogeneity of variances as assumption is violated,
#we can use an alternative procedure (Welch one-way test, that doesn't require that assumption have been implementd in the function oneway.test()
#Anova test with no assumption of equal variances,
oneway.test(longdata_RSB$RF_distance~longdata_RSB$sample_type)
# or pairwise t-test with no assumption of equal variances
#pairwise.t.test(longdata_RSB$RF_distance~longdata_RSB$sample_type),p.adjust.method = "BH",pool.sd = F)
#Check the homogeneiety of variance assumption

#Check the norma;lity assumption with norality plot of residuals
#In the given plot, the quantiles of the residuals are plotted against the quantiles of the normal distribution.
#A 45 degree refernce ine is also plotted
#THe normal probability plot of residuals is used to check the assumption tat the residuals are normally distributed.
#It should  approximately follow a straight line
#As all the points fall approximately along this reference line we can assume normality.

#the conclusion above is supported by the shapiro-WIlk test on the ANOVA residuals (W = 0.95117, p-value = 0.001979) normality violated or not?violated??
#extract the residuals
RF_aov.residuals <- residuals(object=RF_anova)
#run shapiro-wilk test
shapiro.test(x=RF_aov.residuals)

#Non-parametric alternative to one-way ANOVA test is Kruskal-Wallis rank sum test, which can be used when ANOVA assumptions are not met

kruskal.test(longdata_RSB$RF_distance~longdata_RSB$sample_type)


#assumptions
#Outliers can be easily identified using box plot methods, implemented in the R function identify_outliers() [rstatix package].

longdata %>%
  group_by(sample_completeness) %>%
  identify_outliers(RF_distance)

longdata_strat %>%
  group_by(sample_completeness) %>%
  identify_outliers(RF_distance)

#visualization
#Create a violin plot
#violin plot with dot plot (geom_dotplot())
violin_plot <- ggplot(longdata, aes(x=sample_completeness, y=RF_distance,color=sample_completeness))+ geom_violin(trim = FALSE)+theme(legend.position="none")+geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5,binwidth = 0.01)+labs(title="The box plot of RF distance for sampling completeness")
violin_plot
violin_plot <- ggplot(longdata_strat, aes(x=sample_completeness, y=RF_distance,color=sample_completeness))+ geom_violin(trim = FALSE)+theme(legend.position="none")+geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5,binwidth = 0.01)+labs(title="The box plot of RF distance for sampling completeness")
violin_plot
#Normality assumption
#The normality assumption can be checked by computing Shapiro-Wilk test for each sample completeness. If the data is normally distributed, the p-value should be greater than 0.05.

longdata %>%
  group_by(sample_completeness) %>%
  shapiro_test(RF_distance)
longdata_strat %>%
  group_by(sample_completeness) %>%
  shapiro_test(RF_distance)

#QQ plot draws the correlation between a given data and the normal distribution. Create QQ plots for each level backbone tree(sample completeness)

ggqqplot(longdata, "RF_distance", facet.by = "sample_completeness",color ="sample_completeness")
ggqqplot(longdata_strat, "RF_distance", facet.by = "sample_completeness",color ="sample_completeness")
#From the resultant plots, as all the points fall approximately along the reference line, we can assume normality

#The assumption of sphericity will be automatically checked during the computation of the ANOVA test using the R function anova_test() [rstatix package]. The Mauchly’s test is internally used to assess the sphericity assumption

#By using the function get_anova_table() [rstatix] to extract the ANOVA table, the Greenhouse-Geisser sphericity correction is automatically applied to factors violating the sphericity assumption

res.aov <- anova_test(data = longdata, dv = RF_distance, wid = id, within = sample_completeness)
res.aov_strat <- anova_test(data = longdata_strat, dv = RF_distance, wid = id, within = sample_completeness)
get_anova_table(res.aov)
get_anova_table(res.aov_strat)
#Post-hoc tests

#perform multiple pairwise paired t-tests between the levels of the within-subjects factor (Sample_completeness). P-values are adjusted using the Bonferroni multiple testing correction method.

# pairwise comparisons

pwc <- longdata %>%
  pairwise_t_test(
    RF_distance ~ sample_completeness, paired = TRUE,
    p.adjust.method = "bonferroni"
  )
pwc

pwc <- longdata_strat %>%
  pairwise_t_test(
    RF_distance ~ sample_completeness, paired = TRUE,
    p.adjust.method = "bonferroni"
  )
pwc


#Correlation tests need numerical variables, so convert more than one column in R data frame to numeric values
longdat <- lapply(longdata,as.numeric)
longdat_strat <- lapply(longdata_strat,as.numeric)
#Calculate the repeated measures correlation coefficient. id is variable giving the subject name/id for each observation. RF_distance is a numeric variable giving the observations for one measure, sample_completeness is a numeric variable giving the observations for the second measure, longdat is the data frame containing the variables.
rmcorr(id, RF_distance, sample_completeness, longdat)
rmcorr(id, RF_distance, sample_completeness, longdat_strat)


#PathDistance


#Now we can use the above function to generate dataframes in different metrics

#Get the dataframe
trans_P_metric_20R <- data.frame(Unlist_P_List_fullyRandomTree20)
trans_P_metric_40R <- data.frame(Unlist_P_List_fullyRandomTree40)
trans_P_metric_60R <- data.frame(Unlist_P_List_fullyRandomTree60)
trans_P_metric_80R <- data.frame(Unlist_P_List_fullyRandomTree80)
trans_P_metric_20S <- data.frame(Unlist_P_List_stratifiedTree20)
trans_P_metric_40S <- data.frame(Unlist_P_List_stratifiedTree40)
trans_P_metric_60S <- data.frame(Unlist_P_List_stratifiedTree60)
trans_P_metric_80S <- data.frame(Unlist_P_List_stratifiedTree80)
#Assigned an id column to the dataframe
id <- c(1:10)
P_metric_Random <- cbind(id,trans_P_metric_20R,trans_P_metric_40R,trans_P_metric_60R,trans_P_metric_80R)
P_metric_Strati <- cbind(id,trans_P_metric_20S,trans_P_metric_40S,trans_P_metric_60S,trans_P_metric_80S)
# #delete the 1st row of the dataframe
# RF_metric_Random  = RF_metric_Random [-1,]
# #adding unique id
# RF_metric_Random$id <- cumsum(!duplicated(RF_metric_Random[1:2]))
colnames(P_metric_Random) <- c("id","0.2","0.4","0.6","0.8")
colnames(P_metric_Strati) <- c("id","0.2","0.4","0.6","0.8")
dim(P_metric_Random)
class(P_metric_Random)
#Melt the dataframe into a form suitable for easy casting
longdata <- melt(P_metric_Random,
                 id="id",
                 measure=c("0.2","0.4","0.6","0.8"))
longdata_Strati <- melt(P_metric_Strati,
                 id="id",
                 measure=c("0.2","0.4","0.6","0.8"))

colnames(longdata)=c("id","sample_completeness","Path_Distance")
colnames(longdata_Strati)=c("id","sample_completeness","Path_Distance")
longdata
longdata_Strati

#The one-way repeated measures ANOVA can be used to determine whether the means of longdata path distances are significantly different among four backbone levels(sample completeness 20%,40%,60%,80%)
ezANOVA(data = longdata, dv = Path_Distance, wid = id, within = sample_completeness)
ezANOVA(data = longdata_Strati, dv = Path_Distance, wid = id, within = sample_completeness)
#check assumptions
#Outliers can be easily identified using box plot methods, implemented in the R function identify_outliers() [rstatix package]
longdata %>%
  group_by(sample_completeness) %>%
  identify_outliers(Path_Distance)

longdata_Strati %>%
  group_by(sample_completeness) %>%
  identify_outliers(Path_Distance)
#There were no extreme outliers

#visualization
#Create a violin plot
violin_plot <- ggplot(longdata, aes(x=sample_completeness, y=Path_Distance,color=sample_completeness))+ geom_violin(trim = FALSE)+theme(legend.position="none")+geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5,binwidth = 0.01)+labs(title="The box plot of Path Distance for sampling completeness")
violin_plot

violin_plot <- ggplot(longdata_Strati, aes(x=sample_completeness, y=Path_Distance,color=sample_completeness))+ geom_violin(trim = FALSE)+theme(legend.position="none")+geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5,binwidth = 0.01)+labs(title="The box plot of Path Distance for sampling completeness")
violin_plot

#Normality assumption
#The normality assumption can be checked by computing Shapiro-Wilk test for each time point. If the data is normally distributed, the p-value should be greater than 0.05.
longdata %>%
  group_by(sample_completeness) %>%  
  shapiro_test(Path_Distance)
longdata_Strati %>%
  group_by(sample_completeness) %>%  
  shapiro_test(Path_Distance)

#The longdata Path_distance was normally distributed at each sample completeness, as assessed by Shapiro-Wilk’s test (p > 0.05).

#QQ plot draws the correlation between a given data and the normal distribution. Create QQ plots for each level backbone tree(sample completeness)
ggqqplot(longdata, "Path_Distance", facet.by = "sample_completeness",color = "sample_completeness")
ggqqplot(longdata_Strati, "Path_Distance", facet.by = "sample_completeness",color = "sample_completeness")
#From the resultant plots, as all the points fall approximately along the reference line, we can assume normality

#The assumption of sphericity will be automatically checked during the computation of the ANOVA test using the R function anova_test() [rstatix package]. The Mauchly’s test is internally used to assess the sphericity assumption
#By using the function get_anova_table() [rstatix] to extract the ANOVA table, the Greenhouse-Geisser sphericity correction is automatically applied to factors violating the sphericity assumption
# res.aov <- anova_test(data = longdata, dv = RF_distance, wid = id,within = sample_completeness)
# get_anova_table(res.aov)

#Post-hoc tests
#perform multiple pairwise paired t-tests between the levels of the within-subjects factor (Sample_completeness). P-values are adjusted using the Bonferroni multiple testing correction method.

# pairwise comparisons
pwc <- longdata %>%
  pairwise_t_test(
    Path_Distance ~ sample_completeness, paired = TRUE,
    p.adjust.method = "bonferroni"
  )
pwc
pwc <- longdata_Strati %>%
  pairwise_t_test(
    Path_Distance ~ sample_completeness, paired = TRUE,
    p.adjust.method = "bonferroni"
  )
pwc
#All the pairwise differences are statistically significant

#Correlation tests need numerical variables, so convert more than one column in R data frame to numeric values
longdat <- lapply(longdata,as.numeric)
longdat_Strati <- lapply(longdata_Strati,as.numeric)
#Calculate the repeated measures correlation coefficient. id is variable giving the subject name/id for each observation. RF_distance is a numeric variable giving the observations for one measure, sample_completeness is a numeric variable giving the observations for the second measure, longdat is the data frame containing the variables.
rmcorr(id, Path_Distance, sample_completeness, longdat)
rmcorr(id, Path_Distance, sample_completeness, longdat_Strati)
###Simplified function abort R, so let's try one by one
library(phangorn)
#Read trees into R generated from RAxML
aF <- read.tree("RAxML_bestTree.FR100_new")
bF <- read.tree("20_1_new")   
cF <- read.tree("20_2_new") 
dF <- read.tree("20_3_new") 
eF <- read.tree("20_4_new") 
fF <- read.tree("20_5_new") 
gF <- read.tree("20_6_new")   
hF <- read.tree("20_7_new")     
iF <- read.tree("20_8_new") 
jF <- read.tree("20_9_new") 
kF <- read.tree("20_10_new") 
lF <- read.tree("40_1_new") 
mF <- read.tree("40_2_new") 
nF <- read.tree("40_3_new") 
oF <- read.tree("40_4_new") 
pF <- read.tree("40_5_new") 
qF <- read.tree("40_6_new") 
rF <- read.tree("40_7_new") 
sF <- read.tree("40_8_new") 
tF <- read.tree("40_9_new") 
uF <- read.tree("40_10_new") 
vF <- read.tree("60_1_new") 
wF <- read.tree("60_2_new") 
xF <- read.tree("60_3_new") 
yF <- read.tree("60_4_new") 
zF <- read.tree("60_5_new") 
aaF <- read.tree("60_6_new") 
bbF <- read.tree("60_7_new") 
ccF <- read.tree("60_8_new")
ddF <- read.tree("60_9_new")
eeF <- read.tree("60_10_new")
ffF <- read.tree("80_1_new")
ggF <- read.tree("80_2_new")
hhF <- read.tree("80_3_new")
iiF <- read.tree("80_4_new")
jjF <- read.tree("80_5_new")
kkF <- read.tree("80_6_new")
llF <- read.tree("80_7_new")
mmF <- read.tree("80_8_new")
nnF <- read.tree("80_9_new")
ooF <- read.tree("80_10_new")

library(ape)
#Newick to distancematrix
DM_100F <- cophenetic(aF)
DM_20_1F <- cophenetic(bF)
DM_20_2F <- cophenetic(cF)
DM_20_3F <- cophenetic(dF)
DM_20_4F <- cophenetic(eF)
DM_20_5F <- cophenetic(fF)
DM_20_6F <- cophenetic(gF)
DM_20_7F <- cophenetic(hF)
DM_20_8F <- cophenetic(iF)
DM_20_9F <- cophenetic(jF)
DM_20_10F <- cophenetic(kF)
DM_40_1F <- cophenetic(lF)
DM_40_2F <- cophenetic(mF)
DM_40_3F <- cophenetic(nF)
DM_40_4F <- cophenetic(oF)
DM_40_5F <- cophenetic(pF)
DM_40_6F <- cophenetic(qF)
DM_40_7F <- cophenetic(rF)
DM_40_8F <- cophenetic(sF)
DM_40_9F <- cophenetic(tF)
DM_40_10F <- cophenetic(uF)
DM_60_1F <- cophenetic(vF)
DM_60_2F <- cophenetic(wF)
DM_60_3F <- cophenetic(xF)
DM_60_4F <- cophenetic(yF)
DM_60_5F <- cophenetic(zF)
DM_60_6F <- cophenetic(aaF)
DM_60_7F <- cophenetic(bbF)
DM_60_8F <- cophenetic(ccF)
DM_60_9F <- cophenetic(ddF)
DM_60_10F <- cophenetic(eeF)
DM_80_1F <- cophenetic(ffF)
DM_80_2F <- cophenetic(ggF)
DM_80_3F <- cophenetic(hhF)
DM_80_4F <- cophenetic(iiF)
DM_80_5F <- cophenetic(jjF)
DM_80_6F <- cophenetic(kkF)
DM_80_7F <- cophenetic(llF)
DM_80_8F <- cophenetic(mmF)
DM_80_9F <- cophenetic(nnF)
DM_80_10F <- cophenetic(ooF)
#function to generate CADM results
CADM <- function(DM){
  #add two matrics together
  AB <- rbind(DM_100F, DM)
  cadm_global <- CADM.global(AB, 2, 6796)$congruence_analysis[1,]
  return(cadm_global)
}
cadm_global <- CADM.global(AB, 2, 6796)
#cadm_AB <- CADM(DM_B)
#cadm_AB$congruence_analysis
#Use CADM function to compare backbone trees with the complete tree
CADM(DM_20_1F)
CADM(DM_20_2)
CADM(DM_20_3)
CADM(DM_20_4)
CADM(DM_20_5)
CADM(DM_20_6)
CADM(DM_20_7)
CADM(DM_20_8)
CADM(DM_20_9)
CADM(DM_20_10)
CADM(DM_40_1)
CADM(DM_40_2)
CADM(DM_40_3)
CADM(DM_40_4)
CADM(DM_40_5)
CADM(DM_40_6)
CADM(DM_40_7)
CADM(DM_40_8)
CADM(DM_40_9)
CADM(DM_40_10)
CADM(DM_60_1)
CADM(DM_60_2)
CADM(DM_60_3)
CADM(DM_60_4)
CADM(DM_60_5)
CADM(DM_60_6)
CADM(DM_60_7)
CADM(DM_60_8)
CADM(DM_60_9)
CADM(DM_60_10)
CADM(DM_80_1)
CADM(DM_80_2)
CADM(DM_80_3)
CADM(DM_80_4)
CADM(DM_80_5)
CADM(DM_80_6)
CADM(DM_80_7)
CADM(DM_80_8)
CADM(DM_80_9)
CADM(DM_80_10)













#CADM
#Newick to distance matrix. here I have used foreach loop to simplify the code 
DM_List <- foreach(i=1:length(treeList)) %do% (cophenetic(treeList[[i]])/max(cophenetic(treeList[[i]]))) 

#Distance matrix of 100% tree
DM_100 <- cophenetic(parseTreeComp)/max(cophenetic(parseTreeComp))  

#function to generate CADM results
CADM <- function(DM){
  #add two matrices together
  AB <- rbind(DM_100F, DM)
  cadm_global <- CADM.global(AB, 2, 4520)$congruence_analysis[1,]
  return(cadm_global)
}

#Use CADM function to compare backbone trees with the complete tree (foreach loop to iterate through CADM list)
CADM_List <- foreach(i=1:length(treeList)) %do% CADM(DM_List[[i]])

#Using the metric(), generate CADM dataframe
CADM_metric <- metric(CADM_List)

#transpose the dataframe
trans_CADM_metric <- data.frame(t(CADM_metric))

#delete the 1st row of the dataframe
trans_CADM_metric = trans_CADM_metric[-1,]

#adding unique id
trans_CADM_metric$id <- cumsum(!duplicated(trans_CADM_metric[1:2]))
colnames(trans_CADM_metric) <- c("0.2","0.4","0.6","0.8","id")

longdata <- melt(trans_CADM_metric,
                 id="id",
                 measure=c("0.2","0.4","0.6","0.8"))
colnames(longdata)=c("id","sample_completeness","CADM_metric")
longdata

ezANOVA(data = longdata, dv = CADM_metric, wid = id,within = sample_completeness)
#The one-way repeated measures ANOVA can be used to determine whether the means of longdata RF distances are significantly different among four backbone levels(sample completeness 20%,40%,60%,80%)

#check assumptions
#Outliers can be easily identified using box plot methods, implemented in the R function identify_outliers() [rstatix package].

longdata %>%
  group_by(sample_completeness) %>%
  identify_outliers(CADM_metric)


#visualization
#Create a violin plot
violin_plot <- ggplot(longdata, aes(x=sample_completeness, y=CADM_metric,color=sample_completeness))+ geom_violin(trim = FALSE)+theme(legend.position="none")+geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5,binwidth = 0.01)+labs(title="The box plot of CADM for sampling completeness")
violin_plot

#Normality assumption
#The normality assumption can be checked by computing Shapiro-Wilk test for each time point. If the data is normally distributed, the p-value should be greater than 0.05.
longdata %>%
  group_by(sample_completeness) %>%  
  shapiro_test(CADM_metric)


#The longdata RF_distance was normally distributed at each sample completeness, as assessed by Shapiro-Wilk’s test (p > 0.05).

#QQ plot draws the correlation between a given data and the normal distribution. Create QQ plots for each level backbone tree(sample completeness)
ggqqplot(longdata, "CADM_metric", facet.by = "sample_completeness", color = "sample_completeness")

#From the resultant plots, as all the points fall approximately along the reference line, we can assume normality
#The assumption of sphericity will be automatically checked during the computation of the ANOVA test using the R function anova_test() [rstatix package]. The Mauchly’s test is internally used to assess the sphericity assumption

#By using the function get_anova_table() [rstatix] to extract the ANOVA table, the Greenhouse-Geisser sphericity correction is automatically applied to factors violating the sphericity assumption
# res.aov <- anova_test(data = longdata, dv = RF_distance, wid = id,within = sample_completeness)
# get_anova_table(res.aov)

#Post-hoc tests
#perform multiple pairwise paired t-tests between the levels of the within-subjects factor (Sample_completeness). P-values are adjusted using the Bonferroni multiple testing correction method.

# pairwise comparisons
pwc <- longdata %>%
  pairwise_t_test(CADM_metric ~ sample_completeness, paired = TRUE, p.adjust.method = "bonferroni")
pwc

#All the pairwise differences are not statistically significant

#Correlation tests need numerical variables, so convert more than one column in R data frame to numeric values
longdat <- lapply(longdata,as.numeric)

#Calculate the repeated measures correlation coefficient. id is variable giving the subject name/id for each observation. RF_distance is a numeric variable giving the observations for one measure, sample_completeness is a numeric variable giving the observations for the second measure, longdat is the data frame containing the variables.
rmcorr(id, CADM_metric, sample_completeness, longdat)                                  





