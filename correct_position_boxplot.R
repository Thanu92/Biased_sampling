#Sep 20, 2023
#RAxML
Correct_position<- read.csv("Sis_correct_position_new.csv")
class(Correct_position)
names(Correct_position)
library(ggpubr)
ggboxplot(Correct_position,x="sample_type",y="correct_placement_percentage", color = "sample_type",palette = c("#6BAED6" ,"#4292C6", "#2171B5", "#084594","#003366","#F1B6DA","#DE77AE","#C51B7D","#8E0152","#800080","#B8E186")
          ,main="Boxplot of correct sister species position percentage of placement sequences based on sample types",add= "point", width=1,outlier.shape = NA)

#Jan 20, 2024
#EPAng
getwd()
Correct_position<- read.csv("EPAng_Sample_type.csv")
class(Correct_position)
names(Correct_position)
library(ggpubr)
ggboxplot(Correct_position,x="Sample_Type",y="Correct_Percentage", color = "Sample_Type",palette = c("#6BAED6" ,"#4292C6", "#2171B5", "#084594","#003366","#F1B6DA","#DE77AE","#C51B7D","#8E0152","#800080")
          ,main="Boxplot of correct sister species position percentage of placement sequences based on sample types for EPAng ",add= "point", width=1,outlier.shape = NA)
#EPAng_RF_Sample_type
RF<- read.csv("EPAng_RF_sample_type.csv")
names(RF)
bp <- ggboxplot(RF,x="Sample_Type",y="RF", color = "Sample_Type",palette = c("#6BAED6" ,"#4292C6", "#2171B5", "#084594","#003366","#F1B6DA","#DE77AE","#C51B7D","#8E0152","#800080")
                ,main="Boxplot of RF distance of placement sequences for EPAng based on different Sample types ",add= "point", width=1,outlier.shape = NA)
ggpar(bp,orientation = "reverse")




#In App-spaM/build folder
Correct_position<- read.csv("species_position_App-spam.csv")
class(Correct_position)
names(Correct_position)
library(ggpubr)
ggboxplot(Correct_position,x="sample_type",y="correct_position_percentage", color = "sample_type",palette = c("#6BAED6" ,"#4292C6", "#2171B5", "#084594")
          ,main="Boxplot of correct species position percentage based on sample types",add= "point", width=0.8,outlier.shape = NA)



#EPAng
Correct_position<- read.csv("EPAng_correct_positionCor.csv")
class(Correct_position)
names(Correct_position)
library(ggpubr)
ggboxplot(Correct_position,x="sample_type",y="correct_position_percentage", color = "sample_type",palette = c("#6BAED6" ,"#4292C6", "#2171B5", "#084594")
          ,main="Boxplot of correct species position percentage based on EPAng sample types",add= "point", width=0.8,outlier.shape = NA)

#datadriven 9 samples(random,strati and bias)
Correct_position<- read.csv("Bias_correct_sp_position.csv")
class(Correct_position)
names(Correct_position)
ggboxplot(Correct_position,x="Sample_type",y="correct_position_percentage", color = "Sample_type",palette = c("#6BAED6" ,"#4292C6", "#2171B5", "#084594","#F1B6DA","#DE77AE","#C51B7D","#8E0152","#B8E186")
          ,main="Boxplot of correct species position percentage based on sample types",add= "point", width=1,outlier.shape = NA)

getwd()
#8E0152" "#C51B7D" "#DE77AE" "#F1B6DA" "#FDE0EF" "#F7F7F7" "#E6F5D0" "#B8E186" "#7FBC41" "#4D9221" "#276419"
library(ggpubr)
#For RF_distance of three software
RF<- read.csv("RF_dist_softwareTool.csv")
bp <- ggboxplot(RF,x="Software_tool",y="RF_distance", color = "Software_tool",palette = c("#6BAED6" ,"#4292C6", "#2171B5", "#084594","#F1B6DA","#DE77AE","#C51B7D","#8E0152","#B8E186","#7FBC41","#4D9221","#276419")
          ,main="Boxplot of RF distance based on software tool",add= "point", width=1,outlier.shape = NA)
ggpar(bp,orientation = "reverse")


#species position for software tools
Correct_position<- read.csv("SW_tool_correct_position.csv")
ggboxplot(Correct_position,x="software_tool",y="correct_position_percentage", color = "software_tool",palette = c("#6BAED6" ,"#4292C6", "#2171B5", "#084594","#F1B6DA","#DE77AE","#C51B7D","#8E0152","#B8E186","#7FBC41","#4D9221","#276419")
                ,main="Boxplot of correct species position percentage based on software tool",add= "point", width=1,outlier.shape = NA)


#For RF_distance of three sample types
RF<- read.csv("RF_dist_sample_type.csv")

library(ggpubr)
bp <- ggboxplot(RF,x="Sample_type",y="RF_distance", color = "Sample_type",palette = c("#6BAED6" ,"#4292C6", "#2171B5", "#084594","#003366","#F1B6DA","#DE77AE","#C51B7D","#8E0152","#800080","#B8E186")
          ,main="Boxplot of RF distance based on sample types",add= "point", width=1,outlier.shape = NA)
ggpar(bp,orientation = "reverse")

#Genus percentage for Random and startified sampling 20%-99%
Correct_position<- read.csv("Genus_percentage_Random_Strati.csv")
class(Correct_position)
names(Correct_position)
library(ggpubr)
ggboxplot(Correct_position,x="sample_type",y="percentatge_of_correct_placement", color = "sample_type",palette = c("#6BAED6" ,"#4292C6", "#2171B5", "#084594","#003366","#F1B6DA","#DE77AE","#C51B7D","#8E0152","#800080")
          ,main="Boxplot of correct genus position percentage of placement sequences based on sample types",add= "point", width=1,outlier.shape = NA)

#Jan 16, 2023
getwd()
#Correct position for software tools (RAxMLEPA and epa-ng)
# Correct_position<- read.csv("Strati_SWtools_sister.csv")
Correct_position<- read.csv("Strati_SWtools_new.csv")
names(Correct_position)
class(Correct_position)

library(ggpubr)
ggboxplot(Correct_position,x="Sample_Type",y="Correct_percentage", color = "Sample_Type",palette = c("#6BAED6" ,"#4292C6", "#2171B5", "#084594","#003366","#F1B6DA","#DE77AE","#C51B7D","#8E0152","#800080")
          ,main="Boxplot of correct sister position percentage of placement sequences for stratified sampling based on different SW tools ",add= "point", width=1,outlier.shape = NA)

#Correct position for software tools (RAxMLEPA and epa-ng)
# RF<- read.csv("RF_Strati_SWtools.csv")
RF<- read.csv("RF_Strati_SWtools_new.csv")
names(RF)
bp <- ggboxplot(RF,x="Sample_complteness",y="RF_distance", color = "Sample_complteness",palette = c("#6BAED6" ,"#4292C6", "#2171B5", "#084594","#003366","#F1B6DA","#DE77AE","#C51B7D","#8E0152","#800080")
          ,main="Boxplot of RF distance of placement sequences for stratified sampling based on different SW tools ",add= "point", width=1,outlier.shape = NA)
ggpar(bp,orientation = "reverse")

#Feb 11,2024
#datadriven 9 samples(random,strati and bias)
Correct_position<- read.csv("RAxML_SAmple_type_extended.csv")
class(Correct_position)
names(Correct_position)
head(Correct_position)
ggboxplot(Correct_position,x="sample_type",y="correct_placement_percentage", color = "sample_type",palette = c("#6BAED6" ,"#4292C6", "#2171B5", "#084594","#003366","#F1B6DA","#DE77AE","#C51B7D","#8E0152","#800080","#B8E186")
          ,main="Boxplot of correct species position percentage based on sample types",add= "point", width=1,outlier.shape = NA)
