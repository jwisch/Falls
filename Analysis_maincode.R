FILEPATH_DATA<-"C:/Users/julie.wisch/Documents/Transition/DraftedMS_Falls/Data/"
FILEPATH_CODE<-"C:/Users/julie.wisch/Documents/Transition/DraftedMS_Falls/Code/"

library(ggplot2)
library(ggpubr)
library(sjPlot)
library(gridExtra)
library(psych)
library(ppcor)

source(paste(FILEPATH_CODE, "CommonFuncs.R", sep = ""))
source(paste(FILEPATH_CODE, "Falls_Funcs.R", sep = ""))

################################################################################################################################
#Data Cleaning
################################################################################################################################
df_PIB <- read.csv(paste(FILEPATH_DATA, "PIB.csv", sep = ""))
dF_fall <- read.csv(paste(FILEPATH_DATA,"adrcpib_fallfinal_deidentified.csv", sep = ""))
df_FC <- read.csv(paste(FILEPATH_DATA,"C:/Users/julie.wisch/Documents/ADRC/Data/Ances_compositeBB298.csv", sep = ""))
df.MRI<-read.csv(paste(FILEPATH_DATA,"HASD_ACS_DR14_3TMR.csv", sep = ""))
demog<-read.csv(paste(FILEPATH_DATA,"DR_demographics_20190122.csv", sep = ""))



#Converting the ID's to factors
df_PIB$MapID <- as.factor(df_PIB$MapID)

#Converting the dates to dates
df_PIB$PET_Date <-as.Date(df_PIB$PET_Date, format = "%m/%d/%Y") #4 digit years get capital Y, 2 digit years get lowercase y



dF_fall$MapID <- as.factor(dF_fall$MapID)
dF_fall$start <- as.Date(dF_fall$start, format = "%d-%b-%y") #use b if it's the month spelled out 

#Dropping all the PIB columns we don't need.
df_PIB<- df_PIB[,c("MapID", "PET_Date",  "PUP_fSUVR_rsf_TOT_CORTMEAN" )]
#To shring a dataframe down to only the columns you need, list them - like above - inside c().  
#Make sure to put your column names in quotes.  Spell them EXACTLY the same.  Seperate with commas.


df_matched<-MatchbyNearestDate(dF_fall, df_PIB, "MapID", "start", "PET_Date")

df_matched<-df_matched[,c(1:4, 29, 30:32, 74:75)]
rm(dF_fall, df_PIB)

df_matched$timefall1<-as.numeric(as.character(df_matched$timefall1))
df_matched$PIBpos<-as.factor(ifelse(df_matched$PUP_fSUVR_rsf_TOT_CORTMEAN > 1.42, 1, 0))
df_matched$BIRTH<-format(as.Date(df_matched$BIRTH, "%d-%b-%y"), "19%y-%m-%d")
df_matched$age<-as.numeric(difftime(df_matched$BIRTH, df_matched$PET_Date, units = "weeks")/(-52))
colnames(df_matched)[4]<-"falls"

df_matched$timegap<-as.numeric(difftime(df_matched$start, df_matched$PET_Date, units = "weeks")/(-52))

df_matched<-df_matched[df_matched$timegap < 2 & df_matched$timegap > -2,] #keeping only people with PIB scans within 2 years of falls enrollment date

#Now getting resting state

df_FC$DATE_SCANNED<-as.Date(df_FC$DATE_SCANNED, format = "%m/%d/%Y")
colnames(df_FC)[1]<-"MapID"

df_FC<-MatchbyNearestDate(df_matched, df_FC, "MapID", "start", "DATE_SCANNED")
df_FC$timegap<-as.numeric(difftime(df_FC$start, df_FC$DATE_SCANNED, units = "weeks")/(-52))

df_FC<-df_FC[df_FC$timegap < 2 & df_FC$timegap > -2,] #keeping only people with PIB scans within 2 years of falls enrollment date

df_FC[,15:105] <- apply(df_FC[,15:105], 2, fisherz)

#Adding in race


colnames(demog)[1]<-"MapID"

df_FC<-merge(df_FC, demog[,c("MapID", "race2", "apoe")], by = "MapID", all.x = TRUE, all.y = FALSE)
df_FC$apoe4<-ifelse(df_FC$apoe == 33 | df_FC$apoe == 23, 0 , 1)

library(tableone)

listVars<-c("PUP_fSUVR_rsf_TOT_CORTMEAN", "timefall1", "EDUC", "age", "GENDER", "race2", "apoe4")
catVars<-c("GENDER", "race2", "apoe4")
table1 <- CreateTableOne(vars = listVars, data = df_FC, factorVars = catVars, strata = c("falls", "PIBpos"))
table1

rm(demog)
#df_FC<-df_FC[ , !(names(df_FC) %in% c("race2", "apoe", "apoe4"))] #dropping some demographics back out so it doesn't jack anything up later....

################################################################################################################################
################################################################################################################################
#Testing for relationship between amyloid level and time to fall
################################################################################################################################
#There are 71 individuals who fall, 50 who do not.  

#Linking amyloid with falls
my_comparisons<-list(c("0", "1"))

p1<-ggboxplot(df_FC, y = "PUP_fSUVR_rsf_TOT_CORTMEAN", x = "falls", color = "falls", palette = "jco", 
          add = "jitter", outlier.shape = NA) + stat_compare_means(comparisons = my_comparisons) + xlab("Did Participant Fall?") + 
  ylab("Cortical Amyloid Accumulation")+scale_color_manual(values=c("#9ecae1","#3182bd"))+
  scale_x_discrete(labels =c("0" = "No", "1" = "Yes")) + theme(legend.position = "none")
df_FC_plot<-df_FC[,1:14]
df_FC_plot$PIBpos<-revalue(df_FC_plot$PIBpos, c("0"="Amyloid Negative", "1"="Amyloid Positive"))
AmyloidPlot<-ggboxplot(df_FC_plot, y = "PUP_fSUVR_rsf_TOT_CORTMEAN", x = "falls", fill = "PIBpos", palette = "jco", 
              add = "jitter", outlier.shape = NA) + stat_compare_means(comparisons = my_comparisons) + xlab("Did Participant Fall?") + 
  ylab("Cortical Amyloid Accumulation")+scale_fill_manual(values=c("#9ecae1","#3182bd"))+
  scale_x_discrete(labels =c("0" = "No", "1" = "Yes")) + labs(fill = "")

pcor.test(df_FC[df_FC$falls == 1, "PUP_fSUVR_rsf_TOT_CORTMEAN"], df_FC[df_FC$falls == 1, "timefall1"], 
          df_FC[df_FC$falls == 1 ,c("EDUC", "age", "GENDER")], method = "pearson")

#basically the same if you treat amyloid as a continuous variable vs categorical



################################################################################################################################
################################################################################################################################
#Doing PCA to get a single intranetwork signature
################################################################################################################################
colnames(df_FC)[c(15, 28, 40, 51, 61, 70, 78, 85, 91, 96, 100, 103, 105)]<-c("Somatomotor", "Lateral Somatomotor",
                                                                             "Cingulate Operculum", "Auditory",
                                                                             "Default Mode", "Memory", "Vision",
                                                                             "Frontoparietal", "Salience",
                                                                             "Subcortical", "Visual Attention",
                                                                             "Default Attention", "Cerebellar")

df_pca<-prcomp(df_FC[df_FC$falls == 1,c(15, 28, 40, 51, 61, 70, 78, 85, 91, 96, 100, 103, 105)])
df_pca$rotation<-df_pca$rotation * -1 #reversing sign so things are more intuitive
#first component explains 37% of variance, second does 15%
df_FC$Signature<-as.numeric(as.matrix(df_FC[,c(15, 28, 40, 51, 61, 70, 78, 85, 91, 96, 100, 103, 105)])%*%df_pca$rotation[,1])

var_coord_func <- function(loadings, comp.sdev){
  loadings*comp.sdev}
loadings <- df_pca$rotation
sdev <- df_pca$sdev
var.coord <- t(apply(loadings, 1, var_coord_func, sdev)) 
var.cos2 <- var.coord^2
# Compute contributions of each brain region to the component
#result is between 0 and 1...should sum to 1
comp.cos2 <- apply(var.cos2, 2, sum)
contrib <- function(var.cos2, comp.cos2){var.cos2/comp.cos2}
var.contrib <- t(apply(var.cos2,1, contrib, comp.cos2))

MakeBarPlots<-function(Component, ComponentTitle){
  theme_set(theme_bw())  
  pcaPlot<-as.data.frame(Component)
  pcaPlot<-cbind(rownames(pcaPlot), pcaPlot[,1])
  colnames(pcaPlot)<-c("region", "contribution")
  pcaPlot<-as.data.frame(pcaPlot)
  pcaPlot$contribution<-as.numeric(as.character(pcaPlot$contribution))
  pcaPlot <- pcaPlot[order(-pcaPlot$contribution), ]  # sort
  # Diverging Barcharts
  p<-ggplot(pcaPlot, aes(x=reorder(region, contribution), y = contribution, label=contribution)) + 
    geom_bar(stat='identity', width=.5)  +
    #coord_cartesian(ylim = c(-0.4, 0.4)) +
    scale_y_continuous()+ ylab("Contribution")+xlab("Intranetwork Connection")+
    labs(subtitle="37% of variance", 
         title= ComponentTitle) + 
    coord_flip() + 
    theme(legend.position = "none")
  return(p)}

MakeBarPlots(subset(var.contrib[,1], var.contrib[,1] > 0), "Global Rs-Fc Signature")
df_FC_plot$Signature<-df_FC$Signature
################################################################################################################################
################################################################################################################################
#Testing for relationship between networks and time to fall
################################################################################################################################
df_FC$cluster<-as.factor(paste(df_FC$PIBpos, df_FC$falls, sep = "-"))
my_comparisons<-list(c("0", "1"))


p2<-ggboxplot(df_FC, x = "falls", y = "Signature", color = "falls", palette = "jco",
          add = "jitter", outlier.shape = NA) + stat_compare_means(comparisons = my_comparisons) + 
  theme(legend.position = "none", axis.title.y = element_text(size = 11))+ 
  ylab("Global Intranetwork Connection") +
  xlab("Did participant fall?") + 
  scale_color_manual(values=c("#9ecae1","#3182bd"))+
  scale_x_discrete(labels =c("0" = "No", "1" = "Yes")) + theme(legend.position = "none")

SignaturePlot<-ggboxplot(df_FC_plot, y = "Signature", x = "falls", fill = "PIBpos", palette = "jco", 
          add = "jitter", outlier.shape = NA) + stat_compare_means(comparisons = my_comparisons) + xlab("Did Participant Fall?") + 
  ylab("Global Rs-Fc Signature")+scale_fill_manual(values=c("#9ecae1","#3182bd"))+
  scale_x_discrete(labels =c("0" = "No", "1" = "Yes")) + labs(fill = "")

grid.arrange(p1, p2, nrow = 1)

df_FC$falls<-as.factor(df_FC$falls)

################################################################################################################################
################################################################################################################################
#Testing for relationship between amyloid level and networks
################################################################################################################################
#Prior lit has found differences in DAN and DMN dependent on amyloid
#https://academic.oup.com/cercor/article/26/2/695/2366897
#also the raichle paper about dmn
pcor.test(df_FC[, "Signature"], df_FC[, "PUP_fSUVR_rsf_TOT_CORTMEAN"], 
          df_FC[,c("EDUC", "age", "GENDER")], method = "pearson")
#p = 0.21
posfall<-pcor.test(df_FC[df_FC$PIBpos == 1 & df_FC$falls == 1, "Signature"], df_FC[df_FC$PIBpos == 1 & df_FC$falls == 1, "PUP_fSUVR_rsf_TOT_CORTMEAN"], 
          df_FC[df_FC$PIBpos == 1 & df_FC$falls == 1,c("EDUC", "age", "GENDER", "apoe4")], method = "pearson")
#p = 0.012, r = -.75
pcor.test(df_FC[df_FC$PIBpos == 0 & df_FC$falls == 1, "Signature"], df_FC[df_FC$PIBpos == 0 & df_FC$falls == 1, "PUP_fSUVR_rsf_TOT_CORTMEAN"], 
          df_FC[df_FC$PIBpos == 0 & df_FC$falls == 1,c("EDUC", "age", "GENDER", "apoe4")], method = "pearson")
#p = 0.34
posno<-pcor.test(df_FC[df_FC$PIBpos == 1 & df_FC$falls == 0, "Signature"], df_FC[df_FC$PIBpos == 1 & df_FC$falls == 0, "PUP_fSUVR_rsf_TOT_CORTMEAN"], 
          df_FC[df_FC$PIBpos == 1 & df_FC$falls == 0,c("EDUC", "age", "GENDER", "apoe4")], method = "pearson")
#p = 0.08, r = 0.70
pcor.test(df_FC[df_FC$PIBpos == 0 & df_FC$falls == 0, "Signature"], df_FC[df_FC$PIBpos == 0 & df_FC$falls == 0, "PUP_fSUVR_rsf_TOT_CORTMEAN"], 
          df_FC[df_FC$PIBpos == 0 & df_FC$falls == 0,c("EDUC", "age", "GENDER", "apoe4")], method = "pearson")
#p = 0.9
ggplot(df_FC[df_FC$PIBpos == 1,], aes(x = Signature, y = PUP_fSUVR_rsf_TOT_CORTMEAN, color = cluster, shape = cluster))+
  geom_point(show.legend = FALSE)+geom_smooth(method = "lm", show.legend = FALSE)+ylab("Cortical Amyloid Accumulation")+
  xlab("Global Rs-Fc Signature")+
  scale_color_manual(values=c("#9ecae1","#3182bd"))+
  scale_shape_manual(values = c(1, 2))+
  annotate("label", x = 0.65, y = 4, size = 3,
           label = paste("Fallers \n R=", round(posfall$estimate, 3), 
                         "\n p=", round(posfall$p.value, 3))) +
  annotate("label", x = 0.75, y = 1.3, size = 3,
           label = paste("Non-Fallers \n R=", round(posno$estimate, 3), 
                         "\n p=", round(posno$p.value, 3)))  + theme_classic()

ggplot(df_FC, aes(x = Signature, y = PUP_fSUVR_rsf_TOT_CORTMEAN, color = falls, shape = falls))+
  geom_point(show.legend = FALSE)+geom_smooth(method = "lm", show.legend = FALSE)+ylab("Cortical Amyloid Accumulation")+
  xlab("Intranetwork Connectivity Signature")+
  scale_color_manual(values=c("#9ecae1","#3182bd"))+
  scale_shape_manual(values = c(1, 2))

#Can't pick anything up when we lump them all together, but for amyloid positive people...
#Amyloid positive individuals who fall tend to have high amyloid and low signatures OR low amyloid and high signatures
#Apos people who do not fall have either low amyloid-low signature or high amyloid-high signature
#so it seems like you can compensate...if you have a lot of amyloid but your brain is still connected, you're not going to tip over
#if you have crappy brain connections, but you don't have too much amyloid, you're still not going to fall over

#Checking to see if there's a relationship between network strength and number of falls
pcor.test(df_FC[ df_FC$falls == 1, "Signature"], df_FC[df_FC$falls == 1, "totfall"], 
          df_FC[df_FC$falls == 1,c("EDUC", "age", "GENDER", "apoe4")], method = "pearson")

#amyloid positives
pcor.test(df_FC[df_FC$PIBpos == 1 & df_FC$falls == 1, "Signature"], df_FC[df_FC$PIBpos == 1 & df_FC$falls == 1, "totfall"], 
                   df_FC[df_FC$PIBpos == 1 & df_FC$falls == 1,c("EDUC", "age", "GENDER", "apoe4")], method = "pearson")
#amyloid negatives
pcor.test(df_FC[df_FC$PIBpos == 0 & df_FC$falls == 1, "Signature"], df_FC[df_FC$PIBpos == 0 & df_FC$falls == 1, "totfall"], 
          df_FC[df_FC$PIBpos == 0 & df_FC$falls == 1,c("EDUC", "age", "GENDER", "apoe4")], method = "pearson")

ggplot(df_FC[df_FC$falls == 1,], aes(x = Signature, y = totfall, color = PIBpos, shape = PIBpos))+
  geom_point(show.legend = FALSE)+geom_smooth(method = "lm", show.legend = FALSE)+ylab("Total Number of Falls")+
  xlab("Intranetwork Connectivity Signature")+
  scale_color_manual(values=c("#9ecae1","#3182bd"))+
  scale_shape_manual(values = c(1, 2))

#################################################################################################################
#################################################################################################################
#################################################################################################################
#Now adding in volumes


#getting z scores on a region by region basis

normalize<-function(PARAM){(PARAM - mean(PARAM))/sd(PARAM)}

#Now looking at Liang Wang's region

df.MRI$ADsig<-(normalize(df.MRI$MR_LV_INFRTMP)+normalize(df.MRI$MR_LV_MIDTMP)+normalize(df.MRI$MR_LV_SUPERTMP)+
                 normalize(df.MRI$MR_RV_INFRTMP)+normalize(df.MRI$MR_RV_MIDTMP)+normalize(df.MRI$MR_RV_SUPERTMP)+
                 normalize(df.MRI$MR_LV_INFRPRTL)+normalize(df.MRI$MR_LV_SUPERPRTL)+
                 normalize(df.MRI$MR_RV_INFRPRTL)+normalize(df.MRI$MR_RV_SUPERPRTL)+
                 normalize(df.MRI$MR_LV_ENTORHINAL)+normalize(df.MRI$MR_RV_ENTORHINAL)+
                 normalize(df.MRI$MR_LV_PRECUNEUS)+normalize(df.MRI$MR_RV_PRECUNEUS))/14

df.MRI$HippoVol<-(normalize(df.MRI$MR_LV_HIPPOCAMPUS)+normalize(df.MRI$MR_RV_HIPPOCAMPUS))/2

df.MRI$MR_Date<-as.Date(df.MRI$MR_Date, format = "%m/%d/%Y")
colnames(df.MRI)[4]<-"MapID"
df_matched<-MatchbyNearestDate(df_matched, df.MRI[,c("MapID", "MR_Date", "ADsig", "HippoVol")], "MapID", "start", "MR_Date")

df_matched$timegap2<-as.numeric(difftime(df_matched$start, df_matched$MR_Date, units = "weeks")/(-52))

df_matched<-df_matched[df_matched$timegap2 < 2 & df_matched$timegap2 > -2,] #keeping only people with mri scans within 2 years of falls enrollment date
df_matched$PIBpos<-revalue(df_matched$PIBpos, c("0"="Amyloid Negative", "1"="Amyloid Positive"))


p1<-ggboxplot(df_matched, x = "falls", y = "ADsig", color = "falls", palette = "jco",
              add = "jitter", outlier.shape = NA) + stat_compare_means(comparisons = my_comparisons) + 
  theme(legend.position = "none", axis.title.y = element_text(size = 11))+ 
  ylab("AD Signature Regions") +
  xlab("Did participant fall?") + 
  scale_color_manual(values=c("#9ecae1","#3182bd"))+
  scale_x_discrete(labels =c("0" = "No", "1" = "Yes")) + theme(legend.position = "none")


p2<-ggboxplot(df_matched, x = "falls", y = "HippoVol", color = "falls", palette = "jco",
              add = "jitter", outlier.shape = NA) + stat_compare_means(comparisons = my_comparisons) + 
  theme(legend.position = "none", axis.title.y = element_text(size = 11))+ 
  ylab("Hippocampal Volume") +
  xlab("Did participant fall?") + 
  scale_color_manual(values=c("#9ecae1","#3182bd"))+
  scale_x_discrete(labels =c("0" = "No", "1" = "Yes")) + theme(legend.position = "none")
HippoVol<-ggboxplot(df_matched, y = "HippoVol", x = "falls", fill = "PIBpos", palette = "jco", 
          add = "jitter", outlier.shape = NA) + stat_compare_means(comparisons = my_comparisons) + xlab("Did Participant Fall?") + 
  ylab("Hippocampus Volume")+scale_fill_manual(values=c("#9ecae1","#3182bd"))+
  scale_x_discrete(labels =c("0" = "No", "1" = "Yes")) + labs(fill = "")

grid.arrange(p1, p2, nrow = 1)

df_matched_FC<-MatchbyNearestDate(df_matched, df_FC[,c("MapID", "start", "Signature", "apoe4")], "MapID", "MR_Date", "start")
df_matched_FC$falls <-  as.factor(df_matched_FC$falls)
p1<-ggplot(df_matched_FC, aes(x = Signature, y = ADsig, color = falls, shape = falls))+
  geom_point(show.legend = FALSE)+geom_smooth(method = "lm", show.legend = FALSE)+ylab("AD Signature Volume")+
  xlab("Intranetwork Connectivity Signature")+xlim(c(0.4, 1.4))+
  scale_color_manual(values=c("#9ecae1","#3182bd"))+
  scale_shape_manual(values = c(1, 2))

ggplot(df_matched_FC, aes(x = Signature, y = HippoVol, color = falls, shape = falls))+
  geom_point(show.legend = FALSE)+geom_smooth(method = "lm", show.legend = TRUE)+ylab("Hippocampal Volume")+
  xlab("Global Rs-Fc Signature")+
  scale_color_manual(values=c("#9ecae1","#3182bd"))+xlim(c(0.4, 1.4))+
  scale_shape_manual(values = c(1, 2)) + 
  annotate("label", x = 1.25, y = -1.42, size = 3,
           label = "Overall Correlation \n R = 0.226; p = 0.046\nParticipants do not differ by Fall Status\n p = 0.31") +
 theme_classic()
grid.arrange(p1, p2, nrow = 2)


p1<-ggplot(df_matched_FC, aes(y = PUP_fSUVR_rsf_TOT_CORTMEAN, x = ADsig, color = falls, shape = falls))+
  geom_point(show.legend = FALSE)+geom_smooth(method = "lm", show.legend = FALSE)+xlab("AD Signature Volume")+
  ylab("Cortical Amyloid Accumulation")+
  scale_color_manual(values=c("#9ecae1","#3182bd"))+xlim(c(-1, 1))+
  scale_shape_manual(values = c(1, 2))

p2<-ggplot(df_matched_FC, aes(y = PUP_fSUVR_rsf_TOT_CORTMEAN, x = HippoVol, color = falls, shape = falls))+
  geom_point(show.legend = FALSE)+geom_smooth(method = "lm", show.legend = TRUE)+xlab("Hippocampal Volume")+
  ylab("Cortical Amyloid Accumulation")+
  scale_color_manual(values=c("#9ecae1","#3182bd"))+xlim(c(-1, 1))+
  scale_shape_manual(values = c(1, 2)) + theme(legend.position = "bottom")
grid.arrange(p1, p2, nrow = 2)

grid.arrange(AmyloidPlot, SignaturePlot, HippoVol, nrow = 1)

model.null<-lm(PUP_fSUVR_rsf_TOT_CORTMEAN ~ falls + PIBpos + falls:PIBpos + GENDER + apoe4 + EDUC + age, data = df_matched_FC)
model.new<-lm(PUP_fSUVR_rsf_TOT_CORTMEAN ~ PIBpos + GENDER + apoe4 + EDUC + age, data = df_matched_FC)
anova(model.null, model.new)
model.new<-lm(PUP_fSUVR_rsf_TOT_CORTMEAN ~ falls + GENDER + apoe4 + EDUC + age, data = df_matched_FC)
anova(model.null, model.new)
model.new<-lm(PUP_fSUVR_rsf_TOT_CORTMEAN ~ falls + PIBpos +  GENDER + apoe4 + EDUC + age, data = df_matched_FC)
anova(model.null, model.new)


model.null<-lm(Signature ~ falls + PIBpos + falls:PIBpos + GENDER + apoe4 + EDUC + age, data = df_matched_FC)
model.new<-lm(Signature ~ PIBpos + GENDER + apoe4 + EDUC + age, data = df_matched_FC)
anova(model.null, model.new)
model.new<-lm(Signature ~ falls + GENDER + apoe4 + EDUC + age, data = df_matched_FC)
anova(model.null, model.new)
model.new<-lm(Signature ~ falls + PIBpos +  GENDER + apoe4 + EDUC + age, data = df_matched_FC)
anova(model.null, model.new)


model.null<-lm(HippoVol ~ falls + PIBpos + falls:PIBpos + GENDER + apoe4 + EDUC + age, data = df_matched_FC)
model.new<-lm(HippoVol ~ PIBpos + GENDER + apoe4 + EDUC + age, data = df_matched_FC)
anova(model.null, model.new)
model.new<-lm(HippoVol ~ falls + GENDER + apoe4 + EDUC + age, data = df_matched_FC)
anova(model.null, model.new)
model.new<-lm(HippoVol ~ falls + PIBpos +  GENDER + apoe4 + EDUC + age, data = df_matched_FC)
anova(model.null, model.new)


#Are people who are PIBpos and have reduced AD sig more likely to fall?
length(df_matched_FC[df_matched_FC$Signature < median(df_matched_FC$Signature), "MapID"])
length(df_matched_FC[df_matched_FC$Signature < median(df_matched_FC$Signature) & df_matched_FC$falls == 1, "MapID"])
length(df_matched_FC[df_matched_FC$Signature < median(df_matched_FC$Signature) & df_matched_FC$falls == 0, "MapID"])
length(df_matched_FC[df_matched_FC$Signature < median(df_matched_FC$Signature) & df_matched_FC$PIBpos == 1, "MapID"])
length(df_matched_FC[df_matched_FC$Signature < median(df_matched_FC$Signature) & df_matched_FC$PIBpos == 0, "MapID"])


length(df_matched_FC[df_matched_FC$Signature > median(df_matched_FC$Signature)  & df_matched_FC$PIBpos == 1, "MapID"])
length(df_matched_FC[df_matched_FC$Signature > median(df_matched_FC$Signature) & df_matched_FC$falls == 1  & df_matched_FC$PIBpos == 1, "MapID"])
length(df_matched_FC[df_matched_FC$Signature < median(df_matched_FC$Signature) & df_matched_FC$falls == 0  & df_matched_FC$PIBpos == 1, "MapID"])

df_matched_FC$LowSig<-ifelse(df_matched_FC$Signature < median(df_matched_FC$Signature), 1, 0)
#There are 22 PIB positive people
#out of the PIB positive participants, 14 have low signature. 9 of them fell, 5 did not.
#out of the PIB positive participants, 8 have high signature. 3 fell. 5 did not.
prop.test(x = c(9, 5), n = c(14, 14)) 
prop.test(x = c(9, 5, 3, 5), n = c(22, 22, 22, 22))


prop.test(table(df_matched_FC[df_matched_FC$PIBpos == 1, "LowSig"], df_matched_FC[df_matched_FC$PIBpos == 1, "falls"]), correct = FALSE)
prop.test(table(df_matched_FC[df_matched_FC$falls == 1, "LowSig"], df_matched_FC[df_matched_FC$falls == 1, "PIBpos"]), correct = FALSE)

#82 people. 22 are PIB positive. 
