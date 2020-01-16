#Function
MatchbyNearestDate<-function(df1, df2, ID, Date1, Date2){
  z <- lapply(intersect(df1[,ID],df2[,ID]),function(id) {
    df1 <- df1[df1[,ID] == id,]
    df2 <- df2[df2[,ID] == id,]
    
    df1[,"indices"] <- sapply(df1[,Date1],function(d) which.min(abs(df2[,Date2] - d)))
    df2[,"indices"] <- 1:nrow(df2)
    
    merge(df1,df2,by=c(ID,'indices'))
  })
  
  df_matched <- do.call(rbind,z)
  df_matched$indices <- NULL
  return(df_matched)
  
}

CleanMRI<-function(MRI, IDLIST = NULL){
MRI<-MRI[!(MRI$Scanner == "CCIR Vida"),]
colnames(MRI)[5]<-"ID"
MRI<-MRI[,c("ID", "MR_Date", "MR_LT_INFRTMP", "MR_LT_MIDTMP", "MR_LT_SUPERTMP", "MR_RT_INFRTMP", "MR_RT_MIDTMP",
            "MR_RT_SUPERTMP", "MR_LT_INFRPRTL", "MR_LT_SUPERPRTL", "MR_RT_INFRPRTL", "MR_RT_SUPERPRTL",
            "MR_LT_ENTORHINAL", "MR_RT_ENTORHINAL", "MR_LT_PRECUNEUS", "MR_RT_PRECUNEUS",
            "MR_LV_HIPPOCAMPUS", "MR_RV_HIPPOCAMPUS", "MR_TOTV_INTRACRANIAL")]
MRI$MR_Date<-as.Date(MRI$MR_Date, format = "%m/%d/%Y")
MRI<-MRI[complete.cases(MRI),]
MRI<-MRI[MRI$ID %in% IDLIST,]

ScaleandNormalize<-function(COLUMN, COLNAME){
  MEAN<-mean(MRI$MR_TOTV_INTRACRANIAL)
  model<-lm(COLUMN ~ MR_TOTV_INTRACRANIAL, data = MRI)
  MRI[,paste("Normalized", COLNAME, sep = "")]<-scale(COLUMN - (coef(model)[2] * (MRI$MR_TOTV_INTRACRANIAL - MEAN)))
  return(MRI)}

for(i in 3:18){
  MRI<-ScaleandNormalize(MRI[,i], names(MRI[i]))
}


MRI$ADsig<-(rowSums(MRI[,20:33]))/14
MRI$HippoVol<-(MRI$NormalizedMR_LV_HIPPOCAMPUS + MRI$NormalizedMR_RV_HIPPOCAMPUS )/2
MRI<-MRI[,c("ID", "MR_Date", "ADsig", "HippoVol")]
return(MRI)}

GetPACC<-function(x){
  x[,"MMSE"] <- as.numeric(ifelse(!is.na(x[,"MMSE"]) & x[,"MMSE"] > 30, "NA", x[,"MMSE"]))
  for(i in 1:length(x[,1])){
    x$nacount[i]<-sum(is.na(c(x[i, "srtfree"], x[i, "digsym"], x[i, "MEMUNITS"], x[i, "MMSE"])))}
  x[,"srtfree"]<-scale(x[,"srtfree"], na.rm = TRUE)
  x[,"digsym"]<-scale(x[,"digsym"], na.rm = TRUE)
  x[,"MEMUNITS"]<-scale(x[,"MEMUNITS"], na.rm = TRUE)
  x[,"MMSE"]<-scale(x[,"MMSE"], na.rm = TRUE)
  x[is.na(x)] <- 0
  PACC<-(x[,"srtfree"]+x[,"digsym"]+x[,"MEMUNITS"]+x[,"MMSE"])/(4-x[,"nacount"])
  return(PACC)}
