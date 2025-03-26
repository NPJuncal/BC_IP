# Iberian Peninsula Blue Carbon inventory #


setwd("C:/Users/npjun/Dropbox/Seagrasses/Inventarios P. Iberica/R Project/BC_IP")

#Libraries#

library(ggplot2) # plots
library(gridExtra) # grid arrange (combining plots)
library(dplyr) # %>%
library(reshape2) # melt function
library(ggpubr)
library(broom) # glance function
library(corrplot)
library(mapproj)
library(grid) # textGrob for plots
library(rbacon)
library(tidyr)
library(stringr) # function str_wrap for plot labels in two rows
library(ggbreak) # scale_y_break function
library(BlueCarbon)

# load functions ----------------------------------------------------------

std.error <- function(x) sd(x, na.rm=TRUE)/sqrt(length(x))
contar<-function(x) length(which(!is.na(x)))
coef.var<-function(x) sd(x, na.rm=TRUE)/mean(x, na.rm=TRUE)



# functions: estimate_h, estimate_stock, test_extrapolation and estimate_flux
# from library "BlueCarbon" (under development)
# https://github.com/EcologyR/BlueCarbon.git

library(BlueCarbon)


estimate_h <- function(df = NULL) {
  
  # class of the dataframe
  if (is.data.frame(df)==FALSE) {stop("The data provided is not class data.frame, please chaeck data and transforme")}
  
  # name of the columns
  if ("Core.ID" %in% colnames(df)==FALSE) {stop("There is not column named Core.ID. Please, check necessary columns in functions documentation")}
  if ("DMin.D" %in% colnames(df)==FALSE) {stop("There is not column named DMin.D. Please, check necessary columns in functions documentation")}
  if ("DMax.D" %in% colnames(df)==FALSE) {stop("There is not column named DMax.D. Please, check necessary columns in functions documentation")}
  
  # class of the columns
  if (is.numeric(df$DMin.D)==FALSE) {stop("Minimum depth data is not class numeric, please check")}
  if (is.numeric(df$DMax.D)==FALSE) {stop("Maximum depth data is not class numeric, please check")}
  
  #check for NAs in depth columns
  if (sum(is.na(df$DMin.D))>0){stop("Samples minimun depth column has NAs, please check")}
  if (sum(is.na(df$DMax.D))>0){stop("Samples maximun depth column has NAs, please check")}
  
  
  # create individual data frames per each core
  
  df$Core.ID <- factor(df$Core.ID, levels=unique(df$Core.ID))
  X<-split(df, df$Core.ID)
  
  
  columns<-c("EMin","EMax","h")
  Fdf2 = data.frame(matrix(nrow = 0, ncol = length(columns)))
  colnames(Fdf2) = columns
  
  
  for(i in 1:length(X)) {
    
    Data<-as.data.frame(X[i])
    colnames(Data)<-colnames(df)
    
    #check if there is spaces between samples (e.g, first sample ends at 5 cm and next starts at 7)
    space<- c()
    
    for (j in 1:(nrow(Data)-1)) {
      
      # if there are no spaces between samples min and maximun depth of samples remain the same
      if (Data[j,which(colnames(Data)=="DMax.D")] == Data[j+1,which(colnames(Data)=="DMin.D")]) {
        space[j]<-FALSE} else {space[j]<-TRUE}}
    
    if (any(space==TRUE)) {
      # if there are spaces between samples it estimate the medium point between the maximum depth of the sample and the minimun depth of the next sample
      # and divide that distance between both samples
      Data <- cbind(Data, EMin=NA, EMax=NA)
      Data[1,"EMin"]<-0
      Data[nrow(Data),"EMax"]<-Data[nrow(Data),"DMax.D"]
      for (j in 1:(nrow(Data)-1)) {
        if(space[j]==TRUE) {
          Data[j,"EMax"]<-Data[j,"DMax.D"]+((Data[j+1,"DMin.D"]-Data[j,"DMax.D"])/2)
          Data[j+1,"EMin"]<-Data[j,"DMax.D"]+((Data[j+1,"DMin.D"]-Data[j,"DMax.D"])/2)} else {
            Data[j,"EMax"]<-Data[j,"DMax.D"]
            Data[j+1,"EMin"]<-Data[j+1,"DMin.D"]}}
      
    }  else{
      Data <- cbind(Data, EMin=NA, EMax=NA)
      Data$EMin<-Data$DMin.D
      Data$EMax<-Data$DMax.D
      
    }
    
    Data <- cbind(Data, h=NA)
    
    #estimation of the thickness of the sample (h) from the new minimun and max depth of the sample
    
    Data<- Data |> dplyr::mutate (h = EMax-EMin)
    
    temp<-cbind(Data$EMin, Data$EMax, Data$h)
    colnames(temp)<-colnames(Fdf2)
    Fdf2<-rbind(Fdf2, temp)
    
  }
  Fdf<-cbind(df, Fdf2)
  
  return(Fdf)
}

estimate_stock <- function(df = NULL, Depth = 100) {
  
  # class of the dataframe
  if (is.data.frame(df)==FALSE) {stop("The data provided is not class data.frame, please chaeck data and transforme")}
  if (is.numeric(Depth)==FALSE) {stop("The Depth provided is not class numeric, please chaeck data and transforme")}
  
  # name of the columns
  if ("Core.ID" %in% colnames(df)==FALSE) {stop("There is not column named Core.ID. Please, check necessary columns in functions documentation")}
  if ("DMin.D" %in% colnames(df)==FALSE) {stop("There is not column named Min.D. Please, check necessary columns in functions documentation")}
  if ("DMax.D" %in% colnames(df)==FALSE) {stop("There is not column named Max.D. Please, check necessary columns in functions documentation")}
  if ("DDBD" %in% colnames(df)==FALSE) {stop("There is not column named DDBD. Please, check necessary columns in functions documentation")}
  if ("POC" %in% colnames(df)==FALSE) {stop("There is not column named POC. Please, check necessary columns in functions documentation")}
  
  # class of the columns
  if (is.numeric(df$DMin.D)==FALSE) {stop("Minimum depth data is not class numeric, please check")}
  if (is.numeric(df$DMax.D)==FALSE) {stop("Maximum depth data is not class numeric, please check")}
  if (is.numeric(df$DDBD)==FALSE) {stop("Dry Bulk Density data is not class numeric, please check")}
  if (is.numeric(df$POC)==FALSE) {stop("Organic carbon data is not class numeric, please check")}
  
  
  
  df<-df[!is.na(df$POC),]
  
  # estimate thickness of the sample
  
  dfh<-estimate_h(df)
  
  # estimate stocks
  
  X<-split(dfh, dfh$Core.ID)
  
  BCS <- data.frame(Core.ID=character(),
                    S.WC=numeric(),
                    D.Max=numeric(),
                    Stock=numeric())
  
  for(i in 1:length(X)) {
    BCS[i,1]<-names(X[i])
    Data<-as.data.frame(X[i])
    colnames(Data)<-colnames(dfh)
    
    if(nrow(Data)<3) next
    
    else{
      
      #estimation of carbon g cm2 per sample, OCgcm2= carbon density (g cm3) by thickness (h)
      Data <-Data |> dplyr::mutate (OCgcm2 = DDBD*(POC/100)*h)
      
      #estimation of the OC stock in the whole core
      BCS[i,2]<-sum(Data[,which(colnames(Data)=="OCgcm2")])
      BCS[i,3]<-max(Data[,which(colnames(Data)=="EMax")])
      
      #if core exactly the standarization depth, we keep the stock of the whole core
      if(max(Data$EMax)==Depth) {BCS[i,4]<-sum(Data[,which(colnames(Data)=="OCgcm2")])}
      
      
      else{
        
        # if the core longer than the standardization depth we estimate the stock until that depth
        if (max(Data$EMax)>=Depth)
          
        {
          Data<-Data[c(1:(length(which(Data$EMax <=Depth))+1)),]
          
          if(nrow(Data)<3) next
          
          else{
            
            BCS[i,4]<-(((sum(Data[c(1:(nrow(Data)-1)),which(colnames(Data)=="OCgcm2")]))+
                          (Data[nrow(Data),which(colnames(Data)=="OCgcm2")]/((max(Data$EMax)-Data[(nrow(Data)-1),which(colnames(Data)=="EMax")]))
                           *(Depth-Data[(nrow(Data)-1),which(colnames(Data)=="EMax")]))))}}
        
        #if core shorter than than the standardization depth we model the OC acumulated mass with depth and predict the stock at that depth
        else {
          
          Data <-Data |> dplyr::mutate (OCM = cumsum(OCgcm2))
          model<-lm(OCM ~ EMax, data=Data)
          BCS[i,4]<-coef(model)[1] + Depth*coef(model)[2]}}
    }}
  return(BCS)
}

test_extrapolation <- function(df = NULL, Depth = 100) {
  
  # class of the dataframe
  if (is.data.frame(df)==FALSE) {stop("The data provided is not class data.frame, please chaeck data and transforme")}
  if (is.numeric(Depth)==FALSE) {stop("The Depth provided is not class numeric, please chaeck data and transforme")}
  
  # name of the columns
  if ("Core.ID" %in% colnames(df)==FALSE) {stop("There is not column named Core.ID. Please, check necessary columns in functions documentation")}
  if ("DMin.D" %in% colnames(df)==FALSE) {stop("There is not column named Min.D. Please, check necessary columns in functions documentation")}
  if ("DMax.D" %in% colnames(df)==FALSE) {stop("There is not column named Max.D. Please, check necessary columns in functions documentation")}
  if ("DDBD" %in% colnames(df)==FALSE) {stop("There is not column named DDBD. Please, check necessary columns in functions documentation")}
  if ("POC" %in% colnames(df)==FALSE) {stop("There is not column named POC. Please, check necessary columns in functions documentation")}
  
  # class of the columns
  if (is.numeric(df$DMin.D)==FALSE) {stop("Minimum depth data is not class numeric, please check")}
  if (is.numeric(df$DMax.D)==FALSE) {stop("Maximum depth data is not class numeric, please check")}
  if (is.numeric(df$DDBD)==FALSE) {stop("Dry Bulk Density data is not class numeric, please check")}
  if (is.numeric(df$POC)==FALSE) {stop("Organic carbon data is not class numeric, please check")}
  
  
  # we select those cores larger than the standard depth
  
  columns<-colnames(df)
  DataE = data.frame(matrix(nrow = 0, ncol = length(columns)))
  colnames(DataE) = columns
  
  
  
  X <- split(df, df$Core.ID)
  
  for (i in 1:length(X)) {
    Data <- as.data.frame(X[i])
    colnames(Data) <- colnames(df)
    
    if (max(Data$DMax.D) < Depth) next
    
    DataE<-rbind(DataE,Data)}
  
  
  # estimate the stock at the standar depth
  
  ExtS <- data.frame(
    Core.ID = character(),
    EXT.100 = numeric(),
    EXT.90 = numeric(),
    EXT.75 = numeric(),
    EXT.50 = numeric(),
    EXT.25 = numeric()
  )
  
  
  hundreth<- Depth #100% of the extrapolation length
  ninety <- Depth * 0.9 #90% of the extrapolation length
  seventy <- Depth * 0.75 #90% of the extrapolation length
  fifhty <- Depth * 0.50 #90% of the extrapolation length
  twentififty <- Depth * 0.25 #90% of the extrapolation length
  
  
  #estimate observed stock
  temp100<-estimate_stock (DataE, Depth=hundreth)
  ExtS[c(1:nrow(temp100)),"Core.ID"]<-temp100$Core.ID
  ExtS$EXT.100<-temp100$Stock
  
  # estimate stock with a 90, 75, 50 and 25 percentage of the standard depth
  
  X <- split(DataE, DataE$Core.ID)
  
  for (i in 1:length(X)) {
    Data <- as.data.frame(X[i])
    colnames(Data) <- colnames(DataE)
    
    Data<-Data[!is.na(Data$POC),]
    Data<-estimate_h(Data)
    
    #estimation of carbon g cm2 per sample, OCgcm2= carbon density (g cm3) by thickness (h)
    Data <-Data |> dplyr::mutate (OCgcm2 = DDBD*(POC/100)*h)
    
    #estimate organic carbon accumulated mass
    Data <-Data |> dplyr::mutate (OCM = cumsum(OCgcm2))
    
    Data<- subset(Data,Data$DMax.D<=ninety)
    
    if (nrow(Data)>3){
      model90<-lm(OCM ~ DMax.D, data=Data)
      ExtS[i,3]<-coef(model90)[1] + 100*coef(model90)[2]}
    
    Data<- subset(Data,Data$DMax.D<=seventy)
    if (nrow(Data)>3){
      model75<-lm(OCM ~ DMax.D, data=Data)
      ExtS[i,4]<-coef(model75)[1] + 100*coef(model75)[2]}
    
    Data<- subset(Data,Data$DMax.D<=fifhty)
    if (nrow(Data)>3){
      model50<-lm(OCM ~ DMax.D, data=Data)
      ExtS[i,5]<-coef(model50)[1] + 100*coef(model50)[2]}
    
    Data<- subset(Data,Data$DMax.D<=twentififty)
    if (nrow(Data)>3){
      model25<-lm(OCM ~ DMax.D, data=Data)
      ExtS[i,6]<-coef(model25)[1] + 100*coef(model25)[2]}}
  
  
  #############
  # we test the correlation between the stock at 1m estimated from real data and the models
  # and the error of the models
  ##############
  
  
  CorMat <- cor(na.omit(ExtS[, c(2:6)]), method = "pearson")
  
  
  corrplot::corrplot(
    CorMat,
    type = "upper",
    order = "hclust",
    tl.col = "black",
    tl.srt = 45
  )
  
  #write.csv(CorMat,file.path(Folder,"Corr.Extrapolation.csv"),sep=";", dec=",")
  
  
  ExtS <- ExtS |> dplyr::mutate (Error.90 = ((EXT.100 - EXT.90) * 100) / EXT.100)
  ExtS <- ExtS |> dplyr::mutate (Error.75 = ((EXT.100 - EXT.75) * 100) / EXT.100)
  ExtS <- ExtS |> dplyr::mutate (Error.50 = ((EXT.100 - EXT.50) * 100) / EXT.100)
  ExtS <- ExtS |> dplyr::mutate (Error.25 = ((EXT.100 - EXT.25) * 100) / EXT.100)
  
  summary(ExtS)
  
  
  #Global Error
  m.ExtS <- ExtS[, c(1, 7:10)]
  m.ExtS <- reshape::melt(m.ExtS, id = c("Core.ID"))
  
  library("ggplot2")
  
  P1<-ggplot2::ggplot(m.ExtS, aes(variable, value)) + ylab("% of deviation from observed value") + xlab("% of standar depth") +
    geom_boxplot() +
    geom_jitter() +
    scale_x_discrete(labels=c("Error.90" = "90%", "Error.75" = "75%", "Error.50" = "50%", "Error.25" = "25%"))+
    theme(
      #axis.title.x = element_blank(),
      #axis.text.x = element_blank(),
      axis.ticks.x = element_blank()
    )
  
  P2<-ggplot2::ggplot(ExtS, aes(ExtS[, 2], ExtS[, 3])) + xlab("Observed Stock") + ylab("Extrapolated stock") +
    geom_point(aes(ExtS[, 2], ExtS[, 3], color = "90%"), size = 2) +
    geom_point(aes(ExtS[, 2], ExtS[, 4], color = "75%"), size = 2) +
    geom_point(aes(ExtS[, 2], ExtS[, 5], color = "50%"), size = 2) +
    geom_point(aes(ExtS[, 2], ExtS[, 6], color = "25%"), size = 2) +
    theme(text = element_text(size = 15)) +
    #geom_text_repel(aes(label=A[,1]), size=4)+
    xlim(0, 5) + ylim(0, 5) +
    geom_abline()
  
  Extrapolation_plot<-gridExtra::grid.arrange(P1,P2, ncol=2)
  
  return(ExtS)
  return(Extrapolation_plot)
  
}# changed from BlueCarbon library to get begative error values

estimate_flux<- function(df=NULL,TimeFrame=100) {
  
  # class of the dataframe
  if (is.data.frame(df)==FALSE) {stop("The data provided is not class data.frame, please check data and transforme")}
  if (is.numeric(TimeFrame)==FALSE) {stop("The TimeFrame provided is not class numeric, please check data and transforme")}
  
  # name of the columns
  if ("Core.ID" %in% colnames(df)==FALSE) {stop("There is not column named Core.ID. Please, check necessary columns in functions documentation")}
  if ("DMin.D" %in% colnames(df)==FALSE) {stop("There is not column named Min.D. Please, check necessary columns in functions documentation")}
  if ("DMax.D" %in% colnames(df)==FALSE) {stop("There is not column named Max.D. Please, check necessary columns in functions documentation")}
  if ("DDBD" %in% colnames(df)==FALSE) {stop("There is not column named DBD. Please, check necessary columns in functions documentation")}
  if ("POC" %in% colnames(df)==FALSE) {stop("There is not column named fOC. Please, check necessary columns in functions documentation")}
  if ("Age" %in% colnames(df)==FALSE) {stop("There is not column named Age. Please, check necessary columns in functions documentation")}
  
  # class of the columns
  if (is.numeric(df$DMin.D)==FALSE) {stop("Minimum depth data is not class numeric, please check")}
  if (is.numeric(df$DMax.D)==FALSE) {stop("Maximum depth data is not class numeric, please check")}
  if (is.numeric(df$DDBD)==FALSE) {stop("Dry Bulk Density data is not class numeric, please check")}
  if (is.numeric(df$POC)==FALSE) {stop("Organic carbon data is not class numeric, please check")}
  if (is.numeric(df$Age)==FALSE) {stop("Age data is not class numeric, please check")}
  
  
  
  #select those cores with chronological models
  df<-df[!is.na(df$Age),]
  df<-df[!is.na(df$POC),]
  
  df<-estimate_h (df)
  
  X<-split(df, df$Core.ID)
  
  BCF <- data.frame(Core.ID=character(),
                    F.WC=numeric(),
                    A.Max=numeric(),
                    Flux=numeric())
  
  
  for(i in 1:length(X)) {
    BCF[i,1]<-names(X[i])
    Data<-as.data.frame(X[i])
    colnames(Data)<-colnames(df)
    
    if(nrow(Data)<3) next
    
    else{
      
      Data <-Data |> dplyr::mutate (OCgcm2 = DDBD*(POC/100)*h)
      
      #estimation of the average carbon flux for the whole core (OC stock/Max Age)
      BCF[i,2]<-(sum(Data[,which(colnames(Data)=="OCgcm2")]))/max(Data$Age)
      BCF[i,3]<-max(Data$Age)
      
      #estimation of the average carbon flux for the selected TimeFrame (OC stock last 100 yrs/TimeFrame)
      
      
      
      
      if (max(Data$Age)==TimeFrame) {
        
        BCF[i,4]<-((sum(Data[c(1:(nrow(Data))),which(colnames(Data)=="OCgcm2")]))/TimeFrame)
        
      } else {
        Data<-Data[c(1:(length(which(Data$Age <=TimeFrame))+1)),]
        
        if (nrow(Data)<=1) {next} else {
          
          BCF[i,4]<-(
            (sum(Data[c(1:(nrow(Data)-1)),which(colnames(Data)=="OCgcm2")])+
               (Data[nrow(Data),which(colnames(Data)=="OCgcm2")]/((max(Data$Age)-Data[(nrow(Data)-1),which(colnames(Data)=="Age")])))
             *(TimeFrame-Data[(nrow(Data)-1),which(colnames(Data)=="Age")]))/TimeFrame)
        }}
      
    }}
  return(BCF)
}


# Loading and clean the soil data set ----------------------------------------------------

File<-"Raw.csv"

A<-read.csv(File, header=T, sep=";", dec=".")
A<-as.data.frame(A)

#delete cores from bare sediments

unique(A$Ecosystem)

A<-subset(A, !Ecosystem=="Bare sand")
A<-subset(A, !Genus=="Unvegetated Seagrass")
A<-subset(A, !Genus=="Unvegetated Salt Marsh")
A<-subset(A, !Genus=="Unvegetated")

unique(A$Ecosystem)
unique(A$Genus)

#create folder to save tables and grafs

Folder="BC_PI_Results"
dir.create(Folder)

# Summary initial data

length(unique(A$Core.ID))# number of cores

shapiro.test(unique(A$Compresion))
median(A$Compresion)
min(unique(A$Compresion)[-1])
max(A$Compresion)

### correct depth and dry bulk density by core compression (linear correction)###

A<- A %>% mutate (DMin.D = Min.D/(1-(Compresion/100)))
A<- A %>% mutate (DMax.D = Max.D/(1-(Compresion/100)))
A<- A %>% mutate (DDBD = DBD*(1-(Compresion/100)))


# Organic matter to organic carbon ----------------------------------------

#### Estimate a linear model to predict OC from OM for each ecosystem, specie and station ###
#skip those models with R2<0.5 or P value>0.05

#create a list of dataframes with data from each ecosystem, specie, and station (site)
X<-split(A, A$Ecosystem)
X2<-split(A, A$Genus)
X3<-split(A, A$Site.ID)
X<-c(X,X2,X3)

#create empty table to log model data
OCEst <- data.frame(ID=character(),
                    R2=numeric(),
                    P=numeric(),
                    f=numeric(),
                    int=numeric(),
                    slope=numeric()
                    )

for(i in 1:length(X)) {
  OCEst[i,1]<-names(X[i])
  Data<-as.data.frame(X[i])
  colnames(Data)<-colnames(A)


  #we only model those ecosystem, genus, and station with more than 5 samples were OC and LOI were mwasured
  if((nrow(Data %>% filter_at(vars(OM,OC),all_vars(!is.na(.)))))<5) next


  else{

  model<-lm(OC ~ OM, data=Data)

  if(summary(model)$r.squared<0.5 | glance(model)$p.value>0.05 ) next

  else{

  OCEst[i,2]<-summary(model)$r.squared
  OCEst[i,3]<-glance(model)$p.value
  OCEst[i,4]<-summary(model)$fstatistic[1]
  OCEst[i,5]<-model$coefficients[1]
  OCEst[i,6]<-model$coefficients[2]

  }}
} #estimate linear correlation between om and oc in our samples

rownames(OCEst)<-OCEst$ID

write.csv(OCEst,file.path(Folder,"OM-OC_lm.csv")) # Table S1 supp material

## Use the models estimated to predict OC from OM content for those samples with no OC data
#If there is OC data for that sample, keep original OC data,
  #else use function estimated for that Site.ID
    #else function for specie
      #else function for ecosystem

A$POC<- NA

for(i in 1:nrow(A)) {

if (is.na(A[i,which( colnames(A)=="OC" )])==FALSE)
  {A[i,which( colnames(A)=="POC" )]<-A[i,which( colnames(A)=="OC" )]}

    else { if (is.na(OCEst[which(rownames(OCEst)==(A[i,which( colnames(A)=="Site.ID" )])),which(colnames(OCEst)=="int")])==FALSE)

            {A[i,which( colnames(A)=="POC" )]<-
              OCEst[which(rownames(OCEst)==(A[i,which( colnames(A)=="Site.ID" )])),which(colnames(OCEst)=="int" )]+
              (OCEst[which(rownames(OCEst)==(A[i,which( colnames(A)=="Site.ID" )])),which(colnames(OCEst)=="slope" )])*
              A[i,which( colnames(A)=="OM" )] }

            else{ if (is.na(OCEst[which(rownames(OCEst)==(A[i,which( colnames(A)=="Genus" )])),which(colnames(OCEst)=="int")])==FALSE)

                   {A[i,which( colnames(A)=="POC" )]<-
                     OCEst[which(rownames(OCEst)==(A[i,which( colnames(A)=="Genus" )])),which(colnames(OCEst)=="int" )]+
                     (OCEst[which(rownames(OCEst)==(A[i,which( colnames(A)=="Genus" )])),which(colnames(OCEst)=="slope" )])*
                     A[i,which( colnames(A)=="OM" )]}

              else {A[i,which( colnames(A)=="POC" )]<-
                OCEst[which(rownames(OCEst)==(A[i,which( colnames(A)=="Ecosystem" )])),which(colnames(OCEst)=="int" )]+
                (OCEst[which(rownames(OCEst)==(A[i,which( colnames(A)=="Ecosystem" )])),which(colnames(OCEst)=="slope" )])*
                A[i,which( colnames(A)=="OM" )]}}}

            } #Apply linear correlations to estimate OC from OM 

## when OM very low, the estimation can give negative values of OC. We change negative values for 0.

A$POC[A$POC < 0] <- 0




# summary of geochemical variables ----------------------------------------
aggregate(A$DDBD, by=list(A$Ecosystem), FUN=mean, na.rm=T)
aggregate(A$DDBD, by=list(A$Ecosystem), FUN=median, na.rm=T)
aggregate(A$DDBD, by=list(A$Ecosystem), FUN=std.error)
aggregate(A$DDBD, by=list(A$Ecosystem), FUN=contar)


aggregate(A$POC, by=list(A$Ecosystem), FUN=mean, na.rm=T)
aggregate(A$POC, by=list(A$Ecosystem), FUN=median, na.rm=T)
aggregate(A$POC, by=list(A$Ecosystem), FUN=std.error)
aggregate(A$POC, by=list(A$Ecosystem), FUN=contar)



# Soil C stocks estimation 1m ---------------------------------------------

# we eliminate those cores <25cm

X<-split(A, A$Core.ID)

names_A<-colnames(A)
A<-as.data.frame (X[1])
colnames(A)<-names_A

for (i in 2:length(X)){
  
  Data<-as.data.frame (X[i])
  colnames(Data)<-names_A

  if(max(Data$Max.D>25)){
    A<-rbind(A,Data)
  }}

# estimate stocks for whole core and at 1m
BCS<-estimate_oc_stock(A, core="Core.ID", mind="DMin.D", maxd="DMax.D", dbd="DDBD", oc="POC")#estimate_stock estandarize the stock at 1m by default 

write.csv(BCS,file.path(Folder,"Stock_core.csv"))

# number and percentage of cores over 1m depth

# create data.frame with core ID, ecosystem, category and max depth
X<-split(A, A$Core.ID)

max_depth<-data.frame(Core.ID=character(),
                      Ecosystem=character(),
                      Category=character(),
                      max=numeric())

for (i in 1:length(X)){
  Data<-as.data.frame(X[i])
  colnames(Data)<-colnames(A)
  
  max_depth[i,1]<-names(X[i])
  max_depth[i,2]<-Data[1,"Ecosystem"]
  
  if (max_depth[i,2]=="Seagrass") {
    max_depth[i,3]<-Data[1,"Genus"]
  }
  
  if (max_depth[i,2]=="Salt Marsh") {
    max_depth[i,3]<-Data[1,"Tidal.R"]
  }
  
  max_depth[i,4]<-max(Data[,"Max.D"])
  
}

sum(max_depth$max<30)
sum(max_depth$max>50)

sum(max_depth$max<30)*100/(nrow(max_depth))
sum(max_depth$max>50)*100/(nrow(max_depth))


# mean length by category

maxp1<-ggplot(max_depth, aes(Ecosystem, max))+
  geom_boxplot(aes(color=Ecosystem))+
  geom_jitter(aes(color=Ecosystem))+
  ylim(0,150)+
  ylab("Core length (cm)")+
  geom_hline(yintercept=100)+
  geom_hline(yintercept=90, color="grey40", linetype="dashed")+
  geom_hline(yintercept=75, color="grey40", linetype="dashed")+
  geom_hline(yintercept=50, color="grey40", linetype="dashed")+
  geom_hline(yintercept=25, color="grey40", linetype="dashed")+
  theme(axis.title.x = element_blank())

max_depth$Category <- factor(max_depth$Category, levels = c("Posidonia oceanica","Cymodocea nodosa","Zostera marina","Zostera noltii", "High", "Medium", "Low", "Microtidal"),ordered = TRUE)


maxp2<-ggplot(max_depth, aes(Category, max))+
  geom_boxplot(aes(color=Category))+
  geom_jitter(aes(color=Category))+
  ylim(0,150)+
  ylab("Core length (cm)")+
  scale_x_discrete(labels = function(x) str_wrap(x, width = 10))+
  geom_hline(yintercept=100)+
  geom_hline(yintercept=90, color="grey40", linetype="dashed")+
  geom_hline(yintercept=75, color="grey40", linetype="dashed")+
  geom_hline(yintercept=50, color="grey40", linetype="dashed")+
  geom_hline(yintercept=25, color="grey40", linetype="dashed")+
  theme(legend.position = "none", 
        axis.title.x = element_blank())



#check extrapolations


ErrorExt<-test_extrapolation(A, core="Core.ID", mind="DMin.D", maxd="DMax.D", dbd="DDBD", oc="POC")

# add columns of ecosystem and cathegory
for (i in 1:nrow(ErrorExt)) {
  
  ErrorExt[i,"Ecosystem"]<-A[which(A$Core.ID==ErrorExt[i,"core"])[1], "Ecosystem"]
  
  if (ErrorExt[i,"Ecosystem"]=="Seagrass") {
    
    ErrorExt[i,"Cathegory"]<-A[which(A$Core.ID==ErrorExt[i,"core"])[1], "Genus"]
  }
  
  if (ErrorExt[i,"Ecosystem"]=="Salt Marsh") {
    
    ErrorExt[i,"Cathegory"]<-A[which(A$Core.ID==ErrorExt[i,"core"])[1], "Tidal.R"]
  }
  
}

shapiro.test(ErrorExt$Error.25) ### normality (>0.05 normal, <0.05 no normal)
#the error distribution is not normal, give median instead of average

median(ErrorExt$error_25, na.rm = T)


#> median(ErrorExt$error_90, na.rm = TRUE)
#[1] 5.199067
#> median(ErrorExt$error_75, na.rm = TRUE)
#[1] 7.834555
#> median(ErrorExt$error_50, na.rm = TRUE)
#[1] 13.38096
#> median(ErrorExt$error_25, na.rm = TRUE)
#[1]  20.47377


m.ErrorExt <- reshape::melt(ErrorExt[,c(1,11:16)], id = c("core", "Ecosystem", "Cathegory"))



maxp3<-ggplot(m.ErrorExt, aes(variable, value)) + ylab("% of deviation") + xlab("Core length used for the extrapolation (cm)") +
  geom_boxplot(aes(color=Ecosystem)) +
  #geom_jitter() +
  scale_x_discrete(labels=c("Error.90" = "90", "Error.75" = "75", "Error.50" = "50", "Error.25" = "25"))+
  ylim(-100,100)+
  geom_hline(yintercept=0)+
  theme(legend.position = "none", 
    #axis.title.x = element_blank(),
    #axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )

sum(m.ErrorExt$Ecosystem=="Seagrass"&m.ErrorExt$variable=="error_25"&!is.na(m.ErrorExt$value))



ext_plot<-grid.arrange(maxp3, maxp1, maxp2, ncol=2, # fig S4 supp material
             layout_matrix = rbind(c(1, 2),
                                   c(3, 3)))
ggsave(path = Folder,"Extrapolation plot.jpg",ext_plot, units="cm", width = 20, height = 20)


# C fluxes estimation AV 100 years ----------------------------------------


BCF<-estimate_seq_rate(A, core="Core.ID", mind="DMin.D", maxd="DMax.D", dbd="DDBD", oc="POC", age= "Age")

write.csv(BCF,file.path(Folder,"Flux_core.csv"))

# retrieve cores younger than 100 but older than 80 years

for (i in 1:nrow(BCF)) {
  
  if (BCF[i,3] > 80) {BCF[i,4]<-BCF[i, 2]}
  
}


# Check dating method


BCF_method<-BCF
BCF_method$method<-NA

for (i in 1:nrow(BCF_method)) {


core<-subset(A, Core.ID == BCF_method[i,"core"])
core<-core[!is.na(core$Raw.Age),]

if (length(unique(core$D.Method))==1) {BCF_method[i,"method"]<-unique(core$D.Method)}
if (length(unique(core$D.Method))==2) {BCF_method[i,"method"]<-"210Pb and 14C"}}

#% of fluxes estimated from cores were only 14C was available

table(BCF_method$method)
7*100/nrow(BCF_method)



## differences between 100, 75 and 50 time frames ####

BCF_50 <-estimate_seq_rate(A, core="Core.ID", mind="DMin.D", maxd="DMax.D", dbd="DDBD", oc="POC", age= "Age", timeframe = 50)
BCF_100 <-estimate_seq_rate(A, core="Core.ID", mind="DMin.D", maxd="DMax.D", dbd="DDBD", oc="POC", age= "Age", timeframe = 100)
BCF_200 <-estimate_seq_rate(A, core="Core.ID", mind="DMin.D", maxd="DMax.D", dbd="DDBD", oc="POC", age= "Age", timeframe = 200)
BCF_500 <-estimate_seq_rate(A, core="Core.ID", mind="DMin.D", maxd="DMax.D", dbd="DDBD", oc="POC", age= "Age", timeframe = 500)
BCF_1000 <-estimate_seq_rate(A, core="Core.ID", mind="DMin.D", maxd="DMax.D", dbd="DDBD", oc="POC", age= "Age", timeframe = 1000)
BCF_2000 <-estimate_seq_rate(A, core="Core.ID", mind="DMin.D", maxd="DMax.D", dbd="DDBD", oc="POC", age= "Age", timeframe = 2000)

comp_flux<-cbind(BCF_50, BCF_100[,4], BCF_200[,4], BCF_500[,4], BCF_1000[,4], BCF_2000[,4])
colnames(comp_flux)<-c("Core.ID", "F.WC", "A.Max", "50", "100", "200", "500", "1000", "2000")

comp_flux_m<-melt(comp_flux[,c(1,4:9)], id=c("Core.ID"))

ggplot(comp_flux_m, aes(x=as.factor(variable), y=as.numeric(value))) +
  ylab(expression(paste("C flux last 100 years (kg"," ", m^-2,yr^-1, ")")))+ xlab("Time frame")+
  geom_boxplot()+
  geom_jitter()
  #theme(axis.title.x=element_blank())
        #axis.text.x=element_blank(),
        #axis.ticks.x=element_blank())

pairwise.wilcox.test(as.numeric(comp_flux_m$value), comp_flux_m$variable,#are significantly different (p < 0.05)
                     p.adjust.method = "BH")

sum_comp_flux<-data.frame(matrix(NA, nrow = 6, ncol = 3))
  
colnames(sum_comp_flux)<-c("time_frame", "mean", "SE")


sum_comp_flux$time_frame<-c(50, 100, 200, 500, 1000, 2000)
sum_comp_flux$mean<-as.vector(colMeans(comp_flux[,c(4:9)], na.rm=TRUE))
sum_comp_flux$SE<-as.numeric(lapply(comp_flux[,c(4:9)], std.error))

comp_flux_m[c(355:360), "variable"]<-sum_comp_flux$time_frame
comp_flux_m[c(355:360), "value"]<-sum_comp_flux$mean
comp_flux_m$SE<-NA
comp_flux_m[c(355:360), "SE"]<-sum_comp_flux$SE
comp_flux_m$ref<-"data"
comp_flux_m[c(355:360), "ref"]<-"mean"


ggplot(sum_comp_flux, aes(time_frame, mean))+
  geom_point()+
  geom_errorbar(aes(ymin=mean-SE, ymax=mean+SE))


plot_S5<-ggplot(comp_flux_m, aes(variable, value))+ # Fig S5 supp material
  ylab(expression(paste("OC flux (kg"," ", m^-2,yr^-1, ")")))+ xlab("Time frame")+
  scale_color_manual(values=c("black", "blue"))+
  geom_point(aes(color=ref))+
  geom_errorbar(aes(ymin=value-SE, ymax=value+SE), color="blue", width=.5, linewidth=1)


ggsave(path = Folder,"Sup S5.jpg",plot_S5, units="cm", width = 12, height = 10)



# Sediment accumulation rates and age of first meter ----------------------

estimate_sed_rate<- function(df=NULL,depth=100) {

        #select those cores with chronological models
        df<-df[!is.na(df$Age),]
        df<-df[!is.na(df$POC),]
        
        X<-split(df, df$Core.ID)
        
        df_out <- data.frame(Core.ID=character(),
                             Ecosystem=character(),
                             sed_rate=numeric(),
                             Age=numeric())
        
        
        for(i in 1:length(X)) {
          df_out[i,1]<-names(X[i])
          Data<-as.data.frame(X[i])
          colnames(Data)<-colnames(df)
          df_out[i,2]<-Data[1,"Ecosystem"]
          
          if(nrow(Data)<3) next
          
          else{
            
            #model age depth to get age at "depth"
            
            model<-lm(Age ~ Max.D, data=Data)
            
            df_out[i,4]<-predict(model, newdata = data.frame(Max.D = depth))
            
            # divide depth by the age
            
            df_out[i,3]<-depth/(predict(model, newdata = data.frame(Max.D = depth)))
         
          }}
        return(df_out)}
        
sed_rate<-estimate_sed_rate(A, depth=100) # sed rate in cm!

aggregate(sed_rate$sed_rate, by=list(sed_rate$Ecosystem), FUN=median, na.rm=T)
aggregate(sed_rate$Age, by=list(sed_rate$Ecosystem), FUN=median, na.rm=T)


# Biomass -----------------------------------------------------------------

File<-"Biomass.csv"

C<-read.csv(File, header=T, sep=";", dec=".")
C<-as.data.frame(C)

length(unique(C$Sample.ID))# number of biomass samples

# summary of biomass per m2
aggregate(C$Above, by=list(C$Ecosystem), FUN=mean, na.rm=T)
aggregate(C$Above, by=list(C$Ecosystem), FUN=median, na.rm=T)
aggregate(C$Above, by=list(C$Ecosystem), FUN=std.error)
aggregate(C$Above, by=list(C$Ecosystem), FUN=contar)


File<-"C.Plant.csv"

D<-read.csv(File, header=T, sep=";", dec=".")
D<-as.data.frame(D)


#average the OC content in different tissues by specie

TOC.Planta <- as.data.frame(tapply(D$TOC, list(D$Specie, D$Tissue), mean))
TOC.Planta.2 <- as.data.frame(tapply(D$TOC, list(D$Ecosystem, D$Tissue), mean))


#Table3 supp material

TOC.Planta_av<-rbind(TOC.Planta, TOC.Planta.2)

temp1 <- as.data.frame(tapply(D$TOC, list(D$Specie, D$Tissue), std.error))
temp2 <- as.data.frame(tapply(D$TOC, list(D$Ecosystem, D$Tissue), std.error))
TOC.Planta_se<-rbind(temp1, temp2)

temp1 <- as.data.frame(tapply(D$TOC, list(D$Specie, D$Tissue), contar))
temp2 <- as.data.frame(tapply(D$TOC, list(D$Ecosystem, D$Tissue), contar))
TOC.Planta_n<-rbind(temp1, temp2)

f_TOC_planta<-cbind(TOC.Planta_av, TOC.Planta_se, TOC.Planta_n)

write.csv(f_TOC_planta,file.path(Folder,"OC plant tissues.csv"))



#estimate OC in biomass

C$Ab.TOC<- NA

species<-rownames(TOC.Planta)

for(i in 1:nrow(C)) {
  #si hay datos pa biomassa abouve de esa especie usa esa especie
  if(C[i,which( colnames(C)=="Specie" )] %in% species==TRUE)

    {if ((is.na(TOC.Planta[which(rownames(TOC.Planta)==C[i,which( colnames(C)=="Specie" )]),which( colnames(TOC.Planta)=="Aboveground" )]))==FALSE)

      {C[i,which( colnames(C)=="Ab.TOC" )]<-C[i,which( colnames(C)=="Above" )]*
                                              ((TOC.Planta[which(rownames(TOC.Planta)==C[i,which( colnames(C)=="Specie" )]),which( colnames(TOC.Planta)=="Aboveground" )]/100))}

      else { if (is.na(TOC.Planta[which(rownames(TOC.Planta)==C[i,which( colnames(C)=="Specie" )]),which( colnames(TOC.Planta)=="Leaves" )])==FALSE)

            {C[i,which( colnames(C)=="Ab.TOC" )]<-C[i,which( colnames(C)=="Above" )]*
            ((TOC.Planta[which(rownames(TOC.Planta)==C[i,which( colnames(C)=="Specie" )]),which( colnames(TOC.Planta)=="Leaves" )]/100))}

            else {next}}}
  #If there is no TOC for that specie use the average of the ecosystem
  else { C[i,which( colnames(C)=="Ab.TOC" )]<-C[i,which( colnames(C)=="Above" )]*
    ((TOC.Planta.2[which(rownames(TOC.Planta.2)==C[i,which( colnames(C)=="Ecosystem" )]),which( colnames(TOC.Planta.2)=="Aboveground" )]/100))}
}



#Abovegroun biomas by month

table(C$X)

C$X <- factor(C$X, levels = c("January","February","March","April", "May", "June", "July", "August", "September", "October", "November", "December"),ordered = TRUE)
plot (table(C$X))

temp<-subset(C, C$Specie=="Zostera marina")
temp<-subset(C, C$Specie=="Zostera marina")

temp$X <- factor(temp$X, levels = c("January","February","March","April", "May", "June", "July", "August", "September", "October", "November", "December"),ordered = TRUE)

plot (table(temp$X))

ggplot(temp, aes( Ab.TOC, X))+
  geom_boxplot()




#Average by site
B_by_Site <-merge(aggregate( Ab.TOC ~ Site.ID, C, mean), aggregate( Ab.TOC ~ Site.ID, C, sd), by = "Site.ID")
colnames(B_by_Site)<-c("Site.ID", "Mean_Biomass", "SD_Biomass")


# Summary BC per station --------------------------------------------------

#Add information about site.ID to stocks and fluxes data.frames

# Stock 1m. Mean and sd by site
BCS$Site.ID<-NA

for (i in 1:nrow(BCS)) {

  Site<- unique(A[c(which(A$Core.ID==BCS[i,which(colnames(BCS)=="core")])),which(colnames(A)=="Site.ID")])
  BCS[i,which(colnames(BCS)=="Site.ID")]<- Site

}

####### Add here the stocks already estimated to BCS


File<-"P.Stocks.csv"

PS<-read.csv(File, header=T, sep=";", dec=".")
PS<-as.data.frame(PS)

BCS<-rbind(BCS,PS)


# Estimate stocks per site

S_by_Site <-merge(aggregate( stock ~ Site.ID, BCS, mean), aggregate( stock ~ Site.ID, BCS, sd), by = "Site.ID")

colnames(S_by_Site)<-c("Site.ID", "Mean_Stock", "SD_Stock")

# Seq rate 100 yr. Mean and sd by site

BCF$Site.ID<-NA

for (i in 1:nrow(BCF)) {

  Site<- unique(A[c(which(A$Core.ID==BCF[i,which(colnames(BCF)=="core")])),which(colnames(A)=="Site.ID")])
  BCF[i,which(colnames(BCF)=="Site.ID")]<- Site

}

####### Add here the stocks already estimated fluxes con cbin a BCF


#File<-"P.Flux.csv"

#PF<-read.csv(File, header=T, sep=";", dec=".")
#PF<-as.data.frame(PF)

#BCF<-rbind(BCF,PF)

F_by_Site <-merge(aggregate( seq_rate ~ Site.ID, BCF, mean), aggregate( seq_rate ~ Site.ID, BCF, sd), by = "Site.ID")

colnames(F_by_Site)<-c("Site.ID", "Mean_Seq_rate", "SD_Seq_rate")

### Summary table per station (Site)

File<-"GInf.csv"

B<-read.csv(File, header=T, sep=";", dec=".")
B<-as.data.frame(B)


### Final data.frame with stocks at 1m and fluxes at 100 yr by station
#with information about region, coast, genus and Ecosystem

library(dplyr)
library(purrr)


BC_PI <-merge(B, S_by_Site, by = "Site.ID", all = TRUE)

BC_PI<-left_join(BC_PI, F_by_Site, by="Site.ID")

BC_PI<-left_join(BC_PI, B_by_Site, by="Site.ID")

### !!!!!!!!!!!! from g cm2 to kg m2

BC_PI [,c(13:16)]<- BC_PI [,c(13:16)]*10

# we eliminate those station with no data for biomass, stocks or fluxes
BC_PI<-BC_PI[!(is.na(BC_PI$Mean_Stock) & is.na(BC_PI$Mean_Seq_rate) & is.na(BC_PI$Mean_Biomass)),]


write.csv(BC_PI,file.path(Folder,"BC_Station.csv"))

#delete row containing information about unvegetated stations
BC_PI<-BC_PI[!(BC_PI$Genus=="Unvegetated Salt Marsh" | BC_PI$Genus=="Unvegetated Seagrass"| BC_PI$Genus=="Unvegetated"),]
BC_PI<-subset(BC_PI, !is.na(BC_PI$Site.ID))


# summary stations

aggregate(BC_PI$Mean_Biomass, by=list(BC_PI$Genus), FUN=contar)


X<-split(BC_PI, BC_PI$Tidal.R)

for (i in 1:length(X)) {
  
  data<-as.data.frame(X[i])
  colnames(data)<-colnames(BC_PI)
  print((X[i]))
  
  print(shapiro.test(data$Mean_Biomass))
  print(shapiro.test(data$Mean_Stock))
  print(shapiro.test(data$Mean_Flux))
}


# Data distribution (Figure 1)-------------------------------------------------------



WM <- map_data("world")


IP <- subset(WM, region %in% c("Portugal","Spain", "Canary Islands", "Azores", "Madeira Islands"))


DBC_PI = filter(BC_PI, !is.na(Mean_Biomass))
BSg<-subset(DBC_PI, Ecosystem=='Seagrass')
BSm<-subset(DBC_PI, Ecosystem=='Salt Marsh')

DBC_PI = filter(BC_PI, !is.na(Mean_Stock))
SSg<-subset(DBC_PI, Ecosystem=='Seagrass')
SSm<-subset(DBC_PI, Ecosystem=='Salt Marsh')

DBC_PI = filter(BC_PI, !is.na(Mean_Seq_rate))
FSg<-subset(DBC_PI, Ecosystem=='Seagrass')
FSm<-subset(DBC_PI, Ecosystem=='Salt Marsh')



BSgp<-ggplot()+
  geom_polygon(data = WM, aes(x=long, y = lat, group = group), fill = "white", color = "black")+
  geom_polygon(data = IP, aes(x=long, y = lat, group = group), fill = "grey", color = "black")+
  coord_map(xlim = c(-30, 5),ylim = c(27, 45))+
  geom_point(aes(BSg$long,BSg$lat), fill="green",pch=21,size=2.5)+
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        #axis.text.x=element_blank(),
        #axis.ticks.x=element_blank()
  )


BSmp<-ggplot()+  ylab("Latitude")+
  geom_polygon(data = WM, aes(x=long, y = lat, group = group), fill = "white", color = "black")+
  geom_polygon(data = IP, aes(x=long, y = lat, group = group), fill = "grey", color = "black")+
  coord_map(xlim = c(-30, 5),ylim = c(27, 45))+
  geom_point(aes(BSm$long,BSm$lat), fill="blue",pch=21,size=2.5)+
  theme(axis.title.x=element_blank(),
        #axis.text.x=element_blank(),
        #axis.ticks.x=element_blank()
  )


SSgp<-ggplot()+ 
  geom_polygon(data = WM, aes(x=long, y = lat, group = group), fill = "white", color = "black")+
  geom_polygon(data = IP, aes(x=long, y = lat, group = group), fill = "grey", color = "black")+
  coord_map(xlim = c(-30, 5),ylim = c(27, 45))+
  geom_point(aes(SSg$long,SSg$lat), fill="green",pch=21,size=2.5)+
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        #axis.text.x=element_blank(),
        #axis.ticks.x=element_blank()
        )

SSmp<-ggplot()+  ylab("Latitude")+
  geom_polygon(data = WM, aes(x=long, y = lat, group = group), fill = "white", color = "black")+
  geom_polygon(data = IP, aes(x=long, y = lat, group = group), fill = "grey", color = "black")+
  coord_map(xlim = c(-30, 5),ylim = c(27, 45))+
  geom_point(aes(SSm$long,SSm$lat), fill="blue",pch=21,size=2.5)+
  theme(axis.title.x=element_blank(),
        #axis.text.x=element_blank(),
        #axis.ticks.x=element_blank()
  )


FSgp<-ggplot()+ xlab("Longitude")+ ylab("Latitude")+
  geom_polygon(data = WM, aes(x=long, y = lat, group = group), fill = "white", color = "black")+
  geom_polygon(data = IP, aes(x=long, y = lat, group = group), fill = "grey", color = "black")+
  coord_map(xlim = c(-30, 5),ylim = c(27, 45))+
  geom_point(aes(FSg$long,FSg$lat), fill="green",pch=21,size=2.5)+
  theme(axis.title.y=element_blank(),
        #axis.text.y=element_blank(),
        #axis.ticks.y=element_blank()
  )


FSmp<-ggplot()+ xlab("Longitude")+ ylab("Latitude")+
  geom_polygon(data = WM, aes(x=long, y = lat, group = group), fill = "white", color = "black")+
  geom_polygon(data = IP, aes(x=long, y = lat, group = group), fill = "grey", color = "black")+
  coord_map(xlim = c(-30, 5),ylim = c(27, 45))+
  geom_point(aes(FSm$long,FSm$lat), fill="blue",pch=21,size=2.5)



grid.arrange(arrangeGrob(BSmp,top="Salt Marsh", left=textGrob("Above Biomass stock",rot = 90,gp=gpar(fontsize=15))),arrangeGrob(BSgp, top="Seagrass"),
             arrangeGrob(SSmp, left=textGrob("OC Soil stock",rot = 90,gp=gpar(fontsize=15))),SSgp,
             arrangeGrob(FSmp, left=textGrob("OC sequestration rate",rot = 90, gp=gpar(fontsize=15))),FSgp,
             ncol=2, nrow=3,
             widths = c(1,0.9))


# Figure 1
final<-arrangeGrob(arrangeGrob(BSmp,top="Salt Marsh", left=textGrob("Above Biomass stock",rot = 90,gp=gpar(fontsize=15))),arrangeGrob(BSgp, top="Seagrass"),
                                arrangeGrob(SSmp, left=textGrob("OC Soil stock",rot = 90,gp=gpar(fontsize=15))),SSgp,
                                arrangeGrob(FSmp, left=textGrob("OC sequestration rate",rot = 90, gp=gpar(fontsize=15))),FSgp,
                                ncol=2, nrow=3,
                                widths = c(1,0.9))


ggsave(path = Folder,"Sampling sites.jpg",final, units="cm", width = 20, height = 20)

#The patchwork package is another option for laying out multiple plots and it also lines up the plot panels
#Unfortunately, patchwork doesn't provide an easy way to add spanning axis titles (like the bottom, left, and right arguments of grid.arrange)
#so we have to manually set the widths for those grobs, relative to the plot grobs.
#https://community.rstudio.com/t/common-axis-title-in-grid-arrange/96353/2



# Summary table -----------------------------------------------------------

# Tabla de areas

File<-"Dis.csv"

Area<-read.csv(File, header=T, sep=";", dec=".")
Area<-as.data.frame(Area)

colnames(Area)<-c("Coast", "Region", "Total", "Posidonia oceanica", "Cymodocea nodosa", "Zostera marina", "Zostera noltii", "Mixed Seagrass", "Low", "Medium", "High", "Microtidal")


Area$Seagrass <- rowSums(Area[ , c(4:8)], na.rm=TRUE)
Area$"Salt Marsh" <- rowSums(Area[ , c(9:12)], na.rm=TRUE)

# Summary table

Summary<- data.frame(Category=character(),
                     Area=numeric(),
                     nS.Biomass=numeric(),
                     SA.Biomass=numeric(),
                     Av.Biomass=numeric(),
                     M.Biomass=numeric(),
                     MAD.Biomass=numeric(),
                     cv.Biomass=numeric(),
                     nS.Stock=numeric(),
                     SA.Stock=numeric(),
                     Av.Stock=numeric(),
                     M.Stock=numeric(),
                     MAD.Stock=numeric(),
                     cv.Stock=numeric(),
                     nS.SeqRate=numeric(),
                     SA.SeqRate=numeric(),
                     Av.SeqRate=numeric(),
                     M.SeqRate=numeric(),
                     MAD.SeqRate=numeric(),
                     cv.SeqRate=numeric())




X<-split(BC_PI, BC_PI$Ecosystem)
Sg<-subset(BC_PI, Ecosystem=='Seagrass')
X2<-split(Sg, Sg$Genus)
Sm<-Sg<-subset(BC_PI, Ecosystem=='Salt Marsh')
X3<-split(Sm, Sm$Tidal.R)

X<-c(X,X2,X3)


# mad

for(i in 1:length(X)) {
  
  
  Data<-as.data.frame(X[i])
  colnames(Data)<-colnames(BC_PI)
  
  Summary[i,1]<-names(X[i])
  Summary[i,2]<-sum(Area[,names(X[i])], na.rm=TRUE)
  
  DataB = subset(Data, !is.na(Mean_Biomass))
  
  Summary[i,3]<-nrow(DataB)
  Summary[i,4]<-Summary[i,3]/Summary[i,2]
  Summary[i,5]<-mean(DataB$Mean_Biomass)
  Summary[i,6]<-median(DataB$Mean_Biomass)
  Summary[i,7]<-mad(DataB$Mean_Biomass,center=median(DataB$Mean_Biomass),na.rm=T,low=F,high=F)
  Summary[i,8]<-coef.var(DataB$Mean_Biomass)
  
  DataS = subset(Data, !is.na(Mean_Stock))
  Summary[i,9]<-nrow(DataS)
  Summary[i,10]<-Summary[i,9]/Summary[i,2]
  Summary[i,11]<-mean(DataS$Mean_Stock)
  Summary[i,12]<-median(DataS$Mean_Stock)
  Summary[i,13]<-mad(DataS$Mean_Stock,center=median(DataS$Mean_Stock),na.rm=T,low=F,high=F)
  Summary[i,14]<-coef.var(DataS$Mean_Stock)
  
  DataF = subset(Data, !is.na(Mean_Seq_rate))
  Summary[i,15]<-nrow(DataF)
  Summary[i,16]<-Summary[i,15]/Summary[i,2]
  Summary[i,17]<-mean(DataF$Mean_Seq_rate)
  Summary[i,18]<-median(DataF$Mean_Seq_rate)
  Summary[i,19]<-mad(DataF$Mean_Seq_rate,center=median(DataF$Mean_Seq_rate),na.rm=T,low=F,high=F)
  Summary[i,20]<-coef.var(DataF$Mean_Seq_rate)

}


write.csv(Summary,file.path(Folder,"Summary table.csv"))



# Comparisons among categories --------------------------------------------

# normal distribution and significat differences among ecosystems

shapiro.test(BC_PI$Mean_Biomass)
shapiro.test(BC_PI$Mean_Stock) ### normality (>0.05 normal, <0.05 no normal)
shapiro.test(BC_PI$Mean_Seq_rate)

ggplot(BC_PI, aes(x=Mean_Biomass)) +
  geom_histogram()



# boxplots and significant diferences
##### WARNING, geom_signif function only tests data from the plot
# if we use ylim to exclude parts of the plot, that will exclude values in those areas from
# the geom_signif test.
# make sure that the level of signif is the same with and without ylim

EB<-ggplot(BC_PI,aes(Ecosystem, Mean_Biomass))+ ylab(expression(paste("Above Biomass OC stock (kg"," ", m^-2, ")")))+
  geom_boxplot()+
  geom_jitter(aes(color=factor(Coast)), alpha = 0.5)+
  scale_colour_manual(values = c("coral1","seagreen", "skyblue", "purple"))+
  ylim(0,2.25)+ 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  
  geom_signif(comparisons=list(c("Salt Marsh", "Seagrass")),
              test = "wilcox.test",
              na.rm = FALSE,
              map_signif_level = TRUE,
              col = "#5b5b5b",
              size = 0.4, #change line size
              y_position = 1,
              tip_length = 0.005)




ESS<-ggplot(BC_PI,aes(Ecosystem, Mean_Stock))+ ylab(expression(paste("Soil OC stock (kg"," ", m^-2, ")")))+
  geom_boxplot()+
  geom_jitter(aes(color=factor(Coast)), alpha = 0.5)+ 
  scale_colour_manual(values = c("coral1","seagreen", "skyblue", "purple"))+
  #ylim(0,55) + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  
  geom_signif(comparisons=list(c("Salt Marsh", "Seagrass")),
              test = "wilcox.test",
              na.rm = FALSE,
              map_signif_level = TRUE,
              col = "#5b5b5b",
              size = 0.4, #change line size
              y_position = 150,
              tip_length = 0.005)


ESF<-ggplot(BC_PI,aes(Ecosystem, Mean_Seq_rate))+ ylab(expression(paste("OC Seq. rates last 100 years (kg"," ", m^-2,yr^-1, ")")))+
  geom_boxplot()+
  geom_jitter(aes(color=factor(Coast)), alpha = 0.5)+
  scale_colour_manual(values = c("coral1","seagreen", "skyblue", "purple"))+
  ylim(0,0.08)+
  theme(axis.title.x=element_blank())+
  
  geom_signif(comparisons=list(c("Salt Marsh", "Seagrass")),
              test = "wilcox.test",
              na.rm = FALSE,
              map_signif_level = TRUE,
              col = "#5b5b5b",
              size = 0.4, #change line size
              y_position = 0.06,
              tip_length = 0.005)



SPE<-grid.arrange(EB,ESS, ESF, nrow=3)
ggsave(path = Folder,"Summary ecosystems.jpg",SPE, units="cm", width = 12.5, height = 20)

# significate diferences between ecosystems
pairwise.wilcox.test(BC_PI$Mean_Biomass, BC_PI$Ecosystem, #### are significantly different (p < 0.05)
                     p.adjust.method = "BH")

pairwise.wilcox.test(BC_PI$Mean_Stock, BC_PI$Ecosystem, #### are significantly different (p < 0.05)
                     p.adjust.method = "BH")

pairwise.wilcox.test(BC_PI$Mean_Seq_rate, BC_PI$Ecosystem, #### are significantly different (p < 0.05)
                     p.adjust.method = "BH")



# Seagrass ----------------------------------------------------------------

#### only pvalues significatives. If not the plot is too overcrowded
##### WARNING, geom_signif function only tests data from the plot
# if we use ylim to exclude parts of the plot, that will exclude values in those areas from
# the geom_signif test.
# make sure that the level of signif is the same with and without ylim


Sg<- subset(BC_PI, Ecosystem=="Seagrass")

Sg<-subset(Sg, !Genus=="Mixed Seagrass")

Sg$Genus <- factor(Sg$Genus, levels = c("Posidonia oceanica","Cymodocea nodosa","Zostera marina","Zostera noltii"),ordered = TRUE)


# normal distribution

shapiro.test(Sg$Mean_Stock) ### normality (>0.05 normal, <0.05 no normal)

ggplot(Sg, aes(x=Mean_Stock)) +
  geom_histogram()

listids <- list()
for (ids in unique(Sg$Genus)){
  subSg <- subset(x=Sg, subset=Genus==ids)
  # apply the rest of your analysis there using subdf, for instance
  listids[[ids]] <- shapiro.test(subSg$Mean_Biomass)

  print(ggplot(subSg, aes(x=Mean_Biomass)) + ggtitle(ids)+
    geom_histogram())

}


# seagrass summary plot

SS<-  ggplot(Sg,aes(Genus, Mean_Stock))+ ylab(expression(paste("Soil C stock (kg"," ", m^-2, ")")))+
  geom_boxplot()+
  geom_jitter(aes(color=factor(Coast)))+
  scale_colour_manual(values = c("coral1","seagreen", "skyblue", "purple"))+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  #ylim(0,55)+
  
  
  geom_signif(comparisons=list(c("Cymodocea nodosa", "Zostera marina"), c("Cymodocea nodosa", "Zostera noltii"), 
                               c("Zostera marina", "Zostera noltii"), c("Posidonia oceanica", "Cymodocea nodosa"),
                               c("Posidonia oceanica", "Zostera marina"), c("Posidonia oceanica", "Zostera noltii")),
              test = "wilcox.test",
              na.rm = FALSE,
              map_signif_level = TRUE,
              col = "#5b5b5b",
              size = 0.4, #change line size
              y_position = 10,
              step_increase = 0.02,
              tip_length = 0.005)
  


SB<-ggplot(Sg,aes(Genus, Mean_Biomass))+ ylab(expression(paste("Biomass C stock (kg"," ", m^-2, ")")))+
  geom_boxplot()+
  geom_jitter(aes(color=factor(Coast)))+
  scale_colour_manual(values = c("coral1","seagreen", "skyblue", "purple"))+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  #ylim(0,2)+
  
  geom_signif(comparisons=list(c("Cymodocea nodosa", "Zostera marina"), c("Cymodocea nodosa", "Zostera noltii"), 
                               c("Zostera marina", "Zostera noltii"), c("Posidonia oceanica", "Cymodocea nodosa"),
                               c("Posidonia oceanica", "Zostera marina"), c("Posidonia oceanica", "Zostera noltii")),
              test = "wilcox.test",
              na.rm = FALSE,
              map_signif_level = TRUE,
              col = "#5b5b5b",
              size = 0.4, #change line size
              y_position = 0.03,
              step_increase = 0.008,
              tip_length = 0.005)


F100<-ggplot(Sg,aes(Genus, Mean_Flux))+ ylab(expression(paste("C flux last 100 years (kg"," ", m^-2,yr^-1, ")")))+
  geom_boxplot()+
  geom_jitter(aes(color=factor(Coast)))+
  scale_colour_manual(values = c("coral1","seagreen", "skyblue", "purple"))+
  scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 10))+
  #ylim(0,0.08)+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_text(face="italic"))+
  
  
  geom_signif(comparisons=list( c("Cymodocea nodosa", "Zostera noltii"), 
                                c("Posidonia oceanica", "Cymodocea nodosa"),
                                c("Posidonia oceanica", "Zostera noltii")),
              test = "wilcox.test",
              na.rm = FALSE,
              map_signif_level = TRUE,
              col = "#5b5b5b",
              size = 0.4, #change line size
              y_position = 0.05,
              step_increase = 0.05,
              tip_length = 0.01)

  

SgPlot<-grid.arrange(SB,SS, F100, nrow=3, top="Seagrass specie")

#ggsave(path = Folder,"Summary seagrass.jpg",SgPlot, units="cm", width = 13, height = 20)


# Salt marshes ------------------------------------------------------------


Sm<-BC_PI[!(BC_PI$Ecosystem=="Seagrass"),]
Sm = filter(Sm, !is.na(Site.ID))


Sm$Tidal.R <- factor(Sm$Tidal.R, levels = c("High","Medium","Low","Microtidal"),ordered = TRUE)


SS<-ggplot(Sm,aes(Tidal.R, Mean_Stock))+ ylab(expression(paste("Soil C stock (kg"," ", m^-2, ")")))+
  geom_boxplot()+
  #ylim(0,50)+
  geom_jitter(aes(color=factor(Coast)))+
  scale_colour_manual(values = c("coral1","seagreen", "skyblue", "purple"))+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+



  geom_signif(comparisons=list(c("Low", "Microtidal"), c("Medium", "Low"), 
                             c("Medium", "Microtidal"), c("High", "Medium"),
                             c("High", "Low"), c("High", "Microtidal")),
            test = "wilcox.test",
            na.rm = FALSE,
            map_signif_level = TRUE,
            col = "#5b5b5b",
            size = 0.4, #change line size
            y_position = 10,
            step_increase = 0.02,
            tip_length = 0.005)



SB<-ggplot(Sm,aes(Tidal.R, Mean_Biomass))+ ylab(expression(paste("Biomass C stock (kg"," ", m^-2, ")")))+
  geom_boxplot()+
  #ylim(0,2)+
  geom_jitter(aes(color=factor(Coast)))+
  scale_colour_manual(values = c("coral1","seagreen", "skyblue", "purple"))+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +


  geom_signif(comparisons=list(c("Medium", "Low"),c("Medium", "Microtidal"), 
                               c("Microtidal", "Low"), c("High", "Medium"),
                               c("High", "Low"), c("High", "Microtidal")),
              test = "wilcox.test",
              na.rm = FALSE,
              map_signif_level = TRUE,
              col = "#5b5b5b",
              size = 0.4, #change line size
              y_position = 0.035,
              step_increase = 0.008,
              tip_length = 0.005)



F100<-ggplot(Sm,aes(Tidal.R, Mean_Flux))+ ylab(expression(paste("C flux last 100 years (kg"," ", m^-2,yr^-1, ")")))+
  geom_boxplot()+
  ylim(0,0.08)+
  geom_jitter(aes(color=factor(Coast)))+
  scale_colour_manual(values = c("coral1","seagreen", "skyblue", "purple"))+
  scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 10))+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_text(face="italic"))+
  
  geom_signif(comparisons=list(c("Medium", "Low"), 
                                c("High", "Medium"),
                               c("High", "Low")),
              test = "wilcox.test",
              na.rm = FALSE,
              map_signif_level = TRUE,
              col = "#5b5b5b",
              size = 0.4, #change line size
              y_position = 0.05,
              step_increase = 0.05,
              tip_length = 0.01)


SmPlot<-grid.arrange(SB,SS, F100, nrow=3, top="Salt Marsh Tidal Range")
#ggsave(path = Folder,"Summary salt marsh.jpg",SmPlot, units="cm", width = 13, height = 20)



# Figure 4 ----------------------------------------------------------------


#### only pvalues significatives. If not the plot is too overcrowded
##### WARNING, geom_signif function only tests data from the plot
# if we use ylim to exclude parts of the plot, that will exclude values in those areas from
# the geom_signif test.
# make sure that the level of signif is the same with and without ylim



ggplot(Sg,aes(Genus, Mean_Stock))+ ylab(expression(paste("Soil OC stock (kg"," ", m^-2, ")")))+
  geom_boxplot()+
  geom_jitter(aes(color=factor(Coast), alpha = 0.5))+
  scale_colour_manual(values = c("coral1","seagreen", "skyblue", "purple"))+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  #ylim(0,55)+
  
  geom_signif(comparisons=list(c("Cymodocea nodosa", "Zostera marina"), c("Cymodocea nodosa", "Zostera noltii"), 
                               c("Zostera marina", "Zostera noltii"), c("Posidonia oceanica", "Cymodocea nodosa"),
                               c("Posidonia oceanica", "Zostera marina"),c("Posidonia oceanica", "Zostera noltii")),
              test = "wilcox.test",
              na.rm = FALSE,
              map_signif_level = TRUE,
              col = "#5b5b5b",
              size = 0.4, #change line size
              y_position = 10,
              step_increase = 0.02,
              tip_length = 0.005)


        SS<-  ggplot(Sg,aes(Genus, Mean_Stock))+ ylab(expression(paste("Soil OC stock (kg"," ", m^-2, ")")))+
          geom_boxplot()+
          geom_jitter(aes(color=factor(Coast), alpha = 0.5))+
          scale_colour_manual(values = c("coral1","seagreen", "skyblue", "purple"))+
          theme(axis.title.x=element_blank(),
                axis.text.x=element_blank(),
                axis.ticks.x=element_blank())+
          ylim(0,55)+
          
          geom_signif(comparisons=list(c("Cymodocea nodosa", "Zostera marina"), c("Cymodocea nodosa", "Zostera noltii"), 
                                       c("Zostera marina", "Zostera noltii"), c("Posidonia oceanica", "Cymodocea nodosa"),
                                       c("Posidonia oceanica", "Zostera marina")),
                      test = "wilcox.test",
                      na.rm = FALSE,
                      map_signif_level = TRUE,
                      col = "#5b5b5b",
                      size = 0.4, #change line size
                      y_position = 10,
                      step_increase = 0.02,
                      tip_length = 0.005)
        





ggplot(Sg,aes(Genus, Mean_Biomass))+ ylab(expression(paste("Biomass OC stock (kg"," ", m^-2, ")")))+
  geom_boxplot()+
  geom_jitter(aes(color=factor(Coast), alpha = 0.5))+
  scale_colour_manual(values = c("coral1","seagreen", "skyblue", "purple"))+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  #ylim(0,2)+
  
  geom_signif(comparisons=list(c("Cymodocea nodosa", "Zostera marina"), c("Cymodocea nodosa", "Zostera noltii"), 
                               c("Zostera marina", "Zostera noltii"), c("Posidonia oceanica", "Cymodocea nodosa"),
                               c("Posidonia oceanica", "Zostera marina"),c("Posidonia oceanica", "Zostera noltii")),
              test = "wilcox.test",
              na.rm = FALSE,
              map_signif_level = TRUE,
              col = "#5b5b5b",
              size = 0.4, #change line size
              y_position = 0.5,
              step_increase = 0.008,
              tip_length = 0.005)
        




        
        SB<-ggplot(Sg,aes(Genus, Mean_Biomass))+ ylab(expression(paste("Biomass OC stock (kg"," ", m^-2, ")")))+
          geom_boxplot()+
          geom_jitter(aes(color=factor(Coast), alpha = 0.5))+
          scale_colour_manual(values = c("coral1","seagreen", "skyblue", "purple"))+
          theme(axis.title.x=element_blank(),
                axis.text.x=element_blank(),
                axis.ticks.x=element_blank())+
          ylim(0,2)+
          
          geom_signif(comparisons=list(c("Posidonia oceanica", "Cymodocea nodosa"),
                                       c("Posidonia oceanica", "Zostera marina"), c("Posidonia oceanica", "Zostera noltii")),
                      test = "wilcox.test",
                      na.rm = FALSE,
                      map_signif_level = TRUE,
                      col = "#5b5b5b",
                      size = 0.4, #change line size
                      y_position = 0.5,
                      step_increase = 0.008,
                      tip_length = 0.005)
        
 ggplot(Sg,aes(Genus, Mean_Seq_rate))+ ylab(expression(paste("OC seq. rate last 100 yr (kg"," ", m^-2,yr^-1, ")")))+
          geom_boxplot()+
          geom_jitter(aes(color=factor(Coast), alpha = 0.5))+
          scale_colour_manual(values = c("coral1","seagreen", "skyblue", "purple"))+
          scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 10))+
          #ylim(0,0.08)+
          theme(axis.title.x=element_blank(),
                axis.text.x=element_text(face="italic"))+
          
   geom_signif(comparisons=list( c("Cymodocea nodosa", "Zostera noltii"), 
                                 c("Posidonia oceanica", "Cymodocea nodosa"),c("Posidonia oceanica", "Zostera noltii")),
                      test = "wilcox.test",
                      na.rm = FALSE,
                      map_signif_level = TRUE,
                      col = "#5b5b5b",
                      size = 0.4, #change line size
                      y_position = 0.05,
                      step_increase = 0.05,
                      tip_length = 0.01)
    
        
        
        F100<-ggplot(Sg,aes(Genus, Mean_Seq_rate))+ ylab(expression(paste("OC seq. rate last 100 yr (kg"," ", m^-2,yr^-1, ")")))+
          geom_boxplot()+
          geom_jitter(aes(color=factor(Coast), alpha = 0.5))+
          scale_colour_manual(values = c("coral1","seagreen", "skyblue", "purple"))+
          scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 10))+
          ylim(0,0.08)+
          theme(axis.title.x=element_blank(),
                axis.text.x=element_text(face="italic"))+
          
          geom_signif(comparisons=list(c("Posidonia oceanica", "Cymodocea nodosa")),
                      test = "wilcox.test",
                      na.rm = FALSE,
                      map_signif_level = TRUE,
                      col = "#5b5b5b",
                      size = 0.4, #change line size
                      y_position = 0.05,
                      step_increase = 0.05,
                      tip_length = 0.01)
        
        SgPlot<-grid.arrange(SB,SS, F100, nrow=3, top="Seagrass species")
        
        ggsave(path = Folder,"Summary seagrass.jpg",SgPlot, units="cm", width = 13, height = 20)
        

        
       
        
        
ggplot(Sm,aes(Tidal.R, Mean_Stock))+ ylab(expression(paste("Soil OC stock (kg"," ", m^-2, ")")))+
          geom_boxplot()+
          #ylim(0,50)+
          geom_jitter(aes(color=factor(Coast), alpha = 0.5))+
          scale_colour_manual(values = c("coral1","seagreen", "skyblue", "purple"))+
          theme(axis.title.x=element_blank(),
                axis.text.x=element_blank(),
                axis.ticks.x=element_blank())+
          
  geom_signif(comparisons=list(c("Medium", "Microtidal"), c("Medium", "Low"),
                               c("Medium", "High"),
                               c("Microtidal", "Low"), 
                               c("High", "Low"),
                               c("High", "Microtidal")),
                      test = "wilcox.test",
                      na.rm = FALSE,
                      map_signif_level = TRUE,
                      col = "#5b5b5b",
                      size = 0.4, #change line size
                      y_position = 20,
                      step_increase = 0.02,
                      tip_length = 0.005)        
        
        
        
        SS<-ggplot(Sm,aes(Tidal.R, Mean_Stock))+ ylab(expression(paste("Soil OC stock (kg"," ", m^-2, ")")))+
          geom_boxplot()+
          ylim(0,50)+
          geom_jitter(aes(color=factor(Coast), alpha = 0.5))+
          scale_colour_manual(values = c("coral1","seagreen", "skyblue", "purple"))+
          theme(axis.title.x=element_blank(),
                axis.text.x=element_blank(),
                axis.ticks.x=element_blank())+
          
          geom_signif(comparisons=list(),
                      test = "wilcox.test",
                      na.rm = FALSE,
                      map_signif_level = TRUE,
                      col = "#5b5b5b",
                      size = 0.4, #change line size
                      y_position = 20,
                      step_increase = 0.02,
                      tip_length = 0.005)
        
ggplot(Sm,aes(Tidal.R, Mean_Biomass))+ ylab(expression(paste("Biomass OC stock (kg"," ", m^-2, ")")))+
          geom_boxplot()+
          #ylim(0,2)+
          geom_jitter(aes(color=factor(Coast), alpha = 0.5))+
          scale_colour_manual(values = c("coral1","seagreen", "skyblue", "purple"))+
          theme(axis.title.x=element_blank(),
                axis.text.x=element_blank(),
                axis.ticks.x=element_blank()) +
          
  geom_signif(comparisons=list(c("Medium", "Microtidal"), 
                               c("Medium", "Low"),
                               c("Medium", "High"),
                               c("Microtidal", "Low"), 
                               c("High", "Low")
                               #,c("High", "Microtidal")
                               ),
                      test = "wilcox.test",
                      na.rm = FALSE,
                      map_signif_level = TRUE,
                      col = "#5b5b5b",
                      size = 0.4, #change line size
                      y_position = 0.6,
                      step_increase = 0.008,
                      tip_length = 0.005)
        
     ##### add significance *** to high and medium and * to medium and low
        
        SB<-ggplot(Sm,aes(Tidal.R, Mean_Biomass))+ ylab(expression(paste("Above Biomass OC stock (kg"," ", m^-2, ")")))+
          geom_boxplot()+
          ylim(0,2)+
          geom_jitter(aes(color=factor(Coast), alpha = 0.5))+
          scale_colour_manual(values = c("coral1","seagreen", "skyblue", "purple"))+
          theme(axis.title.x=element_blank(),
                axis.text.x=element_blank(),
                axis.ticks.x=element_blank()) +
          
          geom_signif(comparisons=list(c("Medium", "Microtidal"),
                                       c("Microtidal", "Low"), 
                                       c("Medium", "Low"),
                                       c("High", "Low")),
                      test = "wilcox.test",
                      na.rm = FALSE,
                      map_signif_level = TRUE,
                      col = "#5b5b5b",
                      size = 0.4, #change line size
                      y_position = 0.6,
                      step_increase = 0.008,
                      tip_length = 0.005)
        

ggplot(Sm,aes(Tidal.R, Mean_Seq_rate))+ ylab(expression(paste("OC seq. rate last 100 yr (kg"," ", m^-2,yr^-1, ")")))+
          geom_boxplot()+
          ylim(0,0.08)+
          geom_jitter(aes(color=factor(Coast), alpha = 0.5))+
          scale_colour_manual(values = c("coral1","seagreen", "skyblue", "purple"))+
          scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 10))+
          theme(axis.title.x=element_blank(),
                axis.text.x=element_text(face="italic"))+
          
  geom_signif(comparisons=list( c("Medium", "Low"),
                               c("Medium", "High"),
                               c("High", "Low")),
                      test = "wilcox.test",
                      na.rm = FALSE,
                      map_signif_level = TRUE,
                      col = "#5b5b5b",
                      size = 0.4, #change line size
                      y_position = 0.05,
                      step_increase = 0.05,
                      tip_length = 0.01)


        F100<-ggplot(Sm,aes(Tidal.R, Mean_Seq_rate))+ ylab(expression(paste("OC seq. rate last 100 yr (kg"," ", m^-2,yr^-1, ")")))+
          geom_boxplot()+
          #ylim(0,0.08)+
          geom_jitter(aes(color=factor(Coast), alpha = 0.5))+
          scale_colour_manual(values = c("coral1","seagreen", "skyblue", "purple"))+
          scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 10))+
          theme(axis.title.x=element_blank(),
                axis.text.x=element_text(face="italic"))+
          
          geom_signif(comparisons=list(),
                      test = "wilcox.test",
                      na.rm = FALSE,
                      map_signif_level = TRUE,
                      col = "#5b5b5b",
                      size = 0.4, #change line size
                      y_position = 0.05,
                      step_increase = 0.05,
                      tip_length = 0.01)
 
        SmPlot<-grid.arrange(SB,SS, F100, nrow=3, top="Salt Marsh Tidal Range")
        ggsave(path = Folder,"Summary salt marsh.jpg",SmPlot, units="cm", width = 13, height = 20)
        

        

# Effect of tidal range in seagrass meadows --------------------------------

temp<-subset(Sg, Sg$Genus=="Posidonia oceanica")      
        
        pairwise.wilcox.test(temp$Mean_Biomass, temp$Tidal.R, #### are significantly different (p < 0.05)
                             p.adjust.method = "BH")
        
        pairwise.wilcox.test(temp$Mean_Stock, temp$Tidal.R, #### are significantly different (p < 0.05)
                             p.adjust.method = "BH")
        
        pairwise.wilcox.test(temp$Mean_Seq_rate, temp$Tidal.R, #### are significantly different (p < 0.05)
                             p.adjust.method = "BH")

        
        temp<-subset(temp, !is.na(temp$Mean_Stock))  
        plyr::count(temp, "Tidal.R")
        
        
        
        
        SgT<-ggplot(Sg,aes(Genus, Mean_Stock, color=factor(Tidal.R)))+ xlab("Seagrass species")+
          ylab(expression(paste("Soil OC stock (kg"," ", m^-2, ")")))+
          geom_boxplot()+
          geom_point(position = position_jitterdodge(), alpha = 0.5)+
          ylim(0,50)+
          scale_x_discrete(labels = function(x) str_wrap(x, width = 10))   
        
        ggsave(path = Folder,"Tidal seagrass.jpg", SgT, units="cm", width = 12, height = 7)      


# Comparison with literature data ----------------------------------------


File<-"literature.csv"

L<-read.csv(File, header=T, sep=";", dec=".")
L<-as.data.frame(L)

# seagrass
LSg<-subset(L, Ecosystem=="Seagrass")
LSg_ib<-subset(L, Specie=="Posidonia oceanica" | Specie=="Zostera noltii" | Specie=="Zostera marina"| Specie=="Cymodocea nodosa")

comp1<-Sg[,c("Ecosystem", "Genus", "lat","long", "Mean_Biomass", "Mean_Stock", "Mean_Seq_rate")]
colnames(comp1)<-c("Ecosystem", "Specie","lat","long", "Biomass",  "Stock", "Seq.rate")
comp1$ref<-"This Study"
comp2<-LSg_ib[,c("Ecosystem", "Specie","lat","long", "Biomass", "Stock", "Flux.100")]
comp2$ref<-"Literature"
colnames(comp2)<-c("Ecosystem", "Specie","lat","long", "Biomass", "Stock", "Seq.rate", "ref")
comp_sg<-rbind(comp1, comp2)

comp_sg2<-comp_sg
comp_sg2$Specie<-"Seagrass"

comp_sg_f<-rbind(comp_sg,comp_sg2)


pl_sg1<-  
  ggplot(comp_sg_f,aes(Specie, Biomass, color=factor(ref)))+ ylab(expression(paste("AG Biomass OC stock (kg"," ", m^-2, ")")))+
  geom_boxplot()+
  geom_point(position = position_jitterdodge(), alpha = 0.1)+
  ylim(0,2)+
  scale_x_discrete(labels = function(x) str_wrap(x, width = 10))

pairwise.wilcox.test(comp_sg$Biomass, comp_sg$ref, #### are significantly different (p < 0.05)
                     p.adjust.method = "BH")



pl_sg2<-  ggplot(comp_sg_f,aes(Specie, Stock, color=factor(ref)))+ ylab(expression(paste("Soil OC stock (kg"," ", m^-2, ")")))+
  geom_boxplot()+
  geom_point(position = position_jitterdodge(), alpha = 0.5) +
  ylim(0,50)+
  scale_x_discrete(labels = function(x) str_wrap(x, width = 10))


pl_sg3<-ggplot(comp_sg_f,aes(Specie, Seq.rate, color=factor(ref)))+ ylab(expression(paste("OC Seq. rate  (kg"," ", m^-2,y^-1, ")")))+
  geom_boxplot()+
  geom_point(position = position_jitterdodge(), alpha = 0.5) +
  ylim(0,0.1)+
  scale_x_discrete(labels = function(x) str_wrap(x, width = 10))

# significant differences between our data and published data

temp<-subset(comp_sg_f, Specie == "Zostera noltii")

pairwise.wilcox.test(temp$Biomass, temp$ref, #### are significantly different (p < 0.05)
                     p.adjust.method = "BH")

pairwise.wilcox.test(temp$Stock, temp$ref, #### are significantly different (p < 0.05)
                     p.adjust.method = "BH")

pairwise.wilcox.test(temp$Mean_Seq_rate, temp$ref, #### are significantly different (p < 0.05)
                     p.adjust.method = "BH")

# salt marshes

LSm<-subset(L, Ecosystem=="Salt Marsh")

comp1<-Sm[,c("Ecosystem", "Tidal.R", "lat","long", "Mean_Biomass", "Mean_Stock", "Mean_Seq_rate")]
colnames(comp1)<-c("Ecosystem", "Tidal.R","lat","long", "Biomass",  "Stock", "Seq.rate")
comp1$ref<-"This Study"
comp2<-LSm[,c("Ecosystem", "Tidal.R", "lat","long","Biomass", "Stock", "Flux.100")]
comp2$ref<-"Literature"
colnames(comp2)<-c("Ecosystem", "Tidal.R","lat","long", "Biomass", "Stock", "Seq.rate", "ref")
comp_sm<-rbind(comp1, comp2)
comp_sm<-comp_sm[!comp_sm$Tidal.R=="",]

comp_sm2<-comp_sm
comp_sm2$Tidal.R<-"Salt Marsh"

comp_sm_f<-rbind(comp_sm,comp_sm2)

pl_sm1<-  ggplot(comp_sm_f,aes(Tidal.R, Biomass, color=factor(ref)))+ ylab(expression(paste("AG Biomass OC stock (kg"," ", m^-2, ")")))+
  geom_boxplot()+
  geom_point(position = position_jitterdodge(), alpha = 0.1)+
  ylim(0,2)


pl_sm2<-ggplot(comp_sm_f,aes(Tidal.R, Stock, color=factor(ref)))+ ylab(expression(paste("Soil OC stock (kg"," ", m^-2, ")")))+
  geom_boxplot()+
  geom_point(position = position_jitterdodge(), alpha = 0.5)+
  ylim(0,40)


pl_sm3<-ggplot(comp_sm_f,aes(Tidal.R, Seq.rate, color=factor(ref)))+ ylab(expression(paste("OC Seq. rate  (kg"," ", m^-2,y^-1, ")")))+
  geom_boxplot()+
  geom_point(position = position_jitterdodge(), alpha = 0.5) 


comp_lit<-grid.arrange(pl_sg1, pl_sg2, pl_sg3, pl_sm1, pl_sm2, pl_sm3, nrow=2)

ggsave(path = Folder,"Comparissom with literature data.jpg",comp_lit, units="cm", width = 40, height = 20)


# significant differences between our data and published data



temp<-subset(comp_sm_f, Tidal.R == "Medium")

pairwise.wilcox.test(temp$Biomass, temp$ref, #### are significantly different (p < 0.05)
                     p.adjust.method = "BH")

pairwise.wilcox.test(temp$Stock, temp$ref, #### are significantly different (p < 0.05)
                     p.adjust.method = "BH")

pairwise.wilcox.test(temp$Seq.rate, temp$ref, #### are significantly different (p < 0.05)
                     p.adjust.method = "BH")


# significant differences between our data and published data whole ecosystem


temp<-comp_sg_f

pairwise.wilcox.test(temp$Biomass, temp$ref, #### are significantly different (p < 0.05)
                     p.adjust.method = "BH")

pairwise.wilcox.test(temp$Stock, temp$ref, #### are significantly different (p < 0.05)
                     p.adjust.method = "BH")

pairwise.wilcox.test(temp$Seq.rate, temp$ref, #### are significantly different (p < 0.05)
                     p.adjust.method = "BH")



# Cymodocea and Posidonia biomass maps (FS1 and S2) -----------------------


# mapa cymo


cymo_bio = filter(comp_sg, !is.na(Biomass) & Specie=="Cymodocea nodosa")
cymo_posi = filter(comp_sg, !is.na(Biomass) & Specie=="Posidonia oceanica")



#Cymo Canarias
median(arrange(cymo_bio, lat)[c(1:18),"Biomass"])

# Inventory no Canary Island
median(arrange((subset(cymo_bio, ref=="This Study")), lat)[-c(1:18),"Biomass"] )

#other published data
median(subset(cymo_bio, ref=="Published")[,"Biomass"])

# map Cymo

cymo1<-ggplot(cymo_bio, aes(long,lat, color=Biomass))+ ylab("Latitude")+ xlab("Longitude")+
  geom_polygon(data = WM, aes(x=long, y = lat, group = group), fill = "white", color = "black")+
  geom_polygon(data = IP, aes(x=long, y = lat, group = group), fill = "grey", color = "black")+
  coord_map(xlim = c(-20, 40),ylim = c(27, 45))+
  labs(tag = "A") +
  scale_colour_gradient(high = "#132B43", low = "#56B1F7")+
  geom_point(aes(shape=cymo_bio$ref),size=2.5)+
  theme(plot.tag = element_text(),
        legend.title=element_blank())

cymo2<-ggplot(cymo_bio, aes(long, Biomass))+ xlab("Longitude")+ 
  ylab(expression(paste("AG Biomass C stock (kg"," ", m^-2, ")")))+
  labs(tag = "B") +
  geom_point(aes(color=ref))+
  theme(legend.position = "none", 
        plot.tag = element_text())

cymo3<-ggplot(cymo_bio, aes(long, Biomass))+ xlab("Longitude")+ 
  geom_point(aes(color=ref))+
  labs(tag = "C") +
  ylim(0, 0.5)+
  theme(axis.title.y=element_blank())


CymoBiomass<-grid.arrange(cymo1, cymo2, cymo3, nrow=2, 
                          layout_matrix = cbind(c(1,2), c(1,3)),
                          top="Cymodocea nodosa biomass")

ggsave(path = Folder,"Cymo biomass.jpg",CymoBiomass, units="cm", width = 13, height = 20)





#Posidonia maps

posi1<-ggplot(cymo_posi, aes(long,lat, color=Biomass))+ ylab("Latitude")+ xlab("Longitude")+
  geom_polygon(data = WM, aes(x=long, y = lat, group = group), fill = "white", color = "black")+
  geom_polygon(data = IP, aes(x=long, y = lat, group = group), fill = "grey", color = "black")+
  coord_map(xlim = c(-20, 40),ylim = c(27, 45))+
  labs(tag = "A") +
  scale_colour_gradient(high = "#132B43", low = "#56B1F7")+
  geom_point(aes(shape=cymo_posi$ref),size=2.5)+
  theme(plot.tag = element_text(),
        legend.title=element_blank())

posi2<-ggplot(cymo_posi, aes(long, Biomass))+ xlab("Longitude")+ 
  ylab(expression(paste("AG Biomass C stock (kg"," ", m^-2, ")")))+
  labs(tag = "B") +
  geom_point(aes(color=ref))+
  theme(legend.position = "none", 
        plot.tag = element_text())

posi3<-ggplot(cymo_posi, aes(long, Biomass))+ xlab("Longitude")+ 
  geom_point(aes(color=ref))+
  labs(tag = "C") +
  ylim(0, 0.5)+
  theme(axis.title.y=element_blank())


PosiBiomass<-grid.arrange(posi1, posi2, posi3, nrow=2, 
                          layout_matrix = cbind(c(1,2), c(1,3)),
                          top="Posidonia oceanica biomass")


ggsave(path = Folder,"Posidonia biomass.jpg",SgPlot, units="cm", width = 13, height = 20)



## median posi this study others


aggregate(cymo_posi[,c("Biomass", "Stock", "Seq.rate")], by=list(cymo_posi$ref), FUN=median, na.rm=TRUE)

aggregate(cymo_posi[,c("Biomass", "Stock", "Seq.rate")], by=list(cymo_posi$ref), FUN=max, na.rm=TRUE)

aggregate(subset(BC_PI, Genus=="Cymodocea nodosa")[,c("Mean_Biomass")], by=list((subset(BC_PI, Genus=="Cymodocea nodosa"))[,"Coast"]), FUN=median, na.rm=TRUE)




# Tabla con media, mediana, se and n para Iberia y para global values ----------

#get n per variable


iberian_peninsula_seagrass<-rbind(
aggregate(subset(comp_sg, ref=="This Study")[,c("Biomass", "Stock", "Seq.rate")], by=list(subset(comp_sg, ref=="This Study")$Specie), FUN=mean, na.rm=TRUE),
aggregate(subset(comp_sg, ref=="This Study")[,c("Biomass", "Stock", "Seq.rate")], by=list(subset(comp_sg, ref=="This Study")$Specie), FUN=median, na.rm=TRUE),
aggregate(subset(comp_sg, ref=="This Study")[,c("Biomass", "Stock", "Seq.rate")], by=list(subset(comp_sg, ref=="This Study")$Specie), FUN=mad, na.rm=T,low=F,high=F),
aggregate(subset(comp_sg, ref=="This Study")[,c("Biomass", "Stock", "Seq.rate")], by=list(subset(comp_sg, ref=="This Study")$Specie), FUN=contar))
iberian_peninsula_seagrass$fun<-c("mean","mean","mean","mean", "median", "median", "median" , "median","mad","mad","mad","mad", "n", "n", "n" , "n")
iberian_peninsula_seagrass$range<-"IB"


global_seagrass<-rbind(
aggregate(comp_sg[,c("Biomass", "Stock", "Seq.rate")], by=list(comp_sg$Specie), FUN=mean, na.rm=TRUE),
aggregate(comp_sg[,c("Biomass", "Stock", "Seq.rate")], by=list(comp_sg$Specie), FUN=median, na.rm=TRUE),
aggregate(comp_sg[,c("Biomass", "Stock", "Seq.rate")], by=list(comp_sg$Specie), FUN=std.error),
aggregate(comp_sg[,c("Biomass", "Stock", "Seq.rate")], by=list(comp_sg$Specie), FUN=contar))
global_seagrass$fun<-c("mean","mean","mean","mean", "median", "median", "median" , "median","mad","mad","mad","mad", "n", "n", "n" , "n" )
global_seagrass$range<-"G"

iberian_peninsula_sm<-rbind(
  aggregate(subset(comp_sm, ref=="This Study")[,c("Biomass", "Stock", "Seq.rate")], by=list(subset(comp_sm, ref=="This Study")$Tidal.R), FUN=mean, na.rm=TRUE),
  aggregate(subset(comp_sm, ref=="This Study")[,c("Biomass", "Stock", "Seq.rate")], by=list(subset(comp_sm, ref=="This Study")$Tidal.R), FUN=median, na.rm=TRUE),
  aggregate(subset(comp_sm, ref=="This Study")[,c("Biomass", "Stock", "Seq.rate")], by=list(subset(comp_sm, ref=="This Study")$Tidal.R), FUN=mad, na.rm=T,low=F,high=F),
  aggregate(subset(comp_sm, ref=="This Study")[,c("Biomass", "Stock", "Seq.rate")], by=list(subset(comp_sm, ref=="This Study")$Tidal.R), FUN=contar))
iberian_peninsula_sm$fun<-c("mean","mean","mean","mean", "median", "median", "median" , "median","mad","mad","mad","mad", "n", "n", "n" , "n")
iberian_peninsula_sm$range<-"IB"

global_sm<-rbind(
  aggregate(comp_sm[,c("Biomass", "Stock", "Seq.rate")], by=list(comp_sm$Tidal.R), FUN=mean, na.rm=TRUE),
  aggregate(comp_sm[,c("Biomass", "Stock", "Seq.rate")], by=list(comp_sm$Tidal.R), FUN=median, na.rm=TRUE),
  aggregate(comp_sm[,c("Biomass", "Stock", "Seq.rate")], by=list(comp_sm$Tidal.R), FUN=std.error),
  aggregate(comp_sm[,c("Biomass", "Stock", "Seq.rate")], by=list(comp_sm$Tidal.R), FUN=contar))
global_sm$fun<-c("mean","mean","mean","mean", "median", "median", "median" , "median","mad","mad","mad","mad", "n", "n", "n" , "n" )
global_sm$range<-"G"


comp_g<-rbind(comp_sg[,c(1,5:8)], comp_sm[,c(1,5:8)])


iberian_peninsula_eco<-rbind(
  aggregate(subset(comp_g, ref=="This Study")[,c("Biomass", "Stock", "Seq.rate")], by=list(subset(comp_g, ref=="This Study")$Ecosystem), FUN=mean, na.rm=TRUE),
  aggregate(subset(comp_g, ref=="This Study")[,c("Biomass", "Stock", "Seq.rate")], by=list(subset(comp_g, ref=="This Study")$Ecosystem), FUN=median, na.rm=TRUE),
  aggregate(subset(comp_g, ref=="This Study")[,c("Biomass", "Stock", "Seq.rate")], by=list(subset(comp_g, ref=="This Study")$Ecosystem), FUN=mad, na.rm=T,low=F,high=F),
  aggregate(subset(comp_g, ref=="This Study")[,c("Biomass", "Stock", "Seq.rate")], by=list(subset(comp_g, ref=="This Study")$Ecosystem), FUN=contar))
iberian_peninsula_eco$fun<-c("mean","mean", "median", "median","mad","mad" , "n", "n")
iberian_peninsula_eco$range<-"IB"


global_eco<-rbind(
  aggregate(comp_g[,c("Biomass", "Stock", "Seq.rate")], by=list(comp_g$Ecosystem), FUN=mean, na.rm=TRUE),
  aggregate(comp_g[,c("Biomass", "Stock", "Seq.rate")], by=list(comp_g$Ecosystem), FUN=median, na.rm=TRUE),
  aggregate(comp_g[,c("Biomass", "Stock", "Seq.rate")], by=list(comp_g$Ecosystem), FUN=mad, na.rm=T,low=F,high=F),
  aggregate(comp_g[,c("Biomass", "Stock", "Seq.rate")], by=list(comp_g$Ecosystem), FUN=contar))
global_eco$fun<-c("mean","mean","median", "median", "mad","mad", "n", "n" )
global_eco$range<-"G"



table_global<-rbind(iberian_peninsula_seagrass,global_seagrass, iberian_peninsula_sm, global_sm, iberian_peninsula_eco, global_eco )

write.csv(table_global,file.path(Folder,"Global table.csv"))


# Final stocks per ecosystem coast and country -------------------------------------


    # we add global median and SE to summary table in Zostera marina and microtidal marshes fluxes (as we have no data)
    
    Summary2<-Summary
    
    #get median and SE Zostera marina flux from global data 
    Summary2[which(Summary2[,"Category"]=="Zostera marina"), which(colnames(Summary2)=="M.SeqRate")]<-
      table_global[which(table_global$Group.1 == "Zostera marina" & table_global$fun == "median" & table_global$range == "G"), "Seq.rate"]
    
    Summary2[which(Summary2[,"Category"]=="Zostera marina"), which(colnames(Summary2)=="MAD.SeqRate")]<-
      table_global[which(table_global$Group.1 == "Zostera marina" & table_global$fun == "mad" & table_global$range == "G"), "Seq.rate"]
    
    
    #get median and SE microtidal flux from global data 
    Summary2[which(Summary2[,"Category"]=="Microtidal"), which(colnames(Summary2)=="M.SeqRate")]<-
      table_global[which(table_global$Group.1 == "Microtidal" & table_global$fun == "median" & table_global$range == "G"), "Seq.rate"]
    
    Summary2[which(Summary2[,"Category"]=="Microtidal"), which(colnames(Summary2)=="MAD.SeqRate")]<-
      table_global[which(table_global$Group.1 == "Microtidal" & table_global$fun == "mad" & table_global$range == "G"), "Seq.rate"]


#built table with final stocks from Areas and Summary data frames (one table with median and one with MAD)
Median_area_stock<-Area[,-c(3, 13, 14)]
Median_area_stock[,c(3:11)]<-NA
MAD_area_stock<-Median_area_stock

Median_area_seq_rate<-Median_area_stock
MAD_area_seq_rate<-Median_area_stock
Median_area_biomass<-Median_area_stock
MAD_area_biomass<-Median_area_stock

Area2<-Area[,-c(3, 13, 14)]

fill_final_table<- function (df, ncolsum) {

for (i in 3:length(df)) {

  for (j in 1:nrow(df)) {
    
    df[j,i]<-Area2[j,i]*Summary2[which(Summary2[,"Category"]==colnames(df[i])), ncolsum]
  
  }
}
  
  # kg m2 by km2. # multiply by 1000000 to get kg km2
  
  df[,c(3:11)]<-df[,c(3:11)]*1000000
  
  # divide by 1000 to get T of OC and by 1000000 to get Mg T OC
  
  df[,c(3:11)]<-df[,c(3:11)]/(1000*1000000)
  
  return(df)
  }


#final values of T OC in each region and coast by category
Median_area_biomass<- fill_final_table(Median_area_biomass, 6)
MAD_area_biomass<- fill_final_table(MAD_area_biomass, 7)

Median_area_stock<- fill_final_table(Median_area_stock, 12)
MAD_area_stock<- fill_final_table(MAD_area_stock, 13)

Median_area_seq_rate<- fill_final_table(Median_area_seq_rate, 18)
MAD_area_seq_rate<- fill_final_table(MAD_area_seq_rate, 19)


#summary table per coasts 

Median_area_biomass$seagrass<-rowSums(Median_area_biomass[ , c(3:7)], na.rm=TRUE)
Median_area_biomass$"salt marsh"<-rowSums(Median_area_biomass[ , c(8:11)], na.rm=TRUE)
Median_area_stock$seagrass<-rowSums(Median_area_stock[ , c(3:7)], na.rm=TRUE)
Median_area_stock$"salt marsh"<-rowSums(Median_area_stock[ , c(8:11)], na.rm=TRUE)
Median_area_seq_rate$seagrass<-rowSums(Median_area_seq_rate[ , c(3:7)], na.rm=TRUE)
Median_area_seq_rate$"salt marsh"<-rowSums(Median_area_seq_rate[ , c(8:11)], na.rm=TRUE)

agg_biomass_coast <- aggregate(Median_area_biomass[,c(3:13)], by=list(Median_area_biomass$Coast), FUN=sum, na.rm=TRUE)
agg_stock_coast <- aggregate(Median_area_stock[,c(3:13)], by=list(Median_area_stock$Coast), FUN=sum, na.rm=TRUE)
agg_seq_rate_coast <- aggregate(Median_area_seq_rate[,c(3:13)], by=list(Median_area_seq_rate$Coast), FUN=sum, na.rm=TRUE)


Median_area_biomass$country<-c("Spain","Spain","Spain","Spain","Spain","Spain","Spain","Spain","Spain","Portugal","Portugal","Portugal","Spain","Spain","Spain","Portugal","Portugal")
Median_area_stock$country<-c("Spain","Spain","Spain","Spain","Spain","Spain","Spain","Spain","Spain","Portugal","Portugal","Portugal","Spain","Spain","Spain","Portugal","Portugal")
Median_area_seq_rate$country<-c("Spain","Spain","Spain","Spain","Spain","Spain","Spain","Spain","Spain","Portugal","Portugal","Portugal","Spain","Spain","Spain","Portugal","Portugal")


agg_biomass_country <- aggregate(Median_area_biomass[,c(3:13)], by=list(Median_area_biomass$country), FUN=sum, na.rm=TRUE)
agg_stock_country <- aggregate(Median_area_stock[,c(3:13)], by=list(Median_area_stock$country), FUN=sum, na.rm=TRUE)
agg_seq_rate_country <- aggregate(Median_area_seq_rate[,c(3:13)], by=list(Median_area_seq_rate$country), FUN=sum, na.rm=TRUE)


final_area_biomass<-rbind(agg_biomass_coast, agg_biomass_country)
final_area_stock<-rbind(agg_stock_coast, agg_stock_country)
final_area_seq_rate<-rbind(agg_seq_rate_coast, agg_seq_rate_country)


final_table4<-rbind(final_area_biomass, final_area_stock, final_area_seq_rate)
write.csv(final_table4,file.path(Folder,"Table4_region.csv"))

# % of stock and fluxes of Posidonia over the total

tot_stock<- rbind(agg_biomass_country, agg_stock_country)
colSums(tot_stock[,c(2,11,12)])

colSums(tot_stock[,c(2,11,12)])[1]*100/(colSums(tot_stock[,c(2,11,12)])[2]+colSums(tot_stock[,c(2,11,12)])[3])



# biomass vs soil stock
#% of seagrass biomass of total stock (biomass + top meter soil)
(sum(agg_biomass_country$seagrass)*100)/(sum(agg_stock_country$seagrass)+sum(agg_biomass_country$seagrass))
#[1] 0.8632568

#% of salt marsh biomass of total stock (biomass + top meter soil)
(sum(agg_biomass_country$`salt marsh`)*100)/(sum(agg_stock_country$`salt marsh`)+sum(agg_biomass_country$`salt marsh`))
#[1] 3.72984

#summary figures 

final<-as.data.frame(melt(final_area_biomass[,c("Group.1", "seagrass", "salt marsh")], ID="Group.1"))
final$stock<-as.data.frame(melt(final_area_stock[,c("Group.1", "seagrass", "salt marsh")], ID="Group.1"))[,3]
final$flux<-as.data.frame(melt(final_area_seq_rate[,c("Group.1", "seagrass", "salt marsh")], ID="Group.1"))[,3]
colnames(final)<-c("Region", "Ecosystem", "Biomass", "Stock", "Flux")


m_final_r<-melt(subset(final, !Region=="Portugal" & !Region=="Spain"), ID=c("Region", "Ecosystem"))
m_final_c<-melt(subset(final, Region=="Portugal" | Region=="Spain"), ID=c("Region", "Ecosystem"))

# stock figures
p_rsg<-ggplot(subset(m_final_r, Ecosystem=="seagrass" & !variable=="Flux"))+ ylab("Tg OC")+
  geom_bar(aes(x=Region, y=value, fill=variable), position = "stack", stat = "identity", width=0.5)+
  scale_fill_manual(values=c("green4","burlywood3"))+
  ylim(0,21)+
  scale_y_break(c(4, 18))+
  theme(legend.position = "none",
        axis.text.y.right = element_blank(),
        axis.line.y.right = element_blank(),
        axis.ticks.y.right = element_blank()
  )


p_rsm<-ggplot(subset(m_final_r, Ecosystem=="salt marsh" & !variable=="Flux"))+ ylab("Tg OC")+
  geom_bar(aes(x=Region, y=value, fill=variable), width=0.5, position = "stack", stat = "identity")+
  scale_fill_manual(values=c("green4","burlywood3"))+
  ylim(0,4)+
  theme(legend.position = "none")
  
p_csg<-  ggplot(subset(m_final_c, Ecosystem=="seagrass" & !variable=="Flux"))+ ylab("Tg OC")+
  geom_bar(aes(x=Region, y=value, fill=variable), position = "stack", stat = "identity")+
  scale_fill_manual(values=c("green4","burlywood3"),labels=c('Biomass', 'Soil'))+
  ylim(0,21)+
  scale_y_break(c(3, 18))+
  theme(axis.text.y.right = element_blank(),
        axis.line.y.right = element_blank(),
        axis.ticks.y.right = element_blank()
  )


p_csm<-ggplot(subset(m_final_c, Ecosystem=="salt marsh" & !variable=="Flux"))+ ylab("Tg OC")+
  geom_bar(aes(x=Region, y=value, fill=variable), position = "stack", stat = "identity")+
  scale_fill_manual(values=c("green4","burlywood3"), labels=c('Biomass', 'Soil'))+
  ylim(0,4)

#flux figure

f_rsg <-   ggplot(subset(m_final_r, variable=="Flux"), aes(x=Region, y=value, fill=Ecosystem))+ ylab(expression(paste("Tg OC ",y^-1,"")))+
  geom_bar(stat = "identity",  position = "dodge")+
  scale_fill_manual(values=c("seagreen", "lightsalmon"), labels=c('Seagrass', 'Salt Marsh'))+
  ylim(0,0.035)+
  theme(legend.position = "none")
f_csg <-ggplot(subset(m_final_c, variable=="Flux"), aes(x=Region, y=value, fill=Ecosystem))+ 
  ylab(expression(paste("Tg OC ",y^-1,"")))+ xlab("Country")+
  geom_bar(stat = "identity",  position = "dodge")+
  scale_fill_manual(values=c("seagreen", "lightsalmon"), labels=c('Seagrass', 'Salt Marsh'))+
  ylim(0,0.035)




#final figure

stocks_by_region<-grid.arrange( print(p_rsg), print(p_csg), print(p_rsm), print(p_csm),  widths=c(0.6, 0.4))

flux_by_region<-grid.arrange(f_rsg, f_csg, ncol=2, widths=c(0.6, 0.4))

ggsave(path = Folder,"stock_Region.jpg",stocks_by_region, units="cm", width = 22, height = 15)
ggsave(path = Folder,"flux_Region.jpg",flux_by_region, units="cm", width = 22, height = 7)








# 50 cm stocks ------------------------------------------------------------

    BCS_50<-estimate_oc_stock(A, depth= 50,core="Core.ID", mind="DMin.D", maxd="DMax.D", dbd="DDBD", oc="POC")
    
    # Stock 50cm. Mean and sd by site
    BCS_50$Site.ID<-NA
    
    for (i in 1:nrow(BCS_50)) {
      
      Site<- unique(A[c(which(A$Core.ID==BCS_50[i,which(colnames(BCS_50)=="core")])),which(colnames(A)=="Site.ID")])
      BCS_50[i,which(colnames(BCS_50)=="Site.ID")]<- Site
      
    }
    
    # Estimate stocks per site
    
    S_by_Site_50 <-merge(aggregate( stock ~ Site.ID, BCS_50, mean), aggregate( stock ~ Site.ID, BCS_50, sd), by = "Site.ID")
    colnames(S_by_Site_50)<-c("Site.ID", "Mean_Stock", "SD_Stock")
    
    
    #get info by site
    BC_PI_50 <-merge(B, S_by_Site_50, by = "Site.ID", all = TRUE)
    BC_PI_50<-subset(BC_PI_50, !is.na(BC_PI_50$Mean_Stock))
    BC_PI_50 [,c(13:14)]<- BC_PI_50 [,c(13:14)]*10 #from g cm2 to kg m2
    
    # Summary table
    
    Summary_50<- data.frame(Category=character(),
                         Area=numeric(),
                         nS.Stock=numeric(),
                         SA.Stock=numeric(),
                         Av.Stock=numeric(),
                         M.Stock=numeric(),
                         MAD.Stock=numeric(),
                         cv.Stock=numeric()
                         )
    
    
    
    
    X<-split(BC_PI_50, BC_PI_50$Ecosystem)
    Sg_50<-subset(BC_PI_50, Ecosystem=='Seagrass')
    X2<-split(Sg_50, Sg_50$Genus)
    Sm_50<-subset(BC_PI_50, Ecosystem=='Salt Marsh')
    X3<-split(Sm_50, Sm_50$Tidal.R)
    
    X<-c(X,X2,X3)
    
    
    # mad
    
    for(i in 1:length(X)) {
      
      
      Data<-as.data.frame(X[i])
      colnames(Data)<-colnames(BC_PI_50)
      
      Summary_50[i,1]<-names(X[i])
      Summary_50[i,2]<-sum(Area[,names(X[i])], na.rm=TRUE)
      
      Summary_50[i,3]<-nrow(Data)
      Summary_50[i,4]<-Summary_50[i,3]/Summary[i,2]
      Summary_50[i,5]<-mean(Data$Mean_Stock)
      Summary_50[i,6]<-median(Data$Mean_Stock)
      Summary_50[i,7]<-mad(Data$Mean_Stock,center=median(Data$Mean_Stock),na.rm=T,low=F,high=F)
      Summary_50[i,8]<-coef.var(Data$Mean_Stock)
    }
    


# study limitations -------------------------------------------------------


p1<-ggplot(Summary[-c(1,2),], aes(Area, nS.Biomass))+ ylab("Number of biomass stations")+
  xlab(expression(paste("Area (k", m^-2, ")")))+
  geom_point()+
  ggrepel::geom_label_repel(aes(label = Category),
                   box.padding   = 0.35, 
                   point.padding = 0.5,
                   segment.color = 'grey50')

p2<-ggplot(Summary[-c(1,2),], aes(Area, nS.Stock))+ ylab("Number of soil stock stations")+
  xlab(expression(paste("Area (k", m^-2, ")")))+
  geom_point()+
  ggrepel::geom_label_repel(aes(label = Category),
                            box.padding   = 0.35, 
                            point.padding = 0.5,
                            segment.color = 'grey50')


p3<-ggplot(Summary[-c(1,2),], aes(Area, nS.SeqRate))+ylab("Number of seq. rate stations")+
  xlab(expression(paste("Area (k", m^-2, ")")))+
  geom_point()+
  ggrepel::geom_label_repel(aes(label = Category),
                            box.padding   = 0.35, 
                            point.padding = 0.5,
                            segment.color = 'grey50')


supp_data<-grid.arrange(p1, p2, p3,nrow=1)

ggsave(path = Folder,"estaciones por area.jpg",supp_data, units="cm", width = 30, height = 9)


