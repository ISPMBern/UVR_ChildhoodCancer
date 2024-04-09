
#==========================================================================================#
#
# Author: Christian Kreis
#
# Created: 2 November 2021
#
# R version 3.6.0 (2019-04-26) -- "Planting of a Tree"
#
# Survival analysis of UV exposure and childhood hematological cancer risk
#
#==========================================================================================#

rm(list=ls(all=TRUE))

library(Hmisc)
library(survival)

Sys.setlocale("LC_TIME","english")

options(scipen=10000)
options(stringsAsFactors=FALSE)

# working directory
.wd <- "Y:/ENVEPI/RESTRICTED/temp/p_UV/R"

# file paths
.tab <- "textres/DisplayElements"

#==========================================================================================#
# Set global parameters and functions
#==========================================================================================#

.t0 <- Sys.time()

period <- "Year"
#period <- "July"

if (period %in% "Year") {
   exp.cat <- "meanUVI_q"
   exp.cat.mean <- "meanUVI_cat_mean"
   exp.lin <- "meanUVI"
   exp.int <- "meanUVIc"
}

if (period %in% "July") {
   exp.cat <- "meanUVI.Jul_q"
   exp.cat.mean <- "meanUVI.Jul_cat_mean"
   exp.lin <- "meanUVI.Jul"
}

# Minimally and fully adjusted conditional logistic regression for categorical UVR and cancer risk

model.cat <- function(data) {
   fmla <- as.formula(paste0("Case ~ ",exp.cat," + sex + strata(Group)"))
   univar <- clogit(fmla, data=data)
   fmla <- as.formula(paste0("Case ~ ",exp.cat," + totalIR + NO2Dx + urban.cat + ssep_q + cantonRegist + sex + strata(Group)"))
   final <-  clogit(fmla, data=data)
   return(list(univar,final))
}

# ANOVA comparing null model with model including categorical and categorical mean term

model.anova <- function(data) {
   model.0          <- clogit(Case ~               totalIR + NO2Dx + urban.cat + ssep_q + cantonRegist + sex + strata(Group), data=data)
   fmla <- as.formula(paste0("Case ~ ",exp.cat," + totalIR + NO2Dx + urban.cat + ssep_q + cantonRegist + sex + strata(Group)"))
   model.cat      <-  clogit(fmla, data=data)
   fmla <- as.formula(paste0("Case ~ ",exp.cat.mean," + totalIR + NO2Dx + urban.cat + ssep_q + cantonRegist + sex + strata(Group)"))
   model.catmean  <-  clogit(fmla, data=data)
   return(anova(model.0,model.cat,model.catmean))
}

# Crude and adjusted conditional logistic regression for linear UVR and cancer risk

model.lin <- function(data) {
   fmla <- as.formula(paste0("Case ~ ",exp.lin," + sex + strata(Group)"))
   univar <- clogit(fmla, data=data)
   fmla <- as.formula(paste0("Case ~ ",exp.lin," + totalIR + NO2Dx + urban.cat + ssep_q + cantonRegist + sex + strata(Group)"))
   final <-  clogit(fmla, data=data)
   return(list(univar,final))
}

# Comparing full model with model including region and interaction between region and UVI

model.int <- function(data) {
   data <- data[!is.na(data$region),]
   fmla <- as.formula(paste0("Case ~ ",exp.int,"          + totalIR + NO2Dx + urban.cat + ssep_q + cantonRegist + sex + strata(Group)"))
   model.lin      <-  clogit(fmla, data=data)
   print(summary(model.lin))
   fmla <- as.formula(paste0("Case ~ ",exp.int," + region + totalIR + NO2Dx + urban.cat + ssep_q + cantonRegist + sex + strata(Group)"))
   model.lin.reg  <-  clogit(fmla, data=data)
   print(summary(model.lin.reg))
   fmla <- as.formula(paste0("Case ~ ",exp.int," * region + totalIR + NO2Dx + urban.cat + ssep_q + cantonRegist + sex + strata(Group)"))
   model.lin.int  <-  clogit(fmla, data=data)
   print(summary(model.lin.int))
   return(anova(model.lin,model.lin.reg,model.lin.int))
}

#==========================================================================================#
# Import
#==========================================================================================#

# SNC Childhood data #
#--------------------#

# import childhood SNC time to event data set
SNCrset <- readRDS(file=file.path(.wd,"data/UVR_CHM_SNCriskset.rds"))

#==========================================================================================#
# Preparation
#==========================================================================================#

SNCrset$totalIR_q <- cut(SNCrset$totalIR, breaks=quantile(SNCrset$totalIR[SNCrset$HM %in% 0],
                         probs=seq(0,1,0.25), na.rm=TRUE), include.lowest=TRUE)

## Select observation periods in time-to-event data set when cases were diagnosed ##

SNCrset <- SNCrset[SNCrset$ageDx >= SNCrset$age.t1 & SNCrset$ageDx <= SNCrset$age.t2,]

# Check if no single individual is included twice in any given risk set
any(duplicated(SNCrset[SNCrset$ageDx >= SNCrset$age.t1 & SNCrset$ageDx <= SNCrset$age.t2,c("sncidNUM","Group")]))
#Must be FALSE

#SNCrset$region <- relevel(SNCrset$region,ref="Central Switzerland")

SNCrset$meanUVIc <- SNCrset$meanUVI - mean(SNCrset$meanUVI)

#==========================================================================================#
# Analysis
#==========================================================================================#

# Survival analysis #
#-------------------#

## all hematological cancers ##

hema <- SNCrset$Group %in% SNCrset$Group[SNCrset$Case %in% 1 & SNCrset$HM %in% 1]
hema.cat   <- model.cat(SNCrset[hema,])
hema.anova <- model.anova(SNCrset[hema,])
hema.lin   <- model.lin(SNCrset[hema,])

model.int(SNCrset[hema,])

## leukemia ##

leuk <- SNCrset$Group %in% SNCrset$Group[SNCrset$Case %in% 1 & SNCrset$Leukemia %in% 1]
leuk.cat   <- model.cat(SNCrset[leuk,])
leuk.anova <- model.anova(SNCrset[leuk,])
leuk.lin   <- model.lin(SNCrset[leuk,])

model.int(SNCrset[leuk,])

## ALL ##

ALL <- SNCrset$Group %in% SNCrset$Group[SNCrset$Case %in% 1 & SNCrset$ALL %in% 1]
ALL.cat   <- model.cat(SNCrset[ALL,])
ALL.anova <- model.anova(SNCrset[ALL,])
ALL.lin   <- model.lin(SNCrset[ALL,])

model.int(SNCrset[ALL,])

## AML ##

AML <- SNCrset$Group %in% SNCrset$Group[SNCrset$Case %in% 1 & SNCrset$AML %in% 1]
AML.cat   <- model.cat(SNCrset[AML,])
AML.anova <- model.anova(SNCrset[AML,])
AML.lin   <- model.lin(SNCrset[AML,])

model.int(SNCrset[AML,])

## lymphoma ##

lymph <- SNCrset$Group %in% SNCrset$Group[SNCrset$Case %in% 1 & SNCrset$Lymphoma %in% 1]
lymph.cat   <- model.cat(SNCrset[lymph,])
lymph.anova <- model.anova(SNCrset[lymph,])
lymph.lin   <- model.lin(SNCrset[lymph,])

model.int(SNCrset[lymph,])

## Hodgkin lymphoma ##

HL <- SNCrset$Group %in% SNCrset$Group[SNCrset$Case %in% 1 & SNCrset$HL %in% 1]
HL.cat   <- model.cat(SNCrset[HL,])
HL.anova <- model.anova(SNCrset[HL,])
HL.lin   <- model.lin(SNCrset[HL,])

model.int(SNCrset[HL,])

## non-Hodgkin lymphoma ##

NHL <- SNCrset$Group %in% SNCrset$Group[SNCrset$Case %in% 1 & SNCrset$NHL %in% 1]
NHL.cat   <- model.cat(SNCrset[NHL,])
NHL.anova <- model.anova(SNCrset[NHL,])
NHL.lin   <- model.lin(SNCrset[NHL,])

model.int(SNCrset[NHL,])

#==========================================================================================#
# Display elements
#==========================================================================================#

### Tables ###
#------------#

# Characteristics table #
#-----------------------#

characteristics <- cbind(

   with(SNCrset[SNCrset$Case %in% 0,],
      rbind(
         table(              get(exp.cat),useNA="always"),
         table(sex,          get(exp.cat),useNA="always"),
         table(yob_c,        get(exp.cat),useNA="always"),
         table(yearEntry.cat,get(exp.cat),useNA="always"),
         table(urban.cat,    get(exp.cat),useNA="always"),
         table(ssep_q,       get(exp.cat),useNA="always"),
         table(cantonRegist, get(exp.cat),useNA="always"),
         table(NO2_q,        get(exp.cat),useNA="always"),
         table(totalIR_q,    get(exp.cat),useNA="always")
      )
   )
,
   with(SNCrset[SNCrset$Case %in% 1,],
      rbind(
         table(              get(exp.cat),useNA="always"),
         table(sex,          get(exp.cat),useNA="always"),
         table(yob_c,        get(exp.cat),useNA="always"),
         table(yearEntry.cat,get(exp.cat),useNA="always"),
         table(urban.cat,    get(exp.cat),useNA="always"),
         table(ssep_q,       get(exp.cat),useNA="always"),
         table(cantonRegist, get(exp.cat),useNA="always"),
         table(NO2_q,        get(exp.cat),useNA="always"),
         table(totalIR_q,    get(exp.cat),useNA="always")
      )
   )
)

write.table(characteristics, file.path(.wd,.tab,paste0("Characteristics.",period,".txt")),
   quote=FALSE, sep="\t", row.names=TRUE)

#------------------------------------------------------------------------------------------#

# Results table: Categorical exposure #
#-------------------------------------#

# Function to extract the output of the regression models to build tables

extract.cat <- function(dg,model) {data.frame(
   outcome   = capitalize(unlist(strsplit(deparse(substitute(model)),split="[.]"))[1]),
   category  = factor(model[[1]]$xlevels[[1]], levels=model[[1]]$xlevels[[1]]),
   cases     = table(SNCrset[dg,"Case"][complete.cases(SNCrset[dg,names(model[[2]]$assign)])],
                  SNCrset[dg,names(model[[2]]$assign)[1]][complete.cases(SNCrset[dg,names(model[[2]]$assign)])])[2,],
   coef      = c(1,exp(model[[1]]$coefficients)[1:4]),
   lower     = c(1,exp(confint(model[[1]], level=0.95))[1:4,1]),
   upper     = c(1,exp(confint(model[[1]], level=0.95))[1:4,2]),
   adj.coef  = c(1,exp(model[[2]]$coefficients)[1:4]),
   adj.lower = c(1,exp(confint(model[[2]], level=0.95))[1:4,1]),
   adj.upper = c(1,exp(confint(model[[2]], level=0.95))[1:4,2])
)}

association.cat <- rbind(
   extract.cat(hema,hema.cat),
   extract.cat(leuk,leuk.cat),
   extract.cat(ALL,ALL.cat),
   extract.cat(AML,AML.cat),
   extract.cat(lymph,lymph.cat),
   extract.cat(HL,HL.cat),
   extract.cat(NHL,NHL.cat)
)

write.table(association.cat, file.path(.wd,.tab,paste0("Association.UVR.CHM.",period,".categorical.txt")),
   quote=FALSE, sep="\t", row.names=FALSE, na="")

# Results table: Linear exposure #
#--------------------------------#

extract.lin <- function(dg,model) {data.frame(
   outcome   = capitalize(unlist(strsplit(deparse(substitute(model)),split="[.]"))[1]),
   cases     = table(SNCrset[dg,"Case"][complete.cases(SNCrset[dg,names(model[[2]]$assign)])])[2],
   coef      = exp(model[[1]]$coefficients)[1],
   lower     = exp(confint(model[[1]], level=0.95))[1,1],
   upper     = exp(confint(model[[1]], level=0.95))[1,2],
   adj.coef  = exp(model[[2]]$coefficients)[1],
   adj.lower = exp(confint(model[[2]], level=0.95))[1,1],
   adj.upper = exp(confint(model[[2]], level=0.95))[1,2]
)}

association.lin <- rbind(
   extract.lin(hema,hema.lin),
   extract.lin(leuk,leuk.lin),
   extract.lin(ALL,ALL.lin),
   extract.lin(AML,AML.lin),
   extract.lin(lymph,lymph.lin),
   extract.lin(HL,HL.lin),
   extract.lin(NHL,NHL.lin)
)

write.table(association.lin, file.path(.wd,.tab,paste0("Association.UVR.CHM.",period,".linear.txt")),
   quote=FALSE, sep="\t", row.names=FALSE, na="")

# Results table: Model comparison #
#---------------------------------#

extract.anova <- function(dg,model) {data.frame(
   outcome = rep(capitalize(unlist(strsplit(deparse(substitute(model)),split="[.]"))[1]),3),
   loglik  = model$loglik,
   Chisq   = model$Chisq,
   Df      = model$Df,
   p.value = model$"P(>|Chi|)"
)}

compmodels <- rbind(
   extract.anova(hema,hema.anova),
   extract.anova(leuk,leuk.anova),
   extract.anova(ALL,ALL.anova),
   extract.anova(AML,AML.anova),
   extract.anova(lymph,lymph.anova),
   extract.anova(HL,HL.anova),
   extract.anova(NHL,NHL.anova)
)

write.table(compmodels, file.path(.wd,.tab,paste0("Compare.models.UVR.CHM.",period,".txt")),
   quote=FALSE, sep="\t", row.names=FALSE, na="")

#==========================================================================================#

### Figures ###
#-------------#

# Histogram of ambient UV radiation #
#-----------------------------------#

df <- data.frame(
   uvi=c(SNCrset$meanUVI,SNCrset$meanUVI.Jul),
   case=rep(SNCrset$HM,2), 
   period=rep(c("Year","July"),each=nrow(SNCrset)))

tiff(filename=file.path(.wd,paste0("graphres/DisplayElements/Histogram_UVIndex_11_15.tif")), 
  compression="lzw", res=600, pointsize=12, width=14, height=10, units="cm")
#windows(width=14/2.54,height=10/2.54)

ggplot(df[df$case %in% 0,], aes(x=uvi))+
   geom_histogram(bins = 61, color="black", fill="white")+
   facet_grid(period ~ .) +
   xlab("UV Index") +
   ylab("Frequency") +
   scale_x_continuous(breaks=2:7)

dev.off()

Sys.time() - .t0

#==========================================================================================#
#==========================================================================================#
#==========================================================================================#
