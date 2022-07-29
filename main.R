########################################################################################################
# This code was created by Wen Li in 2019 and was revised to adapt to the current format in 2022
# This code provides the implementation of the Conditional Concordance assisted Learning (CCAF)
#   proposed in the paper "Conditional Concordance-assisted Learning for Combining Biomarkers for 
#   Cancer Population Screening" by Li W, Li R, Yan Q, Feng Z, Ning J. It utilizes the helper functions
#   in ....R and .CPP, and applies the proposed method to the example data "sub.csv", and validated the 
#   performance on "val_sub.csv.
########################################################################################################
rm(list=ls())


library(pROC)
library(MASS)
library(Rcpp)

# go to the directory which saves the codes and example data
setwd("")
sourceCpp("LogLikeExactC.cpp")
source("HelperFunctions.R")


#----------------------------     read in data ----------------------------------------#
# Note that in the simulation, we first generated a data in which each subject 
#  was assigned a matching group membership (or matching strata). We named it the "full"
#  set in the code. Then a matched case-control data was created by selecting all cases 
#  and selecting control with 1:1 matching. We named it the "sub" data in the code. Note 
#  that the biomarkers were only measured on subjects in the matched case-control data.

full <- read.csv("full.csv")
sub <- read.csv("sub.csv")

#--------------------------- Set prespecified specificity -----------------------------#
# say the investigators decide to control specificity the level of 0.70
sp <- 0.70 

#--------------------------------- Prepare data ---------------------------------------#  

### append weight based on outcome Y, matching strata Z in the full data set.
tab <- table(full$Y, full$Z)
weight <- data.frame(Z=as.numeric(colnames(tab)), weight=as.numeric(tab[1,]/tab[2,]))   
weight[weight==Inf] = NA

sub <- merge(sub, weight, by="Z")
sub$weight[sub$Y==1] <- 1


### order data sets before fitting model
sub <- sub[order(sub$groupID, -sub$Y),] 

#---------------------------------- run proposed CCAF  --------------------------------#
### prepare starting value
fit.cl <- try(clogit(Y ~ X1 + X2 + strata(groupID), data=sub), silent=FALSE)
coef.cl <- ucli.norm(fit.cl$coefficients)

### fit CCAF
coef.pro <- ProposedExact(ini=seq(-0.5, 0.5, by=0.05) + coef.cl[-1], data=sub, sp=sp)$coef
  
pro.res <- get.train.res(markers=sub[,grep("X[0-9]",colnames(sub)),drop=F],
                           Y=sub$Y,
                           weight=sub$weight,
                           coef=coef.pro,
                           sp=sp,
                           get.full=FALSE,
                           markers.full=NULL,
                           Y.full = NULL
)
  
#------------------ Use coef and cutoff estimated from CCAF to validate -------------#
# read in vaidation set and append weight based on Y, Z.
val.full <- read.csv("val_full.csv")
val.sub <- read.csv("val_sub.csv")
tab <- table(val.full$Y, val.full$Z)
weight <- data.frame(Z=as.numeric(colnames(tab)), weight=as.numeric(tab[1,]/tab[2,])) 
val.sub <- merge(val.sub, weight, by="Z")
val.sub$weight[val.sub$Y==1] <- 1
  
pro.val <- get.val.res(markers=val.sub[,grep("X[0-9]",colnames(val.sub)),drop=F],
                         Y=val.sub$Y,
                         weight=val.sub$weight,
                         coef=coef.pro,
                         c=pro.res$c,
                         get.full=FALSE,
                         markers.full=NULL,
                         Y.full=NULL)



# sensitivity from validation
pro.val$sens
# [1] 0.759

# specificity from validation
pro.val$sp
# [1] 0.5998399

# sensitivity from training
pro.res$sens
# [1] 0.74
