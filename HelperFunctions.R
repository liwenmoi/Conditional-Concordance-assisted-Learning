LoglikeExactRC <- function(coef1less, sub, sp, idIndex){
  ### Caution: assume dataset sorted by groupID, Y (descending)
  
  if (abs(coef1less)>1) {
    negloglike = 1e6
  } else {
    
    coef <- c(sqrt(1 - sum(coef1less^2)), coef1less)     # bc norm = 1 constraint
    
    sub$score <- as.matrix(sub[,grep("X[0-9]",colnames(sub))])%*%coef
    
    c <- find.c(score=sub$score[sub$Y==0], weight=sub$weight[sub$Y==0], sp=sp)
    
    negloglike <- loglike_exact_C(names(idIndex), idIndex=idIndex, score=sub$score, c=c)
    
  }
  
  return(negloglike)
  
}


ProposedExact <- function(ini, data, sp){
  ### given a grid of stating values, data (ordered, with id, Y, X1, X2,.., weight, groupID), sp, 
  ### return normed coefficients,
  data <- data[order(data$groupID, -data$Y),]
  
  groupIDs <- data$groupID
  uniqueIDs <- unique(groupIDs)
  idx <- list()
  for (i in 1:length(uniqueIDs)){
    idx[[paste(uniqueIDs[i])]] <- which(groupIDs==uniqueIDs[i])
  }
  
  ### start calculation
  reses <- c()
  for (start in ini){
    # print(start)
    res = optim(par=start, LoglikeExactRC, sub=data, sp=sp, idIndex=idx)
    reses <- rbind(reses, c(res$convergence, res$value, res$par))
  }
  reses <- as.data.frame(reses); 
  colnames(reses) <- c("conv", "value", "par")
  reses <- subset(reses, conv==0)
  coef <- reses$par[which.min(reses$value)]
  coef <- inv.ucli(coef)
  
  return(list(coef=coef))
  
} 


find.c <- function(score, weight, sp){ 
  ### finding cutoff only needs information in the controls
  ### given the controls' score, disease status, weight and pre-specified specificity
  od <- order(score)                                      
  score.o <- score[od]
  w.o <- weight[od]
  wsum <- cumsum(w.o)
  scaled.w <- wsum/wsum[length(wsum)]
  idx <- which(scaled.w>=sp)[1]
  c = score.o[idx]
}   


specificity.matched <- function(score, weight, tau){
  #### given controls' score and weight
  od <- order(score)                                      
  score.o <- score[od]
  w.o <- weight[od]
  wsum <- cumsum(w.o)
  scaled.w <- wsum/wsum[length(wsum)]
  
  idx <- which(score.o <= tau)
  idx <- idx[length(idx)]
  
  sp = scaled.w[idx]
} 


get.train.res <- function(markers,
                          Y,
                          weight,
                          coef,
                          sp, 
                          get.full=TRUE,
                          markers.full=NULL,
                          Y.full=NULL){
  ### given markers in matrix form, disease status, weight, coefficients, and sp
  ### return c, sens
  
  cmb <- as.matrix(markers) %*% coef
  
  
  ### given the risk index and the true result, return auc ignoring the weight
  roc_emp <- roc(Y, as.numeric(cmb), smooth=FALSE, auc=TRUE,
                 direction="<", levels = c(0,1))
  
  
  c <- find.c(score=cmb[Y==0], weight=weight[Y==0], sp=sp)
  
  ### check sensitivity
  sens <- sum(cmb[Y==1]>c)/sum(Y==1)
  
  
  sp.full = NULL
  if (get.full==TRUE){
    ### check specificity on the full data; no weight
    cmb.full <- as.matrix(markers.full) %*% coef
    sp.full <- sum(cmb.full[Y.full==0]<=c)/sum(Y.full==0)
  }
  
  return(list(cmb=cmb, c=c, sens=sens, sp.full=sp.full, auc=roc_emp$auc))
  
}


get.val.res <- function(markers,
                        Y,
                        weight,
                        coef,
                        c,
                        get.full=TRUE,
                        markers.full=NULL,
                        Y.full=NULL){ 
  ### given coef, c, markers
  ### return sens, sp

  cmb <- as.matrix(markers) %*% coef
  
  ### given the risk index and the true result, return auc ignoring the weight
  roc_emp <- roc(Y, as.numeric(cmb), smooth=FALSE, auc=TRUE,
                 direction="<", levels = c(0,1))
  
  
  sens <- sum(cmb[Y==1]>c)/sum(Y==1)
  sp <- specificity.matched(score=cmb[Y==0], weight=weight[Y==0], tau=c)
  
  
  sp.full = NULL
  if (get.full==TRUE){
    cmb.full <- as.matrix(markers.full) %*% coef
    sp.full <- sum(cmb.full[Y.full==0]<=c)/sum(Y.full==0)
  }
  
  return(list(cmb=cmb, sens=sens, sp=sp, sp.full=sp.full, auc=roc_emp$auc))
  
}



ucli.norm <- function(coef){coef/sqrt(sum(coef^2))}
inv.ucli <- function(coef1less){c(sqrt(1-sum(coef1less^2)), coef1less)}
