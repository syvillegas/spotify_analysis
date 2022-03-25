#########################################################################################
##
## R ROUTINES FOR FINAL PROJECT- STATISTICAL MODELLING AND INFERENCE.
##
## AUTHORS: LORENZO LESKOVAR AND SERGIO VILLEGAS,
## MODIFYING SOME OF THE ROUTINES FROM DAVID ROSSELL AND PAUL ROGNON, UNIVERSITAT POMPEU FABRA
##
#########################################################################################

# INDEX
# 1. SETTING PENALIZATION PARAMETER VIA BIC
# 2. CROSS-VALIDATION

#########################################################################################
## 1. SETTING PENALIZATION PARAMETER VIA BIC
#########################################################################################


lasso.bic <- function(y,x) {
  #Select model in LASSO path with best BIC (using LASSO regression estimates)
  #Input
  # - y: vector with response variable
  # - x: design matrix
  # - extended: whether to use EBIC (Chen and Chen 2008) instead of BIC
  #
  #Output: list with the following elements
  # - coef: LASSO-estimated regression coefficient with lambda set via BIC
  # - ypred: predicted y
  # - lambda.opt: optimal value of lambda
  # - lambda: data.frame with bic and number of selected variables for each value of lambda
  require(glmnet)
  fit <- glmnet(x=x,y=y,family='gaussian',alpha=1)
  
  pred <- cbind(1,x) %*% rbind(fit$a0,fit$beta)
  n <- length(y)
  p <- colSums(fit$beta!=0) + 1
  
  bic <- n * log(colSums((y-pred)^2)/length(y)) + n*(log(2*pi)+1) + log(n)*p 

  sel <- which.min(bic)
  beta <- c(fit$a0[sel],fit$beta[,sel]); names(beta)[1]= 'Intercept'
  ypred <- pred[,sel]
  ans <- list(coef=beta,ypred=ypred,lambda.opt=fit$lambda[sel],lambda=data.frame(lambda=fit$lambda,bic=bic,nvars=p))
  return(ans)
}



#########################################################################################
## 2. CROSS-VALIDATION
#########################################################################################

kfoldCV.mle <- function(y,x,K=10,seed) {
## Perform K-fold cross-validation for least-squares regression estimate
## Input
## - y: response
## - x: data.frame with predictors, intercept should not be present
## - K: number of folds in K-fold cross-validation
## - seed: random number generator seed (optional)
## Output
## - pred: cross-validated predictions for y
## - ssr: residual sum of squares, sum((y-pred)^2)
  if (!missing(seed)) set.seed(seed)
  subset <- rep(1:K,ceiling(nrow(x)/K))[1:nrow(x)]
  subset <- sample(subset,size=nrow(x),replace=FALSE)
  pred <- double(nrow(x))
  num_betas <- c()
  if (ncol(x)>0) {
    for (k in 1:K) {
      sel <- subset==k
      fit <- lm(y[!sel] ~ ., data=x[!sel,,drop=FALSE])
      pred[sel] <- predict(fit, newdata=x[sel,,drop=FALSE])
      b= as.vector(coef(fit))
      num_betas <- c(num_betas, length(b[b!=0]))
    }
  } else {
    for (k in 1:K) {
      sel <- subset==k
      pred[sel] <- mean(y[!sel],na.rm=TRUE)
      num_betas <- c(num_betas, dim(x)[2])
    }
  }
  return(list(pred=pred,ssr=sum((pred-y)^2,na.rm=TRUE),num_betas=num_betas))
}

kfoldCV.lasso <- function(y,x,K=10,seed,criterion='cv') {
## Perform K-fold cross-validation for LASSO regression estimate (lambda set either via cross-val or BIC or EBIC)
## Input
## - y: response
## - x: data.frame with predictors, intercept should not be present
## - K: number of folds in K-fold cross-validation
## - seed: random number generator seed (optional)
## - criterion: the criterion to select the penalization parameter, either cross-val or BIC or EBIC
## Output
## - pred: cross-validated predictions for y
## - ssr: residual sum of squares, sum((y-pred)^2)
  require(glmnet)
  if (!missing(seed)) set.seed(seed)
  subset <- rep(1:K,ceiling(nrow(x)/K))[1:nrow(x)]
  subset <- sample(subset,size=nrow(x),replace=FALSE)
  pred <- double(nrow(x))
  num_betas <- c()
  cat("Starting cross-validation")
  if (ncol(x)>0) {  #if there are some covariates
    for (k in 1:K) {
        sel <- subset==k
        if (criterion=='cv') {
            fit <- cv.glmnet(x=x[!sel,,drop=FALSE], y=y[!sel], alpha = 1, nfolds=10)
            b= as.vector(coef(fit,s='lambda.min'))
            num_betas <- c(num_betas, sum(b!= 0))
            pred[sel] <- b[1] + x[sel,,drop=FALSE] %*% as.matrix(b[-1])
        } else if (criterion=='bic'){
            fit <- lasso.bic(y=y[!sel],x=x[!sel,,drop=FALSE])
            pred[sel] <- fit$coef[1] + x[sel,,drop=FALSE] %*% matrix(fit$coef[-1],ncol=1)
            num_betas <- c(num_betas, sum(fit$coef!=0))
        } else { stop("method.lambda not implemented") }
        cat(".")
    }
  } else { #if there are no covariates, just use the intercept
    for (k in 1:K) {
      sel <- subset==k
      pred[sel] <- mean(y[!sel],na.rm=TRUE)
    }
  }
  cat("\n")
  return(list(pred=pred,ssr=sum((pred-y)^2,na.rm=TRUE),num_betas=num_betas))
}

kfoldCV.bms <- function(y,x,K=10,seed) {
  ## Perform K-fold cross-validation for least-squares regression estimate
  ## Input
  ## - y: response
  ## - x: data.frame with predictors, intercept should not be present
  ## - K: number of folds in K-fold cross-validation
  ## - seed: random number generator seed (optional)
  ## Output
  ## - pred: cross-validated predictions for y
  ## - ssr: residual sum of squares, sum((y-pred)^2)
  if (!missing(seed)) set.seed(seed)
  subset <- rep(1:K,ceiling(nrow(x)/K))[1:nrow(x)]
  subset <- sample(subset,size=nrow(x),replace=FALSE)
  pred <- double(nrow(x))
  num_betas <- c()
  for (k in 1:K) {
    sel <- subset==k
    
    non_constant <- c()
    for (column in colnames(x[!sel,,drop=FALSE])){
      non_constant <- c(non_constant,length(unique(x[!sel,column,drop=FALSE])) > 1)
    }
    
    fit <- modelSelection(y=y[!sel],x=x[!sel,non_constant,drop=FALSE], 
                          priorCoef=zellnerprior(taustd=1),
                          priorDelta=modelbbprior(1,1),
                          priorVar=igprior(alpha=.01, lambda=.01),
                          family = 'normal',verbose=FALSE)
    
    #Here we are doing a point estimate for every fold
    b <- coef(fit)[-c(1,nrow(coef(fit))),1]
    pred[sel] <- x[sel,non_constant,drop=FALSE] %*% b + coef(fit)[1,1]
    num_betas <- c(num_betas, length(b[b != 0]))
  }
  return(list(pred=pred,ssr=sum((pred-y)^2,na.rm=TRUE),num_betas=num_betas))
}


