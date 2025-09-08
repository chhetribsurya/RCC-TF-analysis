get.scores <- function(beta, se, nsamp,mu.y = 0,mu.x = 0){
  
  var.x <- 1/(nsamp)/(se^2)
  cov.xy <- beta*var.x
  sc <- nsamp*(cov.xy + mu.x*mu.y)
  
  return(sc)
}

############

kernel_stats <- function(beta,se,nsamp,kernel = NULL, mu.y = 0,mu.x = 0,weights = NULL){
  
  if(length(beta) ==0 || length(se) == 0){
    stop("non zero beta and se vectors needed")
    
  }
  n.snp <- length(beta)
  if(length(se) != n.snp)
  {
    stop("length of beta and se vectors don't match !!")
  }
  
  if(is.null(kernel)){
    kernel <- diag(1,n.snp);
    print("no kernel provided. defaulting to independent effects");
  }
  
  qx <- get.scores(beta = beta, se = se, nsamp = nsamp,mu.y = mu.y,mu.x = mu.x);
  
  if(is.null(weights)){
    
    weights <- rep(1,n.snp)
  }
  
  qx <- qx*weights
  
  tstat <- t(qx)%*%kernel%*%qx
  
  return(tstat)  
  
}


#############

get.weights.beta <- function (MAF, weights.beta = c(1,1)){
  
  n <- length(MAF)
  weights <- rep(0, n)
  idx <- which(MAF == 0)
  if (length(idx) == n) {
    stop("No polymorphic SNPs")
  }
  else if (length(idx) == 0) {
    weights <- dbeta(MAF, weights.beta[1], weights.beta[2])
  }
  else {
    weights[-idx] <- dbeta(MAF[-idx], weights.beta[1], 
                           weights.beta[2])
  }
  return(weights)
}

##############

get.w.mat <- function(maf,n.samp, kernel = NULL, weights = NULL,is.burden = FALSE){
  
  n.snp <- length(maf)

  if(is.null(weights)){
    
    weights <- rep(1,n.snp)
  }
  
  if(is.null(kernel)){
    kernel <- diag(1,n.snp);
    print("no kernel provided. defaulting to independent effects");
  }
  
  ZTZ <- n.samp*diag(2*maf*(1-maf))
  WZ <- diag(weights)%*%ZTZ%*%diag(weights)  
  
  if (is.burden) {
    A1 <- sum(WZ) - mean(WZ)
  }else{
    
    
    L <- chol(kernel, pivot = TRUE)
    D <- L%*%WZ%*%t(L)
    A1 <- D - mean(D)
  }

  return(A1)
  
}

################


