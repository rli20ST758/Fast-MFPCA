
## The function implements the face algorithm 
face.Cov <- function(Y, argvals, A0, B, Anew, Bnew, G_invhalf, s, Cov=FALSE, pve=0.99, npc=NULL, lambda=NULL, alpha=0.7, 
                     search.grid=TRUE, search.length=100, lower=-20, upper=20){
  
  ######## precalculation for missing data ########
  imputation <- FALSE
  Niter.miss <- 1
  L <- ncol(Y)
  n <- nrow(Y)
  
  Index.miss <- is.na(Y)
  if(sum(Index.miss)>0){
    num.miss <- rowSums(is.na(Y))
    for(i in 1:n){
      if(num.miss[i]>0){
        y <- Y[i,]
        seq <- (1:L)[!is.na(y)]
        seq2 <-(1:L)[is.na(y)]
        t1 <- argvals[seq]
        t2 <- argvals[seq2]
        fit <- smooth.spline(t1,y[seq])
        temp <- predict(fit,t2,all.knots=TRUE)$y
        if(max(t2)>max(t1)) temp[t2>max(t1)] <- mean(y[seq])
        if(min(t2)<min(t1)) temp[t2<min(t1)] <- mean(y[seq])
        Y[i,seq2] <- temp
      }
    }
    imputation <- TRUE
    Niter.miss <- 100
  }
  convergence.vector <- rep(0,Niter.miss)
  iter.miss <- 1
  lambda.input <- lambda
  totalmiss <- mean(Index.miss)
  
  
  while(iter.miss <= Niter.miss&&convergence.vector[iter.miss]==0) {
    ###################################################
    ######## Transform the Data           #############
    ###################################################
    Ytilde <- t(as.matrix(Y%*%B) %*% A0)
    C_diag <- rowSums(Ytilde^2)
    
    ###################################################
    ########  Select Smoothing Parameters #############
    ###################################################
    Y_square <- sum(Y^2)
    Ytilde_square <- sum(Ytilde^2)
    face_gcv <- function(x) {
      lambda <- exp(x)
      lambda_s <- (lambda*s)^2/(1 + lambda*s)^2
      gcv <- sum(C_diag*lambda_s) - Ytilde_square + Y_square
      trace <- sum(1/(1+lambda*s))
      gcv <- gcv/(1-alpha*trace/L/(1-totalmiss))^2
      return(gcv)
    }
    
    
    if(is.null(lambda.input) && iter.miss<=2) {
      if(!search.grid){
        fit <- optim(0,face_gcv,method=method,lower=lower,upper=upper,control=control)
        if(fit$convergence>0) {
          expression <- paste("Smoothing failed! The code is:",fit$convergence)
          print(expression)
        }
        lambda <- exp(fit$par)
      } else {
        Lambda <- seq(lower,upper,length=search.length)
        Length <- length(Lambda)
        Gcv <- rep(0,Length)
        for(i in 1:Length)
          Gcv[i] <- face_gcv(Lambda[i])
        i0 <- which.min(Gcv)
        lambda <- exp(Lambda[i0])
      }
    }
    YS <- matrix.multiply(Ytilde,1/(1+lambda*s),2)
    
    ###################################################
    ####  Eigendecomposition of Smoothed Data #########
    ###################################################
    temp0 <- YS%*%t(YS)/n
    temp <- as.matrix(Anew%*%as.matrix(temp0%*%t(Anew)))
    Eigen <- eigen(temp,symmetric=TRUE)
    A = Eigen$vectors
    Phi = Bnew %*% A
    Sigma = Eigen$values
   
    if(iter.miss>1&&iter.miss< Niter.miss) {
      diff <- norm(YS-YS.temp,"F")/norm(YS,"F")
      if(diff <= 0.02)
        convergence.vector[iter.miss+1] <- 1
    }
    
    YS.temp <- YS
    iter.miss <- iter.miss + 1
    N <- min(n, ncol(B))
    d <- Sigma[1:N]
    d <- d[d>0]
    per <- cumsum(d)/sum(d)
    N <- ifelse (is.null(npc), min(which(per>pve)), min(npc, length(d)))
    
    #########################################
    #######     Principal  Scores   #########
    ########   data imputation      #########
    #########################################
    if(imputation) {
      Phi.N <- Phi[,1:N]
      A.N <- G_invhalf %*% A[,1:N]
      d <- Sigma[1:N]
      sigmahat2  <-  max(mean(Y[!Index.miss]^2) -sum(Sigma),0)
      if(N>1){
        Xi <- solve(t(Phi.N)%*%Phi.N + diag(sigmahat2/d)) %*% t(as.matrix(Y%*%B) %*% A.N)
      } else{
        Xi <- solve(t(Phi.N)%*%Phi.N + sigmahat2/d) %*% t(as.matrix(Y%*%B) %*% A.N)
      }
      Yhat <- t(Phi.N %*% Xi)
      Y <- Y*(1-Index.miss) + Yhat*Index.miss
      if(sum(is.na(Y))>0) print("error")
    }
    
  } ## end of while loop
  
  Phi.N <- Phi[,1:N]
  evalues <- Sigma[1:N]
  Ktilde <- NULL
  if(Cov) {
    Ktilde <- Phi.N %*%  matrix.multiply(t(Phi.N),evalues,2)
  }
  
  return(list(Yhat=Y, decom=temp, Ktilde=Ktilde, evalues=evalues, efunctions=Phi.N))
}




## a function supposed to be in the refund package
quadWeights<- function(argvals, method = "trapezoidal")
{
  ret <- switch(method,
                trapezoidal = {D <- length(argvals)
                1/2*c(argvals[2] - argvals[1], argvals[3:D] -argvals[1:(D-2)], argvals[D] - argvals[D-1])},
                midpoint = c(0,diff(argvals)),
                stop("function quadWeights: choose either trapezoidal or midpoint quadrature rule"))
  
  return(ret)  
}



pspline.setting <- function(x,knots=select.knots(x,35),p=3,m=2,weight=NULL,type="full",
                            knots.option="equally-spaced"){
  
  # x: the marginal data points
  # knots: the list of interior knots or the numbers of interior knots
  # p: degrees for B-splines, with defaults values 3
  # m: orders of difference penalty, with default values 2
  # knots.option: type of knots placement, with default values "equally-spaced"
  
  #require(splines)
  #require(Matrix)
  
  ### design matrix 
  K = length(knots)-2*p-1
  B = spline.des(knots=knots, x=x, ord=p+1, outer.ok=TRUE,sparse=TRUE)$design
  
  bs = "ps"
  
  if(knots.option == "quantile"){
    bs = "bs"
  }
  
  s.object = s(x=x, bs=bs, k=K+p, m=c(p-1,2), sp=NULL)
  object  = smooth.construct(s.object,data = data.frame(x=x),knots=list(x=knots))
  P =  object$S[[1]]
  if(knots.option == "quantile") P = P / max(abs(P))*10 # rescaling
  
  if(is.null(weight)) weight <- rep(1,length(x))
  
  if(type=="full"){
    
    Sig = crossprod(matrix.multiply(B,weight,option=2),B)
    eSig = eigen(Sig)
    V = eSig$vectors
    E = eSig$values
    if(min(E)<=0.0000001) {#cat("Warning! t(B)%*%B is singular!\n");
      #cat("A small identity matrix is added!\n");
      E <- E + 0.000001;
    }
    Sigi_sqrt = matrix.multiply(V,1/sqrt(E))%*%t(V)
    
    tUPU = Sigi_sqrt%*%(P%*%Sigi_sqrt)
    Esig = eigen(tUPU,symmetric=TRUE)
    U = Esig$vectors
    s = Esig$values
    s[(K+p-m+1):(K+p)]=0
    A = B%*%(Sigi_sqrt%*%U)
  }
  
  if(type=="simple"){
    A = NULL
    s = NULL
    Sigi_sqrt = NULL
    U = NULL
  }
  
  List = list("A" = A, "B" = B, "s" = s, "Sigi.sqrt" = Sigi_sqrt, "U" = U, "P" = P)
  
  return(List)
}





matrix.multiply <- function(A,s,option=1){
  if(option==2)
    return(A*(s%*%t(rep(1,dim(A)[2]))))
  if(option==1)
    return(A*(rep(1,dim(A)[1])%*%t(s)))
}
