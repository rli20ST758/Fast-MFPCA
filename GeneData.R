
################################################################################
# I: the number of subjects.
# J: the mean number of visits for each subject.
# L: the number of observations in each curve.
# design: generate regular or irregular spaced data, default as "regular". 
# level: the level of sparse
# sigma: the standard deviation of random errors
# balanced: indicate whether to generate balanced data
################################################################################
GeneData <- function(I = 100, J = 3, L = 100, design = "regular", level = 0.1, 
                     sigma = 1, balanced = FALSE){
  
  K1 <- 4
  K2 <- 4
  K <- K1 + K2
  lambda1 <- 0.5^(0:(K1-1))
  lambda2 <- 0.5^(0:(K2-1))
  tlength <- 1
  t <- seq(0, tlength, length = L)
  tt <- t/tlength
  
  # Eigenfunctions
  f1 <- matrix(0, nrow=K1, ncol=L)
  for ( i in 1:(K1/2) ) {
    f1[2*i-1,] <- sqrt(2/tlength)*sin(i*tt*2*pi)
    f1[2*i,] <- sqrt(2/tlength)*cos(i*tt*2*pi)
  }

  f2 <- matrix(0, nrow=K2, ncol=L)
  f2[1,] <- rep(1, L)*sqrt(1/tlength)
  f2[2,] <- sqrt(3/tlength) * (2*tt - 1)
  f2[3,] <- sqrt(5/tlength) * (6*tt^2 - 6 * tt + 1)
  f2[4,] <- sqrt(7/tlength) * (20*tt^3 - 30*tt^2 + 12 * tt -1)

  # Generate scores
  ## generate number of visits for each subject from poisson distribution
  if(balanced == FALSE){
    J_subj <- pmax(rpois(I, J), 1)
  }else{
    J_subj <- rep(J, I)
  }
  n <- sum(J_subj)
  si1 <- matrix(0, nrow=I, ncol=K1)
  si2 <- matrix(0, nrow=n, ncol=K2)
  for(k in 1:K1) {
    si1[,k] <- rnorm(I, sd=sqrt(lambda1[k]))
  }
  for(k in 1:K2) {
    si2[,k] <- rnorm(n, sd=sqrt(lambda2[k]))
  }
  # Generate errors
  epsilon <- matrix(rnorm(n*L,sd=sigma),nc=L)
  
  
  # Generate dense data
  Y0 <- matrix(0,nrow=n,ncol=L)
  J_ind <- c(0, cumsum(J_subj))
  for(m in 1:I) {
    temp <- apply( ( si1[m,] %*% t(rep(1,L)) ) * f1, 2,sum)
    for(j in 1:J_subj[m]) {
      Y0[J_ind[m]+j ,] <- temp + apply( ( si2[J_ind[m]+j ,] %*% t(rep(1,L)) ) * f2, 2,sum) + epsilon[J_ind[m]+j,]
    }
  }
  
  # Generate sparse data
  if (design == "regular") {
    Y <- Y0
  } else {
    nobs <- floor(level*L)
    Y <- matrix(NA,nrow=n,ncol=L)
    for (i in 1:n) {
      idx <- sample(1:L, nobs)
      Y[i,idx] <- Y0[i,idx]
    }
  }
  
  # return values
  evalues <- list(level1=lambda1, level2=lambda2)
  eigenfunctions <- list(level1=t(f1), level2=t(f2))
  id <- rep(1:I, J_subj)
  
  return(list(Y = Y, evalues = evalues, eigenfunctions = eigenfunctions, id = id))
}



