#Juergen Laeuter (1996) Exact t and F Tests for Analyzing Studies with Multiple Endpoints. BIOMETRICS 52, 964-970 
#theorem 2
lsd.test <- function (resp, alternative = rep(1,dim(as.matrix(resp))[2]), null = NULL, D = function(resp) t(as.matrix(diag(t(resp)%*%resp)))){
#resp is a nxp matrix (not a pxn as in the original article)
#null is the covarriates nxc matrix NOT IMPLEMENTED YET
#alternative vector with labels of the 2 samples
# is the coaviariates to be tested. for the moment it is not used. itdeally it is a matrix(1,n,1)
#D is a qxp transformation matrix

  n <- dim(resp)[1]
  p <- dim(resp)[2]
  
  alternative=as.matrix(alternative)
  
  if(!is.null(null))  {
    null=as.matrix(null)
    H0 = diag(n) - null%*%solve(t(null)%*%null)%*%t(null)
    CH = eigen(H0)
    V = CH$vectors
    L = diag(round(CH$values,4))
    W = V%*%L 
	n <- n-dim(null)[2]
    alternative = (t(W)%*%H0%*%alternative)[1:n,,drop=F]
    resp = (t(W)%*%H0%*%resp)[1:n,,drop=F]
  }
  
  if(is.function(D)) D <- D(resp)
  q <- dim(D)[1]
  
  HH = alternative%*%solve(t(alternative)%*%alternative)%*%t(alternative) 
  yHy <- t(resp)%*%(HH)%*%resp
  yGy <- t(resp)%*%(diag(n) - HH)%*%resp
  DyHyD <- D%*% yHy %*% t(D)
  DyGyD <- D%*% yGy %*% t(D)
  
  F <- (n-q+1-dim(alternative)[2])/(q-1+dim(alternative)[2])*sum(diag(DyHyD%*%solve(DyGyD)))
  p <- 1-pf(F,q-1+dim(alternative)[2],n-q+1-dim(alternative)[2])
  list(p=p,F=F)
}
