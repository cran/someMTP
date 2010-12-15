lsd.test <- function (resp, alternative = 1, null = NULL, D = NULL){

  resp=as.matrix(resp)
  n <- dim(resp)[1]
  p <- dim(resp)[2]
  
  if(length(alternative)==1) alternative <- rep(alternative,dim(resp)[1]) 
  alternative=as.matrix(alternative)
  k <- dim(alternative)[2]

  
  if(!is.null(null))  {
    if(length(null)==1) null <- rep(null,dim(resp)[1]) 
    null <- as.matrix(null)
	h <- dim(null)[2]
    IP0 <- diag(n) - null%*%solve(t(null)%*%null)%*%t(null)
	} else{
    h <- 0
	IP0 <- diag(n) 
   }

  if (is.null(D)) {
   D=diag(t(resp)%*%IP0%*%resp)  
  }else   if(is.function(D)) D <- D(resp=resp,null=null)
  D <- as.matrix(D)  
  q <- dim(D)[2]
   
  H = t(resp)%*%IP0%*%alternative%*%solve(t(alternative)%*%IP0%*%alternative)%*%t(alternative)%*%IP0%*%resp
  G <- t(resp)%*%IP0%*%resp - H
  DHD	 <- t(D) %*% H %*% D
  DGD <- t(D) %*% G %*% D
  df=c(q-1+k,n-h-k+1-q)
  F <- df[2]/df[1]*sum(diag(DHD))/sum(diag(DGD))
  p <- 1-pf(F,df[1],df[2])
  list(F=F,df=df,p=p,D=D)
}