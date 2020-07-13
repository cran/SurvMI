#Imports:
#  survival
#'
#' This function calculates the weighted Cox model from Snappin (1998) when event uncertainty presents
#' The variance was updated to Cook (2000)
#'
#' @param data_list transformed data list
#' @param covariates list of covariates on the RHS of Cox model
#' @param init list of inital value pass to nlm
#' @return A class object with weighted Cox estimation for given data
#' @export
#' @examples


mlogpl_snapinn=function(beta,p1=p1,n=n,x=x,e=e,Z=Z,s=s,w=w) {
  #calculate minus the log of the partial likelihood on page 210
  #result is a scalar
  #beta is the parameter vector of length p1, the goal is to find the value of beta that minimizes the objective function
  #Z is a matrix of covariates. the number of rows is n=number of subjects. the number of columns is p1=the number of covariates in the model
  #e is a vector of length n containing the e_i=number of potential event times for subject i
  #x is a list with number of elements equal to n. the ith element is a vector of length e_i+1 and the elements are the ordered potential event/censoring times x_i,j
  #w is a list with number of elements equal to n. the ith element is a vector of length e_i+1 and the elements are the weights w_i,j that are between 0 and 1 and add up to 1
  #c1 is a scalar that is called c in the paper and defined at the top of p. 211
  num=0
  den=0
  for (i in 1:n) if (e[[i]]>1) for (j in 1:(e[[i]]-1)) {
    num=num+w[[i]][j]
    den=den+w[[i]][j]^2
  }
  c1=num/den

  res=0;
  grad=rep(0,p1)
  for (i in c(1:n)[e>1]) for (j in 1:(e[[i]]-1)) {
    #    for (i in 1:n) if (e[i]>0) for (j in 1:e[i]) {
    xij=x[[i]][j]
    den=0;num=rep(0,p1)
    for (k in 1:n) for (l in c(1:(e[[k]]))[x[[k]]>=xij]) {
      num=num+w[[k]][l]*Z[k,]*exp(sum(beta*Z[k,]))
      den=den+w[[k]][l]*exp(sum(beta*Z[k,]))
    }
    res=res+w[[i]][j]*(sum(beta*Z[i,])-log(den))
    grad=grad+w[[i]][j]*(Z[i,]-num/den)
    #print(grad)

  }
  res=-c1*res
 # attr(res,"gradient")=-c1*grad
  return(res)

}

#mlogpl_snapinn(beta=1)
#nlm2=nlm(mlogpl_snapinn,p=-.5,hessian=T)

##this is the old function - too slow
mlogpl_cook=function(beta,p1=p1,n=n,x=x,e=e,Z=Z,s=s,w=w) {
  res=0
  grad=rep(0,p1)
  for (i in c(1:n)[e>1])  for (j in 1:(e[[i]]-1)) {
    ##xij is the distinct event time
    xij=x[[i]][j]
    den=0
    for (k in 1:n) {
      ##this actually coincidents with snapinn's method except for the constant c1
      ##sum for subjects with censored time smaller than or equal to xij
      if (sum(as.numeric(x[[k]]>=xij))>0){
      l=min(which(x[[k]]>=xij))

      den=den+s[[k]][l]*exp(sum(beta*Z[k,]))
  }
    }

    res=res+w[[i]][j]*(sum(beta*Z[i,])-log(den))
   # print(res)
    #grad=grad+w[[i]][j]*(Z[i,]-num/den)
  }
  res=-res
  return(res)
}

#mlogpl_cook(beta=1,p1=p1,n=n,x=x,e=e,Z=Z,s=s,w=w)
#nlm1=nlm(mlogpl_cook,p=-0.5,hessian=T)

#A new vector-based log partial likelihood function - which is currently used in the nlm
log.parlik <- function(beta,Y,weights,rs,delta){
  Xbeta <- as.vector(Y%*%beta)
  num <- sum(delta*(weights*Xbeta))
  temp <- vector()
  for(i in 1:length(rs)) {
    temp[i] <- delta[i]*weights[i]*log(sum(weights[rs[[i]]]*exp(Xbeta[rs[[i]]])))
    ##last event is censored
    if (is.na(temp[i])){
      temp[i]=0
    }
  }
  den <- sum(temp)
  return(-num + den)
}

#log.parlik(beta=c(1))

#nlm(log.parlik,p=init,hessian=T)

Coxwt<-function(data_list,covariates,init=NULL,BS=FALSE,nBS=1000){

    if (  inherits(data_list, "list") ==FALSE|!exists('time', where=data_list)
          |!exists('prob', where=data_list)|!exists('e', where=data_list)
          |!exists('s', where=data_list)|!exists('weights', where=data_list)
          |!exists('data_uc', where=data_list)) {
      stop("Check the input data list")
    }


  Z=data.frame(data_list$data_uc[,c(covariates)])
  names(Z)=covariates
  n=nrow(Z)
  ##expand covariates to multiple level for MLE - no intercept
  form2=as.formula(paste("~",paste(covariates, collapse="+")))
  x_beta<-model.matrix(form2,Z)[,-1]

  if (!is.null(dimnames(x_beta)[[2]])){
  p1=length(dimnames(x_beta)[[2]])
  }

   if (is.null(dimnames(x_beta)[[2]])){
    p1=length(covariates)
    x_beta<-matrix(x_beta)
  }

  if (!missing(init) && !is.null(init)) {
    if (length(init) != p1) {
      stop("Wrong length for inital values")
    }
  }
  else {init <- rep(0,p1)}

  #beta=rep(1,p1)
  e=data_list$e
  s=data_list$s
  x=data_list$time
  w=data_list$weights

  times<-vector()
  weights<-vector()
  Y<-vector()
  delta<-vector()
  for (i in (1:n) ) {
    for (j in 1:(e[[i]])) {
      xij=x[[i]][j]
      times<-c(times,xij)
      wij=w[[i]][j]
      weights<-c(weights,wij)
      if(j==e[[i]]){
        dij=0
      }
      else{dij=1}
      delta<-c(delta,dij)
    }
    Y<-c(Y,rep(x_beta[i,],(e[[i]])))
  }
  times <- as.vector(times)
  weights <- as.vector(weights)
  Y<-matrix(Y,ncol=p1,byrow=T)

  # Risk set function
  risk.set <- function(t) {which(times >= t)}
  # Risk set at each event time
  rs <- apply(as.matrix(times), 1, risk.set)
  nlm1=nlm(log.parlik,Y=Y,weights=weights,rs=rs,delta=delta,p=init,hessian=T)
  #nlm2=nlm(mlogpl_snapinn,p=init,hessian=T)
  est<-rep(0,p1)
  var<-rep(0,p1)
  pvalue<-rep(0,p1)
  z_stat<-rep(0,p1)
  column_name<-rep(NA,p1)
  for (i in (1:p1)){
  est[i]<-nlm1$est[i]
  var[i]<-solve(nlm1$hess)[i,i]

  pvalue[i]=2*(1-pnorm(abs(est[i]/sqrt(var[i]))))
  }
  z_stat=est/sqrt(var)
  if (!is.null(dimnames(x_beta)[[2]])){
  column_name<-dimnames(x_beta)[[2]]
  }
  else if (is.null(dimnames(x_beta)[[2]])){
    column_name<-covariates
  }
  if (BS == TRUE){
    var_bs<-rep(0,p1)
    beta_bs<-rep(0,p1)
  ##a BS method to calculate estimation variance
  nr=nrow(data_list$data_uc)
  #ind0=c(1:nr)[data_list$data_uc$trt==0]
  #ind1=c(1:nr)[data_list$data_uc$trt==1]
  #nr0=length(ind0)
  #nr1=length(ind1)

  beta_res<-matrix(rep(0,nBS*p1),ncol=p1,byrow=T)
  for (z in 1:nBS) {
    sample<-c(sample(c(1:nr),nr,replace = T))
    e<-data_list$e[sample]
    x<-data_list$time[sample]
    w<-data_list$weights[sample]
    Z<-data.frame(data_list$data_uc[sample,c(covariates)])
    names(Z)=covariates
    x_beta<-model.matrix(form2,Z)[,-1]
    if (is.null(dimnames(x_beta)[[2]])){
      x_beta<-matrix(x_beta)
    }
    #table(Z)
    times<-vector()
    weights<-vector()
    Y<-vector()
    delta<-vector()
    for (i in (1:n) ) {
      for (j in 1:(e[[i]])) {
        xij=x[[i]][j]
        times<-c(times,xij)
        wij=w[[i]][j]
        weights<-c(weights,wij)
        if(j==e[[i]]){
          dij=0
        }
        else{dij=1}
        delta<-c(delta,dij)
      }
      Y<-c(Y,rep(x_beta[i,],(e[[i]])))
    }
    times <- as.vector(times)
    weights <- as.vector(weights)
    Y<-matrix(Y,ncol=p1,byrow=T)

    rs <- apply(as.matrix(times), 1, risk.set)
    nlmbs=nlm(log.parlik,Y=Y,weights=weights,rs=rs,delta=delta,p=init,hessian=T)
    estbs<-rep(0,p1)
    varbs<-rep(0,p1)

    for (i in (1:p1)){
      estbs[i]<-nlmbs$est[i]
      varbs[i]<-solve(nlmbs$hess)[i,i]

    }
    beta_res[z,]=estbs
  }

  for (i in (1:p1)){
    beta_bs[i]<-mean(beta_res[,i])
    var_bs[i]<-var(beta_res[,i])
  }

    fit<-list(coefficients  = est,
                              hr=exp(est),
                              var    = var,
              coefficients_bs=beta_bs,
              var_bs=var_bs,
              z = z_stat,
                              column_name = column_name,
                              pvalue    = pvalue,
              nBS=nBS)
  }
  if (BS == FALSE){
  fit<-list(coefficients  = est,
            hr=exp(est),
            var    = var,
            z = z_stat,
            column_name = column_name,
            pvalue    = pvalue)
  }
  class(fit) <- 'Coxwt'
  return(fit)
}

