#Imports:
#  survival
#' zoo
#'
#' @param data_list the dataset list which has been transformed
#' @param nMI number of imputation
#' @param covariates categorical variable we would like to conduct the Log-rank test
#' @param strata strata variable used in Log-rank test
#' @return A class object with Log-Rank MI estimation for given data
#' @export
#' @examples
#'
##main function
##time and status variables are derived from long formatted data by using data transform function
LRMI<-function(data_list,nMI,covariates,strata=NULL,...){
  form1="Surv(time,cens)~"
  ##only one categorical varaible allowed in the Log-rank test
  form2=NULL
  if (missing(strata) |is.null(strata)) {
    form<-as.formula(paste("survival::",paste(form1,covariates)))
  }
  if (!missing(strata) && !is.null(strata)) {
    form2=paste(paste("+strata(",paste(strata, collapse="+")),")")
    form<-as.formula(paste(paste("survival::",paste(form1,covariates)),form2))
  }
  data1<-data_list$data_uc
  #data1$trt[1:250]<-2
  #data1$strata<-c(rep(1,600),rep(2,200),rep(3,200))

  #strata=c("strata")

  cate_level<-names(table(data1[,covariates]))
  cate_n<-length(cate_level)

  ##if nMI is null we output the weighted Log-rank estimation instead
  if (!is.null(nMI) && nMI>0){
    n<-nrow(data1)
  stat=matrix(0,nrow=nMI,ncol=cate_n)
  obs=matrix(0,nrow=nMI,ncol=cate_n)
  exp=matrix(0,nrow=nMI,ncol=cate_n)
  ngroup=as.numeric(table(data1[,covariates]))
  vararray=array(0,dim=c(nMI,cate_n,cate_n))
  censor=rep(0,n)
  data1$time=rep(0,n)
  data1$cens=rep(0,n)
  varmean=0
  for (isim in 1:nMI) {
    for (i in 1:n) {
      j=which(cumsum(data_list$weights[[i]])>runif(1))[1]
      censor[i]=j
      data1$time[i]=data_list$time[[i]][j]

      if (censor[i]<data_list$e[i] ) {data1$cens[i]=1}
      else if (censor[i]>=data_list$e[i]) {data1$cens[i]=0}
    }

    lr1=survival::survdiff(form,data=data1,...)
    stat[isim,]=c(lr1$obs-lr1$exp)
    obs[isim,]=c(lr1$obs)
    exp[isim,]=c(lr1$exp)

    vararray[isim,,]=lr1$var
    varmean=varmean+lr1$var
  }

  statmean=colMeans(stat)
  varmean=varmean/nMI
  obsmean=colMeans(obs)
  expmean=colMeans(exp)
  x=statmean
  B=matrix(0,cate_n,cate_n)
  for (isim in 1:nMI) {
    Bi=matrix(stat[isim,]-x,nrow=cate_n,ncol=1)
    B=B+Bi%*%t(Bi)
  }
  B=B/(nMI-1)
  var_MI<-varmean+((nMI+1)/nMI)*B
  stat1<-statmean[1:(cate_n-1)]
  chisq<-stat1%*%solve(var_MI[1:(cate_n-1),1:(cate_n-1)])%*%matrix(stat1)
  p=as.numeric(1-pchisq(chisq,df=(cate_n-1)))


  res<-list(est=x,pvalue=p,var=var_MI,est_mat=stat,Var_mat=vararray,chisq=chisq,
            between_var=B, within_var=varmean,nMI=nMI,df=cate_n-1,obsmean=obsmean,expmean=expmean,
            covariates=covariates,ngroup=ngroup,cate_level=dimnames(lr1$n)[[1]])
  attr(res,"est")=x
  attr(res,"var")=varmean+((nMI+1)/nMI)*B
  attr(res,"Est_mat")=stat
  attr(res,"chisq")=chisq
  attr(res,"Var_mat")=vararray
  attr(res,"Between Var")=B
  attr(res,"Within Var")=varmean
  attr(res,"nMI")=nMI
  attr(res,"pvalue")=p
  attr(res,"df")=df
  attr(res,"covariates")=covariates
  attr(res,"ngroup")=ngroup
  attr(res,"Mean observed events")=obsmean
  attr(res,"Mean expected events")=expmean

  }


  if (is.null(nMI)){
  ##another set of estimation based on Cook (2002) (4) and (5)
  LR_cook<-data.frame()

  for (z in (0:cate_n)){
    Z=data.frame(data1[,c(covariates)])
    names(Z)<-c(covariates)

    if (z>0){
    index<-Z[,c(covariates)]==cate_level[z]
    n=sum(as.numeric(index))
    }
    if (z==0){
      index<-c(1:nrow(Z))
      n=length(index)
    }


    e=data_list$e[index]
    s=data_list$s[index]
    x=data_list$time[index]
    w=data_list$weights[index]


    Z <- as.matrix(Z)
    delta<-vector()
    for (i in (1:n) ) {
      for (j in 1:(e[[i]])) {
        if(j==e[[i]]){
          dij=0
        }
        else{dij=1}
        delta<-c(delta,dij)
      }
    }

    times<-unlist(x)
    weights<-unlist(w)

    Y<-c(rep(cate_level[z],length(times)))
    Y<-matrix(Y,ncol=1,byrow=T)

    # Risk set function
    risk.set <- function(t) {which(times >= t)}
    event.set <- function(t) {which(times == t)}
    # Risk set at each event time
    rs <- apply(as.matrix(times), 1, risk.set)
    es <- apply(as.matrix(times), 1, event.set)

    d <- vector()
    r<-vector()
    for(i in 1:length(rs)) {
      r[i] <- sum(weights[rs[[i]]])
      d[i] <- sum(weights[es[[i]]]*delta[es[[i]]])
    }
    if (z==0){
      LR_cook<-data.frame(times,d,r)
      LR_cook<-LR_cook[order(LR_cook$times),]
    }
    if (z>0){
    temp<-data.frame(times,d,r)
    temp<-temp[order(temp$times),]
    names(temp)<-c("times",paste("d", z, sep = ""), paste("r", z, sep = ""))


    LR_cook<-merge(LR_cook,temp,all.x=TRUE,by = c("times"))

    }

  }
  name_index<-names(LR_cook)
  ##fill in missing values
  for (i in (1:ncol(LR_cook))){
    index<-min(which(!is.na(LR_cook[,i])==TRUE))
    if (index>1 & substr(name_index[i],1,1)=="r"){
      LR_cook[(1:(index-1)),i]=rep(LR_cook[index,i],(index-1))
    }
    if (index>1 & substr(name_index[i],1,1)=="d"){
      LR_cook[(1:(index-1)),i]=rep(0,(index-1))
    }
    if (substr(name_index[i],1,1)=="d"){
      LR_cook[is.na(LR_cook[,i]),i] <- 0
    }
    if (substr(name_index[i],1,1)=="r"){
      for (j in (2:nrow(LR_cook))){
        LR_cook[j,i] <- LR_cook[(j-1),i]-LR_cook[(j-1),(i-1)]
      }

    }
  }

  v<-rep(0,cate_n)
  w_d<-matrix(rep(0,cate_n*cate_n),nrow=cate_n)
  for (z in (1:cate_n)){
    v[z]=sum(LR_cook[,(4+(z-1)*2)]-LR_cook[,(5+(z-1)*2)]*(LR_cook[,2])/(LR_cook[,3]))
  }
  for (i in (1:cate_n)){
    for (j in (1:cate_n)){
      if (i==j){
      w_d[i,j]=sum((LR_cook[,2])*(LR_cook[,(5+(i-1)*2)]/(LR_cook[,3]))*(1-LR_cook[,(5+(i-1)*2)]/LR_cook[,3])*
                     ((LR_cook[,3]-LR_cook[,2])/(LR_cook[,3]-1)))
      }
      if (i!=j){
        w_d[i,j]=-sum(LR_cook[,(5+(j-1)*2)]*LR_cook[,(5+(i-1)*2)]*(LR_cook[,2])*(LR_cook[,3]-LR_cook[,2])/
                        ((LR_cook[,3])^2*(LR_cook[,3]-1)))
      }
    }
  }
  lr_c<-t(v[1:(cate_n-1)])%*%solve(w_d[1:(cate_n-1),1:(cate_n-1)])%*%v[1:(cate_n-1)]
  p=as.numeric(1-pchisq(lr_c,df=(cate_n-1)))

  res<-list(est=v,pvalue=p,Var_matrix=w_d,chisq=lr_c,
            df=cate_n-1,nMI=nMI,
            column_name=covariates,cate_n=cate_n,cate_level=cate_level)

  }
  class(res) <- c("LRMI", class(res))
  return(res)
}


