#Imports:
#  survival
#  zoo
#'
#' This function calculates the KM estimator if we would like to incorporate more potential events
#' into analysis
#'
#'
#' @param data_list the data list which has been transformed
#' @param covariates The categorical grouping varaible, if missing then overall KM is returned. Also,
#' if user provide the orig data, the covariates in orig should be the same as it in the transformed data list.
#' @param nMI number of MI
#' @param plot T/F whether output a KM plot, comparing original results and MI results
#' @param data_orig an optional original data set without any uncertainty, this can be supplied by user, this
#' data set must contain time variable and event indicator variable, which can be specified later
#' @param time_var time variable in data_orig, if user provide orig dataset and would like to specify the
#' time and event indicator variable
#' @param event_var event indicator variable in data_orig
#' @return A class object with KM estimation, if plot == TRUE, then KM plot as well
#' @export
#' @examples KMMI(data_list=data_intrim,nMI=100,covariates=c("trt"),plot=TRUE,data_orig=data_orig,
#' time_var=c("time"),event_var=c("cens"))
#'
#'
KMMI<-function(data_list,nMI,covariates=NULL,data_orig=NULL,plot=TRUE,time_var=NULL,event_var=NULL)
  {
  data1<-data_list$data_uc
  n<-nrow(data1)
  if (!is.null(covariates)){
  ngroup<-nlevels(as.factor(data1[,covariates]))
  cate_level<-names(table(data1[,covariates]))
  }
  if (is.null(covariates)){
    ngroup=1
  }
  ##count how many time points in each group
  num_event<-data.frame(data1,count=unlist(data_list$e))

  ##if data_orig is not missing then we need variable for time and event
  if (!is.null(data_orig) && !missing(time_var) && !is.null(time_var)) {
    names(data_orig)[names(data_orig) == time_var] <- "time"
  }

  if (!is.null(data_orig) && !missing(event_var) && !is.null(event_var)) {
    names(data_orig)[names(data_orig) == event_var] <- "cens"
  }
  form1="Surv(time,cens)~"

  if (!missing(covariates) && !is.null(covariates)) {
  form2=paste(covariates, collapse="+")
  }
  if (missing(covariates) | is.null(covariates)) {
    form2=1
  }
  form<-as.formula(paste("survival::",paste(form1,form2)))

  ##resultant KM estimation is for all potential event timepoints
  timelist<-NULL
  trtlist<-NULL
  for (i in (1:n)){
    temp<-c(unlist(data_list$time[i]))
    #timelist<-c(timelist,temp[-(length(temp))])
    timelist<-c(timelist,temp)
  }
  ntime<-length(timelist)
  if (ngroup==1){
    trtlist<-rep(1,ntime)
  }
  if (ngroup>1){
  for (i in (1:ngroup)){
    temp<-rep(cate_level[i],sum(num_event[which(num_event[,c(covariates)]==cate_level[i]),]$count))
    trtlist<-c(trtlist,temp)
    }
  }


  frame<-data.frame(time=timelist,covariates=trtlist)
  frame<-frame[order(frame$covariates,frame$time),]
  timemin<-aggregate(frame$time, by=list(frame$covariates), FUN=min)


  ##if nMI is null we output the weighted KM estimation instead
  if (!is.null(nMI) && nMI>0){
    km_res<-matrix(rep(0, nMI*ntime),ncol = nMI)
    km_sd<-matrix(rep(0, nMI*ntime),ncol=nMI)
    for (z in 1:nMI)
    {

      data_temp<-data_list$data_uc
      data_temp$time=rep(0,n)
      data_temp$cens=rep(0,n)
      censor<-rep(0,n)
      for (i in 1:n) {

        #generate censor indicator, if cum prob > runif[1] then event
        j=which(cumsum(data_list$weights[[i]])> runif(1))[1]
        censor[i]=j
        data_temp$time[i]=data_list$time[[i]][j]

        if (censor[i]<data_list$e[i] ) {data_temp$cens[i]=1}
        else if (censor[i]>=data_list$e[i]) {data_temp$cens[i]=0}
      }


      KM <- survival::survfit(form, type="kaplan-meier", conf.type="log",data=data_temp)
      cov<-NULL
      for (i in (1:ngroup)){
        if (ngroup>1){
        temp<-c(rep(cate_level[i],KM$strata[i]))
        cov<-c(cov,temp)
        }
        if (ngroup==1){
          cov=rep(1,length(KM$time))
        }
      }
      km_res_temp<-data.frame(time=KM$time,surv=KM$surv,sd=KM$std.err*KM$surv,covariates=cov)

      frame<-merge(frame,km_res_temp, by = c("covariates","time"),all.x=TRUE)
      frame<-frame[order(frame$covariates,frame$time),]
      ##the na.locf should by group
      for (i in (1:ngroup)){
        frame$surv[is.na(frame$surv)&frame$time==timemin$x[i]&frame$covariates==timemin$Group.1[i]]<-1
      }

      frame$surv<-zoo::na.locf(frame$surv)
      frame$sd[frame$surv==1]<-0
      frame$sd<-zoo::na.locf(frame$sd)
      km_res[,z]<-c(frame$surv)
      km_sd[,z]<-c(frame$sd)

      frame<-frame[,c("time","covariates")]

    }
    #between MI estimators variances
    b<-apply(km_res,1,var,na.rm=TRUE)
    #sum of within MI estimator variance
    w<-apply(km_sd^2,1,sum,na.rm=TRUE)
    #final KM estimaor and variance

    km_mi<-data.frame(frame, est=apply(km_res,1,mean,na.rm=TRUE),
                      sd_total=sqrt(w/nMI+(1+1/nMI)*b),sd_b=sqrt(b),sd_w=sqrt(w/nMI))

  }
 # test<-data.frame(time=km_mi$time,km_res,est=km_mi$est,trt=km_mi$covariates)
 # write.csv(df_x, file = "C:/Users/Yiming.Chen/Documents/adjudication/example.csv",row.names=FALSE)


  if (!is.null(data_orig) ){
fit_orig<-survival::survfit(form, data = data_orig)
names(data_orig)[names(data_orig) == covariates] <- "covariates"

}

##another set of generalized KM estimator based on Cook (2000) fromula (1) and (3)
KM_cook<-data.frame()


for (z in (1:ngroup)){
if (ngroup >1){
  Z=data.frame(data1[,c(covariates)])
  names(Z)<-c(covariates)
index<-Z[,c(covariates)]==cate_level[z]
Z <- as.matrix(Z)
}
if (ngroup ==1){
  index<-data_list$e>=1
}
n=sum(as.numeric(index))

e=data_list$e[index]
s=data_list$s[index]
x=data_list$time[index]
w=data_list$weights[index]



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

if (ngroup >1){
Y<-c(rep(cate_level[z],length(times)))
}
if (ngroup ==1){
  Y<-rep(1,length(times))
}
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

  ##variance
 w_var<-matrix(rep(0,n*length(times)),ncol=length(times))
 s_var<-matrix(rep(0,n*length(times)),ncol=length(times))
    for (j in (1:n)){
      for (i in (1:length(times))){
        if (sum(as.numeric((x[[j]]==times[i])))>0){
      w_var[j,i]<-delta[i]*w[[j]][which((x[[j]]==times[i]))]
        }
        else  {w_var[j,i]=0}
        if (sum(which((x[[j]]<=times[i])))>0){
          s_var[j,i]<-delta[i]*s[[j]][max(which((x[[j]]<=times[i])))]
        }
        else {
          s_var[j,i]<-1
        }
    }

    }

  ##row stands for subjects, column for event time
  d_matrix<-matrix(rep(d,n),nrow=n,byrow=T)
  #d_matrix[1,]
  r_matrix<-matrix(rep(r,n),nrow=n,byrow=T)
  #r_matrix[,1]
  ##this matrix is subjects' variance at each time point, if non-event time, the column will be all 0
  Cook_var<-(w_var-s_var*d_matrix/r_matrix)/(r_matrix-d_matrix)
  ##order to match event time
  index<-order(times)
  Cook_var<-Cook_var[,index]

  ##for tk < t, row cumsum
  Cook_sd<-sqrt(colSums(t((apply(Cook_var,1,cumsum))^2)))

  temp<-data.frame(times,d,Y,r,prod=(1-(d/r)))
  temp<-temp[order(temp$Y,temp$times),]
  temp$est<-cumprod(temp$prod)
  temp$sd<-Cook_sd

  KM_cook<-rbind(temp,KM_cook)

}


##depends on data_orig, nMI number, we output corresponding plots, should we add CI?
if (plot == TRUE && !is.null(nMI)){

  if (is.null(data_orig)){
  plot(1, type="n", main=expression(paste("MI Kaplan-Meier Plot")),
       xlab="Time", ylab="Survival Probability", xlim=c(0, max(km_mi$time)),
       ylim=c(0, 1))
  }
  if (!is.null(data_orig)){
    plot(fit_orig, main=expression(paste("MI Kaplan-Meier Plot")),
         xlab="Time", ylab="Survival Probability", lwd=2, col=1:ngroup)
  }
  if (ngroup>1){
  legend(x="topright", col=1:ngroup, lwd=2, legend=cate_level)

  for (i in (1:ngroup)){
    lines(est~time,data=km_mi[km_mi$cov==cate_level[i],],type="S", col = i, lty = 5)
   #lines(est~times,data=KM_cook[KM_cook$Y==cate_level[i],],type="s", col = i, lty = 1)
  }
  }
  if (ngroup==1){
    lines(est~time,data=km_mi,type="S", col = 1, lty = 5)
  }

}
if (plot == TRUE && is.null(nMI)){

  plot(1, type="n", main=expression(paste("Weighted Kaplan-Meier Plot")),
                           xlab="Time", ylab="Survival Probability", xlim=c(0, max(KM_cook$times)),
                           ylim=c(0, 1))
  if (ngroup>1){
  legend(x="topright", col=1:ngroup, lwd=2, legend=cate_level)

  for (i in (1:ngroup)){
    lines(est~times,data=KM_cook[KM_cook$Y==cate_level[i],],type="S", col =i, lty = 6)
  }
  }
  if (ngroup==1){
    lines(est~times,data=KM_cook,type="S", col = 1, lty = 6)
  }

}

if (!is.null(nMI) && ngroup>1){
  return(list(KM_mi=km_mi,
              #fraction_km=new,
              KM_cook=KM_cook,
              ngroup=ngroup,cate_level=cate_level,nMI=nMI,covariates=covariates))
}
if (!is.null(nMI) && ngroup==1){
  return(list(KM_mi=km_mi,
              KM_cook=KM_cook,
              nMI=nMI))
}
if (is.null(nMI)){
return(list(ngroup=ngroup,KM_cook=KM_cook,cate_level=cate_level))
}
}

