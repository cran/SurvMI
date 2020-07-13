#Imports:
#'
#' This function simulate two types of events with uncertainty based on exponential distribution,
#' can be extended to multiple groups case or use different distribution later
#'
#' @param n number of total subjects we would like
#' @param true_hr true hr between trt and control
#' @param haz_c event rates in control group
#' @return A dataset with one potential event per row
#' @export
#' @examples df_x<-data_sim(n=1000,0.8,haz_c=0.5/365)
#' df_y<-data_sim2(data_list=data_intrim,covariates=c("trt"),percentage=0.2)
#'
#'


data_sim<-function(n=200,true_hr=0.8,haz_c=1/365){

id = NA
x=w=p=list()


for (i in 1:(n/2)) {
  starttime=365*runif(1) #time from start of study till this subject is randomized - one year accrual
  maxtime=3*365-starttime #follow up time - 3 years
  xi=rexp(1,haz_c) #use exponential distribution to generate time to event variable
  if (xi>maxtime) { #if event time larger than follow up time then censored
    pi=0
    wi=1
    xi=maxtime
  } else
  { #if first random generated event time does not exceed follow up time, then multiple event time will be generated
    pi=runif(1)
    wi=pi
    sumwi=wi
    xij=xi+rexp(1,haz_c)
    xi=c(xi,xij)
    while (xij<maxtime) {
      pij=runif(1)
      wij=pij*(1-sumwi)
      wi=c(wi,wij)
      pi=c(pi,pij)
      sumwi=sumwi+wij
      xij=xij+rexp(1,haz_c)
      xi=c(xi,xij)
    }
    xi[length(xi)]=maxtime
    wi=c(wi,1-sumwi)
  }
  x=c(x,list(xi))
  w=c(w,list(wi))
  p=c(p,list(pi))
  id[i] = i
}



for (i in (n/2+1):n) {
  starttime=365*runif(1) #time from start of study till this subject is randomized
  maxtime=3*365-starttime
  xi=rexp(1,true_hr*haz_c) # a different exponential distribution will be used for event time simulation
  if (xi>maxtime) {
    wi=1
    pi=0
    xi=maxtime
  } else
  {
    pi=runif(1)
    wi=pi
    sumwi=wi
    xij=xi+rexp(1,true_hr*haz_c)
    xi=c(xi,xij)
    while (xij<maxtime) {
      pij=runif(1)
      wij=pij*(1-sumwi)
      wi=c(wi,wij)
      sumwi=sumwi+wij
      xij=xij+rexp(1,true_hr*haz_c)
      xi=c(xi,xij)
      pi=c(pi,pij)
    }
    xi[length(xi)]=maxtime
    wi=c(wi,1-sumwi)
  }
  x=c(x,list(xi))
  w=c(w,list(wi))
  p=c(p,list(pi))
  id[i] = i
}

trt=c(rep(0,n/2),rep(1,n/2))

trt_long<-c(rep(0,length(unlist(x))))
time_long<-c(rep(0,length(unlist(x))))
prob_long<-c(rep(0,length(unlist(x))))
id_long<-c(rep(0,length(unlist(x))))
cum<-0
for (i in (1:length(x))){

  index<-cum+length(unlist(x[[i]]))
  trt_long[(cum+1):index]=trt[i]
  id_long[(cum+1):index]=id[i]
  time_long[(cum+1):index]=matrix(unlist(x[[i]]))
  if (length(unlist(x[[i]]))>1){
  prob_long[(cum+1):index]=matrix(c(unlist(p[[i]]),0))
  }
  if (length(unlist(x[[i]]))==1){
    prob_long[(cum+1):index]=0
  }
  cum=index
}

df_x <- data.frame(id_long,trt_long, time_long, prob_long)

return(df_x)
}

##a slightly different simulation function based on previous simulation
##we asked user to input a data in long format with event prob for each record,
##we can simplify it into this format with more potential events only,
##or subject can supply a simliar long format with most subjects only have one record with certainty
##some censored subjects would have one or more potential events
##this may only work for one covariate case, but covariate can have multiple levels
##covariate here is the covariate we would like to put restriction on constant hazard assumption
##ie. the variable we specify hazard ratio on in the data_sim function
data_sim2<-function(data_list,covariates,percentage){
    n<-nrow(data_list$data_uc)
    ##for those have potential events before censoring, we pick an event with largest prob in MI
    ##makes that event the "certain" event in this simulation
    km_orig_index<-unlist(lapply(data_list$weights,which.max))
    censor_index<-unlist(data_list$e)
    ##true censored subjects has km_orig_index==censor_index and we should not touch them
    ##data1 is non-censored subjects
    data1<-data.frame(data_list$data_uc,test=as.numeric(km_orig_index<censor_index))
    ##data2 is censored subjects with potential events - the set we are going to manipulate
    ##however, this set is comparable small
    data2<-data.frame(data_list$data_uc,test=as.numeric(censor_index==km_orig_index&censor_index>1))

    #data2<-data.frame(data_list$data_uc,test=as.numeric(censor_index==km_orig_index))
    ##calculate event rates per group
    r<-prop.table(table(data1[,c(covariates,"test")]),1)[,2]

    ngroup<-nlevels(as.factor(data1[,covariates]))
    cate_level<-names(table(data1[,covariates]))
    sgroup<-cumsum(as.numeric(table(data1[,covariates])))
    ##simulation prob
    ##for censored subjects with potential events, how much of them are simulated to have events
    np<-table(data1[,c(covariates,"test")])[,2]/table(data2[,c(covariates,"test")])[,2]
      km_orig_time<-NULL
      censor_time<-NULL
      p<-rep(NA,n)
      km_time_p<-rep(NA,n)
      np1<-rep(NA,ngroup)
    for (j in 1:ngroup){
      if (j==1){
        start=1
        end=sgroup[j]
        np1[1]<-percentage*runif(1)
      }

      else if (j>1) {start=(sgroup[j-1]+1)
      end=sgroup[j]
      np1[j]<-np1[1]*(np[j]/np[1])}
      if (np1[j]>1) {
        stop("Not enough potential events for simulation, specify smaller percentage or change event rate")
      }
    for (i in (start:end)){
      temp1<-data_list$time[[i]][km_orig_index[i]]
      temp2<-data_list$time[[i]][censor_index[i]]
      km_orig_time<-c(km_orig_time,temp1)
      censor_time<-c(censor_time,temp2)

      ##the imputed potential events ratio are expected to be the same as certain events
      if (censor_index[i]==km_orig_index[i]&censor_index[i]>1&(runif(1)<=np1[j])){
      #all censored subjects with potential events now have event probability
        #if (censor_index[i]==km_orig_index[i]&runif(1)<=np1[j]){
        p[i]=max(data_list$prob[[i]][1:km_orig_index[i]])
        #p[i]=runif(1)
        km_time_p[i]<-data_list$time[[i]][which.max(data_list$prob[[i]][1:km_orig_index[i]])]
      }
    }
    }

    temp<-data.frame(time=km_orig_time,cens=(1-as.numeric(censor_index==km_orig_index)),
                        prob=p,data_list$data_uc,km_time_p,censor_time)

    var_list<-colnames(data_list$data_uc)
    more<-temp[!is.na(temp$km_time_p),]
    #table(more$trt)
    more1<-data.frame(time=more$km_time_p,prob=more$p,more[,var_list])
    ##for those with certain events, add records for censoring
    temp2<-temp[temp$cens==1,]
    more2<-data.frame(time=temp2$censor_time,prob=(1-temp2$cens),temp2[,var_list])
    temp$prob<-temp$cens
    temp<-temp[,names(more1)]
    final<-rbind(temp,more1,more2)
  return(final)
}

