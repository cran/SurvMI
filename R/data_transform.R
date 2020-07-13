#Imports:
#
#'
#' This function transforms data from long format to list. Computing corresponding weights and survival
#' probability. This function should be called if we would like to calculate the weighted Cox model
#'
#' @param data the dataset in long format with a row for each potential event
#' @param var_list is the list of identification variables, such as: c("id_long","trt_long")
#' @param time is the time variable need to be transformed, e.g. time_long
#' @param prob is the prob variable need to be transformed, e.g. prob_long
#' @param var_list_new is the new names for the id variables, if missing, previous names would be used
#' @return A tranformed list object
#' @export
#' @examples data_intrim<-uc_data_transform(data=df_x, var_list=c("id_long","trt_long"),
#' var_list_new=c("id","trt"),time="time_long", prob="prob_long")
#'

uc_data_transform<- function (data,var_list,var_list_new,time,prob){
##how many unique row we should create
id_w=unique(data[,c(var_list)])

n=nrow(id_w)

if (sum(data[,prob]==0)!= n) {
  stop("There must be one censoring record for each subject with event prob = 0")
}

if (sum(!complete.cases(data))>0) {
  stop("No missing data allowed")
}

if (sum(data[,prob]<0)>0|sum(data[,prob]>1)>0) {
  stop("Event probability must be within 0 and 1")
}

if (!is.null(var_list_new)&&length(var_list)!=length(var_list_new)) {
  stop("Variable name length should be the same for identification variables")
}

if (is.null(var_list)) {
  stop("Identification variables cannot be missing")
}

if (is.null(time)) {
  stop("Time variable cannot be missing")
}

if (is.null(prob)) {
  stop("Prob variable cannot be missing")
}

time_p=weights=prob_e=e_list=surv=list()

form1=paste(paste("cbind(Freq =",time),")~")
form2=paste(var_list, collapse="+")
form<-as.formula(paste(form1,form2))
e<-aggregate(form, data = data, FUN = function(x){NROW(x)})

##it's important for e to presever same order as id_w
id_w$order  <- 1:nrow(id_w)
e<-merge(id_w,e)
e<-e[order(e$order), ]
for (i in 1:n) {
  temp<-merge(data,e[i,var_list])
  temp<-temp[order(temp[,c(time)]),]
  if (e$Freq[i]==1) { #if there is only one event
    prob_e[i]=list(0)
    time_p[i]=list(temp[,c(time)])
  } else{

    time_p[i]=list(temp[,c(time)])
    prob_e[i]=list(temp[,c(prob)])
    ##this is inefficient
   #time_p[i]=list(data[interaction(data[,var_list]) %in% interaction(e[i,var_list]),c(time)])
   #prob_e[i]=list(data[interaction(data[,var_list]) %in% interaction(e[i,var_list]),c(prob)])
  }
}

##calculate weights and survival prob
for (i in 1:n) {
  prod=0
  sumwij=0
  w=NULL
  prod=NULL
  s=NULL
  if (e$Freq[i]==1) { #if there is only one time then censored
   w=1
   s=1
   prod=1
  } else{
    for (j in (1:e$Freq[i])){
      if (j == 1) {wij = prob_e[[i]][1]
      prod=(1-prob_e[[i]][1])
      sumwij=prob_e[[i]][1]
      s=1}
      else if (j > 1 & j < e$Freq[i]){
        wij = prod*prob_e[[i]][j]
        s=c(s,prod)
        prod=prod*(1-prob_e[[i]][j])
        sumwij=sumwij+wij
        }
      else if (j == e$Freq[i]){wij = 1-sumwij
      s=c(s,prod)}
      w=c(w,wij)
    }
  }
  weights=c(weights,list(w))
  e_list=c(e_list,list(e$Freq[i]))
  surv=c(surv,list(s))
}

##output this wide format maynot be necessary, the intrim dataset will only be used in calculation
##
#id_w$time<-NA
#id_w$prob<-NA
#id_w$weights<-NA
#for (i in 1:n) {
#  id_w$time[i]<-time_p[i]
#  id_w$prob[i]<-prob_e[i]
#  id_w$weights[i]<-weights[i]
#}
id_w<-id_w[,c(var_list)]
if (!is.null(var_list_new)){
names(id_w)<-var_list_new
}

uc_list_all<-list(time=time_p,prob=prob_e,weights=weights,e = e_list,s=surv,data_uc=id_w,data_long=data)
return(uc_list_all)
}

