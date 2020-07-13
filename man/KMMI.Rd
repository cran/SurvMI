\name{KMMI}
\alias{KMMI}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Kaplan-Meier estimation with event uncertainty
%%  ~~function to do ... ~~
}
\description{
KM estimation for survival data when event uncertainty presents. KM plot will be output if plot=TRUE specfied.
}
\usage{
KMMI(data_list,nMI,covariates,data_orig = NULL,plot = TRUE,
time_var=NULL,event_var=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data_list}{
 The data list which has been transformed from the long format by uc_data_transform function.
}
  \item{nMI}{
Number of imputations (>1). If missing, weighted statistics would be output instead.
}
  \item{covariates}{
The grouping varaible, no need to be factorized. If missing then the overall KM is returned.
}
  \item{plot}{
T/F, whether output a KM plot, the plot potentially contains KM curves from original dataset and imputed/weighted dataset.
}
  \item{data_orig}{
The original data without any uncertain events. If supplies then user can compare results from certain events only and all possible events.
}
  \item{time_var}{
Time variable in data_orig. If user provides the orig dataset then user need to specify the time and event indicator variable in the orignal dataset.
}
  \item{event_var}{
Event indicator variable in the original data set.
}
}

\value{
\item{KM_mi}{A dataset contains MI estimation and variance at all potential event time}
\item{KM_cook}{A dataset contains weighted KM estimation and variance at all potential event time}
\item{ngroup}{Number of groups}
\item{cate_level}{Values of the categorical variable}
\item{nMI}{Number of imputed datasets}
}
\references{
[1]Cook TD. Adjusting survival analysis for the presence of unadjudicated study events. Controlled clinical trials. 2000;21(3):208-222.

[2]Cook TD, Kosorok MR. Analysis of time-to-event data with incomplete event adjudication. Journal of the american statistical association. 2004;99(468):1140-1152.

[3]Klein JP, Moeschberger ML. Survival Analysis : Techniques for Censored and Truncated Data. New York: Springer; 1997.

[4]Rubin DB. Multiple Imputation for Nonresponse in Surveys. New York: Wiley; 1987
}

\author{
Yiming Chen
}

\seealso{
\code{\link{uc_data_transform}}
}
\examples{
##an example with more potential event case
##data_orig was created as keeping the event with largest weights for individuals
df_x<-data_sim(n=500,0.8,haz_c=0.5/365)
data_intrim<-uc_data_transform(data=df_x,
                               var_list=c("id_long","trt_long"),
                               var_list_new=c("id","trt"),
                               time="time_long",
                               prob="prob_long")
df_y<-data_sim2(data_list=data_intrim,covariates=c("trt"),percentage=1)
data_orig<-df_y[df_y$prob==0|df_y$prob==1,]
data_orig<-data_orig[!duplicated(data_orig$id),]
data_orig$cens<-data_orig$prob

\donttest{
##weighted estimation
KM_res<-KMMI(data_list=data_intrim,nMI=NULL,covariates=c("trt"),plot=TRUE,data_orig=NULL)

##MI estimation
KMMI(data_list=data_intrim,nMI=1000,covariates=c("trt"),plot=TRUE,data_orig=NULL)

data_intrim2<-uc_data_transform(data=df_y, var_list=c("id","trt"),
                               var_list_new=NULL,time="time", prob="prob")

KMMI(data_list=data_intrim2,nMI=1000,covariates=c("trt"),plot=TRUE,data_orig=data_orig,
time_var=c("time"),event_var=c("cens"))
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ survival }% use one of  RShowDoc("KEYWORDS")
\keyword{ nonparametric}% __ONLY ONE__ keyword per line
