\name{CoxMI}
\alias{CoxMI}
\title{Cox PH model with MI method
}
\description{
CoxMI function estimated Cox model with uncertain endpoints by using MI method. Users have to provide survival data in a long format with rows for all potential events, together with corresponding event probabilities. The long format data should be transformed by the uc_data_transform function into a data list before feed into the function.
}
\usage{
CoxMI(data_list,nMI=1000,covariates=NULL,id=NULL,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data_list}{
 The data list which has been transformed from the long format by the uc_data_transform function.
}
  \item{nMI}{
Number of imputations (>1).
}
  \item{covariates}{
Vector of covariates on the RHS of Cox model. Categorical variables need to be encoded as factor variables before entering the model. This encoding has to be done before the data transform step.
}
  \item{id}{
Vector of id variable if Andersen-Gill model is required.
}
  \item{\dots}{
Other arguments passed on to coxph().
}
}
\details{
Calculates the estimated parameters as in the usual Cox proportional hazards model when event uncertainties present. The data are assumed to consist of potential event times with probabilities or weights between 0 and 1 corresponding to the probability that an event occurred at each time.

}
\value{

\item{est}{Estimated vector of coefficients in the model}
\item{var}{Estimated variance of the coefficients}
\item{betamat}{Matrix containing estimate of coefficient from each imputed dataset}
\item{Var_mat}{Array containing variances for each imputed dataset}
\item{Between Var}{Between imputation variance}
\item{Within Var}{Mean within imputed dataset variance}
\item{nMI}{Number of imputed datasets}
\item{pvalue}{Estimated two-sided p-value}
\item{en}{Expected events count - mean event count of imputed datasets}
}
\references{
[1] Rubin DB. Multiple Imputation for Nonresponse in Surveys. New York: Wiley; 1987
}
\author{
Yiming Chen, John Lawrence
}

\seealso{
  \code{\link{Coxwt}}, \code{\link{CoxMI.summ}}.
}
\examples{
set.seed(128)
df_x<-data_sim(n=500,true_hr=0.8,haz_c=0.5/365)
df_x$f.trt<-as.factor(df_x$trt_long)
data_intrim<-uc_data_transform(data=df_x,
                               var_list=c("id_long","f.trt"),
                               var_list_new=c("id","trt"),
                               time="time_long",
                               prob="prob_long")
#nMI=10 used in the example below to reduce the time needed
#but a large number as nMI=1000 is recommended in practice

fit<-CoxMI(data_list=data_intrim,nMI=10,covariates=c("trt"))
CoxMI.summ(fit)

\donttest{
fit<-CoxMI(data_list=data_intrim,nMI=1000,covariates=c("trt"),id=c("id"))
CoxMI.summ(fit)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{survival }% use one of  RShowDoc("KEYWORDS")
