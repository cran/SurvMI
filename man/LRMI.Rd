\name{LRMI}
\alias{LRMI}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Log-rank test with events uncertainty
%%  ~~function to do ... ~~
}
\description{
This function conducts the Log-rank test with respect to uncertain endpoints, by MI or weighted
method.
}
\usage{
LRMI(data_list, nMI, covariates, strata = NULL,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data_list}{
The data list which has been transformed from the long format by uc_data_transform function.
}
  \item{nMI}{
Number of imputation (>1). If missing, weighted statistics would be output instead.
}
  \item{covariates}{
The categorical variable used in the Log-rank test. No need to factorlize numeric variables.
}
  \item{strata}{
Strata variable may required by the Log-rank test
}
  \item{\dots}{
Other arguments passed on to survdiff().
}
}

\value{
\item{est}{Estimated LR statistics, either from the MI method or weighted method}
\item{var}{Estimated variance matrix}
\item{est_mat}{Matrix containing estimate of statistics from each imputed dataset}
\item{Var_mat}{Array containing variances for each imputed dataset}
\item{Between Var}{Between imputation variance}
\item{Within Var}{Mean within imputed dataset variance}
\item{nMI}{Number of imputed datasets}
\item{pvalue}{Estimated two-sided Chi-square test p-value}
\item{df}{Degree of freedom}
\item{covariates}{covariates}
\item{ngroup}{Number of groups}
\item{obsmean}{Mean of observed events count across imputations}
\item{expmean}{Mean of expected events count across imputations}
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
\code{\link{uc_data_transform}}, \code{\link{LRMI.summ}}
}
\examples{
df_x<-data_sim(n=500,0.8,haz_c=0.5/365)
data_intrim<-uc_data_transform(data=df_x,
                               var_list=c("id_long","trt_long"),
                               var_list_new=c("id","trt"),
                               time="time_long",
                               prob="prob_long")

#nMI=10 used in the example below to reduce the time needed
#but a large number as nMI=1000 is recommended in practice
fit<-LRMI(data_list=data_intrim,nMI=10,covariates=c("trt"),strata=NULL)
LRMI.summ(fit)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ survival }% use one of  RShowDoc("KEYWORDS")

