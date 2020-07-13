\name{Coxwt}
\alias{Coxwt}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{Weighted Cox PH model estimation
}
\description{
Estimate the Cox PH model by weighted partial likelihood. Event weights are calcualted with respect to event probabilities.
}
\usage{
Coxwt(data_list,covariates,init=NULL,BS=FALSE,nBS=1000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data_list}{
 The data list which has been transformed from the long format by the uc_data_transform function.
}
  \item{covariates}{
The vector of varaible on the RHS of the Cox model.
}
  \item{init}{
The initial value of covariates vector in the likelihood, length matches the length of covariates.
}
  \item{BS}{
T/F, whether conduct estimation via the Bootstrap method.
}
  \item{nBS}{
Number of BS, only effective if BS=TRUE.
}
}

\value{
\item{coefficients}{Estimated vector of coefficients in the model}
\item{var}{Estimated variance of the coefficients}
\item{hr}{Estimated hazard ratios in the model}
\item{z}{Wald test statistics}
\item{pvalue}{Estimated two-sided p-value}
\item{coefficients_bs}{Bootstrapped coefficient estimation}
\item{var_bs}{Bootstrapped variance estimation}
\item{column_name}{Column name}
}
\references{
[1]Cook TD. Adjusting survival analysis for the presence of unadjudicated study events. Controlled clinical trials. 2000;21(3):208-222.

[2]Cook TD, Kosorok MR. Analysis of time-to-event data with incomplete event adjudication. Journal of the american statistical association. 2004;99(468):1140-1152.

[3]Snapinn SM. Survival analysis with uncertain endpoints. Biometrics. 1998;54(1):209-218.
}

\author{
Yiming Chen, John Lawrence
}

\seealso{
  \code{\link{CoxMI}}, \code{\link{Coxwt.summ}}.
}
\examples{
df_x<-data_sim(n=500,0.8,haz_c=0.5/365)
data_intrim<-uc_data_transform(data=df_x,
                               var_list=c("id_long","trt_long"),
                               var_list_new=c("id","trt"),
                               time="time_long",
                               prob="prob_long")
 fit<-Coxwt(data_list=data_intrim,covariates=c("trt"),init=c(1),BS=FALSE)
 Coxwt.summ(fit)

##an example if we would like to check the BS variance
\donttest{
fit2<-Coxwt(data_list=data_intrim,covariates=c("trt"),init=c(1),BS=TRUE, nBS = 100)
Coxwt.summ(fit2)
}
}
\keyword{ survival }% use one of  RShowDoc("KEYWORDS")
