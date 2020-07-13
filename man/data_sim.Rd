\name{data_sim}
\alias{data_sim}
\alias{data_sim2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Simulated survival data with uncertain endpoints from exponential distribution.
}
\description{
data_sim function simulates data from a hypothetic 1:1 two-arms clinical trial, with one year uniform accrual period and three years follow-up.

data_sim2 function simplifies data list generated from above function to a more events only case. Note this function is only used for demonstration purpose.
}
\usage{
data_sim(n=200,true_hr=0.8,haz_c=1/365)
data_sim2(data_list,covariates,percentage)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{Total number of subject.
}
  \item{true_hr}{True hazard ratio between trt and control.
}
  \item{haz_c}{True event rate in the control arm.
}
  \item{data_list}{
The data list which has been transformed from the long format by uc_data_transform function.
}

  \item{covariates}{
The covariate we pose the true HR.
}

  \item{percentage}{
The percentage of censored subjects with potential events we would like to ultilize in the analysis. Ideally, with more potential events added, more power gain of imputation.
}
}

\value{
Dataframe. Simulated datasets with event probabilities and potential event date.
}

\author{
Yiming Chen, John Lawrence
}

\examples{
df_x<-data_sim(n=500,true_hr=0.8,haz_c=1/365)
data_intrim<-uc_data_transform(data=df_x,
                               var_list=c("id_long","trt_long"),
                               var_list_new=c("id","trt"),
                               time="time_long",
                               prob="prob_long")
df_y<-data_sim2(data_list=data_intrim,covariates=c("trt"),percentage=0.2)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{datagen}% use one of  RShowDoc("KEYWORDS")
\keyword{manip}% __ONLY ONE__ keyword per line
