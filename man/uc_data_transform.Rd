\name{uc_data_transform}
\alias{uc_data_transform}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Transform long formatted time-to-event data into a data list
}
\description{
This function transforms data from long format (one record per event) to a datalist with length as unique subject number. The transformation is required before fitting other models from the package.
}
\usage{
uc_data_transform(data,var_list,var_list_new,time,prob)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
The dataset in long format with a row for each potential event. For ceonsoring record, the event prob should be 0. It should include id, time and prob variables at a minimum. If any covariates are included in the call to the function, then these variables should also be included. A censoring record is required for each subject. Categorical variables need to be encoded as factor varaible before transformationif they are expected to be in the Cox model.
}
  \item{var_list}{
The list of identification variables, such as: c("id_long","trt_long").
}
  \item{time}{
The time variable need to be transofirmed, e.g. time_long.
}
  \item{prob}{
The prob variable need to be transformed, e.g. prob_long.
}
  \item{var_list_new}{
The character vector contains the new names for the id variables defined in the var_list, if missing, previous variable names would be used.
}
}

\value{
\item{time}{The list of all potential event time }
\item{prob}{The list of all potential event probabilities}
\item{weights}{The list of all potential event weights }
\item{e}{The list of individual potential event count}
\item{s}{The list of all survival probabilities }
\item{data_uc}{The dataset contains unique information of each subject}
\item{data_long}{The dataset contains the original data in long format}
}

\author{
Yiming Chen
}
\examples{
df_x<-data_sim(n=1000,true_hr=0.8,haz_c=0.5/365)
df_x$f.trt<-as.factor(df_x$trt_long)
data_intrim<-uc_data_transform(data=df_x,
                               var_list=c("id_long","f.trt"),
                               var_list_new=c("id","trt"),
                               time="time_long",
                               prob="prob_long")

}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{manip}% use one of  RShowDoc("KEYWORDS")
\keyword{list}
