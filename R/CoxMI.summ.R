CoxMI.summ <-
  function(x, digits=3)
  {
    if (!is.null(cl<- x$call)) {
      cat("Call:\n")
      dput(cl)
      cat("\n")
    }

    coef <- x$est
    se <- sqrt(diag(x$var))
    pvalue <-x$pvalue
    hr<-exp(x$est)
    within_se<-sqrt(diag(x$within_var))
    between_se<-sqrt(diag(x$between_var))
    column_name<-x$column_name
    tmp <- cbind(coef, hr,se, within_se, between_se,coef/se,pvalue)
    dimnames(tmp) <- list(column_name, c("coef","HR" ,"se(coef)","within","between", "z", "p"))
    cat("\n")
    cat("MI Cox model estimation with ", x$nMI, "imputations",  sep=" ")
    cat("\n")

    prmatrix(tmp)
    invisible(x)
    s1=paste("n= ", x$n, sep=" ")
    s2=paste("number of expected events=", x$en, sep=" ")
    cat(s1,s2,sep=", ")
    cat("\n")
  }
