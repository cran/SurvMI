LRMI.summ <-
  function(x, digits=3)
  {
    if (!is.null(cl<- x$call)) {
      cat("Call:\n")
      dput(cl)
      cat("\n")
    }

    if (!is.null(x$nMI)){
      tmp <- cbind(x$ngroup, x$obsmean,x$expmean)

      dimnames(tmp) <- list(x$cate_level,c("N","Observed" ,"Expected"))

      cat("\n")
      cat("MI Log-rank test with ", x$nMI, "imputations",  sep=" ")
      cat("\n")

      prmatrix(tmp)
      invisible(x)
    }

    if (is.null(x$nMI)){

      cat("\n")
      cat("Weighted Log-rank test")
      cat("\n")
    }
    s1=paste("Chisq=", format(round(x$chisq, 3)), sep=" ")
    s2=paste(paste("on ", x$df, sep=" "),"degrees of freedom,")
    s3=paste("p= ", format(round(x$pvalue, 3)), sep=" ")
    cat(s1,s2,s3,sep=" ")
    cat("\n")
  }
