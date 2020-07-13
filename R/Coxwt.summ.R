
Coxwt.summ <-function(x, digits=3)
{
  if (!is.null(cl<- x$call)) {
    cat("Call:\n")
    dput(cl)
    cat("\n")
  }

  if (is.null( x$nBS)) {

  coef <- x$coefficients
  se <- sqrt(x$var)
  pvalue <-x$pvalue
  hr<-x$hr
  z<-x$z
  column_name<-x$column_name
  tmp <- cbind(coef, hr,se, z,pvalue)
  dimnames(tmp) <- list(column_name, c("coef","HR" ,"se(coef)", "z", "p"))
  cat("\n")
  cat("Weighted Cox model estimation")
  cat("\n")
  prmatrix(tmp)
  invisible(x)
  }
  if (!is.null( x$nBS)) {


    coef <- x$coefficients
    coef_bs <- x$coefficients_bs
    se <- sqrt(x$var)
    se_bs <- sqrt(x$var_bs)
    pvalue <-x$pvalue
    hr<-x$hr
    z<-x$z
    column_name<-x$column_name
    tmp <- cbind(coef, hr,se,coef_bs,se_bs, z,pvalue)
    dimnames(tmp) <- list(column_name, c("Coef","HR" ,"se(Coef)", "Bootstrapped coef",
                                         "Bootstrapped se(Coef)","z", "p"))
    cat("\n")
    cat("Weighted Cox model estimation, with number of Bootstrap =", x$nBS,  sep=" ")
    cat("\n")
    prmatrix(tmp)
    invisible(x)
  }
}
