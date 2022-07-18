`%nin%` <- purrr::negate(`%in%`)

#' Harmonic mean
#'
#' @param a Vector of values
#'
#' @return The harmonic mean of the values: 1/mean(1/a).
#' @export
#'
#' @examples
mean.harmonic <- function(a){
  1/mean(1/a) #compute the harmonic mean
}

#' Aggregate influence curves to the cluster level
#'
#' Aggregates individual-level influence curve (IC) estimates to the cluster
#' level. This allows for appropriate statistical inference that respects the
#' cluster as the independent unit, whether the target parameter is an
#' individual-level or cluster-level effect.
#'
#' @param recordIC The individual-level IC estimates
#' @param id Vector of cluster membership IDs
#'
#' @return A vector of aggregated IC values with length equal to the number of
#'   unique cluster IDs.
#' @export
#'
#' @examples
aggregate_IC <- function(recordIC, id) {
  if (is.null(id)) return(recordIC)
  aggregatedIC <- as.matrix(stats::aggregate(recordIC, list(id=id), sum)[, -1, drop = FALSE])
  num.records <- nrow(recordIC)
  num.clusters <- nrow(aggregatedIC)
  aggregatedIC <- aggregatedIC * num.clusters / num.records
  return(aggregatedIC)
}



#' Calculate confidence intervals and p-values on relative or absolute scale
#'
#' @param goal String specifying the scale of the target parameter. Default is
#'   \code{RD}, risk/rate difference. Any other values assume that input values
#'   are given on the log scale, and the function will exponentiate the
#'   estimated target parameter and confidence interval bounds to output an
#'   artihmetic risk/rate ratio.
#' @param psi True value (if known) of the target parameter, for example, in a
#'   simulation study.
#' @param psi.hat Estimated value of the target parameter.
#' @param se Standard error of the estimated target parameter.
#' @param df Degrees of freedom for the Student's \emph{t} distribution as an
#'   approximation of the asymptotic normal distribution. If \code{df > 40}, the
#'   value is ignored and a normal distribution is used for inference.
#' @param sig.level Desired significance (alpha) level. Defaults to 0.05.
#' @param one.sided Logical indicating that a one-sided test is desired.
#'   Defaults to \code{FALSE}.
#' @param alt.smaller If one-sided test is desired, is the alternative
#'   hypothesis that the intervention arm will have a smaller value that the
#'   control arm? For example, if you expect a public health intervention to
#'   reduce the mean of a disease outcome, use alt.smaller = TRUE (this
#'   corresponds to a null hypothesis that the intervention did not reduce the
#'   mean disease outcome).
#'   
#' @return A one-row data frome with the estimated target parameter value
#'   (\code{est}), the (two-sided) confidence interval \code{CI.lo},
#'   \code{CI/hi}, the standard error of the estimate, the (possibly one-sided)
#'   p-value, and bias/coverage/rejection indicators (if true value or target
#'   parameter is supplied). NOTE: If \code{goal != "RD"}, the output standard
#'   error will be on the log scale.
#' @export
#'
#' @examples
get.inference <- function(goal = 'RD', psi = NA, psi.hat, se, df = 99, sig.level = 0.05, 
                          one.sided = F, alt.smaller = NULL){
  
  # test statistic
  # (on the log-transformed scale if goal is arithmetic RR or odds ratio)
  tstat <- psi.hat / se
  
  if(df > 40){
    # assume normal distribution
    cutoff <- stats::qnorm(sig.level/2, lower.tail=F)
    # one.sided hypothesis test 
    if(one.sided){
      pval <- stats::pnorm(tstat, lower.tail=alt.smaller) 
    } else{
      pval<- 2*stats::pnorm(abs(tstat), lower.tail=F) 
    }
  } else {
    # use Student's t-distribution
    # print('Using t-distribution')
    cutoff <- stats::qt(sig.level/2, df=df, lower.tail=F)
    # one.sided hypothesis test if specified
    if(one.sided){
      pval <- stats::pt(tstat, df=df, lower.tail= alt.smaller ) 
    } else{
      pval <- 2*stats::pt(abs(tstat), df=df, lower.tail=F)
    }
  }
  # 95% confidence interval 
  CI.lo <- (psi.hat - cutoff*se)
  CI.hi <- (psi.hat + cutoff*se)
  
  # If on relative (log) scale, transform back 
  if(goal != 'RD'){
    psi.hat <- exp(psi.hat)
    CI.lo <- exp(CI.lo)
    CI.hi <- exp(CI.hi)
  }  
  
  # bias
  bias <- (psi.hat - psi)
  # confidence interval coverage?
  cover <-  CI.lo <= psi & psi <= CI.hi
  # reject the null?
  reject <- as.numeric(pval < sig.level)
  return(data.frame(est = psi.hat, CI.lo, CI.hi, se = se, pval, bias, cover, reject))
}