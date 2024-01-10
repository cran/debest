#' Estimate parameters of the distribution of a*Beta(1, b) .
#' @param x - vector of real positive numbers.
#' @return estimates of a and b. For a, it is the range of x. For b, two estimates are provided: the MLE, and the method of moment (matching mean).
#' @references Hong Zhang, Jie Pu, Shibing Deng, Satrajit Roychoudhury, Haitao Chu and Douglas Robinson. "Study Duration Prediction for Clinical Trials with Time-to-Event Endpoints Using Mixture Distributions Accounting for Heterogeneous Population", arXiv:2401.00540.
#' @examples
#' getBeta(as.numeric(dat_udca$entry.dt))
#' @export
#' @import stats
getBeta = function(x){
  a = max(x) - min(x)
  x_std = (x - min(x))/(a + 1e-10)
  # MLE
  mean_log = mean(log(1-x_std))
  out_MLE = -1/mean_log
  # Matching mean
  out_MoM = (1-mean(x_std))*(mean(x_std)*(1-mean(x_std))/var(x_std)-1)

  return(list(a=a, b_MLE=out_MLE, b_Mean=out_MoM))
}
