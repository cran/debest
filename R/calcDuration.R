#' Calculate the study duration based on Weibull distributions.
#' @param d - scalar, target number of events.
#' @param n - scalar, sample size.
#' @param proportion - vector of percentages of the subgroup.
#' @param SHAPEV - vector of shape parameters of Weibull survival distributions of the subgroups.
#' @param SCALEV - vector of scale parameters of Weibull survival distributions of the subgroups.
#' @param LAMW - vector of exponential drop-out distribution parameters of the subgroups.
#' @param A - vector of enrollment durations of the subgroups.
#' @param BETA - vector of beta distribution parameters of the subgroups.
#' @param conf.level - scalar, confidence level, default 0.9.
#' @param nsim - scalar, number of repetitions, default 1e4.
#' @return study duration estimate, d_med, and the confidence interval (d_lower, d_upper), as well as all the realizations, Z_d, of the study duration from the simulation.
#' @references Hong Zhang, Jie Pu, Shibing Deng, Satrajit Roychoudhury, Haitao Chu and Douglas Robinson. "Study Duration Prediction for Clinical Trials with Time-to-Event Endpoints Using Mixture Distributions Accounting for Heterogeneous Population", arXiv:2401.00540.
#' @examples
#' res_weibull = getWeilbull(dat_udca)
#' res_beta = getBeta(as.numeric(dat_udca$entry.dt))
#' prop = c(table(dat_udca$group)/length(dat_udca$group))
#' SHAPEV = res_weibull$shape
#' SCALEV = res_weibull$scale
#' LAMW = rep(-log(1 - 0.1)/6, 4)
#' A = rep(res_beta$a/30.416, 4) # convert days to months
#' BETA = rep(res_beta$b_Mean, 4)
#' myres1 = calcDuration(d=50, n=169, proportion=prop, SHAPEV, SCALEV, LAMW=LAMW, A, BETA)
#' c(myres1$d_lower, myres1$d_med, myres1$d_upper)
#' # drop-out will make the target number of events not achievable
#' myres2 = calcDuration(d=80, n=169, proportion=prop, SHAPEV, SCALEV, LAMW=LAMW, A, BETA)
#' c(myres2$d_lower, myres2$d_med, myres2$d_upper)
#' # If there is no drop-out
#' myres3 = calcDuration(d=80, n=169, proportion=prop, SHAPEV, SCALEV, LAMW=rep(0, 4), A, BETA)
#' c(myres3$d_lower, myres3$d_med, myres3$d_upper)
#' @export
#' @import stats
calcDuration = function(d, n, proportion, SHAPEV, SCALEV, LAMW, A, BETA, conf.level = 0.9, nsim = 1e4){
  rF_list = vector(mode='list', length=length(proportion))
  n_group = round(n*proportion)
  hi = cumsum(n_group)
  lo = c(1, hi[1:(length(hi)-1)]+1)
  for(i in 1:length(rF_list)){
    rF_list[[i]] = function(x)rFkl_weibull(nsim=x, shapeV=SHAPEV[i], scaleV=SCALEV[i],
                                           lamW=LAMW[i], a=A[i], beta=BETA[i])
  }
  Z_d = rep(NA, nsim)
  for(j in 1:nsim){
    Z_out = rep(NA, sum(n_group))
    for(i in 1:length(hi)){
      Z_out[lo[i]:hi[i]] = c(rF_list[[i]](n_group[i]))
    }
    Z_d[j] = sort(Z_out)[d]
  }
  d_med = median(Z_d)
  d_upper = quantile(Z_d, probs = 0.5 + conf.level/2)
  d_lower = quantile(Z_d, probs = 0.5 - conf.level/2)
  return(list(Z_d = Z_d,
              d_med= d_med,
              d_upper = d_upper,
              d_lower = d_lower))
}

