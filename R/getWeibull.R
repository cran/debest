#' Estimate parameters of the Weibull survival distributions for each subgroup.
#' @param dat - data frame of three columns: 1) time: follow up time; 2) status: indicator, 0=no event, 1=event; 3) group, integer 1,...,n, for each subgroup of patients.
#' @return shape and scale parameters of the Weibull distribution for each subgroup.
#' @references Hong Zhang, Jie Pu, Shibing Deng, Satrajit Roychoudhury, Haitao Chu and Douglas Robinson. "Study Duration Prediction for Clinical Trials with Time-to-Event Endpoints Using Mixture Distributions Accounting for Heterogeneous Population", arXiv:2401.00540.
#' @examples
#' # dat_udca already has time, status and group columns defined,
#' getWeilbull(dat_udca)
#' @export
#' @import stats survival flexsurv
getWeilbull = function(dat){
  ngroup = max(dat$group)
  fit_weilbull = sapply(1:ngroup,
                        function(x)exp(flexsurvreg(Surv(time, status) ~ 1,
                                                   data=dat[dat$group==x,],
                                                   dist="weibull")$coef))
  return(list(shape = fit_weilbull[1,],
         scale = fit_weilbull[2,]))
}

