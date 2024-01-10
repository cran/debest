#' @importFrom stats median pgamma quantile rbeta rexp rweibull var


# Formula (3) from Zhang et al. (2024). arXiv:2401.00540
# F_list is a list of CDFs that are corresponding to each subgroup of patients.
F_mix = function(t, proportion, F_list){
  sapply(t, function(x)sum(sapply(F_list, mapply, x) * proportion))
}

# Formula (7) from Zhang et al. (2024). arXiv:2401.00540
Fkl = function(t, lamV, lamW, a, beta){
  lamSum = lamV+lamW
  lamV/lamSum*(1-(max_vec(0,1-t/a))^beta-exp(-lamSum*(t-a))*beta*gamma(beta)/(a*lamSum)^beta*(pgamma(a,shape=beta,rate=lamSum)-pgamma(max_vec(0,a-t),shape=beta,rate=lamSum)))
}

# Equation below formula (2) from Machida et al. (2021). doi:10.1002/sim.8911
F_machida = function(t, lamV, lamW, a){
  lamSum = lamV+lamW
  lamV/lamSum * (1 - exp(-lamSum*t)/
                   (a*lamSum)*(exp(lamSum*a) - 1))
}

# Formula (8) from Zhang et al. (2024). arXiv:2401.00540
Fkl_unif = function(t, lamV, lamW, a){
  lamSum = lamV+lamW
  lamV/lamSum*(1-max_vec(0,1-t/a)-exp(-lamSum*t)/(a*lamSum)*(exp(lamSum*min_vec(a, t))-1))
}

# Formula (9) from Zhang et al. (2024). arXiv:2401.00540
FA = function(t, a, r1, lam1, lam2){
  r2 = 1-r1
  part1 = max_vec(0, 1-t/a)
  part2 = r1*exp(-lam1*t)/(a*lam1)*(exp(lam1*min_vec(a, t))-1)
  part3 = r2*exp(-lam2*t)/(a*lam2)*(exp(lam2*min_vec(a, t))-1)
  1 - part1 - part2 - part3
}

# Formula (10) from Zhang et al. (2024). arXiv:2401.00540
FE = function(t, a, r1, lam1){
  1 - max_vec(0, 1-r1*t/a) - r1*exp(-lam1*t)/(a*lam1)*(exp(lam1*min_vec(a/r1, t))-1)
}

# Helper function for the main calcDuration function
rFkl_weibull = function(nsim=1e4, shapeV, scaleV, lamW, a, beta){
  Z = rep(Inf, nsim);
  V = rweibull(nsim, shape=shapeV, scale=scaleV);
  if(lamW>0){
    W = rexp(nsim, rate=lamW);
    id_event = which(V<W)
  }else{
    id_event = 1:nsim
  }
  U = a*rbeta(nsim,shape1=1, shape2=beta);
  Z[id_event] = V[id_event] + U[id_event]
  Z
}



max_vec = function(lower, vec){
  # return max(lower, vec[i]) for every i
  sapply(vec, function(x)max(lower, x))
}

min_vec = function(upper, vec){
  # return min(upper, vec[i]) for every i
  sapply(vec, function(x)min(upper, x))
}

# Fkl_exp_sim = function(t, lamV, lamW, a, beta, nsim=1e4){
#   Z = rFkl_exp(nsim=nsim, lamV=lamV, lamW=lamW, a=a, beta=beta)
#   sapply(t, function(x)mean(Z<x))
# }
#
# Fkl_weibull_sim = function(t, shapeV, scaleV, lamW, a, beta, nsim=1e4){
#   Z = rFkl_weibull(nsim=nsim, shapeV=shapeV, scaleV=scaleV, lamW=lamW, a=a, beta=beta)
#   sapply(t, function(x)mean(Z<x))
# }


# rFkl_exp = function(nsim=1e4, lamV, lamW, a, beta){
#   Z = rep(Inf, nsim);
#   V = rexp(nsim, rate=lamV);
#   if(lamW>0){
#     W = rexp(nsim, rate=lamW);
#     id_event = which(V<W)
#   }else{
#     id_event = 1:nsim
#   }
#   U = a*rbeta(nsim,shape1=1, shape2=beta);
#   Z[id_event] = V[id_event] + U[id_event]
#   Z
# }
