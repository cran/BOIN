#'
#' Determine the dose combination for the next cohort of new patients for drug-combination trials that aim to find a MTD
#'
#' Determine the dose combination for the next cohort of new patients for drug-combination trials that aim to find a MTD
#'
#' @usage next.comb(target, npts, ntox, dose.curr, n.earlystop=100,
#'                  p.saf=0.6*target, p.tox=1.4*target, cutoff.eli=0.95,
#'                  extrasafe=FALSE, offset=0.05)
#'
#' @param target the target DLT rate
#' @param npts a \code{J*K} matrix \code{(J<=K)} containing the number of patients treated at each dose combination
#' @param ntox a \code{J*K} matrix \code{(J<=K)} containing the number of patients experienced
#'             dose-limiting toxicity at each dose combination
#' @param dose.curr the current dose combination
#' @param n.earlystop the early stopping parameter. If the number of patients
#'                    treated at the current dose reaches \code{n.earlystop},
#'                    stop the trial and select the MTD based on the observed data.
#'                    The default value \code{n.earlystop=100} essentially turns
#'                    off this type of early stopping.
#' @param p.saf the highest toxicity probability that is deemed subtherapeutic
#'              (i.e. below the MTD) such that dose escalation should be undertaken.
#'              The default value is \code{p.saf=0.6*target}.
#' @param p.tox the lowest toxicity probability that is deemed overly toxic such
#'              that deescalation is required. The default value is
#'              \code{p.tox=1.4*target}.
#' @param cutoff.eli the cutoff to eliminate an overly toxic dose for safety.
#'                   We recommend the default value of (\code{cutoff.eli=0.95})
#'                   for general use.
#' @param extrasafe set \code{extrasafe=TRUE} to impose a more stringent stopping rule
#' @param offset a small positive number (between \code{0} and \code{0.5}) to control how strict the
#'               stopping rule is when \code{extrasafe=TRUE}. A larger value leads to a more
#'               strict stopping rule. The default value \code{offset=0.05} generally works well.
#'
#' @details This function is used to determine dose combination for conducting combination trials.
#'          Given the currently observed data, \code{next.comb()} determines dose combination for
#'          treating the next cohort of new patients. The currently observed data include: the
#'          number of patients treated at each dose combination (i.e., \code{npts}),
#'          the number of patients who experienced dose-limiting toxicities at each dose
#'          combination (i.e., \code{ntox}), and the level of current dose (i.e., \code{dose.curr}).
#'
#' @return the recommended dose for treating the next cohort of patients (\code{$next_dc}).
#'
#' @author Suyu Liu and Ying Yuan
#'
#' @references Liu S. and Yuan, Y. (2015). Bayesian Optimal Interval Designs for Phase I Clinical
#'             Trials, \emph{Journal of the Royal Statistical Society: Series C}, 64, 507-523.
#'
#'            Lin R. and Yin, G. (2017). Bayesian Optimal Interval Designs for Dose Finding in
#'            Drug-combination Trials, \emph{Statistical Methods in Medical Research}, 26, 2155-2167.
#'
#'           Yan, F., Zhang, L., Zhou, Y., Pan, H., Liu, S. and Yuan, Y. (2020).BOIN: An R Package
#'            for Designing Single-Agent and Drug-Combination Dose-Finding Trials Using Bayesian Optimal
#'            Interval Designs. \emph{Journal of Statistical Software}, 94(13),1-32.<doi:10.18637/jss.v094.i13>.
#'
#'
#'
#' @seealso  Tutorial: \url{http://odin.mdacc.tmc.edu/~yyuan/Software/BOIN/BOIN2.6_tutorial.pdf}
#'
#'           Paper: \url{http://odin.mdacc.tmc.edu/~yyuan/Software/BOIN/paper.pdf}
#'
#' @examples
#'
#' ## determine the dose combination for the next cohort of new patients
#' n <- matrix(c(3, 0, 0, 0, 0, 7, 6, 0, 0, 0, 0, 0, 0, 0, 0), ncol=5, byrow=TRUE)
#' y <- matrix(c(0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0), ncol=5, byrow=TRUE)
#' nxt.comb <- next.comb(target=0.3, npts=n, ntox=y, dose.curr=c(2, 2))
#' summary(nxt.comb)
#'
#' @export
next.comb <- function (target, npts, ntox, dose.curr, n.earlystop = 100, p.saf = 0.6 *
                         target, p.tox = 1.4 * target, cutoff.eli = 0.95, extrasafe = FALSE,
                       offset = 0.05)
{
  if (npts[dose.curr[1], dose.curr[2]] == 0) {
    stop("dose entered is not the current dose")

  }
  if (target < 0.05) {
   stop("the target is too low! ")

  }
  if (target > 0.6) {
   stop("the target is too high! ")

  }
  if ((target - p.saf) < (0.1 * target)) {
   stop("the probability deemed safe cannot be higher than or too close to the target! ")

  }
  if ((p.tox - target) < (0.1 * target)) {
   stop("the probability deemed toxic cannot be lower than or too close to the target! ")

  }
  if (offset >= 0.5) {
   stop("the offset is too large! ")

  }
  if (n.earlystop <= 6) {
    warning("the value of n.earlystop is too low to ensure good operating characteristics. Recommend n.earlystop = 9 to 18. ")

  }
  temp = get.boundary(target, ncohort = 150, cohortsize = 1,
                      n.earlystop, p.saf, p.tox, cutoff.eli, extrasafe, offset)$boundary_tab
  b.e = temp[2, ]
  b.d = temp[3, ]
  b.elim = temp[4, ]
  lambda1 = log((1 - p.saf)/(1 - target))/log(target * (1 -
                                                          p.saf)/(p.saf * (1 - target)))
  lambda2 = log((1 - target)/(1 - p.tox))/log(p.tox * (1 -
                                                         target)/(target * (1 - p.tox)))
  n = npts
  y = ntox
  earlystop = 0
  d = dose.curr
  nc = n[d[1], d[2]]
  ndose = length(npts)
  elimi = matrix(rep(0, ndose), dim(n)[1], dim(n)[2])
  if (n[d[1], d[2]] >= n.earlystop) {
    cat("Terminate the trial because the number of patients treated at (",
        d[1], ", ", d[2], ") has reached ", n.earlystop, ".")
    d = c(99, 99)
    earlystop = 1
  }
  if (!is.na(b.elim[nc])) {
    if (d[1] == 1 && d[2] == 1 && y[d[1], d[2]] >= b.elim[nc]) {
      d = c(99, 99)
      earlystop = 1
      cat("Terminate the trial because the lowest dose is overly toxic ")
    }
    if (extrasafe) {
      if (d[1] == 1 && d[2] == 1 && n[1, 1] >= 3) {
        if (1 - pbeta(target, y[1, 1] + 1, n[1, 1] -
                      y[1, 1] + 1) > cutoff.eli - offset) {
          d = c(99, 99)
          earlystop = 1
          cat("Terminate the trial because the lowest dose is overly toxic ")
        }
      }
    }
  }
  for (i in 1:dim(n)[1]) {
    for (j in 1:dim(n)[2]) {
      if (n[i, j] > 0 && (!is.na(b.elim[n[i, j]]))) {
        if (y[i, j] >= b.elim[n[i, j]]) {
          elimi[i:dim(n)[1], j:dim(n)[2]] = 1
        }
      }
    }
  }
  out = list(next_dc = c(NA, NA))
  if (earlystop == 0) {
    if (y[d[1], d[2]] <= b.e[nc]) {
      n.temp=n; n.temp[n.temp==0]=1
      phat.mat=y/n.temp #p.hat for doses

      elevel = matrix(c(1, 0, 0, 1), 2)
      pr_H0 = rep(0, length(elevel)/2)
      nn = pr_H0
      for (i in seq(1, length(elevel)/2, by = 1)) {
        if (d[1] + elevel[1, i] <= dim(n)[1] && d[2] +
            elevel[2, i] <= dim(n)[2]) {
          if (elimi[d[1] + elevel[1, i], d[2] + elevel[2,i]] == 0) {
            if(i==1 & d[1]+1<=dim(n)[1]){
              if(any(phat.mat[d[1]+1,1:d[2]]>=lambda2)){
                pr_H0[i]=0
              }else{
                yn = y[d[1] + elevel[1, i], d[2] + elevel[2,i]]
                nn[i] = n[d[1] + elevel[1, i], d[2] + elevel[2,i]]
                pr_H0[i] <- pbeta(lambda2, yn + 0.5, nn[i] -yn + 0.5) -
                  pbeta(lambda1, yn + 0.5, nn[i] -yn + 0.5)
              }

            }
            if(i==2 & d[2]+1<=dim(n)[2]){
              if(any(phat.mat[1:d[1],d[2]+1]>=lambda2)){
                pr_H0[i]=0
              }else{
                yn = y[d[1] + elevel[1, i], d[2] + elevel[2,i]]
                nn[i] = n[d[1] + elevel[1, i], d[2] + elevel[2,i]]
                pr_H0[i] <- pbeta(lambda2, yn + 0.5, nn[i] -yn + 0.5) -
                  pbeta(lambda1, yn + 0.5, nn[i] -yn + 0.5)
              }
            }

          }
        }
      }

      pr_H0 = pr_H0 + nn * 5e-04
      if (max(pr_H0) == 0) {
        d = d
      }else{

        k = which(pr_H0 == max(pr_H0))[as.integer(runif(1) *length(which(pr_H0 == max(pr_H0))) + 1)]
        d = d + c(elevel[1, k], elevel[2, k])
      }

    }else if (y[d[1], d[2]] >= b.d[nc]) {
      delevel = matrix(c(-1, 0, 0, -1), 2)
      pr_H0 = rep(0, length(delevel)/2)
      nn = pr_H0
      for (i in seq(1, length(delevel)/2, by = 1)) {
        if (d[1] + delevel[1, i] > 0 && d[2] + delevel[2,
                                                       i] > 0) {
          yn = y[d[1] + delevel[1, i], d[2] + delevel[2,
                                                      i]]
          nn[i] = n[d[1] + delevel[1, i], d[2] + delevel[2,
                                                         i]]
          pr_H0[i] = pbeta(lambda2, yn + 0.5, nn[i] -
                             yn + 0.5) - pbeta(lambda1, yn + 0.5, nn[i] -
                                                 yn + 0.5)
        }
      }
      pr_H0 = pr_H0 + nn * 5e-04
      if (max(pr_H0) == 0) {
        d = d
      }
      else {
        k = which(pr_H0 == max(pr_H0))[as.integer(runif(1) *
                                                    length(which(pr_H0 == max(pr_H0))) + 1)]
        d = d + c(delevel[1, k], delevel[2, k])
      }
    }
    else {
      d = d
    }
    out = list(next_dc = d)
  }
  class(out)<-"boin"
  return(out)
}
