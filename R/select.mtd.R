#'
#' Select the maximum tolerated dose (MTD) for single agent trials
#'
#' Select the maximum tolerated dose (MTD) when the single-agent trial is completed
#'
#' @usage select.mtd(target, npts, ntox, cutoff.eli=0.95, extrasafe=FALSE, offset=0.05,
#'                  boundMTD=FALSE,p.tox=1.4*target)
#'
#' @param target the target DLT rate
#' @param npts a vector containing the number of patients treated at each dose level
#' @param ntox a vector containing the number of patients who experienced dose-limiting
#'              toxicity at each dose level
#' @param cutoff.eli the cutoff to eliminate overly toxic doses for safety. We recommend
#'                    the default value of (\code{cutoff.eli=0.95}) for general use.
#' @param extrasafe set \code{extrasafe=TRUE} to impose a more strict stopping rule for
#'                   extra safety
#' @param offset a small positive number (between \code{0} and \code{0.5}) to control how strict the
#'                stopping rule is when \code{extrasafe=TRUE}. A larger value leads to
#'                a more strict stopping rule. The default value \code{offset=0.05}
#'                generally works well.
#' @param boundMTD set \code{boundMTD=TRUE} to impose the condition: the isotonic estimate of toxicity
#'                 probability for the selected MTD must be less than de-escalation boundary.
#' @param p.tox the lowest toxicity probability that is deemed overly toxic such
#'              that deescalation is required. The default value is
#'                \code{p.tox=1.4*target}.
#'
#' @details \code{select.mtd()} selects the MTD based on isotonic estimates of toxicity
#'          probabilities. \code{select.mtd()} selects as the MTD dose \eqn{j^*}, for which the
#'          isotonic estimate of the DLT rate is closest to the target. If there
#'          are ties, we select from the ties the highest dose level when the estimate
#'          of the DLT rate is smaller than the target, or the lowest dose level
#'          when the estimate of the DLT rate is greater than the target. The
#'          isotonic estimates are obtained by the pooled-adjacent-violators algorithm
#'          (PAVA) (Barlow, 1972).
#'
#' @return  \code{select.mtd()} returns (1) target toxicity probability (\code{$target}), (2) selected MTD (\code{$MTD}),
#' (3) isotonic estimate of the DLT probablity at each dose and associated \eqn{95\%} credible interval (\code{$p_est}),
#' and (4) the probability of overdosing defined as \eqn{Pr(toxicity>\code{target}|data)} (\code{$p_overdose})
#'
#'
#' @note  The MTD selection and dose escalation/deescalation rule are two independent
#'        components of the trial design. When appropriate, another dose selection
#'        procedure (e.g., based on a fitted logistic model) can be used to select
#'        the MTD after the completion of the trial using the BOIN design.
#'
#'
#' @author Suyu Liu, Yanhong Zhou, and Ying Yuan
#'
#' @references Liu S. and Yuan, Y. (2015). Bayesian Optimal Interval Designs for
#'            Phase I Clinical Trials, \emph{Journal of the Royal Statistical Society:
#'            Series C}, 64, 507-523.
#'
#'            Yan, F., Zhang, L., Zhou, Y., Pan, H., Liu, S. and Yuan, Y. (2020).BOIN: An R Package
#'            for Designing Single-Agent and Drug-Combination Dose-Finding Trials Using Bayesian Optimal
#'            Interval Designs. \emph{Journal of Statistical Software}, 94(13),1-32.<doi:10.18637/jss.v094.i13>.
#'
#'
#'             Yuan Y., Hess K.R., Hilsenbeck S.G. and Gilbert M.R. (2016). Bayesian Optimal Interval Design: A
#'        Simple and Well-performing Design for Phase I Oncology Trials, \emph{Clinical Cancer Research}, 22, 4291-4301.
#'
#'
#' @seealso Tutorial: \url{http://odin.mdacc.tmc.edu/~yyuan/Software/BOIN/BOIN2.6_tutorial.pdf}
#'
#' Paper: \url{http://odin.mdacc.tmc.edu/~yyuan/Software/BOIN/paper.pdf}
#'
#' @examples
#'
#' ### select the MTD for BOIN single agent trial
#' n <- c(3, 3, 15, 9, 0)
#' y <- c(0, 0, 4, 4, 0)
#' selmtd <- select.mtd(target=0.3, npts=n, ntox=y)
#' summary(selmtd)
#' plot(selmtd)
#'
#' @export
select.mtd <- function (target, npts, ntox, cutoff.eli = 0.95, extrasafe = FALSE,
                        offset = 0.05, boundMTD = FALSE, p.tox=1.4*target)
{
  pava <- function(x, wt = rep(1, length(x))) {
    n <- length(x)
    if (n <= 1)
      return(x)
    if (any(is.na(x)) || any(is.na(wt))) {
      stop("Missing values in 'x' or 'wt' not allowed")
    }
    lvlsets <- (1:n)
    repeat {
      viol <- (as.vector(diff(x)) < 0)
      if (!(any(viol)))
        break
      i <- min((1:(n - 1))[viol])
      lvl1 <- lvlsets[i]
      lvl2 <- lvlsets[i + 1]
      ilvl <- (lvlsets == lvl1 | lvlsets == lvl2)
      x[ilvl] <- sum(x[ilvl] * wt[ilvl])/sum(wt[ilvl])
      lvlsets[ilvl] <- lvl1
    }
    x
  }

  lambda_d = log((1 - target)/(1 - p.tox))/log(p.tox * (1 -target)/(target * (1 - p.tox)))

  y = ntox
  n = npts
  ndose = length(n)
  elimi = rep(0, ndose)
  for (i in 1:ndose) {
    if (n[i] >= 3) {
      if (1 - pbeta(target, y[i] + 1, n[i] - y[i] + 1) >
          cutoff.eli) {
        elimi[i:ndose] = 1
        break
      }
    }
  }
  if (extrasafe) {
    if (n[1] >= 3) {
      if (1 - pbeta(target, y[1] + 1, n[1] - y[1] + 1) >
          cutoff.eli - offset) {
        elimi[1:ndose] = 1
      }
    }
  }
  if (elimi[1] == 1 || sum(n[elimi == 0]) == 0) {
    selectdose = 99
  }
  else {
    adm.set = (n != 0) & (elimi == 0)
    adm.index = which(adm.set == T)
    y.adm = y[adm.set]
    n.adm = n[adm.set]
    phat = (y.adm + 0.05)/(n.adm + 0.1)
    phat.var = (y.adm + 0.05) * (n.adm - y.adm + 0.05)/((n.adm +
                                                           0.1)^2 * (n.adm + 0.1 + 1))
    phat = pava(phat, wt = 1/phat.var)
    phat = phat + (1:length(phat)) * 1e-10
    if(boundMTD){
      if(all(phat>lambda_d)){selectdose=99}else{
        phat=phat[phat<=lambda_d]
        selectd = sort(abs(phat - target), index.return = T)$ix[1]
        selectdose = adm.index[selectd]
      }

    }else{
      selectd = sort(abs(phat - target), index.return = T)$ix[1]
      selectdose = adm.index[selectd]
      }


  }

    trtd = (n != 0)
    poverdose = pava(1 - pbeta(target, y[trtd] + 0.05, n[trtd] -
                                 y[trtd] + 0.05))
    phat.all = pava((y[trtd] + 0.05)/(n[trtd] + 0.1), wt = 1/((y[trtd] +
                                                                 0.05) * (n[trtd] - y[trtd] + 0.05)/((n[trtd] + 0.1)^2 *
                                                                                                       (n[trtd] + 0.1 + 1))))
    lowerCIs=pava(qbeta(0.025, y[trtd] + 0.05,n[trtd] - y[trtd] + 0.05),wt = 1/((y[trtd] +
                                                                 0.05) * (n[trtd] - y[trtd] + 0.05)/((n[trtd] + 0.1)^2 *
                                                                                                       (n[trtd] + 0.1 + 1))))
    upperCIs=pava(qbeta(0.975, y[trtd] + 0.05,n[trtd] - y[trtd] + 0.05),wt = 1/((y[trtd] +
                                                                 0.05) * (n[trtd] - y[trtd] + 0.05)/((n[trtd] + 0.1)^2 *
                                                                                                       (n[trtd] + 0.1 + 1))))

    A1 = A2 = A3 = A4 = NULL
    k = 1
    for (i in 1:ndose) {
      if (n[i] > 0) {
        A1 = append(A1, formatC(phat.all[k], digits = 2,
                                format = "f"))
        A2 = append(A2, formatC(lowerCIs[i], digits = 2, format = "f"))
        A3 = append(A3, formatC(upperCIs[i], digits = 2, format = "f"))
        A4 = append(A4, formatC(poverdose[k], digits = 2,
                                format = "f"))
        k = k + 1
      }
      else {
        A1 = append(A1, "----")
        A2 = append(A2, "----")
        A3 = append(A3, "----")
        A4 = append(A4, "----")
      }
    }
    p_est = data.frame(cbind(dose = 1:length(npts), phat = A1,
                             CI = paste("(", A2, ",", A3, ")", sep = "")))
    out = list(target = target, MTD = selectdose, p_est = p_est,
               p_overdose = A4)

  class(out)<-"boin"
  return(out)
}
