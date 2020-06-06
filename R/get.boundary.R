#'
#' Generate the optimal dose escalation and deescalation boundaries for conducting the trial.
#'
#' Use this function to generate the optimal dose escalation and deescalation boundaries for conducting the trial.
#'
#'
#' @param target the target DLT rate
#' @param ncohort the total number of cohorts
#' @param cohortsize the cohort size
#' @param n.earlystop the early stopping parameter. If the number of patients treated at
#'                    the current dose reaches \code{n.earlystop}, stop the trial
#'                    and select the MTD based on the observed data. The default
#'                    value \code{n.earlystop=100} essentially turns off the type
#'                    of early stopping.
#' @param p.saf the highest toxicity probability that is deemed subtherapeutic
#'              (i.e., below the MTD) such that dose escalation should be made.
#'              The default value is \code{p.saf = 0.6 * target}.
#' @param p.tox the lowest toxicity probability that is deemed overly toxic such
#'              that deescalation is required. The default value is
#'              \code{p.tox=1.4*target}.
#' @param cutoff.eli the cutoff to eliminate an overly toxic dose for safety.
#'                   We recommend the default value (\code{cutoff.eli=0.95}) for general use.
#' @param extrasafe set \code{extrasafe=TRUE} to impose a more strict stopping rule for extra safety,
#'               expressed as the stopping boundary value in the result .
#' @param offset a small positive number (between 0 and 0.5) to control how strict
#'               the stopping rule is when \code{extrasafe=TRUE}. A larger value leads
#'               to a more strict stopping rule. The default value
#'               (\code{offset=0.05}) generally works well.
#'
#' @details The dose escalation and deescalation boundaries are all we need to run a
#'          phase I trial when using the BOIN design. The decision of which dose to
#'          administer to the next cohort of patients does not require complicated
#'          computations, but only a simple comparison of the observed DLT rate
#'          at the current dose with the dose escalation and deescalation boundaries.
#'          If the observed DLT rate at the current dose is smaller than or equal
#'          to the escalation boundary, we escalate the dose; if the observed toxicity
#'          rate at the current dose is greater than or equal to the deescalation boundary,
#'          we deescalate the dose; otherwise, we retain the current dose. The dose
#'          escalation and deescalation boundaries are chosen to minimize the probability
#'          of assigning patients to subtherapeutic or overly toxic doses, thereby
#'          optimizing patient ethics. \code{get.boundary()} also outputs the elimination
#'          boundary, which is used to avoid treating patients at overly toxic doses based
#'          on the following Bayesian safety rule: if \eqn{Pr(p_j > \phi | m_j , n_j ) > 0.95} and
#'          \eqn{n_j \ge 3}, dose levels \eqn{j} and higher are eliminated from the trial, where \eqn{p_j} is
#'          the toxicity probability of dose level \eqn{j}, \eqn{\phi} is the target DLT rate,
#'          and \eqn{m_j} and \eqn{n_j} are the number of toxicities and patients treated at dose level \eqn{j}.
#'          The trial is terminated if the lowest dose is eliminated.
#'
#'
#'          The BOIN design has two built-in stopping rules: (1) stop the trial if the lowest dose is eliminated
#'          due to toxicity, and no dose should be selected as the MTD; and (2) stop the trial
#'          and select the MTD if the number of patients treated at the current dose reaches
#'          \code{n.earlystop}. The first stopping rule is a safety rule to protect patients
#'          from the case in which all doses are overly toxic. The rationale for the second
#'          stopping rule is that when there is a large number (i.e., \code{n.earlystop})
#'          of patients assigned to a dose, it means that the dose-finding algorithm has
#'          approximately converged. Thus, we can stop the trial early and select the MTD
#'          to save the sample size and reduce the trial duration. For some applications,
#'          investigators may prefer a more strict safety stopping rule than rule (1) for
#'          extra safety when the lowest dose is overly toxic. This can be achieved by
#'          setting \code{extrasafe=TRUE}, which imposes the following more strict safety
#'          stopping rule: stop the trial if (i) the number of patients treated at the
#'          lowest dose >=3, and (ii) \eqn{Pr(toxicity\ rate\ of\ the\ lowest\ dose > \code{target} | data)
#'          > \code{cutoff.eli}-\code{offset}}. As a tradeoff, the strong stopping rule will decrease the
#'          MTD selection percentage when the lowest dose actually is the MTD.
#'
#' @return  \code{get.boundary()} returns a list object, including the dose escalation and de-escalation
#'          boundaries \code{$lambda_e} and \code{$lambda_d} and the corresponding decision tables
#'          \code{$boundary_tab} and \code{$full_boundary_tab}. If \code{extrasafe=TRUE}, the function also returns
#'          a (more strict) safety stopping boundary \code{$stop_boundary}.
#'
#'
#' @note We should avoid setting the values of \code{p.saf} and \code{p.tox} very close to the
#'       \code{target}. This is because the small sample sizes of typical phase I trials prevent us from
#'       differentiating the target DLT rate from the rates close to it. In addition,
#'       in most clinical applications, the target DLT rate is often a rough guess,
#'       and finding a dose level with a DLT rate reasonably close to the target rate
#'       will still be of interest to the investigator. The default values provided by
#'       \code{get.boundary()} are generally reasonable for most clinical applications.
#'
#' @references Liu S. and Yuan, Y. (2015). Bayesian Optimal Interval Designs for Phase I
#'             Clinical Trials, \emph{Journal of the Royal Statistical Society: Series C}, 64, 507-523.
#'
#'             Yuan Y., Hess K.R., Hilsenbeck S.G. and Gilbert M.R. (2016). Bayesian Optimal Interval Design: A
#'        Simple and Well-performing Design for Phase I Oncology Trials, \emph{Clinical Cancer Research}, 22, 4291-4301.
#'
#' @seealso Tutorial: \url{http://odin.mdacc.tmc.edu/~yyuan/Software/BOIN/BOIN2.6_tutorial.pdf}
#'
#'          Paper: \url{http://odin.mdacc.tmc.edu/~yyuan/Software/BOIN/paper.pdf}
#'
#' @author Suyu Liu and Ying Yuan
#'
#' @examples
#'
#' ## get the dose escalation and deescalation boundaries for BOIN design with
#' ## the target DLT rate of 0.3, maximum sample size of 30, and cohort size of 3
#' bound <- get.boundary(target=0.3, ncohort=10, cohortsize=3)
#' summary(bound) # get the descriptive summary of the boundary
#' plot(bound)    # plot the flowchart of the design with boundaries
#'
#' @import stats
#' @export
get.boundary <- function (target, ncohort, cohortsize, n.earlystop = 100, p.saf = 0.6 *
                            target, p.tox = 1.4 * target, cutoff.eli = 0.95, extrasafe = FALSE,
                          offset = 0.05)
{
  density1 <- function(p, n, m1, m2) {
    pbinom(m1, n, p) + 1 - pbinom(m2 - 1, n, p)
  }
  density2 <- function(p, n, m1) {
    1 - pbinom(m1, n, p)
  }
  density3 <- function(p, n, m2) {
    pbinom(m2 - 1, n, p)
  }
  if (target < 0.05) {
    stop("the target is too low! ")

  }
  if (target > 0.6) {
    stop("the target is too high!")

  }
  if ((target - p.saf) < (0.1 * target)) {
    stop("the probability deemed safe cannot be higher than or too close to the target!")
  }
  if ((p.tox - target) < (0.1 * target)) {
    stop("the probability deemed toxic cannot be lower than or too close to the target!")
  }
  if (offset >= 0.5) {
    stop("the offset is too large!")
  }
  if (n.earlystop <= 6) {
    warning("the value of n.earlystop is too low to ensure good operating characteristics. Recommend n.earlystop = 9 to 18.")

  }
  npts = ncohort * cohortsize
  ntrt = NULL
  b.e = NULL
  b.d = NULL
  elim = NULL
  for (n in 1:npts) {
    lambda1 = log((1 - p.saf)/(1 - target))/log(target *
                                                  (1 - p.saf)/(p.saf * (1 - target)))
    lambda2 = log((1 - target)/(1 - p.tox))/log(p.tox * (1 -
                                                           target)/(target * (1 - p.tox)))
    cutoff1 = floor(lambda1 * n)
    cutoff2 = ceiling(lambda2 * n)
    ntrt = c(ntrt, n)
    b.e = c(b.e, cutoff1)
    b.d = c(b.d, cutoff2)
    elimineed = 0
    if (n < 3) {
      elim = c(elim, NA)
    }
    else {
      for (ntox in 1:n) {
        if (1 - pbeta(target, ntox + 1, n - ntox + 1) >
            cutoff.eli) {
          elimineed = 1
          break
        }
      }
      if (elimineed == 1) {
        elim = c(elim, ntox)
      }
      else {
        elim = c(elim, NA)
      }
    }
  }
  for (i in 1:length(b.d)) {
    if (!is.na(elim[i]) && (b.d[i] > elim[i]))
      b.d[i] = elim[i]
  }
  boundaries0 = rbind(ntrt, b.e, b.d, elim)[, 1:min(npts, n.earlystop)]
  rownames(boundaries0) = c("Number of patients treated", "Escalate if # of DLT <=",
                            "Deescalate if # of DLT >=", "Eliminate if # of DLT >=")
  colnames(boundaries0) = rep("", min(npts, n.earlystop))
  out = list()
  if (cohortsize > 1) {
    out = list(lambda_e = lambda1, lambda_d = lambda2,
               boundary_tab = boundaries0[,(1:floor(min(npts, n.earlystop)/cohortsize)) * cohortsize],
               full_boundary_tab = boundaries0)
  }
  else out = list(lambda_e = lambda1, lambda_d = lambda2, boundary_tab = boundaries0[,
                                                                                     (1:floor(min(npts, n.earlystop)/cohortsize)) * cohortsize])
  if (extrasafe) {
    stopbd = NULL
    ntrt = NULL
    for (n in 1:npts) {
      ntrt = c(ntrt, n)
      if (n < 3) {
        stopbd = c(stopbd, NA)
      }
      else {
        for (ntox in 1:n) {
          if (1 - pbeta(target, ntox + 1, n - ntox +
                        1) > cutoff.eli - offset) {
            stopneed = 1
            break
          }
        }
        if (stopneed == 1) {
          stopbd = c(stopbd, ntox)
        }
        else {
          stopbd = c(stopbd, NA)
        }
      }
    }
    stopboundary = rbind(ntrt, stopbd)[, 1:min(npts, n.earlystop)]
    rownames(stopboundary) = c("The number of patients treated at the lowest dose  ",
                               "Stop the trial if # of DLT >=        ")
    colnames(stopboundary) = rep("", min(npts, n.earlystop))
    out = c(out, list(target = target, cutoff = cutoff.eli - offset, stop_boundary = stopboundary))
  }
 class(out)<-"boin"
  return(out)
}
