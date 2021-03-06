#'
#' Select the maximum tolerated dose (MTD) or MTD contour for drug combination trials
#'
#' Select the maximum tolerated dose (MTD) or MTD contour after the drug combination trial is
#' completed using the BOIN design or waterfall design
#'
#'
#' @param target the target DLT rate
#' @param npts a \code{J*K} matrix \code{(J<=K)} containing the number of patients treated at each dose combination
#' @param ntox a \code{J*K} matrix \code{(J<=K)} containing the number of patients experienced
#'             dose-limiting toxicity at each dose combination
#' @param cutoff.eli the cutoff to eliminate an overly toxic dose for safety.
#'                   We recommend the default value of (\code{cutoff.eli=0.95})
#'                   for general use.
#' @param extrasafe set \code{extrasafe=TRUE} to impose a more strict stopping
#'                  rule for extra safety
#' @param offset a small positive number (between \code{0} and \code{0.5}) to control how
#'               strict the stopping rule is when \code{extrasafe=TRUE}. A
#'               larger value leads to a more strict stopping rule. The
#'               default value \code{offset=0.05} generally works well.
#' @param boundMTD set \code{boundMTD=TRUE} to impose the condition: the isotonic estimate of toxicity
#'                 probability for the selected MTD must be less than de-escalation boundary.
#' @param p.tox the lowest toxicity probability that is deemed overly toxic such
#'              that deescalation is required. The default value is
#'                \code{p.tox=1.4*target}.
#' @param mtd.contour set \code{mtd.contour=TRUE} to select the MTD contour,
#'                    otherwise select a single MTD. The value of \code{mtd.contour}
#'                    should be consistent with that in \code{get.oc.comb()}.
#'
#'
#' @return \code{select.mtd.comb()} returns returns (1) target toxicity probability (\code{$target}),
#' (2) selected MTD or MTD contour (\code{$MTD}),
#' (3) isotonic estimate of the DLT probablity at each dose (\code{$p_est}).
#'
#'
#' @details \code{select.mtd.comb()} selects a MTD or the MTD contour based
#'          on matrix isotonic estimates of toxicity probabilities, depending on
#'          \code{mtd.contour} is set as \code{TRUE} or \code{FALSE}. The (matrix)
#'          isotonic estimates are obtained by the R package (Iso::biviso).
#'
#' @note The MTD selection and dose escalation/deescalation rule are two independent
#'       components of the trial design. When appropriate, another dose selection
#'       procedure (e.g., based on a fitted logistic model) can be used to select
#'       the MTD after the completion of the trial using the BOIN or waterfall design.
#'
#' @author Suyu Liu, Liangcai Zhang, Yanhong Zhou, and Ying Yuan
#'
#' @references Liu S. and Yuan, Y. (2015). Bayesian Optimal Interval Designs for Phase I Clinical
#'             Trials, \emph{Journal of the Royal Statistical Society: Series C}, 64, 507-523.
#'
#'            Lin R. and Yin, G. (2017). Bayesian Optimal Interval Designs for Dose Finding in
#'            Drug-combination Trials, \emph{Statistical Methods in Medical Research}, 26, 2155-2167.
#'
#'            Yan, F., Zhang, L., Zhou, Y., Pan, H., Liu, S. and Yuan, Y. (2020).BOIN: An R Package
#'            for Designing Single-Agent and Drug-Combination Dose-Finding Trials Using Bayesian Optimal
#'            Interval Designs. \emph{Journal of Statistical Software}, 94(13),1-32.<doi:10.18637/jss.v094.i13>.
#'
#'
#'            Zhang L. and Yuan, Y. (2016). A Simple Bayesian Design to Identify the Maximum
#'            Tolerated Dose Contour for Drug Combination Trials, \emph{Statistics in Medicine}, 35, 4924-4936.
#'
#' @seealso  Tutorial: \url{http://odin.mdacc.tmc.edu/~yyuan/Software/BOIN/BOIN2.6_tutorial.pdf}
#'
#'           Paper: \url{http://odin.mdacc.tmc.edu/~yyuan/Software/BOIN/paper.pdf}
#'
#' @examples
#'
#' ### drug-combination trial to find a single MTD
#'
#' ## Select the MTD based on the data from a 3x5 combination trial
#' ## matrix n contains the number of patients treated at each dose combination
#' ## matrix y contains the number of patients experienced toxicity at each dose combination
#' n <- matrix(c(3, 5, 0, 0, 0, 7, 6, 15, 0, 0, 0, 0, 4, 0, 0), ncol=5, byrow=TRUE)
#' y <- matrix(c(0, 1, 0, 0, 0, 1, 1, 4, 0, 0, 0, 0, 2, 0, 0), ncol=5, byrow=TRUE)
#' sel.comb <- select.mtd.comb(target=0.3, npts=n, ntox=y)
#' summary(sel.comb)
#' plot(sel.comb)
#'
#'
#' ### drug-combination trial to find the MTD contour
#'
#' ## Select the MTD contour based on the data from a 3x4 combination trial
#' ## matrix n contains the number of patients treated at each dose combination
#' ## matrix y contains the number of patients experienced toxicity at each dose combination
#' n <- matrix(c(6, 9, 24, 0,  6, 24, 9, 0,  12, 18, 0, 0), ncol=4, byrow=TRUE)
#' y <- matrix(c(0, 1,  5, 0,  1,  5, 4, 0,  1, 5, 0, 0), ncol=4, byrow=TRUE)
#' sel.comb2 <- select.mtd.comb(target=0.3, npts=n, ntox=y, mtd.contour=TRUE)
#' summary(sel.comb2)
#' plot(sel.comb2)
#'
#' @export
select.mtd.comb <- function (target, npts, ntox, cutoff.eli = 0.95, extrasafe = FALSE,
                             offset = 0.05, boundMTD=FALSE, p.tox=1.4*target,mtd.contour = FALSE)
{
  lambda_d = log((1 - target)/(1 - p.tox))/log(p.tox * (1 -target)/(target * (1 - p.tox)))
  y = ntox
  n = npts
  if (nrow(n) > ncol(n) | nrow(y) > ncol(y)) {
   stop("npts and ntox should be arranged in a way (i.e., rotated) such that for each of them, the number of rows is less than or equal to the number of columns.")

  }
  elimi = matrix(0, dim(n)[1], dim(n)[2])
  if (extrasafe) {
    if (n[1, 1] >= 3) {
      if (1 - pbeta(target, y[1, 1] + 1, n[1, 1] - y[1,
                                                     1] + 1) > cutoff.eli - offset) {
        elimi[, ] = 1
      }
    }
  }
  for (i in 1:dim(n)[1]) {
    for (j in 1:dim(n)[2]) {
      if (n[i, j] >= 3) {
        if (1 - pbeta(target, y[i, j] + 1, n[i, j] -
                      y[i, j] + 1) > cutoff.eli) {
          elimi[i:dim(n)[1], j] = 1
          elimi[i, j:dim(n)[2]] = 1
          break
        }
      }
    }
  }

  selectdose=NULL

  if (elimi[1] == 1) {
    selectdose = c(99, 99)
    selectdoses = matrix(selectdose, nrow = 1)
  }else {
    phat = (y + 0.05)/(n + 0.1)
    phat = round(Iso::biviso(phat, n + 0.1, warn = TRUE)[, ],2)
   # phat.out = phat
    lower.mat=qbeta(0.025,y+0.05,n-y+0.05)
    lower.mat=round(Iso::biviso(lower.mat),2)

    upper.mat=qbeta(0.975,y+0.05,n-y+0.05)
    upper.mat=round(Iso::biviso(upper.mat),2)
    phat.out<-matrix(paste0(format(phat,digits=1),"(",lower.mat,", ",upper.mat,")"),byrow=FALSE,nrow=dim(phat)[1])
    colnames(phat.out)=paste0("B",1:dim(n)[2])
    rownames(phat.out)=paste0("A",1:dim(n)[1])
    phat.out.noCI=round(phat,2)
    phat.out[n == 0] = "NA"
    phat[elimi == 1] = 1.1
    phat = phat * (n != 0) + (1e-05) * (matrix(rep(1:dim(n)[1],
                                                   each = dim(n)[2], len = length(n)), dim(n)[1], byrow = T) +
                                          matrix(rep(1:dim(n)[2], each = dim(n)[1], len = length(n)),
                                                 dim(n)[1]))

    if(boundMTD){
      if(all(phat[n!=0]>lambda_d)){
        selectdose = c(99, 99)
        selectdoses = matrix(selectdose, nrow = 1)
      }else{
        phat[phat>lambda_d]=10}}

     if(is.null(selectdose)){
      phat[n == 0] = 10
      selectdose = which(abs(phat - target) == min(abs(phat -
                                                         target)), arr.ind = TRUE)


    if (length(selectdose) > 2)
      selectdose = selectdose[1, ]
    aa = function(x) as.numeric(as.character(x))
    if (mtd.contour == TRUE) {
      selectdoses = cbind(row = 1:dim(n)[1], col = rep(99,
                                                       dim(n)[1]))
      for (k in dim(n)[1]:1) {
        kn = n[k, ]
        ky = y[k, ]
        kelimi = elimi[k, ]
        kphat = phat[k, ]
        if (kelimi[1] == 1 || sum(n[kelimi == 0]) ==
            0) {
          kseldose = 99
        }else {
          adm.set = (kn != 0) & (kelimi == 0)
          adm.index = which(adm.set == T)
          y.adm = ky[adm.set]
          n.adm = kn[adm.set]
          selectd = sort(abs(kphat[adm.set] - target),
                         index.return = T)$ix[1]
          kseldose = adm.index[selectd]
        }
        selectdoses[k, 2] = ifelse(is.na(kseldose), 99,
                                   kseldose)
        if (k < dim(n)[1])
          if (selectdoses[k + 1, 2] == dim(n)[2])
            selectdoses[k, 2] = dim(n)[2]
        if (k < dim(n)[1])
          if (aa(selectdoses[k + 1, 2]) == dim(n)[2] &
              aa(selectdoses[k + 1, 2]) == aa(selectdoses[k,
                                                          2]))
            selectdoses[k, 2] = 99
      }
    }else {
      selectdoses = matrix(99, nrow = 1, ncol = 2)
      selectdoses[1, ] = matrix(selectdose, nrow = 1)
    }

      selectdoses = matrix(selectdoses[selectdoses[, 2] !=
                                         99, ], ncol = 2)
    }

    colnames(selectdoses) = c("DoseA", "DoseB")

  }
  if (mtd.contour == FALSE) {
    if (selectdoses[1, 1] == 99 && selectdoses[1, 2] == 99) {
      cat("All tested doses are overly toxic. No MTD is selected! \n")}
    #  out=list(target = target, MTD = 99, p_est = matrix(NA,nrow = dim(npts)[1], ncol = dim(npts)[2]))
   # }
   # else {

      out=list(target = target, MTD = selectdoses, p_est=phat.out.noCI,p_est_CI = phat.out)
   # }

    class(out)<-"boin"
    return(out)
  }

  else {
    if (length(selectdoses) == 0) {
      cat("All tested doses are overly toxic. No MTD is selected! \n")
      out=list(target = target, MTD = 99, p_est = matrix(NA,nrow = dim(npts)[1], ncol = dim(npts)[2]))
    }
    else {

      out=list(target = target, MTD = selectdoses,  p_est=phat.out.noCI,p_est_CI = phat.out)
    }

    class(out)<-"boin"
    return(out)
  }
}
