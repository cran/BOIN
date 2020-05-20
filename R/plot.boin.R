#'
#' Plot the flowchart and simulation results for BOIN designs
#'
#' Plot the objects returned by other functions, including (1) flowchart of BOIN design;
#' (2) operating characteristics of the design, including selesction percentage and the
#' number of patients treated at each dose;
#' (3) the estimate of toxicity probability for each dose and corresponding 95\% credible interval
#'
#'
#' @param x the object returned by other functions
#' @param ... ignored arguments
#' @param name the name of the object to be plotted.
#'             User doesn't need to input this parameter.
#'
#' @return \code{plot()} returns a figure or a series of figures depending on the object entered
#'
#' @author Suyu Liu, Liangcai Zhang, Yanhong Zhou,  and Ying Yuan
#'
#' @examples
#'
#' ###### single-agent trial ######
#'
#' ## get dose escalation and deescalation boundaries for conducting the trial
#' bound <- get.boundary(target=0.3, ncohort=10, cohortsize=3)
#' plot(bound)
#'
#'
#' ## get the operating characteristics for BOIN single agent trial
#' oc <- get.oc(target=0.3, p.true=c(0.05,0.15,0.3,0.45,0.6),
#'    ncohort=10, cohortsize=3, ntrial=1000)
#' summary(oc)
#' plot(oc)
#'
#'
#' ## select the MTD based on the trial data
#' n <- c(3, 3, 15, 9, 0)
#' y <- c(0, 0, 4, 4, 0)
#' selmtd <- select.mtd(target=0.3, npts=n, ntox=y)
#' summary(selmtd)
#' plot(selmtd)
#'
#'
#' ###### drug-combination trial ######
#'
#' ##### combination trial to find a single MTD  ######
#'
#' ## get the operating characteristics for BOIN combination trial
#' p.true <- matrix(c(0.01,0.03,0.10,0.20,0.30,
#'                 0.03,0.05,0.15,0.30,0.60,
#'                 0.08,0.10,0.30,0.60,0.75), byrow=TRUE, ncol=5)
#'
#' oc.comb <- get.oc.comb(target=0.3, p.true, ncohort=20, cohortsize=3, n.earlystop=12,
#'      startdose=c(1,1),ntrial=100)
#' summary(oc.comb)
#' plot(oc.comb)
#'
#'
#' ## select a MTD based on the trial data
#' n <- matrix(c(3, 5, 0, 0, 0, 7, 6, 15, 0, 0, 0, 0, 4, 0, 0), ncol=5, byrow=TRUE)
#' y <- matrix(c(0, 1, 0, 0, 0, 1, 1, 4, 0, 0, 0, 0, 2, 0, 0), ncol=5, byrow=TRUE)
#' sel.comb <- select.mtd.comb(target=0.3, npts=n, ntox=y)
#' summary(sel.comb)
#' plot(sel.comb)
#'
#'
#' ##### combination trial to find a MTD contour (e.g., multiple MTDs)  #####
#'
#' ## get the operating characteristics for BOIN waterfall design
#' p.true <- matrix(c(0.01, 0.10, 0.20, 0.30,
#'                 0.03, 0.15, 0.30, 0.60,
#'                 0.08, 0.30, 0.60, 0.75), byrow=TRUE, ncol=4)
#'
#' oc.comb2 <- get.oc.comb(target=0.3, p.true, ncohort=c(8,6,6), cohortsize=3, n.earlystop=12,
#'        startdose=c(1,1), ntrial=100, mtd.contour=TRUE)
#' summary(oc.comb2)
#' plot(oc.comb2)
#'
#'
#' ## select the MTD contour based on the trial data
#' n <- matrix(c(6, 9, 24, 0,  6, 24, 9, 0,  12, 18, 0, 0), ncol=4, byrow=TRUE)
#' y <- matrix(c(0, 1,  5, 0,  1,  5, 4, 0,  1, 5, 0, 0), ncol=4, byrow=TRUE)
#' sel.comb2 <- select.mtd.comb(target=0.3, npts=n, ntox=y, mtd.contour=TRUE)
#' summary(sel.comb2)
#' plot(sel.comb2)
#'
#' @importFrom grDevices dev.flush dev.hold devAskNewPage
#' @importFrom graphics abline arrows arrows axis barplot legend mtext par plot points rect segments text
#' @export
plot.boin<- function (x,..., name = deparse(substitute(x)))
{
  new.obj = unlist(strsplit(name, split = "\\$"))
  strpattern = "none"
  if (length(new.obj) >= 2) {
    strpattern = new.obj[2]
  }
  assign("objectPlot", get(new.obj[1]))
  get.flowchart<-function(){
    lambda2 = round(objectPlot$lambda_e, 3)
    lambda1 = round(objectPlot$lambda_d, 3)
    # if (.Platform$OS.type == "windows") {
    #   dev.new(height = 7.36, width = 6.27, rescale = "fit")
    # }
    # else {
    #   dev.new(height = 7.367442, width = 6.580645)
    # }
    par(lwd = 1.5, mar = c(1, 1, 1, 1))

    plot(0, type = "n", xlim = c(1, 10), ylim = c(-3,
                                                  10.2), xaxt = "n", yaxt = "n", bty = "n", pch = "",
         ylab = "", xlab = "")
    theta = seq(0, 2 * pi, len = 100)
    r = 0.56
    x = 5 + 2 * r * cos(theta)
    y = 10 + r * sin(theta)
    points(x, y, type = "l")
    arrows(5, 10 - 0.56, 5, 8.5, length = 0.15)
    rect(4, 7.5, 6, 8.5)
    arrows(5, 7.5, 5, 6.5, length = 0.15)
    arrows(10, 8, 6, 8, length = 0.15)
    r = 0.5
    x = 2 + 2 * r * cos(theta)
    y = 5.5 + r * sin(theta)
    points(x, y, type = "l")
    arrows(4, 5.5, 3, 5.5, length = 0.15)
    segments(4, 5.5, 5, 6.5)
    segments(4, 5.5, 5, 4.5)
    segments(5, 6.5, 6, 5.5)
    segments(5, 4.5, 6, 5.5)
    arrows(5, 4.5, 5, 3.5, length = 0.15)
    segments(5, 3.5, 4, 2.5)
    segments(4, 2.5, 5, 1.5)
    segments(5, 1.5, 6, 2.5)
    segments(6, 2.5, 5, 3.5)
    segments(2, 2.5, 4, 2.5)
    arrows(2, 2.5, 2, 0.5, length = 0.15)
    segments(5, 1.5, 5, 1.3)
    arrows(5, 0.9, 5, 0.5, length = 0.15)
    segments(6, 2.5, 8, 2.5)
    arrows(8, 2.5, 8, 0.5, length = 0.15)
    rect(1, -0.5, 3, 0.5)
    rect(4, -0.5, 6, 0.5)
    rect(7, -0.5, 9, 0.5)
    segments(2, -0.5, 2, -1.25)
    segments(5, -0.5, 5, -2)
    segments(8, -0.5, 8, -1.25)
    segments(2, -1.25, 8, -1.25)
    segments(5, -2, 10, -2)
    segments(10, 8, 10, -2)
    text(5, 10, labels = "Start \n at the prespecified \n starting dose",
         cex = 0.8)
    text(5, 8, labels = "Treat a patient or a \n cohort of patients",
         cex = 0.8)
    text(2, 5.5, labels = "Stop the trial and \n select the MTD",
         cex = 0.8)
    text(3.5, 5.8, labels = "Yes", cex = 0.8)
    text(5, 5.6, labels = "Reach \n the maximum \n sample size",
         cex = 0.8)
    text(5.2, 4.2, labels = "No", cex = 0.8)
    text(3, 2.8, labels = expression("" <= ""), cex = 0.8)
    text(3.4, 2.8, labels = lambda2, cex = 0.8)
    text(7, 2.8, labels = expression("" >= ""), cex = 0.8)
    text(7.4, 2.8, labels = lambda1, cex = 0.8)
    text(5, 2.4, labels = "Compute \n the DLT rate* \n at the current \n dose",
         cex = 0.8)
    text(5, 1.06, labels = paste("Within (", lambda2,
                                 ", ", lambda1, ")", sep = ""), cex = 0.8)
    text(2, 0, labels = "Escalate the dose", cex = 0.8)
    text(5, 0, labels = "Retain the current \n dose",
         cex = 0.8)
    text(8, 0, labels = "De-escalate the \n dose", cex = 0.8)
    text(par("usr")[2]/2, -3, expression(paste("* DLT rate = ",
                                               frac("Total number of patients who experienced DLT at the current dose",
                                                    "Total number of patients treated at the current dose"),
                                               sep = "")), cex = 0.8, adj = c(0.5, NA))
  }

  if (!is.element(strpattern, c("none", names(objectPlot)))) {
    warning("Please double check and specify the variable to be plotted...\n")
  }
  else {
    #determine if flowchart is plotted, which argument is ignored if flowchart is plotted
    if (!is.null(objectPlot$boundary_tab) | (!is.null(objectPlot$percentstop) &
                                             strpattern == "flowchart")) {
     get.flowchart()
    }
    else if (!is.null(objectPlot$percentstop)) { #oc for single-agent trial is entered

      get.flowchart()
      oask <- devAskNewPage(TRUE)
      on.exit(devAskNewPage(oask))
      dev.flush()
      dev.hold()
      par(mar = c(5, 6, 4, 2))
        bplot = barplot(objectPlot$selpercent, ylab = "selection percentage (%)",
                        ylim = c(0, 100), cex.names = 1, xaxt = "n",
                        cex.lab = 1.3)
        mtext("Selection percentage", 3, line = 0, cex = 1.3)
        axis(1, at = bplot, labels = seq(1, length(objectPlot$selpercent)))
        mtext("Dose level", 1, line = 2, cex = 1)
        dev.flush()
        dev.hold()
        bplot = barplot(objectPlot$npatients, ylab = "Number of patients",
                        ylim = c(0, sum(objectPlot$npatients)), cex.names = 1,
                        beside = FALSE, xaxt = "n", cex.lab = 1.3)
        axis(1, at = bplot, labels = seq(1, length(objectPlot$npatients)))
        mtext("Patient allocation", 3, line = 0, cex = 1.3)
        mtext("Dose level", 1, line = 2, cex = 1)
        dev.flush()
        dev.hold()
        bplot = barplot(objectPlot$ntox, ylab = "Number of toxicities",
                        ylim = c(0, sum(objectPlot$ntox)), cex.names = 1,
                        beside = FALSE, xaxt = "n", cex.lab = 1.3)
        axis(1, at = bplot, labels = seq(1, length(objectPlot$ntox)))
        mtext("Observed toxicity", 3, line = 0, cex = 1.3)
        mtext("Dose level", 1, line = 2, cex = 1)


    }
    else if (!is.null(objectPlot$pcs) | !is.null(objectPlot$pcs.contour)) {

      if (is.null(objectPlot$pcs.contour)) { ##Diagram for combination trial when finding single-agent
        J = nrow(objectPlot$p.true)
        K = ncol(objectPlot$p.true)
        xlab = "Drug B"
        ylab = "Drug A"
        if (J > K) {
          S = J
          J = K
          K = S
          xlab = "Drug A"
          ylab = "Drug B"
        }
        xmax = K * 2 - 1
        xmin = 1
        ymax = 2 * J
        ymin = 1
        ptcex = 1.5
        par(mar = c(5, 5, 2, 2))
        plot(1:xmax, xlim = c(xmin, xmax), ylim = c(ymin,
                                                    ymax + 0.5), pch = "", axes = F, xlab = xlab,
             ylab = ylab, cex.axis = 1, cex.lab = 1)
        for (i in seq(1, xmax, by = 2)) for (j in seq(1,
                                                      ymax - 1, by = 2)) points(i, j, pch = 1,
                                                                                cex = ptcex)
        if (J > 2 & K > 2) {
          arrows(3 + 0.1 * J/K, 3, 3 + 1, 3, col = 3,
                 length = 0.06, lty = 1, lwd = 2)
          arrows(3 - 0.1 * J/K, 3, 3 - 1, 3, col = 2,
                 length = 0.06, lty = 1, lwd = 2)
          arrows(3, 3 - 0.1 * J/K, 3, 3 - 1, col = 2,
                 length = 0.06, lty = 1, lwd = 2)
          arrows(3, 3 + 0.1 * J/K, 3, 3 + 1, col = 3,
                 length = 0.06, lty = 1, lwd = 2)
          points(3, 3, pch = 19, cex = ptcex)
        }
        if (J == 2 | K == 2) {
          arrows(1 + 0.2 * J/K, 1, 1 + 1, 1, col = 3,
                 length = 0.06, lty = 1, lwd = 2)
          arrows(1, 1 + 0.2 * J/K, 1, 1 + 1, col = 3,
                 length = 0.06, lty = 1, lwd = 2)
          points(1, 1, pch = 19, cex = ptcex)
        }
        xx = quantile(1:xmax, c(1/4, 3/4))
        text(x = xx[1] + 0.1 * ymax/(xmax + 1), y = ymax,
             labels = "escalation", pos = 4)
        arrows(xx[1] - 1, ymax, xx[1], ymax, col = 3,
               length = 0.06, lty = 1, lwd = 2)
        text(x = xx[2] + 0.1 * ymax/(xmax + 1), y = ymax,
             labels = "de-escalation", pos = 4)
        arrows(xx[2] - 1, ymax, xx[2], ymax, col = 2,
               length = 0.06, lty = 1, lwd = 2)
      }
      else {##Waterfall design
        J = nrow(objectPlot$p.true)
        K = ncol(objectPlot$p.true)
        xlab = "Drug B"
        ylab = "Drug A"
        if (J > K) {
          S = J
          J = K
          K = S
          xlab = "Drug A"
          ylab = "Drug B"
        }
        xmax = K * 2
        xmin = -0.8
        ymax = J * (J + 1) + 2 * J
        ymin = 1 - 0.5
        ptcex = 1.5
        par(mar = c(5, 5, 4, 2))
        plot(1:xmax, xlim = c(xmin, xmax), ylim = c(ymin,
                                                    ymax + 0.5), pch = "", axes = F, xlab = xlab,
             ylab = ylab, cex.axis = 1, cex.lab = 1)
        active.rows = NULL
        for (j in 1:(J + 1)) active.rows = c(active.rows,
                                             1:J + (J + 2) * (j - 1))
        for (i in seq(1, xmax, by = 2)) for (j in 1:ymax) if (is.element(j,
                                                                         active.rows))
          points(i, j, pch = 1, cex = ptcex)
        text(-0.5, sort(seq(ymax - (J - 1)/2, (J +
                                                 1)/2, len = J + 1), decreasing = FALSE),
             paste("(", letters[seq(J + 1, 1)], ")", sep = ""),
             cex = 1)
        segments(0.5, ymax + 0.5, xmax, ymax + 0.5)
        segments(0.5, ymax - (J - 1) - 0.5, 0.5, ymax +
                   0.5)
        segments(0.5, ymax - (J - 1) - 0.5, 1.5, ymax -
                   (J - 1) - 0.5)
        segments(1.5, ymax - (J - 1) - 0.5, 1.5, ymax -
                   0.5)
        segments(1.5, ymax - 0.5, xmax, ymax - 0.5)
        segments(xmax, ymax - 0.5, xmax, ymax + 0.5)
        mtds = NULL
        crows = NULL
        arrows(1, ymax - J + 1 + 0.15 * J/K, 1, ymax -
                 J + 1 + 0.6, col = 1, length = 0.06, lty = 1,
               lwd = 2)
        tmpx = sort(sample(seq(1, xmax - 2, 2), J -
                             1), decreasing = FALSE)
        for (j in 1:(J - 1)) {
          crow = ymax - J * j - 2 * j - j
          crows = c(crows, crow)
          rect(2, crow - 0.5, xmax, crow + 0.5)
          points(tmpx[j], crow + 1, pch = 8, cex = ptcex)
          mtds = rbind(mtds, c(tmpx[j], J - j + 1))
          if (tmpx[j] + 2 < xmax - 1) {
            arrows(tmpx[j] + 2 + 0.2, crow, tmpx[j] +
                     2 + 0.5, crow, col = 1, length = 0.08,
                   lty = 1, lwd = 2)
          }
          else {
            arrows(tmpx[j] + 2 - 0.2, crow, tmpx[j] +
                     2 - 0.5, crow, col = 1, length = 0.08,
                   lty = 1, lwd = 2)
          }
        }
        mtds = rbind(mtds, c(xmax - 1, 1))
        for (j in 1:J) {
          x = mtds[j, 1]
          y = mtds[j, 2]
          points(rep(x, y), y + (J + 2) * 0:(y - 1),
                 pch = 8, cex = ptcex)
        }
      }

      oask <- devAskNewPage(TRUE)
      on.exit(devAskNewPage(oask))
      dev.flush()
      dev.hold()
      par(mar = c(5, 6, 4, 2))

      sel.comb=objectPlot$selpercent
      rownames(sel.comb)=1:dim(sel.comb)[1]
      colnames(sel.comb)=1:dim(sel.comb)[2]
      barplot(sel.comb,beside=TRUE,ylab="Selection percentage (%)",
              xlab="Drug B",ylim = c(0, 100),legend.text=rownames(sel.comb),
              args.legend=list(title="Drug A",horiz=TRUE,x="top"))
      mtext("Selection percentage", 3, line = 1, cex = 1.3)


      dev.flush()
      dev.hold()
      npts.comb=objectPlot$npatients
      rownames(npts.comb)=1:dim(npts.comb)[1]
      colnames(npts.comb)=1:dim(npts.comb)[2]
      barplot(npts.comb,beside=TRUE,ylab="Number of patients",
              xlab="Drug B",ylim=c(0,sum(npts.comb,na.rm=TRUE)),
              legend.text=rownames(npts.comb),
              args.legend=list(title="Drug A",horiz=TRUE,x="top"))
      mtext("Patient allocation", 3, line = 1, cex = 1.3)
      dev.flush()
      dev.hold()
      ntox.comb=objectPlot$ntox
      rownames(ntox.comb)=1:dim(ntox.comb)[1]
      colnames(ntox.comb)=1:dim(ntox.comb)[2]
      barplot(ntox.comb,beside=TRUE,ylab="Number of toxicities",
              xlab="Drug B", ylim=c(0,sum(ntox.comb,na.rm=TRUE)),
              legend.text=rownames(ntox.comb),
              args.legend=list(title="Drug A",horiz=TRUE,x="top"))
      mtext("Observed toxicity", 3, line = 1, cex = 1.3)

    }
    else if (!is.null(objectPlot$MTD)) {##select mtd
      if (objectPlot$MTD[1] == 99) {
        warning("All tested doses are overly toxic. No MTD is selected!\n")
      }
      else {
        if (!is.null(objectPlot$p_est)) {

          par(mfrow = c(1, 1), mar = c(5, 5, 4, 2))
          if (length(objectPlot$MTD) >= 2) {
            p_est.comb=objectPlot$p_est
            rownames(p_est.comb)=1:dim(p_est.comb)[1]
            colnames(p_est.comb)=1:dim(p_est.comb)[2]
            barplot(p_est.comb,beside=TRUE,ylab="DLT rate",
                    ylim=c(0,round(max(p_est.comb,na.rm=TRUE)*1.5,1)),xlab="Drug B",legend.text=rownames(p_est.comb),
                    args.legend=list(title="Drug A",horiz=TRUE,x="top"))
          }
          else {
            p_est = objectPlot$p_est
            p_hat = p_est[, 2]
            ci = p_est[, 3]
            ci = gsub("[\\(\\)]", "", ci)
            conf.intv = matrix(unlist(strsplit(ci, ",")),
                               byrow = TRUE, ncol = 2)
            if (p_est[1, 2] == "----") {
              warning("The trial is stopped since the lowest dose is too toxic.\n")
            }
            else {
              numbs = ifelse(sum(p_hat == "----") ==
                               0, length(p_hat), min(which(p_hat ==
                                                             "----")) - 1)
              numbs2 = length(p_hat)
              phatx = as.numeric(as.character(p_hat[1:numbs]))
              lwr = as.numeric(as.character(conf.intv[1:numbs,
                                                      1]))
              upr = as.numeric(as.character(conf.intv[1:numbs,
                                                      2]))
              par(mar = c(5, 5, 4, 2))
              plot(1:numbs2, ylim = c(0, 1), xlab = "Dose level",
                   ylab = "DLT rate", pch = "", xaxt = "n",
                   cex.lab = 1.3)
              axis(1, at = 1:numbs2, labels = 1:numbs2)
              abline(h = objectPlot$target, lty = 2,
                     col = 2)
              points(1:numbs, phatx, pch = 19)
              arrows(x0 = 1:numbs, x1 = 1:numbs, y0 = lwr,
                     y1 = upr, code = 3, angle = 90, length = 0.1)
              if (numbs < numbs2) {
                points((numbs + 1):numbs2, seq(min(1,
                                                   max(phatx, na.rm = T) + 0.05), min(max(phatx,
                                                                                          na.rm = T) + 0.2, 1), length = numbs2 -
                                                 numbs), pch = "*", cex = 1.5)
                legend("topleft", "*   no patient treated")
              }
            }
          }
        }
        else {
          warning("Please set verbose=TRUE to get more details of the results.\n")
        }
      }
    }
    else {
      warning("Please double check and specify the variable to be plotted...\n")
    }
  }
}



###Original code
# plot.boin <- function (x,type, ..., name = deparse(substitute(x)))
# {
#   new.obj = unlist(strsplit(name, split = "\\$"))
#   strpattern = "none"
#   if (length(new.obj) >= 2) {
#     strpattern = new.obj[2]
#   }
#   assign("objectPlot", get(new.obj[1]))
#   if (!is.element(strpattern, c("none", names(objectPlot)))) {
#     warning("Please double check and specify the variable to be plotted...\n")
#   }
#   else {
#     if (!is.null(objectPlot$boundary_tab) | (!is.null(objectPlot$percentstop) &
#                                              strpattern == "flowchart")) {
#       lambda2 = round(objectPlot$lambda_e, 3)
#       lambda1 = round(objectPlot$lambda_d, 3)
#       if (.Platform$OS.type == "windows") {
#         dev.new(height = 7.36, width = 6.27, rescale = "fit")
#       }
#       else {
#         dev.new(height = 7.367442, width = 6.580645)
#       }
#       par(lwd = 1.5, mar = c(1, 1, 1, 1))
#       plot(0, type = "n", xlim = c(1, 10), ylim = c(-3,
#                                                     10.2), xaxt = "n", yaxt = "n", bty = "n", pch = "",
#            ylab = "", xlab = "")
#       theta = seq(0, 2 * pi, len = 100)
#       r = 0.56
#       x = 5 + 2 * r * cos(theta)
#       y = 10 + r * sin(theta)
#       points(x, y, type = "l")
#       arrows(5, 10 - 0.56, 5, 8.5, length = 0.15)
#       rect(4, 7.5, 6, 8.5)
#       arrows(5, 7.5, 5, 6.5, length = 0.15)
#       arrows(10, 8, 6, 8, length = 0.15)
#       r = 0.5
#       x = 2 + 2 * r * cos(theta)
#       y = 5.5 + r * sin(theta)
#       points(x, y, type = "l")
#       arrows(4, 5.5, 3, 5.5, length = 0.15)
#       segments(4, 5.5, 5, 6.5)
#       segments(4, 5.5, 5, 4.5)
#       segments(5, 6.5, 6, 5.5)
#       segments(5, 4.5, 6, 5.5)
#       arrows(5, 4.5, 5, 3.5, length = 0.15)
#       segments(5, 3.5, 4, 2.5)
#       segments(4, 2.5, 5, 1.5)
#       segments(5, 1.5, 6, 2.5)
#       segments(6, 2.5, 5, 3.5)
#       segments(2, 2.5, 4, 2.5)
#       arrows(2, 2.5, 2, 0.5, length = 0.15)
#       segments(5, 1.5, 5, 1.3)
#       arrows(5, 0.9, 5, 0.5, length = 0.15)
#       segments(6, 2.5, 8, 2.5)
#       arrows(8, 2.5, 8, 0.5, length = 0.15)
#       rect(1, -0.5, 3, 0.5)
#       rect(4, -0.5, 6, 0.5)
#       rect(7, -0.5, 9, 0.5)
#       segments(2, -0.5, 2, -1.25)
#       segments(5, -0.5, 5, -2)
#       segments(8, -0.5, 8, -1.25)
#       segments(2, -1.25, 8, -1.25)
#       segments(5, -2, 10, -2)
#       segments(10, 8, 10, -2)
#       text(5, 10, labels = "Start \n at the prespecified \n starting dose",
#            cex = 0.8)
#       text(5, 8, labels = "Treat a patient or a \n cohort of patients",
#            cex = 0.8)
#       text(2, 5.5, labels = "Stop the trial and \n select the MTD",
#            cex = 0.8)
#       text(3.5, 5.8, labels = "Yes", cex = 0.8)
#       text(5, 5.6, labels = "Reach \n the maximum \n sample size",
#            cex = 0.8)
#       text(5.2, 4.2, labels = "No", cex = 0.8)
#       text(3, 2.8, labels = expression("" <= ""), cex = 0.8)
#       text(3.4, 2.8, labels = lambda2, cex = 0.8)
#       text(7, 2.8, labels = expression("" >= ""), cex = 0.8)
#       text(7.4, 2.8, labels = lambda1, cex = 0.8)
#       text(5, 2.4, labels = "Compute \n the DLT rate* \n at the current \n dose",
#            cex = 0.8)
#       text(5, 1.06, labels = paste("Within (", lambda2,
#                                    ", ", lambda1, ")", sep = ""), cex = 0.8)
#       text(2, 0, labels = "Escalate the dose", cex = 0.8)
#       text(5, 0, labels = "Retain the current \n dose",
#            cex = 0.8)
#       text(8, 0, labels = "De-escalate the \n dose", cex = 0.8)
#       text(par("usr")[2]/2, -3, expression(paste("* DLT rate = ",
#                                                  frac("Total number of patients who experienced DLT at the current dose",
#                                                       "Total number of patients treated at the current dose"),
#                                                  sep = "")), cex = 0.8, adj = c(0.5, NA))
#     }
#     else if (!is.null(objectPlot$percentstop)) {
#       dev.new()
#       op <- par(no.readonly = TRUE)
#       dev.off()
#       par(op)
#       par(mfrow = c(1, 1), mar = c(5, 5, 4, 2))
#       if (length(new.obj) == 2) {
#         if (strpattern == "selpercent") {
#           bplot = barplot(objectPlot$selpercent, ylab = "selection percentage (%)",
#                           ylim = c(0, 100), cex.names = 1, xaxt = "n",
#                           cex.lab = 1.3)
#           axis(1, at = bplot, labels = seq(1, length(objectPlot$selpercent)))
#           mtext("Dose level", 1, line = 0, cex = 1.3)
#         }
#         else if (strpattern == "npatients") {
#           bplot = barplot(objectPlot$npatients, ylab = "Number of patients",
#                           ylim = c(0, max(objectPlot$npatients)), cex.names = 1,
#                           beside = FALSE, xaxt = "n", cex.lab = 1.3)
#           axis(1, at = bplot, labels = seq(1, length(objectPlot$npatients)))
#           mtext("Dose level", 1, line = 0, cex = 1.3)
#         }
#         else if (strpattern == "ntox") {
#           bplot = barplot(objectPlot$ntox, ylab = "number of toxicities",
#                           ylim = c(0, max(objectPlot$ntox)), cex.names = 1,
#                           beside = FALSE, xaxt = "n", cex.lab = 1.3)
#           axis(1, at = bplot, labels = seq(1, length(objectPlot$ntox)))
#           mtext("Dose level", 1, line = 0, cex = 1.3)
#         }
#       }
#       else {
#         par(mfrow = c(2, 1), mar = c(5, 5, 4, 2))
#         bplot = barplot(objectPlot$selpercent, ylab = "selection percentage (%)",
#                         ylim = c(0, 100), cex.names = 1, xaxt = "n",
#                         cex.lab = 1.3)
#         mtext("selection percentage", 3, line = 0, cex = 1.3)
#         axis(1, at = bplot, labels = seq(1, length(objectPlot$selpercent)))
#         mtext("Dose level", 1, line = 0, cex = 1.3)
#         bplot = barplot(objectPlot$npatients, ylab = "Number of patients",
#                         ylim = c(0, max(objectPlot$npatients)), cex.names = 1,
#                         beside = FALSE, xaxt = "n", cex.lab = 1.3)
#         axis(1, at = bplot, labels = seq(1, length(objectPlot$npatients)))
#         mtext("Patient allocation", 3, line = 0, cex = 1.3)
#         mtext("Dose level", 1, line = 0, cex = 1.3)
#       }
#     }
#     else if (!is.null(objectPlot$pcs) | !is.null(objectPlot$pcs.contour)) {
#       # if (!requireNamespace("epade", quietly = TRUE)) {
#       #   install.packages("epade", repos = "http://cran.us.r-project.org",
#       #                    dependencies = TRUE)
#       #   if (!requireNamespace("epade", quietly = TRUE))
#       #     stop(paste("Package: ", "epade", " not found!!!",
#       #                sep = ""))
#       # }
#       dev.new()
#       op <- par(no.readonly = TRUE)
#       dev.off()
#       par(op)
#       par(mfrow = c(1, 1), mar = c(5, 5, 2, 2))
#       if (strpattern == "selpercent") {
#         # requireNamespace("epade", quietly = TRUE)
#         # epade::bar3d.ade(objectPlot$selpercent, wall = 6,
#         #                  zw = 1, xw = 0.2, yticks = seq(0, 100, 20),
#         #                  zticks = seq(1:nrow(objectPlot$selpercent)),
#         #                  xticks = seq(1:ncol(objectPlot$selpercent)),
#         #                  xlab = "Drug B", ylab = "selection percentage (%)",
#         #                  zlab = "Drug A", col = "lavender", alpha = 0.5)
#         sel.comb=objectPlot$selpercent
#         rownames(sel.comb)=1:dim(sel.comb)[1]
#         colnames(sel.comb)=1:dim(sel.comb)[2]
#         barplot(sel.comb,beside=TRUE,ylab="Selection percentage (%)",
#                 xlab="Drug B",ylim = c(0, 100),legend.text=rownames(sel.comb),
#                 args.legend=list(title="Drug A",horiz=TRUE,x="top"))
#
#       }
#       else if (strpattern == "npatients") {
#         requireNamespace("epade", quietly = TRUE)
#         # epade::bar3d.ade(objectPlot$npatients, wall = 6,
#         #                  zw = 1, xw = 0.2, yticks = round(seq(1, max(objectPlot$npatients),
#         #                                                       max(objectPlot$npatients)/7), 0), zticks = seq(1:nrow(objectPlot$npatients)),
#         #                  xticks = seq(1:ncol(objectPlot$npatients)),
#         #                  xlab = "Drug B", ylab = "Number of patients",
#         #                  zlab = "Drug A", col = "lavender", alpha = 0.5)
#         npts.comb=objectPlot$npatients
#         rownames(npts.comb)=1:dim(npts.comb)[1]
#         colnames(npts.comb)=1:dim(npts.comb)[2]
#         barplot(npts.comb,beside=TRUE,ylab="Number of patients",
#                 xlab="Drug B",ylim=c(0,ceiling(max(npts.comb,na.rm=TRUE)*1.3)),
#                 legend.text=rownames(npts.comb),
#                 args.legend=list(title="Drug A",horiz=TRUE,x="top"))
#       }
#       else if (strpattern == "ntox") {
#         # requireNamespace("epade", quietly = TRUE)
#         # epade::bar3d.ade(objectPlot$ntox, wall = 6, zw = 1,
#         #                  xw = 0.2, yticks = round(seq(1, max(objectPlot$ntox),
#         #                                               max(objectPlot$ntox)/7), 0), zticks = seq(1:nrow(objectPlot$ntox)),
#         #                  xticks = seq(1:ncol(objectPlot$ntox)), xlab = "Drug B",
#         #                  ylab = "number of toxicities", zlab = "Drug A",
#         #                  col = "lavender", alpha = 0.5)
#         ntox.comb=objectPlot$ntox
#         rownames(ntox.comb)=1:dim(ntox.comb)[1]
#         colnames(ntox.comb)=1:dim(ntox.comb)[2]
#         barplot(ntox.comb,beside=TRUE,ylab="Number of toxicities",
#                 xlab="Drug B", ylim=c(0,ceiling(max(ntox.comb,na.rm=TRUE)*1.3)),
#                 legend.text=rownames(ntox.comb),
#                 args.legend=list(title="Drug A",horiz=TRUE,x="top"))
#
#       }
#       else if (strpattern == "flowchart") {
#         if (is.null(objectPlot$pcs.contour)) {
#           J = nrow(objectPlot$p.true)
#           K = ncol(objectPlot$p.true)
#           xlab = "Drug B"
#           ylab = "Drug A"
#           if (J > K) {
#             S = J
#             J = K
#             K = S
#             xlab = "Drug A"
#             ylab = "Drug B"
#           }
#           xmax = K * 2 - 1
#           xmin = 1
#           ymax = 2 * J
#           ymin = 1
#           ptcex = 1.5
#           par(mar = c(5, 5, 2, 2))
#           plot(1:xmax, xlim = c(xmin, xmax), ylim = c(ymin,
#                                                       ymax + 0.5), pch = "", axes = F, xlab = xlab,
#                ylab = ylab, cex.axis = 1, cex.lab = 1)
#           for (i in seq(1, xmax, by = 2)) for (j in seq(1,
#                                                         ymax - 1, by = 2)) points(i, j, pch = 1,
#                                                                                   cex = ptcex)
#           if (J > 2 & K > 2) {
#             arrows(3 + 0.1 * J/K, 3, 3 + 1, 3, col = 3,
#                    length = 0.06, lty = 1, lwd = 2)
#             arrows(3 - 0.1 * J/K, 3, 3 - 1, 3, col = 2,
#                    length = 0.06, lty = 1, lwd = 2)
#             arrows(3, 3 - 0.1 * J/K, 3, 3 - 1, col = 2,
#                    length = 0.06, lty = 1, lwd = 2)
#             arrows(3, 3 + 0.1 * J/K, 3, 3 + 1, col = 3,
#                    length = 0.06, lty = 1, lwd = 2)
#             points(3, 3, pch = 19, cex = ptcex)
#           }
#           if (J == 2 | K == 2) {
#             arrows(1 + 0.2 * J/K, 1, 1 + 1, 1, col = 3,
#                    length = 0.06, lty = 1, lwd = 2)
#             arrows(1, 1 + 0.2 * J/K, 1, 1 + 1, col = 3,
#                    length = 0.06, lty = 1, lwd = 2)
#             points(1, 1, pch = 19, cex = ptcex)
#           }
#           xx = quantile(1:xmax, c(1/4, 3/4))
#           text(x = xx[1] + 0.1 * ymax/(xmax + 1), y = ymax,
#                labels = "escalation", pos = 4)
#           arrows(xx[1] - 1, ymax, xx[1], ymax, col = 3,
#                  length = 0.06, lty = 1, lwd = 2)
#           text(x = xx[2] + 0.1 * ymax/(xmax + 1), y = ymax,
#                labels = "de-escalation", pos = 4)
#           arrows(xx[2] - 1, ymax, xx[2], ymax, col = 2,
#                  length = 0.06, lty = 1, lwd = 2)
#         }
#         else {
#           J = nrow(objectPlot$p.true)
#           K = ncol(objectPlot$p.true)
#           xlab = "Drug B"
#           ylab = "Drug A"
#           if (J > K) {
#             S = J
#             J = K
#             K = S
#             xlab = "Drug A"
#             ylab = "Drug B"
#           }
#           xmax = K * 2
#           xmin = -0.8
#           ymax = J * (J + 1) + 2 * J
#           ymin = 1 - 0.5
#           ptcex = 1.5
#           par(mar = c(5, 5, 4, 2))
#           plot(1:xmax, xlim = c(xmin, xmax), ylim = c(ymin,
#                                                       ymax + 0.5), pch = "", axes = F, xlab = xlab,
#                ylab = ylab, cex.axis = 1, cex.lab = 1)
#           active.rows = NULL
#           for (j in 1:(J + 1)) active.rows = c(active.rows,
#                                                1:J + (J + 2) * (j - 1))
#           for (i in seq(1, xmax, by = 2)) for (j in 1:ymax) if (is.element(j,
#                                                                            active.rows))
#             points(i, j, pch = 1, cex = ptcex)
#           text(-0.5, sort(seq(ymax - (J - 1)/2, (J +
#                                                    1)/2, len = J + 1), decreasing = FALSE),
#                paste("(", letters[seq(J + 1, 1)], ")", sep = ""),
#                cex = 1)
#           segments(0.5, ymax + 0.5, xmax, ymax + 0.5)
#           segments(0.5, ymax - (J - 1) - 0.5, 0.5, ymax +
#                      0.5)
#           segments(0.5, ymax - (J - 1) - 0.5, 1.5, ymax -
#                      (J - 1) - 0.5)
#           segments(1.5, ymax - (J - 1) - 0.5, 1.5, ymax -
#                      0.5)
#           segments(1.5, ymax - 0.5, xmax, ymax - 0.5)
#           segments(xmax, ymax - 0.5, xmax, ymax + 0.5)
#           mtds = NULL
#           crows = NULL
#           arrows(1, ymax - J + 1 + 0.15 * J/K, 1, ymax -
#                    J + 1 + 0.6, col = 1, length = 0.06, lty = 1,
#                  lwd = 2)
#           tmpx = sort(sample(seq(1, xmax - 2, 2), J -
#                                1), decreasing = FALSE)
#           for (j in 1:(J - 1)) {
#             crow = ymax - J * j - 2 * j - j
#             crows = c(crows, crow)
#             rect(2, crow - 0.5, xmax, crow + 0.5)
#             points(tmpx[j], crow + 1, pch = 8, cex = ptcex)
#             mtds = rbind(mtds, c(tmpx[j], J - j + 1))
#             if (tmpx[j] + 2 < xmax - 1) {
#               arrows(tmpx[j] + 2 + 0.2, crow, tmpx[j] +
#                        2 + 0.5, crow, col = 1, length = 0.08,
#                      lty = 1, lwd = 2)
#             }
#             else {
#               arrows(tmpx[j] + 2 - 0.2, crow, tmpx[j] +
#                        2 - 0.5, crow, col = 1, length = 0.08,
#                      lty = 1, lwd = 2)
#             }
#           }
#           mtds = rbind(mtds, c(xmax - 1, 1))
#           for (j in 1:J) {
#             x = mtds[j, 1]
#             y = mtds[j, 2]
#             points(rep(x, y), y + (J + 2) * 0:(y - 1),
#                    pch = 8, cex = ptcex)
#           }
#         }
#       }
#       else {
#         if (strpattern == "none") {
#           dev.new()
#           op <- par(no.readonly = TRUE)
#           dev.off()
#           par(op)
#           par(mfrow = c(2, 1), mar = c(4, 3, 2, 2))
#           # requireNamespace("epade")
#           # epade::bar3d.ade(objectPlot$selpercent, wall = 6,
#           #                  zw = 1, xw = 0.2, main = "selection percentage",
#           #                  yticks = seq(0, 100, 20), zticks = seq(1:nrow(objectPlot$selpercent)),
#           #                  xticks = seq(1:ncol(objectPlot$selpercent)),
#           #                  xlab = "Drug B", ylab = "selection percentage (%)",
#           #                  zlab = "Drug A", col = "lavender", alpha = 0.5)
#           # epade::bar3d.ade(objectPlot$npatients, wall = 6,
#           #                  zw = 1, xw = 0.2, main = "Patient allocation",
#           #                  yticks = round(seq(1, max(objectPlot$npatients),
#           #                                     max(objectPlot$npatients)/7), 0), zticks = seq(1:nrow(objectPlot$npatients)),
#           #                  xticks = seq(1:ncol(objectPlot$npatients)),
#           #                  xlab = "Drug B", ylab = "Number of patients",
#           #                  zlab = "Drug A", col = "lavender", alpha = 0.5)
#
#           sel.comb=objectPlot$selpercent
#           rownames(sel.comb)=1:dim(sel.comb)[1]
#           colnames(sel.comb)=1:dim(sel.comb)[2]
#           barplot(sel.comb,beside=TRUE,ylab="Selection percentage (%)",
#                   xlab="Drug B", ylim=c(0,100),
#                   legend.text=rownames(sel.comb),
#                   args.legend=list(title="Drug A",horiz=TRUE,x="top"))
#
#
#           npts.comb=objectPlot$npatients
#           rownames(npts.comb)=1:dim(npts.comb)[1]
#           colnames(npts.comb)=1:dim(npts.comb)[2]
#           barplot(npts.comb,beside=TRUE,ylab="Selection percentage (%)",
#                   xlab="Drug B", ylim=c(0,ceiling(max(npts.comb,na.rm=TRUE)*1.3)),
#                   legend.text=rownames(npts.comb),
#                   args.legend=list(title="Drug A",horiz=TRUE,x="top"))
#
#
#
#
#         }
#         else {
#           warning(paste("The variable (", strpattern,
#                         ") to be plotted cannot be found in data object: ",
#                         new.obj[1], ".\n", sep = ""))
#         }
#       }
#     }
#     else if (!is.null(objectPlot$MTD)) {
#       if (objectPlot$MTD[1] == 99) {
#         warning("All tested doses are overly toxic. No MTD is selected!\n")
#       }
#       else {
#         if (!is.null(objectPlot$p_est)) {
#           dev.new()
#           op <- par(no.readonly = TRUE)
#           dev.off()
#           par(op)
#           par(mfrow = c(1, 1), mar = c(5, 5, 4, 2))
#           if (length(objectPlot$MTD) >= 2) {
#             # if (!requireNamespace("epade", quietly = TRUE)) {
#             #   install.packages("epade", repos = "http://cran.us.r-project.org",
#             #                    dependencies = TRUE)
#             #   requireNamespace("epade")
#             #   if (!requireNamespace("epade", quietly = TRUE))
#             #     stop(paste("Package: ", "epade", " not found!!!",
#             #                sep = ""))
#             # }
#             # requireNamespace("epade", quietly = TRUE)
#             # epade::bar3d.ade(objectPlot$p_est, wall = 6,
#             #                  zw = 1, xw = 0.2, yticks = round(seq(0,
#             #                                                       1, 0.1), 1), zticks = seq(1:nrow(objectPlot$p_est)),
#             #                  xticks = seq(1:ncol(objectPlot$p_est)),
#             #                  xlab = "Drug B", ylab = "DLT rate", zlab = "Drug A",
#             #                  col = "lavender", alpha = 0.5)
#             p_est.comb=objectPlot$p_est
#             rownames(p_est.comb)=1:dim(p_est.comb)[1]
#             colnames(p_est.comb)=1:dim(p_est.comb)[2]
#             barplot(p_est.comb,beside=TRUE,ylab="DLT rate",
#                     ylim=c(0,round(max(p_est.comb,na.rm=TRUE)*1.5,1)),xlab="Drug B",legend.text=rownames(p_est.comb),
#                     args.legend=list(title="Drug A",horiz=TRUE,x="top"))
#           }
#           else {
#             p_est = objectPlot$p_est
#             p_hat = p_est[, 2]
#             ci = p_est[, 3]
#             ci = gsub("[\\(\\)]", "", ci)
#             conf.intv = matrix(unlist(strsplit(ci, ",")),
#                                byrow = TRUE, ncol = 2)
#             if (p_est[1, 2] == "----") {
#               warning("The trial is stopped since the lowest dose is too toxic.\n")
#             }
#             else {
#               numbs = ifelse(sum(p_hat == "----") ==
#                                0, length(p_hat), min(which(p_hat ==
#                                                              "----")) - 1)
#               numbs2 = length(p_hat)
#               phatx = as.numeric(as.character(p_hat[1:numbs]))
#               lwr = as.numeric(as.character(conf.intv[1:numbs,
#                                                       1]))
#               upr = as.numeric(as.character(conf.intv[1:numbs,
#                                                       2]))
#               par(mar = c(5, 5, 4, 2))
#               plot(1:numbs2, ylim = c(0, 1), xlab = "Dose level",
#                    ylab = "DLT rate", pch = "", xaxt = "n",
#                    cex.lab = 1.3)
#               axis(1, at = 1:numbs2, labels = 1:numbs2)
#               abline(h = objectPlot$target, lty = 2,
#                      col = 2)
#               points(1:numbs, phatx, pch = 19)
#               arrows(x0 = 1:numbs, x1 = 1:numbs, y0 = lwr,
#                      y1 = upr, code = 3, angle = 90, length = 0.1)
#               if (numbs < numbs2) {
#                 points((numbs + 1):numbs2, seq(min(1,
#                                                    max(phatx, na.rm = T) + 0.05), min(max(phatx,
#                                                                                           na.rm = T) + 0.2, 1), length = numbs2 -
#                                                  numbs), pch = "*", cex = 1.5)
#                 legend("topleft", "*   no patient treated")
#               }
#             }
#           }
#         }
#         else {
#           warning("Please set verbose=TRUE to get more details of the results.\n")
#         }
#       }
#     }
#     else {
#       warning("Please double check and specify the variable to be plotted...\n")
#     }
#   }
# }
