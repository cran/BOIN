\name{BOIN-package}
\alias{BOIN-package}
\alias{BOIN}
\docType{package}
\title{
Bayesian Optimum INterval (BOIN) Design for Finding the Maximum Tolerated Dose (MTD) for Single-Agent and Drug-Combination Trials
}
\description{
The package is used to design single-agent or drug-combination phase I clinical trials using the BOIN design. The BOIN design is motivated by the top priority and concern of clinicians, which is to effectively treat patients and minimize the chance of exposing them to subtherapeutic or overly toxic doses. The prominent advantage of the BOIN design is that it achieves simplicity and a superior performance at the same time. The BOIN design is algorithm-based and can be implemented in a simple way similar to the traditional ``3+3" design. The BOIN design yields an average performance that is comparable to that of the continual reassessment method (CRM, one of the best model-based designs) in terms of selecting the MTD, but has a substantially lower risk of assigning patients to subtherapeutic or overly toxic doses. 
}
\details{
\tabular{ll}{
Package: \tab BOIN\cr
Type: \tab Package\cr
Version: \tab 2.0\cr
Date: \tab 2015-08-20\cr
License: \tab GPL-2\cr
}
%~~ An overview of how to use the package, including the ~~
%~~ most important functions ~~
}
\author{
Suyu Liu < syliu@mdanderson.org > and Ying Yuan < yyuan@mdanderson.org >
%Maintainer: %Who to complain to <yourfault@somewhere.net>
%~~ The author and/or maintainer of the package ~~
}
\references{
%% ~put references to the literature/web site here ~
Liu S. and Yuan, Y. (2015). Bayesian Optimal Interval Designs for Phase I Clinical Trials, \emph{Journal of the Royal Statistical Society: Series C}, 64, 507-523.

Lin R. and Yin, G. (2015). Bayesian Optimal Interval Designs for Dose Finding in Drug-combination Trials, \emph{Statistical Methods in Medical Research}, to appear.
}

%~~ Optionally other standard keywords, one per line, from ~~
%~~ file KEYWORDS in the R documentation directory ~~
%~~ \keyword{ package }
%~~ \keyword{ decision }
%~~ \keyword{ error }
%~~ \keyword{ adaptive}
%~~ \keyword{ design }
%~~ \keyword{ dose }
%~~ \keyword{ finding }
%~~ \keyword{ maximum }
%~~ \keyword{ tolerance }
\seealso{
%~~ Optional links to other man pages, e.g. ~~
%~~ \code{\link[<pkg>:<pkg>-package]{<pkg>}} ~~
Tutorial: \url{http://odin.mdacc.tmc.edu/~yyuan/Software/BOIN/BOIN2.0_tutorial.pdf}

Paper: \url{http://odin.mdacc.tmc.edu/~yyuan/Software/BOIN/paper.pdf}
}
\examples{
###############################################
##   Example 1. Single-agent Phase I Trial   ##	
###############################################

# Obtain operating characteristics
get.oc(target=0.3, p.true=c(0.05, 0.15, 0.3, 0.45, 0.6), ncohort=10, 
cohortsize=3, ntrial=1000)

# Obtain dose escalation boundaries for trial conduct
get.boundary(target=0.3, ncohort=10, cohortsize=3)	
	
# Select the MTD when the trial is completed
n<-c(3, 3, 15, 9, 0)
y<-c(0, 0, 4, 4, 0)
select.mtd(target=0.3, npts=n, ntox=y)

#########################################################
##   Example 2. A 3x5 Drug-combination Phase I Trial   ##
#########################################################

# Obtain operating characteristics
p.true = matrix(c(0.05,0.10,0.15,0.30,0.45, 0.10,0.15,0.30,0.45,0.55, 
0.15,0.30,0.45,0.50,0.60), ncol=5, byrow=TRUE)
get.oc.comb(target=0.3, p.true, ncohort=40, cohortsize=1, ntrial=1000)

# Obtain dose escalation boundaries for trial conduct
get.boundary(target=0.3, ncohort=40, cohortsize=1)	

# make the decision of dose escalation/deescalation during the course of trial conduct
# matrix n contains the number of patients treated at each dose combination
# matrix y contains the number of patients experienced toxicity at each dose combination
n<-matrix(c(3, 0, 0, 0, 0,   7, 6, 0, 0, 0,  0, 0, 0, 0, 0), ncol=5, byrow=TRUE)     
y<-matrix(c(0, 0, 0, 0, 0,   1, 1, 0, 0, 0,  0, 0, 0, 0, 0), ncol=5, byrow=TRUE)  
conduct.comb(target=0.3, npts=n, ntox=y, dose.curr=c(2, 2))

# Select the MTD when the combination trial is completed
n<-matrix(c(3, 5, 0, 0, 0,   7, 6, 15, 0, 0,  0, 0, 4, 0, 0), ncol=5, byrow=TRUE)    
y<-matrix(c(0, 1, 0, 0, 0,   1, 1, 4, 0, 0,   0, 0, 2, 0, 0), ncol=5, byrow=TRUE)    
select.mtd.comb(target=0.3, npts=n, ntox=y)
}
