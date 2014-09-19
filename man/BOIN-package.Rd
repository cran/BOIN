\name{BOIN-package}
\alias{BOIN-package}
\alias{BOIN}
\docType{package}
\title{
Bayesian Optimum Interval (BOIN) design for finding the maximum tolerated dose (MTD)
}
\description{
Package is used to design phase I clinical trials using the BOIN design. The BOIN design is motivated by top priority and concern of clinicians, which is to effectively treat patients and minimize the chance of exposing them to subtherapeutic or overly toxic doses. The prominent advantage of the BOIN design is that it achieves simplicity and superior performance at the same time. The BOIN design is algorithm-based and can be implemented in a simple way similar to the traditional ``3+3" design. The BOIN design yields average performance comparable to the continual reassessment method (CRM, one of the best model-based designs) in terms of selecting the MTD, but has a substantially lower risk of assigning patients to subtherapeutic or overly toxic doses. 
}
\details{
\tabular{ll}{
Package: \tab BOIN\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2014-06-05\cr
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
Liu S. and Yuan, Y. (2014). Bayesian Optimal Interval Designs for Phase I Clinical Trials, \emph{Journal of the Royal Statistical Society: Series C}, to appear
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
Tutorial: \url{http://odin.mdacc.tmc.edu/~yyuan/Software_release/OptInterval/tutorial.pdf}

Paper: \url{http://odin.mdacc.tmc.edu/~yyuan/Software_release/OptInterval/paper.pdf}
}
\examples{

	
###### Trial Design #####
get.boundary(target=0.3, ncohort=10, cohortsize=3)	
	
###### Obtain operating characteristics #####
get.oc(target=0.3, p.true=c(0.05, 0.15, 0.3, 0.45, 0.6), ncohort=10, cohortsize=3, ntrial=1000)

##### Select the MTD when the trial is completed #####
n<-c(3, 3, 15, 9, 0)
y<-c(0, 0, 4, 4, 0)
select.mtd(target=0.3, ntox=y, npts=n)


}
