\name{conduct.comb}
\alias{conduct.comb}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Determine the dose for the next cohort of new patients during the trial conduct for drug combination trials 
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Determine the dose for the next cohort of new patients during the trial conduct for drug combination trials 
}
\usage{
conduct.comb(target, npts, ntox, dose.curr, n.earlystop=100,  p.saf="default", 
p.tox="default", cutoff.eli=0.95, extrasafe=FALSE, offset=0.05)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{target}{
%%     ~~Describe \code{target} here~~
target toxicity rate
}
  \item{npts}{
%%     ~~Describe \code{npts} here~~
a matrix containing the number of patients treated at each dose combination
}
  \item{ntox}{
%%     ~~Describe \code{ntox} here~~
a matrix containing the number of patients who experienced dose-limiting toxicities at each dose combination
}
  \item{dose.curr}{
%%     ~~Describe \code{cohortsize} here~~
the current dose combination
}
 \item{n.earlystop}{
%%     ~~Describe \code{ndmax} here~~
Early stopping parameter. If the number of patients treated at the current dose reaches \code{n.earlystop}, stop the trial and select the MTD based on the observed data. The default value \code{n.earlystop=100} essentially turns off this type of early stopping. 
}
  \item{p.saf}{
%%     ~~Describe \code{p.saf} here~~
the highest toxicity probability that is deemed subtherapeutic (i.e. below the MTD) such that dose escalation should be undertaken. The default value is \code{p.saf=0.6 x target}.}
  \item{p.tox}{
%%     ~~Describe \code{p.tox} here~~
the lowest toxicity probability that is deemed overly toxic such that deescalation is required. The default value is \code{p.tox=1.4 x target}.
}
  \item{cutoff.eli}{
%%     ~~Describe \code{cutoff.eli} here~~
the cutoff to eliminate an overly toxic dose for safety. We recommend the default value of \code{(cutoff.eli=0.95)} for general use 
}
  \item{extrasafe}{
%%     ~~Describe \code{ntrial} here~~
set \code{extrasafe=TRUE} to impose a more stringent stopping rule
}
  \item{offset}{
%%     ~~Describe \code{ntrial} here~~
a small positive number (between 0 and 0.5) to control how strict the stopping rule is when \code{extrasafe=TRUE}. A larger value leads to a more strict stopping rule. The default value \code{offset=0.05} generally works well.
}
}
\details{
This function is used for conducting combination trials. Given the currently observed data, 
\code{conduct.comb()} determines dose combination for treating the next cohort of new patients. The currently observed data include: the number of patients treated at each dose combination (i.e., \code{npts}), the number of patients who experienced dose-limiting toxicities at each dose combination (i.e.,\code{ntox}), and the level of current dose (i.e., \code{dose}).
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
\code{conduct.comb()} returns the dose for treating the next cohort of new patients. 
}

\references{
%% ~put references to the literature/web site here ~
Liu S. and Yuan, Y. (2015). Bayesian Optimal Interval Designs for Phase I Clinical Trials, \emph{Journal of the Royal Statistical Society: Series C}, 64, 507-523.

Lin R. and Yin, G. (2015). Bayesian Optimal Interval Designs for Dose Finding in Drug-combination Trials, \emph{Statistical Methods in Medical Research}, to appear.
}

\author{
Suyu Liu and Ying Yuan (yyuan@mdanderson.org)
}


\seealso{
Tutorial: \url{http://odin.mdacc.tmc.edu/~yyuan/Software/BOIN/BOIN2.0_tutorial.pdf}

Paper: \url{http://odin.mdacc.tmc.edu/~yyuan/Software/BOIN/paper.pdf}
}

\examples{
## Consider a 3x5 drug combination phase I trial aiming to find the MTD with 
## a target toxicity rate of 0.3. Assume that the current dose is (2, 2), 
## matrix n contains the number of patients treated at each combinations, 
## and matrix y the number of patients experienced dose-limiting toxicity 
## at each combinations.

n<-matrix(c(3, 0, 0, 0, 0,   7, 6, 0, 0, 0,  0, 0, 0, 0, 0), ncol=5, byrow=TRUE)    
y<-matrix(c(0, 0, 0, 0, 0,   1, 1, 0, 0, 0,  0, 0, 0, 0, 0), ncol=5, byrow=TRUE)    
conduct.comb(target=0.3, npts=n, ntox=y, dose.curr=c(2, 2))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ operating characteristics }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
