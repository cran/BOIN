\name{select.mtd.comb}
\alias{select.mtd.comb}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Select the maximum tolerated dose (MTD) for drug combination trials
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
\code{select.mtd.comb} is used to select the maximum tolerated dose (MTD) when the drug combination trial is completed. 
}
\usage{
select.mtd.comb(target, npts, ntox, cutoff.eli=0.95, extrasafe=FALSE, 
offset=0.05, print=TRUE)
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
a matrix containing the number of patients experienced dose-limiting toxicity at each dose combination
}
  \item{cutoff.eli}{
%%     ~~Describe \code{cutoff.eli} here~~
the cutoff to eliminate overly toxic doses for safety. We recommend the default value of \code{(cutoff.eli=0.95)} for general use. 
}
  \item{extrasafe}{
%%     ~~Describe \code{ntrial} here~~
set \code{extrasafe=TRUE} to impose a more strict stopping rule for extra safety
}
  \item{offset}{
%%     ~~Describe \code{ntrial} here~~
a small positive number (between 0 and 0.5) to control how strict the stopping rule is when \code{extrasafe=TRUE}. A larger value leads to a more strict stopping rule. The default value \code{offset=0.05} generally works well.
}
  \item{print}{
%%     ~~Describe \code{cutoff.eli} here~~
prints out the dose selection result. 
}
}
\details{
\code{select.mtd.comb()} selects the MTD based on isotonic estimates of toxicity probabilities. 
\code{select.mtd.comb()} selects as the MTD dose \emph{j*}, for which the isotonic estimate of the toxicity rate is closest to the \code{target}. If there are ties, we select from the ties the highest dose level when the estimate of the toxicity rate is smaller than the \code{target}, or the lowest dose level when the estimate of the toxicity rate is greater than the \code{target}. The (matrix) isotonic estimates are obtained by the pooled-adjacent-violators algorithm (PAVA) (Barlow, 1972).
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
\code{select.mtd.comb()} returns the MTD based on the trial data.
}

\references{
%% ~put references to the literature/web site here ~
Liu S. and Yuan, Y. (2015). Bayesian Optimal Interval Designs for Phase I Clinical Trials, \emph{Journal of the Royal Statistical Society: Series C}, 64, 507-523.

Lin R. and Yin, G. (2015). Bayesian Optimal Interval Designs for Dose Finding in Drug-combination Trials, \emph{Statistical Methods in Medical Research}, to appear.
}

\author{
Suyu Liu and Ying Yuan (yyuan@mdanderson.org)
}
\note{
The MTD selection and dose escalation/deescalation rule are two independent components of the trial design. When appropriate, another dose selection procedure (e.g., based on a fitted logistic model) can be used to select the MTD after the completion of the trial using the BOIN design.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
Tutorial: \url{http://odin.mdacc.tmc.edu/~yyuan/Software/BOIN/BOIN2.0_tutorial.pdf}

Paper: \url{http://odin.mdacc.tmc.edu/~yyuan/Software/BOIN/paper.pdf}
}

\examples{
## Select the MTD based on the data from a 3x5 combination trial
## matrix n contains the number of patients treated at each dose combination
## matrix y contains the number of patients experienced toxicity at each dose combination

n<-matrix(c(3, 5, 0, 0, 0, 7, 6, 15, 0, 0, 0, 0, 4, 0, 0), ncol=5, byrow=TRUE)    
y<-matrix(c(0, 1, 0, 0, 0, 1, 1, 4, 0, 0, 0, 0, 2, 0, 0), ncol=5, byrow=TRUE)     
select.mtd.comb(target=0.3, npts=n, ntox=y)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{selection of the MTD}
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
