\name{select.mtd}
\alias{select.mtd}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Select the maximum tolerated dose (MTD)
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
\code{select.mtd} is used to select the maximum tolerated dose (MTD) when the trial is completed. 
}
\usage{
select.mtd(target, ntox, npts, cutoff.eli = 0.95, print=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{target}{
%%     ~~Describe \code{target} here~~
target toxicity rate
}
  \item{ntox}{
%%     ~~Describe \code{ntox} here~~
a vector containing the number of patients experienced dose-limiting toxicity at each dose level
}
  \item{npts}{
%%     ~~Describe \code{npts} here~~
a vector containing the number of patients treated at each dose level
}
  \item{cutoff.eli}{
%%     ~~Describe \code{cutoff.eli} here~~
the cutoff to eliminate overly toxic doses for safety. We recommend the default value of \code{(cutoff.eli=0.95)} for general use. 
}
  \item{print}{
%%     ~~Describe \code{cutoff.eli} here~~
prints out the dose selection result. 
}
}
\details{
\code{select.mtd} selects the MTD based on isotonic estimates of toxicity probabilities. \code{select.mtd} selects as the MTD dose \emph{j*}, for which the isotonic estimate of toxicity rate is closet to the \code{target}. If there are ties, we select from ties the highest dose level when the estimate of toxicity rate is smaller than \code{target}, or the lowest dose level when the estimate of toxicity rate is greater than \code{target}. The isotonic estimates are obtained by the pooled-adjacent-violators algorithm (PAVA) (Barlow, 1972).
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
\code{select.mtd} returns the MTD based on the trial data.
}

\references{
%% ~put references to the literature/web site here ~
Liu S. and Yuan, Y. (2014). Bayesian Optimal Interval Designs for Phase I Clinical Trials, \emph{Journal of the Royal Statistical Society: Series C}, to appear
}

\author{
Suyu Liu and Ying Yuan (yyuan@mdanderson.org)
}
\note{
The MTD selection and dose escalation/deescalation rule are two independent components of the trial design. When appropriate, other dose selection procedure (e.g., based on a fitted logistic model) can also be used to select the MTD after the completion of the trial using the BOIN design.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
Tutorial: \url{http://odin.mdacc.tmc.edu/~yyuan/Software_release/OptInterval/tutorial.pdf}

Paper: \url{http://odin.mdacc.tmc.edu/~yyuan/Software_release/OptInterval/paper.pdf}
}

\examples{
##### Select the MTD based on the trial data #####
n<-c(3, 3, 15, 9, 0)    #the number of patients treated at 5 investigational doses
y<-c(0, 0, 4, 4, 0)     #the number of patients experienced toxicity at 5 doses
select.mtd(target=0.3, ntox=y, npts=n)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{selection of the MTD}
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
