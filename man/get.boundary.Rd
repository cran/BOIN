\name{get.boundary}
\alias{get.boundary}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Generate dose escalation and deescalation boundaries
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Generate the optimal dose escalation and deescalation boundaries for conducting the trial.
%The output is presented in two forms: based on the observed toxicity rate and based on the observed dose-limiting toxicity (DLT). The latter is recommended for easier application by clinical researchers.
}
\usage{
get.boundary(target, ncohort, cohortsize = 3, p.saf = "default", p.tox = "default", 
design = 1, cutoff.eli = 0.95, extrasafe = FALSE, offset = 0.05, print = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{target}{
%%     ~~Describe \code{target} here~~
target toxicity rate
}
  \item{ncohort}{
%%     ~~Describe \code{ncohort} here~~
the total number of cohort
}
  \item{cohortsize}{
%%     ~~Describe \code{cohortsize} here~~
the cohort size
}
  \item{p.saf}{
%%     ~~Describe \code{p.saf} here~~
the highest toxicity probability that is deemed subtherapeutic (i.e. below the MTD) such that dose escalation should be made.The default value is \code{p.saf=0.6xtarget}. 
}
  \item{p.tox}{
%%     ~~Describe \code{p.tox} here~~
the loweset toxicity probability that is deemed overly toxic such that deescalation is required. The default value is \code{p.tox=1.4xtarget.}
}
  \item{design}{
%%     ~~Describe \code{design} here~~
to use the local optimal design (the default, \code{design=1}) or the global optimal design (\code{design=2}). We generally recommend the local optimal design.
}
  \item{cutoff.eli}{
%%     ~~Describe \code{cutoff.eli} here~~
the cutoff to eliminate the overly toxic dose for safety. We recommend the default value of \code{(cutoff.eli=0.95)} for general use. 
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
%%     ~~Describe \code{print} here~~
prints out the boundary results.
}
}

\details{
%%  ~~ If necessary, more details than the description above ~~
The dose escalation and deescalation boundaries are all we need to run a phase I trial when using the BOIN design. The decision of which dose to administer to the next cohort of patients does not require complicated computations, but only a simple comparison of the observed toxicity rate at the current dose with the dose escalation and deescalation boundaries. If the observed toxicity rate at the current dose is smaller than or equal to the escalation boundary, we escalate the dose; if the observed toxicity rate at the current dose is greater than or equal to the deescalation boundary, we deescalate the dose; otherwise,
we retain the current dose. The dose escalation and deescalation boundaries are chosen to minimize the probability of assigning patients to subtherapeutic or overly toxic doses, thereby optimizing patient ethics. 

\code{get.boundary} also outputs elimination boundary, which is used to prevent treating patients at overly toxic doses based on the following Bayesian safety rule: 

if \eqn{pr(p_j>\phi | m_j,n_j)>0.95} and \eqn{n_j \geq 3}, dose levels \emph{j} and higher are eliminated from the trial,

where \eqn{p_j} is the toxicity probability of dose level \emph{j}, \eqn{\phi} is the target toxicity rate, \eqn{m_j} and \eqn{n_j} are the number of toxicities and patients treated at dose level \emph{j}. The trial is terminated if the lowest dose is eliminated.

The BOIN design has a built-in stopping rule, i.e., stop the trial if the lowest dose is eliminated due to toxicity. For some applications, investigators may prefer a more strict stopping rule for extra safety when the lowest dose is overly toxic. Setting \code{extrasafe=TRUE} imposes the following more strict stopping rule: stop the trial if (1) the number of patients treated at the lowest dose >=3, and (2) Pr(toxicity rate of the lowest dose > \code{target} | data) > \code{cutoff.eli-offset}. As a tradeoff, the strong stopping rule will decrease the selection percentage of the MTD when the lowest dose actually is the MTD.
}

\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
\code{get.boundary} returns the optimal dose escalation and deescalation boundaries for running the trial. The dose elimination boundary is also returned for preventing continuously exposing patients to overly toxic doses.
}

\references{
%% ~put references to the literature/web site here ~
Liu S. and Yuan, Y. (2014). Bayesian Optimal Interval Designs for Phase I Clinical Trials, \emph{Journal of the Royal Statistical Society: Series C}, to appear
}

\author{
Suyu Liu and Ying Yuan (yyuan@mdanderson.org)
}

\note{
%%  ~~further notes~~
We recommend to use the local BOIN design (i.e., the default option) because of its better operating characteristics and simplicity. We should avoid setting the values of \code{p.saf} and \code{p.tox} very close to \code{target}. This is because the small sample sizes of typical phase I trials prevent us from differentiating the target toxicity rate from the rates close to it. In addition, in most clinical applications, the target toxicity rate is often a rough guess, and finding a dose level with a toxicity rate reasonably close to the target rate will still prove to be of interest to the investigator. The default values provided by \code{get.boundary} are generally reasonable for most clinical applications.  
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
Tutorial: \url{http://odin.mdacc.tmc.edu/~yyuan/Software_release/OptInterval/tutorial.pdf}

Paper: \url{http://odin.mdacc.tmc.edu/~yyuan/Software_release/OptInterval/paper.pdf}
}

\examples{
## Consider a phase I trial aiming to find the MTD with a target toxicity rate of 0.3 
## the maximum sample size is 30 patients in cohort size of 3

get.boundary(target=0.3, ncohort=10, cohortsize=3)
%##---- Should be DIRECTLY executable !! ----
%##-- ==>  Define data, use random,
%##--	or do  help(data=index)  for the standard data sets.
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ dose escalation/deescalation boundaries }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
