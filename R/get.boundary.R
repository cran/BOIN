get.boundary <- function(target, ncohort, cohortsize, n.earlystop=100, p.saf="default", p.tox="default", design=1, cutoff.eli=0.95, extrasafe=FALSE, offset=0.05, print=TRUE)
{
	density1 <- function(p, n, m1, m2) {pbinom(m1, n, p)+1-pbinom(m2-1, n, p);}
	density2 <- function(p, n, m1) {1-pbinom(m1, n, p);}
	density3 <- function(p, n, m2) {pbinom(m2-1, n, p);}
	
# if the user does not provide p.saf and p.tox, use the default values
	if(p.saf=="default") p.saf=0.6*target;
	if(p.tox=="default") p.tox=1.4*target;
	
### simple error checking
	if(target<0.05) {cat("Error: the target is too low! \n"); return();}
	if(target>0.6)  {cat("Error: the target is too high! \n"); return();}
	if((target-p.saf)<(0.1*target)) {cat("Error: the probability deemed safe cannot be higher than or too close to the target! \n"); return();}
	if((p.tox-target)<(0.1*target)) {cat("Error: the probability deemed toxic cannot be lower than or too close to the target! \n"); return();}
	if(offset>=0.5) {cat("Error: the offset is too large! \n"); return();}
	if(n.earlystop<=6) {cat("Warning: the value of n.earlystop is too low to ensure good operating characteristics. Recommend n.earlystop = 9 to 18 \n"); return();}
    
### numerical search for the boundaries that minimize decision errors of dose escalation/deescalation
	npts = ncohort*cohortsize;
	ntrt=NULL; b.e=NULL; b.d=NULL; elim=NULL;
	for(n in 1:npts)
	{
		error.min=3;
		for(m1 in 0:(n-1))
		{
			for(m2 in (m1+1):n)
			{
				if(design==2)
				{
					error1 = integrate(density1, lower=p.saf, upper=p.tox, n, m1, m2)$value/(p.tox-p.saf);
					error2 = integrate(density2, lower=0, upper=p.saf, n, m1)$value/p.saf;
					error3 = integrate(density3, lower=p.tox, upper=1, n, m2)$value/(1-p.tox);
				}
				else
				{
					error1 = pbinom(m1, n, target)+1-pbinom(m2-1, n, target);
					error2 = 1-pbinom(m1, n, p.saf);
					error3 = pbinom(m2-1, n, p.tox);
				}
				error=error1+error2+error3;
				if(error<error.min) {error.min=error; cutoff1=m1; cutoff2=m2;}
			}
		}
		ntrt = c(ntrt, n);
		b.e = c(b.e, cutoff1);
		b.d = c(b.d, cutoff2);
		
		elimineed=0; # indicating whether elimination is needed
		if(n<3) { elim = c(elim, NA); }  # require treating at least 3 patients before eliminating a dose
		else
		{
			for(ntox in 1:n) #determine elimination boundary, prior beta(1,1) is used in beta-binomial model
			{
				if(1-pbeta(target, ntox+1, n-ntox+1)>cutoff.eli) {elimineed=1; break;}
			}
			if(elimineed==1) { elim = c(elim, ntox); }
			else { elim = c(elim, NA); } # set the elimination boundary large such that no elimination will actually occurs
		}
	}
	for(i in 1:length(b.d)) { if(!is.na(elim[i]) && (b.d[i]>elim[i])) b.d[i]=elim[i]; }
	boundaries = rbind(ntrt, b.e, b.d, elim)[,1:min(npts, n.earlystop)];
	rownames(boundaries) = c("Number of patients treated", "Escalate if # of DLT <=",
							 "Deescalate if # of DLT >=", "Eliminate if # of DLT >=" );
	colnames(boundaries) = rep("", min(npts, n.earlystop));
	
	if(print)
    {
        if(design==1)
        {
            lambda1  = log((1-p.saf)/(1-target))/log(target*(1-p.saf)/(p.saf*(1-target)));
            lambda2  = log((1-target)/(1-p.tox))/log(p.tox*(1-target)/(target*(1-p.tox)));
			cat("Escalate dose if the observed toxicity rate at the current dose <= ", lambda1, "\n");
			cat("Deescalate dose if the observed toxicity rate at the current dose >= ", lambda2, "\n\n");
			cat("This is equivalent to the following decision boundaries\n");
        }
        if(design==2) { cat("The decision boundaries for the global BOIN design are given by \n"); }
        print(boundaries[, (1:floor(min(npts, n.earlystop)/cohortsize))*cohortsize]);
        
        if(cohortsize>1)
        {
            cat("\n"); cat("A more completed version of the decision boundaries is given by\n");
            print(boundaries);
        }
        
        cat("\n");
        if(!extrasafe) cat("Default stopping rule: stop the trial if the lowest dose is eliminated.\n");
	}
	
	if(extrasafe)
    {
        stopbd=NULL;
		ntrt=NULL;
        for(n in 1:npts)
        {
            ntrt = c(ntrt, n);
            if(n<3) { stopbd = c(stopbd, NA); }  # require treating at least 3 patients before stop a trial
			else
			{
				for(ntox in 1:n) #determine stopping boundary, prior beta(1,1) is used in beta-binomial model
				{
					if(1-pbeta(target, ntox+1, n-ntox+1)>cutoff.eli-offset) {stopneed=1; break;}
				}
				if(stopneed==1) { stopbd = c(stopbd, ntox); }
				else { stopbd = c(stopbd, NA); } 
			}
        }
		
		stopboundary = rbind(ntrt, stopbd)[, 1:min(npts, n.earlystop)];
		rownames(stopboundary) = c("Number of patients treated at the lowest dose  ", "stop the trial if # of DLT >=        ");
		colnames(stopboundary) = rep("", min(npts, n.earlystop));
		
		if(print)
		{
			cat("\n");
			cat("In addition to the default stopping rule (i.e., stop the trial if the lowest dose is eliminated), \n"); 
            cat("the following more strict stopping safety rule will be used for extra safety: \n")
			cat(" stop the trial if (1) the number of patients treated at the lowest dose >= 3 AND", "\n", 
				"(2) Pr(the toxicity rate of the lowest dose >", target, "| data) > ", cutoff.eli-offset, ",\n",
				"which corresponds to the following stopping boundaries:\n")
			print(stopboundary);
		}
	}																														
	if(!print) return(boundaries);
}
