get.oc <- function(target, p.true, ncohort, cohortsize=3, startdose=1, p.saf="default", p.tox="default", design=1, cutoff.eli=0.95, extrasafe=FALSE, offset=0.05, ntrial=1000)
{   	
## if the user does not provide p.saf and p.tox, set them to the default values
	if(p.saf=="default") p.saf=0.6*target;
	if(p.tox=="default") p.tox=1.4*target;
	
## simple error checking
	if(target<0.05) {cat("Error: the target is too low! \n"); return(1);}
	if(target>0.6)  {cat("Error: the target is too high! \n"); return(1);}
	if((target-p.saf)<(0.1*target)) {cat("Error: the probability deemed safe cannot be higher than or too close to the target! \n"); return(1);}
	if((p.tox-target)<(0.1*target)) {cat("Error: the probability deemed toxic cannot be lower than or too close to the target! \n"); return(1);}
	if(offset>=0.5) {cat("Error: the offset is too large! \n"); return(1);}
	
	set.seed(6);
	ndose=length(p.true);	
	npts = ncohort*cohortsize;
	Y=matrix(rep(0, ndose*ntrial), ncol=ndose); # store toxicity outcome
	N=matrix(rep(0, ndose*ntrial), ncol=ndose); # store the number of patients
	dselect = rep(0, ntrial); # store the selected dose level
	
## obtain dose escalation and deescalation boundaries
	temp=get.boundary(target, ncohort, cohortsize, p.saf, p.tox, design, cutoff.eli, extrasafe, print=FALSE); 	
	b.e=temp[2,];   # escalation boundary
	b.d=temp[3,];   # deescalation boundary
	b.elim=temp[4,];  # elimination boundary
	
	
################## simulate trials ###################
	for(trial in 1:ntrial)
	{
		y<-rep(0, ndose);    ## the number of DLT at each dose level
		n<-rep(0, ndose);    ## the number of patients treated at each dose level
		earlystop=0;         ## indiate whether the trial terminates early
		d=startdose;         ## starting dose level
		elimi = rep(0, ndose);  ## indicate whether doses are eliminated
		
		for(i in 1:ncohort)  
		{  			
### generate toxicity outcome
			y[d] = y[d] + sum(runif(cohortsize)<p.true[d]);
			n[d] = n[d] + cohortsize;
			nc = n[d]/cohortsize;
			
## determine if the current dose should be eliminated
			if(!is.na(b.elim[nc]))
			{
				if(y[d]>=b.elim[nc]) 
				{      
					elimi[d:ndose]=1;
					if(d==1) {earlystop=1; break;} 
				}
## implement the extra safe rule by decreasing the elimination cutoff for the lowest dose
				if(extrasafe)  
				{
					if(d==1 && y[1]>=3)
					{
						if(1-pbeta(target, y[1]+1, n[1]-y[1]+1)>cutoff.eli-offset) {earlystop=1; break;}
					}
				}
			}
			
## dose escalation/de-escalation
			if(y[d]<=b.e[nc] && d!=ndose) { if(elimi[d+1]==0) d=d+1; }
			else if(y[d]>=b.d[nc] && d!=1) { d=d-1; }
			else { d=d; }
		}
		Y[trial,]=y;
		N[trial,]=n;
		if(earlystop==1) { dselect[trial]=99; }
		else  { dselect[trial]=select.mtd(target, y, n, cutoff.eli, print=FALSE); }
	}
  	
# output results 
	selpercent=rep(0, ndose);
	for(i in 1:ndose) { selpercent[i]=sum(dselect==i)/ntrial*100; }
	cat("selection percentage at each dose level (%):\n");
	cat(formatC(selpercent, digits=1, format="f"), sep="  ", "\n");
	cat("number of patients treated at each dose level:\n");
	cat(formatC(apply(N,2,mean), digits=1, format="f"), sep ="  ", "\n");
	cat("number of toxicity observed at each dose level:\n");
	cat(formatC(apply(Y,2,mean), digits=1, format="f"), sep ="  ", "\n");
	cat("average number of toxicities:", formatC(sum(Y)/ntrial, digits=1, format="f"), "\n");
	cat("average number of patients:", formatC(sum(N)/ntrial, digits=1, format="f"), "\n");
	cat("percentage of early stopping due to toxicity:", formatC(sum(dselect==99)/ntrial*100, digits=1, format="f"), "% \n");
	
	if(length(which(p.true==target))>0) # if MTD exist
	{
		cat("risk of poor allocation:", formatC(mean(N[, p.true==target]<npts/ndose)*100, digits=1, format="f"), "% \n");
		cat("risk of high toxicity:", formatC(mean(rowSums(Y)>target*npts)*100, digits=1, format="f"), "% \n");
	}
	else
	{
		cat("risk of poor allocation: NA (no dose has the toxicity rate =",  target, ")\n");
		cat("risk of high toxicity:: NA (no dose has the toxicity rate =", target, ")\n");
	}
}