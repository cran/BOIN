get.oc <- function(target, p.true, ncohort, cohortsize, n.earlystop=100, startdose=1, p.saf="default", p.tox="default", cutoff.eli=0.95, extrasafe=FALSE, offset=0.05, ntrial=1000)
{   	
## if the user does not provide p.saf and p.tox, set them to the default values
	if(p.saf=="default") p.saf=0.6*target;
	if(p.tox=="default") p.tox=1.4*target;
	
## simple error checking
	if(target<0.05) {cat("Error: the target is too low! \n"); return();}
	if(target>0.6)  {cat("Error: the target is too high! \n"); return();}
	if((target-p.saf)<(0.1*target)) {cat("Error: the probability deemed safe cannot be higher than or too close to the target! \n"); return();}
	if((p.tox-target)<(0.1*target)) {cat("Error: the probability deemed toxic cannot be lower than or too close to the target! \n"); return();}
	if(offset>=0.5) {cat("Error: the offset is too large! \n"); return();}
    if(n.earlystop<=6) {cat("Warning: the value of n.earlystop is too low to ensure good operating characteristics. Recommend n.earlystop = 9 to 18 \n"); return();}
	
	set.seed(6);
	ndose=length(p.true);	
	npts = ncohort*cohortsize;
	Y=matrix(rep(0, ndose*ntrial), ncol=ndose); # store toxicity outcome
	N=matrix(rep(0, ndose*ntrial), ncol=ndose); # store the number of patients
	dselect = rep(0, ntrial); # store the selected dose level
	
## obtain dose escalation and deescalation boundaries
	temp=get.boundary(target, ncohort, cohortsize, n.earlystop, p.saf, p.tox, cutoff.eli, extrasafe, print=FALSE);
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
            
            if(n[d]>=n.earlystop) break;
			
## determine if the current dose should be eliminated
			if(!is.na(b.elim[n[d]]))
			{
				if(y[d]>=b.elim[n[d]]) 
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
			if(y[d]<=b.e[n[d]] && d!=ndose) { if(elimi[d+1]==0) d=d+1; }
			else if(y[d]>=b.d[n[d]] && d!=1) { d=d-1; }
			else { d=d; }
		}
		Y[trial,]=y;
		N[trial,]=n;
		if(earlystop==1) { dselect[trial]=99; }
		else  { dselect[trial]=select.mtd(target, n, y, cutoff.eli, extrasafe, offset, print=FALSE); }
	}
  	
# output results 
	selpercent=rep(0, ndose);
    nptsdose = apply(N,2,mean);
    ntoxdose = apply(Y,2,mean);
    
	for(i in 1:ndose) { selpercent[i]=sum(dselect==i)/ntrial*100; }
	cat("selection percentage at each dose level (%):\n");
	cat(formatC(selpercent, digits=1, format="f"), sep="  ", "\n");
	cat("number of patients treated at each dose level:\n");
	cat(formatC(nptsdose, digits=1, format="f"), sep ="  ", "\n");
	cat("number of toxicity observed at each dose level:\n");
	cat(formatC(ntoxdose, digits=1, format="f"), sep ="  ", "\n");
	cat("average number of toxicities:", formatC(sum(Y)/ntrial, digits=1, format="f"), "\n");
	cat("average number of patients:", formatC(sum(N)/ntrial, digits=1, format="f"), "\n");
	cat("percentage of early stopping due to toxicity:", formatC(sum(dselect==99)/ntrial*100, digits=1, format="f"), "% \n");
	
	if(length(which(p.true==target))>0) # if MTD exist
	{
		cat("risk of poor allocation:", formatC(mean(N[, p.true==target]<npts/ndose)*100, digits=1, format="f"), "% \n");
		cat("risk of high toxicity:", formatC(mean(rowSums(Y)>target*npts)*100, digits=1, format="f"), "% \n");
	}
    
    invisible(data.frame(target=target, p.true=p.true, ncohort=ncohort, cohortsize = cohortsize,
    ntotal=ncohort*cohortsize, startdose = startdose, p.saf = p.saf, p.tox = p.tox,  
    cutoff.eli = cutoff.eli, extrasafe = extrasafe, offset = offset, ntrial = ntrial,
    dose=1:ndose, selpercent=selpercent, nptsdose=nptsdose, ntoxdose=ntoxdose, totaltox=sum(Y)/ntrial, totaln=sum(N)/ntrial, pctearlystop=sum(dselect== 99)/ntrial * 100))
}