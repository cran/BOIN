####### tested.  8/20/2015

conduct.comb <- function(target, npts, ntox, dose.curr, n.earlystop=100, p.saf="default", p.tox="default", cutoff.eli=0.95, extrasafe=FALSE, offset=0.05)
{
## if the user does not provide p.saf and p.tox, set them to the default values
    if(p.saf=="default") p.saf=0.6*target;
    if(p.tox=="default") p.tox=1.4*target;
    
## simple error checking
    if(npts[dose.curr[1], dose.curr[2]]==0)  {cat("Error: dose entered is not the current dose \n"); return(1);}
    if(target<0.05) {cat("Error: the target is too low! \n"); return(1);}
    if(target>0.6)  {cat("Error: the target is too high! \n"); return(1);}
    if((target-p.saf)<(0.1*target)) {cat("Error: the probability deemed safe cannot be higher than or too close to the target! \n"); return(1);}
    if((p.tox-target)<(0.1*target)) {cat("Error: the probability deemed toxic cannot be lower than or too close to the target! \n"); return(1);}
    if(offset>=0.5) {cat("Error: the offset is too large! \n"); return();}
    if(n.earlystop<=6) {cat("Warning: the value of n.earlystop is too low to ensure good operating characteristics. Recommend n.earlystop = 9 to 18 \n"); return();}
    
## obtain dose escalation and de-escalation boundaries
    temp=get.boundary(target, ncohort=150, cohortsize=1, n.earlystop, p.saf, p.tox, cutoff.eli, extrasafe, offset, print=FALSE);
    b.e=temp[2,];   # escalation boundary
    b.d=temp[3,];   # deescalation boundary
    b.elim=temp[4,];  # elimination boundary    
    lambda1  = log((1-p.saf)/(1-target))/log(target*(1-p.saf)/(p.saf*(1-target)));
    lambda2  = log((1-target)/(1-p.tox))/log(p.tox*(1-target)/(target*(1-p.tox)));
    
    n=npts;
    y=ntox;
	earlystop=0;
    d=dose.curr;
    nc = n[d[1],d[2]];
	ndose=length(npts);
	elimi = matrix(rep(0, ndose),dim(n)[1],dim(n)[2]);  ## indicate whether doses are eliminated
    
	
## determine if early termination is needed
	
	if(n[d[1],d[2]]>=n.earlystop)        
	{
		cat("Terminate the trial because the number of patients treated at (", d[1], ", ", d[2], ") has reached", n.earlystop,  "\n");
		d=c(99, 99); earlystop=1; 
	}
	
	if(!is.na(b.elim[nc]))
    {
        if(d[1]==1 && d[2]==1 && y[d[1],d[2]]>=b.elim[nc])
        {
			d=c(99, 99); earlystop=1; 
			cat("Terminate the trial because the lowest dose is overly toxic \n");
        }
        
## implement the extra safe rule by decreasing the elimination cutoff for the lowest dose
        if(extrasafe)
        {
            if(d[1]==1 && d[2]==1 && y[1,1]>=3)
            {
                if(1-pbeta(target, y[1,1]+1, n[1,1]-y[1,1]+1)>cutoff.eli-offset)
                {
					d=c(99, 99); earlystop=1; 
					cat("Terminate the trial because the lowest dose is overly toxic \n");
				}
            }
        }
    }
	
	
## determine elimination status for combinations
	for(i in 1:dim(n)[1])
	{
		for(j in 1:dim(n)[2])
		{
			if(n[i,j]>0 && (!is.na(b.elim[n[i,j]])))
			{
				if(y[i,j]>=b.elim[n[i,j]])
				{
					elimi[i:dim(n)[1], j:dim(n)[2]]=1;
				}
			}
		}
	}
    
	if(earlystop==0)
	{
## dose escalation/de-escalation
		if(y[d[1],d[2]]<=b.e[nc])
		{
			elevel=matrix(c(1,0,0,1),2);
			pr_H0=rep(0,length(elevel)/2)
			nn=pr_H0;
			for ( i in seq(1,length(elevel)/2,by=1))
			{ if (d[1]+elevel[1,i]<=dim(n)[1] && d[2]+elevel[2,i]<=dim(n)[2])
				{
					if (elimi[d[1]+elevel[1,i],d[2]+elevel[2,i]]==0)
					{
						yn=y[d[1]+elevel[1,i],d[2]+elevel[2,i]];
						nn[i]=n[d[1]+elevel[1,i],d[2]+elevel[2,i]];
						pr_H0[i]<-pbeta(lambda2,yn+0.5,nn[i]-yn+0.5)-pbeta(lambda1,yn+0.5,nn[i]-yn+0.5)
					}
				}
			}
			pr_H0=pr_H0+nn*0.0005;  ## break ties
			
			if (max(pr_H0)==0) {d=d} else
			{
				k=which(pr_H0==max(pr_H0))[as.integer(runif(1)*length(which(pr_H0==max(pr_H0)))+1)];
				d=d+c(elevel[1,k],elevel[2,k]);
			}
			
		}
		else if(y[d[1],d[2]]>=b.d[nc])
		{
			delevel=matrix(c(-1,0,0,-1),2)
			pr_H0=rep(0,length(delevel)/2)
			nn=pr_H0;
			for ( i in seq(1,length(delevel)/2,by=1))
			{
				if (d[1]+delevel[1,i]>0 && d[2]+delevel[2,i]>0)
				{
					yn=y[d[1]+delevel[1,i],d[2]+delevel[2,i]];
					nn[i]=n[d[1]+delevel[1,i],d[2]+delevel[2,i]];
					pr_H0[i]=pbeta(lambda2,yn+0.5,nn[i]-yn+0.5)-pbeta(lambda1,yn+0.5,nn[i]-yn+0.5)
				}
			}
			pr_H0=pr_H0+nn*0.0005; ## break ties
			
			if (max(pr_H0)==0) {d=d}  else
			{
				k=which(pr_H0==max(pr_H0))[as.integer(runif(1)*length(which(pr_H0==max(pr_H0)))+1)];
				d=d+c(delevel[1,k],delevel[2,k]);
			}
		}
		else { d=d; }
		cat("The recommended dose combination for the next cohort of patients is (", d[1], ", ", d[2], ")", "\n");
	}
	invisible(d);
}






