## tested 8/28/2015

select.mtd.comb <- function(target, npts, ntox, cutoff.eli=0.95, extrasafe=FALSE, offset=0.05, print=TRUE)
{
    ## isotonic transformation using the pool adjacent violator algorithm (PAVA)
    ## determine whether the dose has been eliminated during the trial
    
    y=ntox;
    n=npts;
    elimi=matrix(0,dim(n)[1],dim(n)[2]);
    for(i in 1:dim(n)[1])
    {
        for (j in 1:dim(n)[2])
        {
            if(n[i,j]>=3) {if(1-pbeta(target, y[i,j]+0.5, n[i,j]-y[i,j]+0.5)>cutoff.eli)
            {elimi[i:dim(n)[1],j]=1;elimi[i,j:dim(n)[2]]=1; break;}}
        }
        
        if(extrasafe)
        {
            if(n[i,j]>=3) {if(1-pbeta(target, y[i,j]+0.5, n[i,j]-y[i,j]+0.5)>cutoff.eli-offset)
            {elimi[i:dim(n)[1],j]=1;elimi[i,j:dim(n)[2]]=1; break;}}
        }
        
    }
    
    if(elimi[1]==1) { selectdose=c(99, 99); } ## no dose should be selected if the first dose is already very toxic
    else
    {
        phat = (y+0.05)/(n+0.1);
        phat[elimi==1]=1.1
        ## perform the isotonic transformation using PAVA
        phat=Iso::biviso(phat,n+0.1,warn=TRUE)[,];
		phat.out=phat; phat.out[n==0]=NA;
        ## break the ties
        phat = phat*(n!=0)+(1E-5)*(matrix(rep(1:dim(n)[1], each = dim(n)[2], 
        len = length(n)),dim(n)[1],byrow=T) +
        matrix(rep(1:dim(n)[2], each = dim(n)[1], len = length(n)),dim(n)[1]))
        ## select dose closest to the target as the MTD
		phat[n==0]=10; ## so that the dose without treating patients will not be selected
        selectdose=which(abs(phat-target) == min(abs(phat-target)), arr.ind = TRUE)
		if(length(selectdose)>2) selectdose=selectdose[1,]  ##if there are still ties, randomly pick the first one.
    }
	
    if(print==TRUE)
    {
        if(selectdose[1]==99 && selectdose[2]==99) { cat("All tested doses are overly toxic. No MTD is selected! \n")}
        else { cat("The MTD is dose combination (", selectdose[1], ", ", selectdose[2], ") \n\n"); 
			cat("Isotonic estimates of toxicity probablities for combinations are \n"); 
			for (i in 1:dim(phat.out)[1]){
				cat(formatC(phat.out[i,], digits=2, format="f", width=5), sep="  ", "\n");
			}
            cat("\n");
            cat("NOTE: no estimate is provided for the doses at which no patient was treated.")
		}
    }
    else { return(selectdose); }
}



