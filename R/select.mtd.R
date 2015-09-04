select.mtd <- function(target, npts, ntox, cutoff.eli=0.95, extrasafe=FALSE, offset=0.05, print=TRUE)
{
## isotonic transformation using the pool adjacent violator algorithm (PAVA)
	pava <- function (x, wt = rep(1, length(x))) 
	{
		n <- length(x)
		if (n <= 1) 
		return(x)
		if (any(is.na(x)) || any(is.na(wt))) {
			stop("Missing values in 'x' or 'wt' not allowed")
		}
		lvlsets <- (1:n)
		repeat {
			viol <- (as.vector(diff(x)) < 0)
			if (!(any(viol))) 
			break
			i <- min((1:(n - 1))[viol])
			lvl1 <- lvlsets[i]
			lvl2 <- lvlsets[i + 1]
			ilvl <- (lvlsets == lvl1 | lvlsets == lvl2)
			x[ilvl] <- sum(x[ilvl] * wt[ilvl])/sum(wt[ilvl])
			lvlsets[ilvl] <- lvl1
		}
		x
	}
## determine whether the dose has been eliminated during the trial
	y=ntox;
	n=npts;
	ndose=length(n);
	elimi=rep(0, ndose);
	for(i in 1:ndose)
	{
		if(n[i]>=3) {if(1-pbeta(target, y[i]+1, n[i]-y[i]+1)>cutoff.eli) {elimi[i:ndose]=1; break;}}
	}
    if(extrasafe)
    {
        if(n[1]>=3) {if(1-pbeta(target, y[1]+1, n[1]-y[1]+1)>cutoff.eli-offset) {elimi[1:ndose]=1;}}
    }
	
## no dose should be selected (i.e., selectdose=99) if the first dose is already very toxic or 
## all uneliminated doses are never used to treat patients
	if(elimi[1]==1 || sum(n[elimi==0])==0) { selectdose=99; } 
	else
	{
		adm.set = (n!=0) & (elimi==0);
		adm.index = which(adm.set==T);
		y.adm = y[adm.set];
		n.adm = n[adm.set];
		
## poster mean and variance of toxicity probabilities using beta(0.05, 0.05) as the prior 
		phat = (y.adm+0.05)/(n.adm+0.1); 
		phat.var = (y.adm+0.05)*(n.adm-y.adm+0.05)/((n.adm+0.1)^2*(n.adm+0.1+1));
		
## perform the isotonic transformation using PAVA
		phat = pava(phat, wt=1/phat.var) 
		phat = phat + (1:length(phat))*1E-10 ## break ties by adding an increasingly small number 
		selectd = sort(abs(phat-target), index.return=T)$ix[1]  ## select dose closest to the target as the MTD
		selectdose = adm.index[selectd];
	}
	
	if(print==TRUE)
	{
		if(selectdose==99) { cat("All tested doses are overly toxic. No MTD is selected! \n")}
		else { cat("The MTD is dose level ", selectdose, "\n\n"); }
		
## output summary statistics
		poverdose = pava(1-pbeta(target, y+0.05, n-y+0.05));
		phat.all = pava((y+0.05)/(n+0.1), wt=1/((y+0.05)*(n-y+0.05)/((n+0.1)^2*(n+0.1+1))));
		cat("Dose    Posterior DLT             95%                  \n", sep="");
		cat("Level     Estimate         Credible Interval   Pr(toxicity>", target, "|data)\n", sep="");
		for(i in 1:ndose)
		{
            if(n[i]>0)
            {
                cat(" ", i, "        ", formatC(phat.all[i], digits=2, format="f"), "         (", formatC(qbeta(0.025, y[i]+0.05, n[i]-y[i]+0.05), digits=2, format="f"), 
				", ", formatC(qbeta(0.975, y[i]+0.05, n[i]-y[i]+0.05), digits=2, format="f"), ")            ", 
				formatC(poverdose[i], digits=2, format="f"), "\n");
            }
            else  # no estimate output for doses never used to treat patients
            {  cat(" ", i, "        ", "----", "         (", "------------", ")            ", 
				"----", "\n");
            }
		}
        cat("NOTE: no estimate is provided for the doses at which no patient was treated.")
		
	}
	else
	{
		return(selectdose);
	}
}

