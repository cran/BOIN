## tested, 8/21/2015

get.oc.comb <- function(target, p.true, ncohort, cohortsize, n.earlystop=100, startdose=c(1,1), p.saf="default", p.tox="default", cutoff.eli=0.95, extrasafe=FALSE, offset=0.05, ntrial=1000)
{
    ## if the user does not provide p.saf and p.tox, set them to the default values
    if(p.saf=="default") p.saf=0.6*target;
    if(p.tox=="default") p.tox=1.4*target;
    
    ## simple error checking
    if(target<0.05) {cat("Error: the target is too low! \n"); return(1);}
    if(target>0.6)  {cat("Error: the target is too high! \n"); return(1);}
    if((target-p.saf)<(0.1*target)) {cat("Error: the probability deemed safe cannot be higher than or too close to the target! \n"); return(1);}
    if((p.tox-target)<(0.1*target)) {cat("Error: the probability deemed toxic cannot be lower than or too close to the target! \n"); return(1);}
    if(offset>=0.5) {cat("Error: the offset is too large! \n"); return();}
    if(n.earlystop<=6) {cat("Warning: the value of n.earlystop is too low to ensure good operating characteristics. Recommend n.earlystop = 9 to 18 \n"); return();}
    
    
    set.seed(6);
    ndose=length(p.true);
    npts = ncohort*cohortsize;
    Y<-array(matrix(rep(0,length(p.true)*ntrial),dim(p.true)[1]),dim=c(dim(p.true),ntrial)) # toxicity outcome
    N <-array(matrix(rep(0,length(p.true)*ntrial),dim(p.true)[1]),dim=c(dim(p.true),ntrial)) # number of patients
    dselect =matrix(rep(0, 2*ntrial), ncol=2); # the selected dose level
    
    ## obtain dose escalation and de-escalation boundaries
    temp=get.boundary(target, ncohort, cohortsize, n.earlystop, p.saf, p.tox, cutoff.eli, extrasafe, print=FALSE);
    b.e=temp[2,];   # escalation boundary
    b.d=temp[3,];   # deescalation boundary
    b.elim=temp[4,];  # elimination boundary
    
    lambda1  = log((1-p.saf)/(1-target))/log(target*(1-p.saf)/(p.saf*(1-target)));
    lambda2  = log((1-target)/(1-p.tox))/log(p.tox*(1-target)/(target*(1-p.tox)));
    
    ################## simulate trials ###################
    for(trial in 1:ntrial)
    {
        y<-matrix(rep(0, ndose),dim(p.true)[1],dim(p.true)[2]);    ## the number of DLT at each dose level
        n<-matrix(rep(0, ndose),dim(p.true)[1],dim(p.true)[2]);    ## the number of patients treated at each dose level
        earlystop=0;         ## indiate whether the trial terminates early
        d=startdose;         ## starting dose level
        elimi = matrix(rep(0, ndose),dim(p.true)[1],dim(p.true)[2]);  ## indicate whether doses are eliminated
        
        for(pp in 1:ncohort)
        {
            ### generate toxicity outcome
            y[d[1],d[2]] = y[d[1],d[2]] + sum(runif(cohortsize)<p.true[d[1],d[2]]);
            n[d[1],d[2]] = n[d[1],d[2]] + cohortsize;
            if(n[d[1],d[2]]>=n.earlystop) break;
            nc = n[d[1],d[2]];
			
            ## determine if the current dose should be eliminated
            if(!is.na(b.elim[nc]))
            {
                if(y[d[1],d[2]]>=b.elim[nc])
                {
                    for (i in min(d[1],dim(p.true)[1]):dim(p.true)[1])
                    { for (j in min(d[2],dim(p.true)[2]):dim(p.true)[2]) { elimi[i,j]=1;}}
                    if(d[1]==1&&d[2]==1) {d=c(99, 99); earlystop=1; break;}
                }
                
                ## implement the extra safe rule by decreasing the elimination cutoff for the lowest dose
                if(extrasafe)
                {
                    if(d[1]==1&&d[2]==1 && y[1,1]>=3)
                    {
                        if(1-pbeta(target, y[1,1]+1, n[1,1]-y[1,1]+1)>cutoff.eli-offset) {d=c(99, 99); earlystop=1; break;}
                    }
                }
            }
            
            ## dose escalation/de-escalation
            if(y[d[1],d[2]]<=b.e[nc])
            {
                elevel=matrix(c(1,0,0,1),2);
                pr_H0=rep(0,length(elevel)/2)
                nn=pr_H0;
                for ( i in seq(1,length(elevel)/2,by=1))
                { if (d[1]+elevel[1,i]<=dim(p.true)[1] && d[2]+elevel[2,i]<=dim(p.true)[2])
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
            } else { d=d; }
        }
        
        Y[,,trial]=y;
        N[,,trial]=n;
        if(earlystop==1) { dselect[trial,]=c(99, 99); }
        else
        {
            selcomb=select.mtd.comb(target, n, y, cutoff.eli, extrasafe, offset, print=FALSE);
			dselect[trial,1]=selcomb[1];
			dselect[trial,2]=selcomb[2];
        }
    }
    
    # output results
    selpercent=matrix(rep(0, ndose),dim(p.true)[1],dim(p.true)[2]);
    nptsdose  =apply(N,c(1,2),mean, digits=2, format="f");
    ntoxdose  =apply(Y,c(1,2),mean, digits=2, format="f");
    
    for(i in 1:dim(p.true)[1])
    for (j in 1:dim(p.true)[2]){
    { selpercent[i,j]=sum(dselect[,1]==i&dselect[,2]==j)/ntrial*100; }}
  
  
    {
        cat("True toxicity rate of dose combinations:\n");
        for (i in 1:dim(p.true)[1]){
            cat(formatC(p.true[i,], digits=2, format="f", width=5), sep="  ", "\n");
        }
        cat("\n");
        cat("selection percentage at each dose combination (%):\n");
        for (i in 1:dim(p.true)[1]){
            cat(formatC(selpercent[i,], digits=2, format="f", width=5), sep="  ", "\n");
        }
        cat("\n");
        cat("number of patients treated at each dose combination:\n");
        for (i in 1:dim(p.true)[1]){
            cat(formatC(apply(N,c(1,2),mean)[i,], digits=2, format="f", width=5), sep ="  ", "\n");}
        cat("\n");
        cat("number of toxicity observed at each dose combination:\n");
        for (i in 1:dim(p.true)[1]){
            cat(formatC(apply(Y,c(1,2),mean)[i,], digits=2, format="f", width=5), sep ="  ", "\n");}
        cat("\n");
        cat("average number of toxicities:", formatC(sum(Y)/ntrial, digits=1, format="f"), "\n");
        cat("average number of patients:", formatC(sum(N)/ntrial, digits=1, format="f"), "\n");
        cat("selection percentage of MTD:", formatC(sum(selpercent[which(p.true==target, arr.ind = TRUE)]), digits=1, format="f"), "\n");
        cat("percentage of patients treated at MTD:", formatC(sum(nptsdose[which(p.true==target, arr.ind = TRUE)])/npts*100, digits=1, format="f"), "\n");
        cat("percentage of early stopping due to toxicity:", formatC(100-sum(selpercent), digits=2, format="f"), "\n");
    }
    mtdpercent=sum(selpercent[which(p.true==target, arr.ind = TRUE)]);
    mtdpts=sum(nptsdose[which(p.true==target, arr.ind = TRUE)])/npts*100;
    
    invisible(list(selpercent=selpercent, nptsdose=nptsdose, ntoxdose=ntoxdose, mtdpercent=mtdpercent, mtdpts=mtdpts));
}