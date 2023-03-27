#' A function for surveys using simple random sampling
#'
#' This function allows you to make estimation or sample size determinations for surveys using simple random sampling.
#' @param y The data vector to be used for estimation.
#' @param B The desired bound on the estimation. Required for sample size determination.
#' @param N The population size. Required.
#' @param p The population proportion.
#' @param sig2 The population variance.
#' @param sample.size.det Logical statement indicating if the action should be sample size determination. If false, estimation occurs.
#' @param estimate Indicate the desired statistic to be estimated. Options include: mean, total, and proportion.
#' @param systematic Logical statement indicating if a systematic sample of a randomly ordered population is to occur. This results in a k value computed in the sample size determination.
#' @param alpha Significance level.
#' @keywords survey sampling
#' @keywords elementary survey sampling
#' @keywords estimation
#' @keywords sample size determination
#' @keywords simple random sample
#' @keywords systematic random sample
#' @export
#' @examples
#' survsamp.srs(B=.1,N=100,p=.75,
#'    sample.size.det = TRUE,estimate="proportion")
#' survsamp.srs(y=rnorm(50,5,2),N=100)#,sig2=1)
#' survsamp.srs(B=75,N=100,sig2=1,estimate="total",
#'    sample.size.det = TRUE,systematic=TRUE)


survsamp.srs=function(y=NULL,B=NULL,N=NULL,
                  p=NULL,sig2=NULL,
                  sample.size.det=FALSE,
                  estimate="mean",systematic=FALSE,
                  alpha=.05){

  if (is.null(N)){print("Please specify a population size.")}

      if (estimate=="mean"){
        if (sample.size.det){
          D=B^2/(2*qnorm(1-alpha/2))
          n=ceiling(N*sig2/((N-1)*D+sig2))
          df=data.frame(Object=c("Estimate:","alpha:","N:","n:",
                          "sigma2:","B:"),
                        Value=c(estimate,alpha,N,n,sig2,B))
          if (systematic==TRUE){
            k=floor(N/n)
            df=data.frame(Object=c("Estimate:","alpha:","N:","n:",
                                   "sigma2:","B:","k:"),
                          Value= c(estimate,alpha,N,n,sig2,
                                   B,k))
            }
        } else {
          n=length(y)
          #FPC check
          FPC=n/N>1/20#if true, then needed
          muhat=mean(y)
          s2=var(y)

          varmuhat=ifelse(is.null(sig2),
                        ifelse(FPC,(1-n/N),1)*(s2/n),
                        ifelse(FPC,((N-n)/(N-1)),1)*(sig2/n))
        B=qt(1-alpha/2,n-1)*sqrt(varmuhat)
        ci=paste("(",round(muhat-B,5),",",round(muhat+B,5),")",sep="")
        df=data.frame(Object=c("Estimate:","alpha:","N:","n:","FPC:",
                        "muhat:","s2:","sig2:","vhatmuhat:",
                        "B:","ci:"),
                      Value=c(estimate,alpha,N,n,FPC,muhat,s2,
                              ifelse(is.null(sig2),NA,sig2),
                              varmuhat,B,ci))
      }} else
        if (estimate=="total"){
          if (sample.size.det){
            D=B^2/(2*qnorm(1-alpha/2)*N^2)
            n=ceiling(N*sig2/((N-1)*D+sig2))
            df=data.frame(Object=c("Estimate:","alpha:","N:","n:",
                            "sigma2:","B:"),
                         Value= c(estimate,alpha,N,n,sig2,B))
            if (systematic==TRUE){
              k=floor(N/n)
              df=data.frame(Object=c("Estimate:","alpha:","N:","n:",
                                     "sigma2:","B:","k:"),
                            Value= c(estimate,alpha,N,n,sig2,B,k))
            }
          } else {
            n=length(y)
            #FPC check
            FPC=n/N>1/20#if true, then needed
            muhat=mean(y)
            s2=var(y)
            tauhat=muhat*N
            vartauhat=N^2*ifelse(is.null(sig2),
                               ifelse(FPC,(1-n/N),1)*(s2/n),
                               ifelse(FPC,((N-n)/(N-1)),1)*(sig2/n))
          B=qt(1-alpha/2,n-1)*sqrt(vartauhat)
          ci=paste("(",round(tauhat-B,5),",",round(tauhat+B,5),")",sep="")
          df=data.frame(Object=c("Estimate:","alpha:","N:","n:","FPC:",
                          "muhat:","tauhat:","s2:","sig2:","vhattauhat:",
                          "B:","ci:"),
                        Value=c(estimate,alpha,N,n,FPC,
                          muhat,tauhat,s2,
                          ifelse(is.null(sig2),NA,sig2),
                          vartauhat,B,ci))
        }} else
          if (estimate=="proportion"){
            if (sample.size.det){
              D=B^2/(2*qnorm(1-alpha/2))
              n=ceiling(N*p*(1-p)/((N-1)*D+p*(1-p)))
              df=data.frame(Object=c("Estimate:","alpha:","N:","n:",
                              "p:","B:"),
                            Value=c(estimate,alpha,N,n,p,B))
              if (systematic==TRUE){
                k=floor(N/n)
                df=data.frame(Object=c("Estimate:","alpha:","N:","n:",
                                       "p:","B:","k:"),
                              Value= c(estimate,alpha,N,n,p,B,k))
              }
            } else {
              n=length(y)
              #FPC check
              FPC=n/N>1/20#if true, then needed
              phat=sum(y)/n#must be formatted as 1/0 where 1 is the level used for calculating the prop
              varphat=ifelse(is.null(p),
                           ifelse(FPC,1-n/N,1)*((phat*(1-phat))/(n-1)),
                           ifelse(FPC,((N-n)/(N-1)),1)*((p*(1-p))/(n-1)))
            B=qnorm(1-alpha/2)*sqrt(varphat)
            ci=paste("(",round(phat-B,5),",",round(phat+B,5),")",sep="")
            df=data.frame(Object=c("Estimate:","alpha:","N:","n:","FPC:",
                            "phat:","vhatphat:","B:","ci:"),
                          Value=c(estimate,alpha,N,n,FPC,phat,varphat,B,ci))
            }} else {print("Please check your estimate input!")}

  print(df,row.names=FALSE)
  }



#' A function for surveys using simple random sampling when ratio estimation is of interest
#'
#' This function allows you to make estimation or sample size determinations for the ratio wtih surveys using simple random sampling.
#' @param x The independent data vector to be used for estimation.
#' @param y The dependent data vector to be used for estimation.
#' @param B The desired bound on the estimation. Required for sample size determination.
#' @param N The population size. Required for sample size determination.
#' @param mux The population mean of x.
#' @param taux The population total of x. This must be supplied if N is NULL.
#' @param sig2 The population variance.
#' @param sample.size.det Logical statement indicating if the action should be sample size determination. If false, estimation occurs.
#' @param estimate Indicate the desired statistic to be estimated. Options include: ratio, mean, and total. This is only needed for sample size determination as all three are reported in estimation.
#' @param alpha Significance level.
#' @keywords survey sampling
#' @keywords elementary survey sampling
#' @keywords estimation
#' @keywords sample size determination
#' @keywords simple random sample
#' @export
#' @examples
#' x=rnorm(50,5,1)
#' survsamp.ratio(x=x,y=x+rnorm(50),
#' estimate="ratio",mux=5,N=100)


survsamp.ratio=function(x=NULL,y=NULL,B=NULL,N=NULL,
                        mux=NULL,taux=NULL,sig2=NULL,
                        sample.size.det=FALSE,
                        estimate="ratio",
                        alpha=.05){


  if (sample.size.det==TRUE) {
    if (is.null(N)){print("Please specify a population size.")}
    sig2=ifelse(is.null(x)==FALSE&is.null(y)==FALSE,var(y-(mean(y)/mean(x))*x),sig2)
    D=ifelse(estimate=="ratio",B^2*mux^2/(2*qnorm(1-alpha/2)),
             ifelse(estimate=="mean",B^2/(2*qnorm(1-alpha/2)),
                    ifelse(estimate=="total",B^2/(2*qnorm(1-alpha/2)*N^2),NA)))
    if (is.na(D)){print("Please check your estimate input!")}
    n=ceiling(N*sig2/(N*D+sig2))
    df=data.frame(Object=c("Estimate:","alpha:","N:","n:",
                           "sigma2:","B:"),
                  Value=c(estimate,alpha,N,n,sig2,B))
  } else {
    n=length(y)
    #FPC check
    FPC=ifelse(is.null(N),FALSE,n/N>1/20)#if true, then needed
    ybar=mean(y);xbar=mean(x)
    mux=ifelse(is.null(mux)==F,mux,ifelse(is.null(taux)==F&is.null(N)==F,taux/N,xbar))
    taux=ifelse(is.null(taux)==F,taux,ifelse(is.null(N),NULL,N*mux))
    r=ybar/xbar
    tauhaty=r*taux
    muhaty=r*mux
    sr2=var(y-r*x)
    vhatr=ifelse(is.null(mux),
                 ifelse(FPC,(1-n/N),1)*(1/xbar^2)*(sr2/n),
                 ifelse(FPC,(1-n/N),1)*(1/mux^2)*(sr2/n))
    vhattauhaty=taux^2*vhatr
    vhatmuhaty=ifelse(FPC,(1-n/N),1)*(sr2/n)
    B.R=qt(1-alpha/2,n-1)*sqrt(vhatr)
    cir=paste("(",round(r-B.R,5),",",round(r+B.R,5),")",sep="")
    B.mu=qt(1-alpha/2,n-1)*sqrt(vhatmuhaty)
    cimu=paste("(",round(muhaty-B.mu,5),",",round(muhaty+B.mu,5),")",sep="")
    B.tau=qt(1-alpha/2,n-1)*sqrt(vhattauhaty)
    citau=paste("(",round(tauhaty-B.tau,5),",",round(tauhaty+B.tau,5),")",sep="")
    df=data.frame(Object=c("Estimate:","alpha:","N:","n:","FPC:",
                           "Rhat:","muhaty:","tauhaty:",
                           "vhatrhat:","vhatmuhaty:","vhattauhaty:",
                           "B for R:","B for muy:","B for tauy:",
                           "R ci:","muy ci:","tauy ci:"),
                  Value=c(estimate,alpha,ifelse(is.null(N),NA,N),n,FPC,r,muhaty,tauhaty,
                          vhatr,vhatmuhaty,vhattauhaty,
                          B.R,B.mu,B.tau,
                          cir,cimu,citau))

    plot(y~x,ylab="Y",xlab="X")
    abline(a=0,b=r,col="green")
    abline(lm(y~x),col=4)
    legend("bottomright",c("Ratio Estimate","Linear Model"),
           lty=1,col=c("green",4))
  }
  print(df,row.names=FALSE)
}



#' A function for surveys using stratified random sampling
#'
#' This function allows you to make estimation or sample size determinations for surveys using stratified random sampling.
#' @param y The data to be used for estimation. Must be a list with length equal to the number of strata.
#' @param B The desired bound on the estimation. Required for sample size determination.
#' @param N The total population size.
#' @param n The total sample size.
#' @param Ni The population size per strata.
#' @param ni The sample size per strata.
#' @param pi The population proportion per strata.
#' @param sig2 The population variance per strata unless equal variance across strata.
#' @param sample.size.det Logical statement indicating if the action should be sample size determination. If false, estimation occurs.
#' @param budget The budget for sample collection. Can be included in sample size determination.
#' @param ss.det.cost A vector pf cost to collect samples in a given strata. Can be included in sample size determination.
#' @param ss.det.a A vector of the allocation which indicates the proportion to be collected for each strats. Can be included in sample size determination.
#' @param estimate Indicate the desired statistic to be estimated. Options include: mean, total, and proportion.
#' @param alpha Significance level.
#' @keywords survey sampling
#' @keywords elementary survey sampling
#' @keywords estimation
#' @keywords sample size determination
#' @keywords stratified random sample
#' @export
#' @examples
#' survsamp.strat(B=50,N=100,
#'    Ni=c(46,54),sig2=c(1,2),
#'    estimate="total",sample.size.det = TRUE)
#' survsamp.strat(y=list(rnorm(25,5,1),rnorm(25,6,2)),
#'    N=100,n=50,ni=c(25,25),Ni=c(46,54),estimate="mean")

survsamp.strat=function(y=NULL,B=NULL,N=NULL,n=NULL,
                        Ni=NULL,ni=NULL,
                        pi=NULL,sig2=NULL,#if know pop. parameters
                        sample.size.det=FALSE,budget=NULL,
                        ss.det.cost=NULL,ss.det.a=NULL,
                        estimate="mean",
                        alpha=.05){

    if (is.null(N)){print("Please specify a population size.")}

    #FPC check
    FPC=any(ni/Ni<1/20)#if true, then needed
    FPC=ifelse(FPC==TRUE,FALSE,TRUE)

      if (estimate=="mean"){
        if (sample.size.det){
          D=B^2/(2*qnorm(1-alpha/2))
#          sig2=ifelse(rep(is.null(sig2),length(sig2)),si2,sig2)
          if (is.null(ss.det.a)==FALSE){#given a
            ai=ss.det.a
          } else
            if (is.null(ss.det.cost)==FALSE){#diff cost, diff var
              ai=(Ni*sqrt(sig2)/sqrt(ss.det.cost))/sum(Ni*sqrt(sig2)/sqrt(ss.det.cost))
            } else
              if (length(sig2)!=1){#same cost, diff var
              ai=Ni*sqrt(sig2)/sum(Ni*sqrt(sig2))
            } else {#same cost, same var
              ai=Ni/N
            }
          n=ifelse(is.null(budget),
                   ceiling((sum(Ni^2*sig2/ai))/(sum(Ni)^2*D+sum(Ni*sig2))),
                   floor(budget/sum(ss.det.cost*ai)))
          ni=ifelse(rep(is.null(budget),length(ai)),ceiling(n*ai),round(n*ai))
          n=sum(ni)
           df=data.frame(Object=c("Estimate:","alpha:","N:","n:","ni:",
                                 "sigma2:","B:","ai:","ci:","Budget:"),
                        Value=c(estimate,alpha,N,n,paste(as.character(ni),collapse=", "),
                                paste(as.character(sig2),collapse=", "),B,
                                paste(as.character(ai),collapse=", "),
                                ifelse(is.null(ss.det.cost),NA,
                                       paste(as.character(ss.det.cost),collapse=", ")),
                                ifelse(is.null(budget),NA,paste(as.character(budget),collapse=", "))))
        } else {
          if (is.list(y)==FALSE){print("Supply the y values as a list!")}
          ybari=rep(NA,length(y))
          for (i in 1:length(y)){ybari[i]=mean(y[[i]])}
          muhat=sum(Ni*ybari)/N
          si2=rep(NA,length(y))
          for (i in 1:length(y)){si2[i]=var(y[[i]])}
          vhatybar=(1/N^2)*sum(Ni^2*ifelse(FPC,(1-ni/Ni),1)*(si2/ni))
      B=qt(1-alpha/2,n-1)*sqrt(vhatybar)
      ci=paste("(",round(muhat-B,5),",",round(muhat+B,5),")",sep="")
      df=data.frame(Object=c("Estimate:","alpha:","N:","Ni:","n:","ni:","FPC:",
                             "muhat:","si2:","vhatmuhat:",
                             "B:","ci:"),
                    Value=c(estimate,alpha,N,paste(as.character(Ni),collapse=", "),
                            n,paste(as.character(ni),collapse=", "),
                            FPC,muhat,paste(si2,collapse=", "),
                            vhatybar,B,ci))
      }} else
      if (estimate=="total"){
        if (sample.size.det){
          D=B^2/(2*qnorm(1-alpha/2)*N^2)
          sig2=ifelse(rep(is.null(sig2),length(sig2)),si2,sig2)
          if (is.null(ss.det.a)==FALSE){#given a
            ai=ss.det.a
          } else
            if (is.null(ss.det.cost)==FALSE){#diff cost, diff var
              ai=(Ni*sqrt(sig2)/sqrt(ss.det.cost))/sum(Ni*sqrt(sig2)/sqrt(ss.det.cost))
            } else
              if (length(sig2)!=1){#same cost, diff var
                ai=Ni*sqrt(sig2)/sum(Ni*sqrt(sig2))
              } else {#same cost, same var
                ai=Ni/N
              }
          n=ifelse(is.null(budget),
                   ceiling((sum(Ni^2*sig2/ai))/(sum(Ni)^2*D+sum(Ni*sig2))),
                   floor(budget/sum(ss.det.cost*ai)))
          ni=ifelse(rep(is.null(budget),length(ai)),ceiling(n*ai),round(n*ai))
          n=sum(ni)
          df=data.frame(Object=c("Estimate:","alpha:","N:","Ni:","n:","ni:",
                                 "sigma2:","B:","ai:","ci:","Budget:"),
                        Value=c(estimate,alpha,N,paste(Ni,collapse=", "),n,
                                paste(as.character(round(n*ai)),collapse=", "),
                                paste(as.character(sig2),collapse=", "),B,paste(as.character(ai),collapse=", "),
                                ifelse(is.null(ss.det.cost),NA,
                                       paste(as.character(ss.det.cost),collapse=", ")),
                                ifelse(is.null(budget),NA,paste(as.character(budget),collapse=", "))))
          } else {
          if (is.list(y)==FALSE){print("Supply the y values as a list!")}
          ybari=rep(NA,length(y))
          for (i in 1:length(y)){ybari[i]=mean(y[[i]])}
          muhat=sum(Ni*ybari)/N
          si2=rep(NA,length(y))
          for (i in 1:length(y)){si2[i]=var(y[[i]])}
          tauhat=N*muhat
          vhattauhat=sum(Ni^2*ifelse(FPC,(1-ni/Ni),1)*(si2/ni))
        B=qt(1-alpha/2,n-1)*sqrt(vhattauhat)
        ci=paste("(",round(tauhat-B,5),",",round(tauhat+B,5),")",sep="")
        df=data.frame(Object=c("Estimate:","alpha:","N:","Ni:","n:","ni:","FPC:",
                               "tauhat:","si2:","vhattauhat:",
                               "B:","ci:"),
                      Value=c(estimate,alpha,N,paste(as.character(Ni),collapse=", "),
                              n,paste(as.character(ni),collapse=", "),
                              FPC,tauhat,paste(si2,collapse=", "),
                              vhattauhat,B,ci))
        }} else
        if (estimate=="proportion"){
          if (sample.size.det){
            D=B^2/(2*qnorm(1-alpha/2))
            if (is.null(ss.det.a)==FALSE){#given a
              ai=ss.det.a
            } else
              if (is.null(ss.det.cost)==FALSE){#diff cost, diff var
                ai=(Ni*sqrt(pi*(1-pi)/ss.det.cost))/sum(Ni*sqrt(pi*(1-pi))/sqrt(ss.det.cost))
              } else
                if (length(pi)!=1){#same cost, diff p
                  ai=Ni*sqrt(pi*(1-pi))/sum(Ni*sqrt(pi*(1-pi)))
                } else {#same cost, same p
                  ai=Ni/N
                }
            n=ifelse(is.null(budget),
                     ceiling((sum(Ni^2*pi*(1-pi))/ai)/(sum(Ni)^2*D+sum(Ni*pi*(1-pi)))),
                     floor(budget/sum(ss.det.cost*ai)))
            ni=ifelse(rep(is.null(budget),length(ai)),ceiling(n*ai),round(n*ai))
            n=sum(ni)
            df=data.frame(Object=c("Estimate:","alpha:","N:","Ni:","n:","ni:",
                                   "pi:","B:","ai:","ci:","Budget:"),
                          Value=c(estimate,alpha,N,paste(as.character(Ni),collapse=", "),n,paste(as.character(ni),collapse=", "),
                                  paste(as.character(pi),collapse=", "),B,paste(as.character(ai),collapse=", "),
                                  ifelse(is.null(ss.det.cost),NA,
                                         paste(as.character(ss.det.cost),collapse=", ")),
                                  ifelse(is.null(budget),NA,paste(as.character(budget),collapse=", "))))
            } else {
              if (is.list(y)==FALSE){print("Supply the y values as a list!")}
              phati=rep(NA,length(y))
            for (i in 1:length(y)){phati[i]=sum(y[[i]])/ni[i]}
            phat=sum(Ni*phati)/N
            vhatphat=(1/N^2)*sum(Ni^2*(1-ni/Ni)*(phati*(1-phati)/(ni-1)))
          B=qnorm(1-alpha/2)*sqrt(vhatphat)
          ci=paste("(",round(phat-B,5),",",round(phat+B,5),")",sep="")
          df=data.frame(Object=c("Estimate:","alpha:","N:","Ni:","n:","ni:","FPC:",
                                 "phat:","vhatphat:",
                                 "B:","ci:"),
                        Value=c(estimate,alpha,N,paste(as.character(Ni),collapse=", "),
                                n,paste(as.character(ni),collapse=", "),
                                FPC,phat,
                                vhatphat,B,ci))
          }} else {print("Please check your estimate input!")}
    print(df,row.names=FALSE)
        }



#' A function for surveys using simple random sampling when estimation of the difference between two groups is of interest
#'
#' This function allows you to make estimation or sample size determinations for the difference between two groups with surveys using simple random sampling.
#' @param y1 The data vector from population 1 to be used for estimation.
#' @param y2 The data vector from population 2 to be used for estimation.
#' @param N1 The population size for population 1. Required.
#' @param N2 The population size for population 2. Required.
#' @param p1 The population 1 proportion.
#' @param p2 The population 2 proportion.
#' @param sig2.1 The population 1 variance.
#' @param sig2.2 The population 2 variance.
#' @param estimate Indicate the desired statistic to be estimated. Options include: mean, total, and proportion.
#' @param alpha Significance level.
#' @param prop.dpnt Logical statement indicating if the two proportions are dependent.
#' @keywords survey sampling
#' @keywords elementary survey sampling
#' @keywords estimation
#' @keywords sample size determination
#' @keywords simple random sample
#' @export
#' @examples
#' survsamp.diff(y1=rnorm(25,5,1),y2=rnorm(25,10,2),
#'      N1=46,N2=54,sig2.1=1,sig2.2=2,
#'      estimate="mean")
#' survsamp.diff(y1=rbinom(25,1,.75),y2=rbinom(25,1,.5),
#'      N1=46,N2=54,
#'      estimate="proportion",prop.dpnt=FALSE)

survsamp.diff=function(y1=NULL,y2=NULL,N1=NULL,N2=NULL,
                       sig2.1=NULL,sig2.2=NULL,
                       p1=NULL,p2=NULL,
                       estimate="mean",
                       alpha=.05,
                       prop.dpnt=FALSE){

  if ((is.null(N1)|is.null(N2))&estimate=="mean"){print("Please specify a population size for both samples.")}

  n1=length(y1)
  muhat1=mean(y1)
  s2.1=var(y1)
  phat1=sum(y1)/n1#must be formatted as 1/0 where 1 is the level used for calculating the prop
  n2=length(y2)
  muhat2=mean(y2)
  s2.2=var(y2)
  phat2=sum(y2)/n2#must be formatted as 1/0 where 1 is the level used for calculating the prop


  if (estimate=="mean"){
    #FPC check
    FPC=(n1/N1>1/20)&(n2/N2>1/20)#if true, then needed

    d=muhat1-muhat2
    varybard=ifelse(is.null(sig2.1),ifelse(FPC,(1-n1/N1),1)*s2.1/n1,
                    ifelse(FPC,(N1-n1)/(N1-1),1)*sig2.1/n1)+
              ifelse(is.null(sig2.2),ifelse(FPC,(1-n2/N2),1)*s2.2/n2,
                     ifelse(FPC,(N2-n2)/(N2-1),1)*sig2.2/n2)
    B=qt(1-alpha/2,min(n1-1,n2-1))*sqrt(varybard)
    ci=paste("(",round(d-B,5),",",round(d+B,5),")",sep="")
    df=data.frame(Object=c("Estimate:","alpha:","N1:","N2:","n1:","n2:","FPC:",
                           "d:","vhatybard:","B:","ci:"),
                  Value=c(estimate,alpha,N1,N2,n1,n2,FPC,d,
                          varybard,B,ci))
  } else
    if (estimate=="proportion"){
      pdhat=phat1-phat2
      varpdhat=phat1*(1-phat1)/n1+phat2*(1-phat2)/n2+
        ifelse(prop.dpnt,2*phat1*phat2/(n1),0)
      B=qnorm(1-alpha/2)*sqrt(varpdhat)
      ci=paste("(",round(pdhat-B,5),",",round(pdhat+B,5),")",sep="")
      df=data.frame(Object=c("Estimate:","alpha:","n1:","n2:",
                             "pdhat:","vhatpdhat:","B:","ci:"),
                    Value=c(estimate,alpha,n1,n2,pdhat,
                            varpdhat,B,ci))
    } else {print("Difference estimates only valid for means and proportions!")}
  print(df,row.names=FALSE)
}



#' A function for surveys using cluster random sampling
#'
#' This function allows you to make estimation or sample size determinations for surveys using cluster random sampling.
#' @param y The total of all observations in the sample needed for mean and total estimation.
#' @param B The desired bound on the estimation. Required for sample size determination.
#' @param a The total number of elements in cluster with the characteristic of interest. Needed for estimating the proportion.
#' @param N The number of clusters in the population. Required.
#' @param m The number of elements per sampled cluster.
#' @param M The number of elements in the population.
#' @param p The population proportion.
#' @param sig2 The population variance.
#' @param sample.size.det Logical statement indicating if the action should be sample size determination. If false, estimation occurs.
#' @param estimate Indicate the desired statistic to be estimated. Options include: mean, total, and proportion.
#' @param alpha Significance level.
#' @keywords survey sampling
#' @keywords elementary survey sampling
#' @keywords estimation
#' @keywords sample size determination
#' @keywords cluster random sample
#' @export
#' @examples
#' survsamp.cluster(N=12,M=100,sig2=35000,estimate="mean",
#'    B=9,sample.size.det=TRUE)
#' set.seed(123019)
#' m=rpois(9,10)
#' survsamp.cluster(y=rpois(9,100),m=m,N=12,M=100,estimate="total")
#' survsamp.cluster(a=m-rpois(9,3),m=m,N=12,M=100,estimate="proportion")

survsamp.cluster=function(y=NULL,B=NULL,a=NULL,N=NULL,
                          m=NULL,M=NULL,
                          p=NULL,sig2=NULL,
                          sample.size.det=FALSE,
                          estimate="mean",
                          alpha=.05){

  if (is.null(N)){print("Please specify a population size.")}

  n=ifelse(is.null(y)&is.null(a),length(m),ifelse(is.null(a),length(y),length(a)))
  ybarcl=sum(y)/sum(m)
  ybart=sum(y)/n
  Mbar=ifelse(is.null(M)==F,M/N,sum(m)/n)#mean cluster size in pop (or sample)
  tauhat=ifelse(is.null(M),N*ybart,M*ybarcl)
  sr2=var(y-ybarcl*m)
  st2=var(y-ybart)
  phat=sum(a)/sum(m)
  sp2=var(a-phat*m)

  #FPC check
  FPC=n/N>1/20#if true, then needed

  if (estimate=="mean"){
    if (sample.size.det){
      sig2=ifelse(is.null(sig2),sr2,sig2)
      D=B^2*Mbar^2/(2*qnorm(1-alpha/2))
      n=ceiling(N*sig2/(N*D+sig2))
      df=data.frame(Object=c("Estimate:","alpha:","N:","n:","M:","Mbar:","sig2:","B:"),
                    Value=c(estimate,alpha,N,n,
                            ifelse(is.null(M),NA,M),Mbar,sig2,B))
    } else {
      vhatybar=ifelse(FPC,(1-n/N),1)*(sr2/(n*Mbar^2))
      B=qt(1-alpha/2,df=n-1)*sqrt(vhatybar)
      ci=paste("(",round(ybarcl-B,5),",",round(ybarcl+B,5),")",sep="")
    df=data.frame(Object=c("Estimate:","alpha:","N:","n:","M:","Mbar:","FPC:",
                           "ybarcl:","sr2:","vhatybar:",
                           "B:","ci:"),
                  Value=c(estimate,alpha,N,n,
                          ifelse(is.null(M),NA,M),Mbar,FPC,ybarcl,
                          sr2,vhatybar,B,ci))
  }} else
    if (estimate=="total"){
      if (sample.size.det){
        sig2=ifelse(is.null(sig2)&is.null(M),st2,ifelse(is.null(sig2),sr2,sig2))
        D=B^2/(2*qnorm(1-alpha/2)*N^2)
        n=ceiling(N*sig2/(N*D+sig2))
        df=data.frame(Object=c("Estimate:","alpha:","N:","n:","M:","Mbar:","sig2:","B:"),
                      Value=c(estimate,alpha,N,n,ifelse(is.null(M),NA,M),Mbar,sig2,B))
      } else {
        vhattauhat=N^2*ifelse(FPC,(1-n/N),1)*(ifelse(is.null(M),st2,sr2)/n)
        B=qt(1-alpha/2,df=n-1)*sqrt(vhattauhat)
        ci=paste("(",round(tauhat-B,5),",",round(tauhat+B,5),")",sep="")
        df=data.frame(Object=c("Estimate:","alpha:","N:","n:","M:","Mbar:","FPC:",
                               "muhat:","ybarcl:","tauhat:","st2:","sr2:","vhattauhat:",
                               "B:","ci:"),
                      Value=c(estimate,alpha,N,n,ifelse(is.null(M),NA,M),
                              Mbar,FPC,ybarcl,ybarcl,tauhat,
                              st2,sr2,vhattauhat,B,ci))
      }
    } else
      if (estimate=="proportion"){
        if (sample.size.det){
          sig2=ifelse(is.null(sig2),sp2,sig2)
          D=B^2*Mbar^2/(2*qnorm(1-alpha/2))
          n=ceiling(N*sig2/(N*D+sig2))
          df=data.frame(Object=c("Estimate:","alpha:","N:","n:","M:","Mbar:",
                                 "sig2:","B:"),
                        Value=c(estimate,alpha,N,n,
                                ifelse(is.null(M),NA,M),Mbar,sig2,B))
        } else {
          vhatphat=ifelse(FPC,(1-n/N),1)*sp2/(n*Mbar^2)
          B=qnorm(1-alpha/2)*sqrt(vhatphat)
          ci=paste("(",round(phat-B,5),",",round(phat+B,5),")",sep="")
          df=data.frame(Object=c("Estimate:","alpha:","N:","n:","M:","Mbar:","FPC:",
                                 "phat:","vhatphat:","B:","ci:"),
                        Value=c(estimate,alpha,N,n,
                                ifelse(is.null(M),NA,M),
                                Mbar,FPC,phat,
                                vhatphat,B,ci))
        }
      } else {print("Please check your estimate input!")}
  print(df,row.names=FALSE)
}



#' A function for population size estimation.
#'
#' This function allows you to make estimation or sample size determinations for surveys using simple random sampling.
#' @param method The method of sampling. Options include: direct, inverse, and quadrant.
#' @param t Size of the first sample (direct and inverse sampling methods).
#' @param n Size of the second sample (direct and inverse sampling methods).
#' @param s The number of retagged individuals (direct and inverse sampling methods).
#' @param mbar The estimated average number of elements per selected quadrant.
#' @param a The area size of the selected quadrants.
#' @param m A vector of the number of elements per selected quadrant.
#' @param n The number of selected quadrants.
#' @param A The total area.
#' @param alpha Significance level.
#' @param evenly.dispersed Logical statement indicating if the elements are evenly dispersed.
#' @keywords survey sampling
#' @keywords elementary survey sampling
#' @keywords estimation
#' @keywords population size estimation
#' @export
#' @examples
#' est.pop.size(method="direct",n=200,t=300,s=62)
#' est.pop.size(method="inverse",n=100,t=150,s=35)
#' est.pop.size(method="quadrant",a=16,
#'             m=c(rep(0,13),rep(1,8),rep(2,12),rep(3,10),rep(4,5),rep(5,2)))
#' #not enough information above to est. M, only lambda (density)
#' est.pop.size(method="quadrant",
#'             mbar=40,a=10,A=60*8,n=20)


est.pop.size=function(method="direct",
                      n=NULL,t=NULL,s=NULL,
                      mbar=NULL,
                      a=NULL,m=NULL,A=NULL,
                      alpha=0.05,
                      evenly.dispersed=TRUE){
  if (method=="direct"){
    Nhat=ceiling(n*t/s)
    VhatNhat=t^2*n*(n-s)/s^3
    B=ceiling(qnorm(1-alpha/2)*sqrt(VhatNhat))
    ci=paste("(",ceiling(Nhat-B),",",ceiling(Nhat+B),")",sep="")
    Nhatc=ceiling((t+1)*(n+1)/(s+1)-1)
    VhatNhatc=(t+1)*(n+1)*(t-s)*(n-s)/((s+1)^2*(s+2))
    Bc=ceiling(qnorm(1-alpha/2)*sqrt(VhatNhatc))
    c.ci=paste("(",ceiling(Nhatc-Bc),",",ceiling(Nhatc+Bc),")",sep="")
    df=data.frame(Object=c("Method:","n:","t:","s:",
                           "Nhat:","VhatNhat:","B:","ci:",
                           "Nhatc:","VhatNhatc:","Bc:","c.ci:"),
                  Value=c(method,n,t,s,
                          Nhat,VhatNhat,B,ci,
                          Nhatc,VhatNhatc,Bc,c.ci))
  } else
    if (method=="inverse"){
    Nhat=ceiling(n*t/s)
    VhatNhat=t^2*n*(n-s)/(s^2*(s+1))
    B=ceiling(sqrt(VhatNhat)*qnorm(1-alpha/2))
    ci=paste("(",ceiling(Nhat-B),",",ceiling(Nhat+B),")",sep="")
    df=data.frame(Object=c("method:","n:","t:","s:",
                           "Nhat:","VhatNhat:","B:","ci:"),
                  Value=c(method,n,t,s,
                          Nhat,VhatNhat,B,ci))
  } else
    if (method=="quadrant"){
    n=ifelse(is.null(m),n,length(m))
    mbar=ifelse(is.null(m),mbar,mean(m))
    lambdahat=mbar/a
    Mhat=ifelse(is.null(A),NA,ceiling(lambdahat*A))
    if (evenly.dispersed){
      Vhatlambdahat=lambdahat/(a*n)
      B.lbd=qnorm(1-alpha/2)*sqrt(Vhatlambdahat)
      VhatMhat=ifelse(is.null(A),NA,(A^2*lambdahat)/(a*n))
      B=ceiling(ifelse(is.null(A),NA,qnorm(1-alpha/2)*sqrt(VhatMhat)))
    } else {
      sm2=var(m)
      Vhatlambdahat=(1/a^2)*(sm2/n)
      B.lbd=qnorm(1-alpha/2)*sqrt(Vhatlambdahat)
      VhatMhat=ifelse(is.null(A),NA,(A^2/a^2)*(sm2/n))
      B=ceiling(ifelse(is.null(A),NA,qnorm(1-alpha/2)*sqrt(VhatMhat)))
    }
    ci.lbd=paste("(",round(lambdahat-B.lbd,5),",",round(lambdahat+B.lbd,5),")",sep="")
    ci=ifelse(is.null(A),NA,paste("(",ceiling(Mhat-B),",",ceiling(Mhat+B),")",sep=""))
    df=data.frame(Object=c("method:","n:","mbar:","lambdahat:",
                           "Vhatlambdahat:","B.lambda:","ci.lambda:",
                           "Mhat:","VhatMhat:","B.M:","ci.M:"),
                  Value=c(method,n,mbar,lambdahat,
                          Vhatlambdahat,
                          B.lbd,ci.lbd,
                          Mhat,VhatMhat,B,ci))
    }
  print(df,row.names=FALSE)
}

