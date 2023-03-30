
#Run functions at the bottom before running 'manipulate' calls above
install.packages('manipulate')
require(manipulate)

manipulate(EnvironmentSD=slider(0.1,0.75),StockUncertainty=slider(0,0.3),StockAssessFrequency=slider(1,4),MSY_Buffer=picker(1,0.9,0.8,0.7,0.6),Sim_Length=slider(60,300),
           fishplot(MSY_Buffer,EnvironmentSD,StockUncertainty,StockFreq=StockAssessFrequency,tplot=Sim_Length))

manipulate(Proportion_MPA=slider(0.25,0.75),pAuto=checkbox(FALSE,"pMPA Auto"),Survival=slider(0,0.85),SU=checkbox(FALSE,"Stock Uncertainty"),EnvironmentSD=slider(0.1,0.75),Sim_Length=slider(60,300),
           mpaplot(Proportion_MPA,pAuto,Survival,0.25*SU,EnvironmentSD,tplot=Sim_Length))





mpaplot=function(p=0.5,pAuto=TRUE,S=0,ssd=0,esd=0.1,r=0.9,tplot=30,trun=5*tplot){
  if(S==0) nA=1 else nA=ceiling(log(0.05)/log(S)); Sp=S-S^nA;
  set.seed(2); errH=rnorm(trun,1,ssd); E=rnorm(trun,1,esd)-1; 
  nl=seq(0,1.1,len=1e3); Y=popdyn(nl,0,r=r)+0*Sp*nl-nl; 
  Bmsy=nl[which.max(Y)]; Cmsy=max(Y)+Sp*Bmsy; Fmsy=Cmsy/(Cmsy+Bmsy);
  Fmsy2=max(Y)/(max(Y)+Bmsy);
  
  if(pAuto) p=(1-Fmsy) - Fmsy*(Sp*Bmsy/popdyn(Bmsy,0,r=r))
  tN=pN=Bmsy*t(t(S^(0:(nA-1))/sum(S^(0:(nA-1))))); hN=0*tN; tH=rH=0;
  for(t in 2:trun){
    tnew=c(popdyn(sum(tN[,t-1]),E[t],r=r), tN[head(1:nA,-1),t-1]*S)
    Ht=errH[t]*Fmsy2; tN=cbind(tN,c((1-Ht)*tnew)); tH=c(tH,Ht*sum(tnew));
    
    rnewRec=popdyn(sum(pN[,t-1],hN[,t-1]),E[t],r=r)
    pN=cbind(pN,c(p*rnewRec,pN[head(1:nA,-1),t-1]*S))
    hnew=c((1-p)*rnewRec,hN[head(1:nA,-1),t-1]*S)
    hMPA=1; hN=cbind(hN,c((1-hMPA)*hnew)); rH=c(rH,hMPA*sum(hnew));
  }
  pNt=tail(colSums(pN),tplot); tNt=tail(colSums(tN),tplot); 
  CEXM=1.2; CEXL=1.5; par(mfcol=c(2,1+(S>0)),mar=c(4,4,4,2)); YLM=c(0,max(1,tNt,pNt));
  plot(tNt,ylim=YLM,type="o",las=1,pch=16,cex.main=CEXM,cex.lab=CEXL,ylab="Biomass",main=paste0(c("Traditional management.  ","")[1+(S>0)]," Proportion not caught annually=",round(1-Fmsy,2),",  Mean catch=",round(mean(tH,na.rm=TRUE),2)),xlab=""); 
  matplot(cbind(pmax(pNt*1,0.05),pNt*0.2),ylim=YLM,col=c(4,2),cex.lab=CEXL,type="o",las=1,pch=16,cex.main=CEXM,ylab="Biomass",main=paste0(c("MPA-based management.  ","")[1+(S>0)]," Proportion system in MPA=",round(p,2),",  Mean catch=",round(mean(rH,na.rm=TRUE),2)),xlab="Years"); 
  if(S>0) barplot(100*rowMeans(tN)/sum(rowMeans(tN)),cex.lab=CEXL,names.arg=as.character(1:nA),cex.main=CEXM,xlab="",las=1,ylab="% of population",main="Traditional management: Age structure")
  if(S>0) barplot(100*rowMeans(pN)/sum(rowMeans(pN)),cex.lab=CEXL,col=4,names.arg=as.character(1:nA),cex.main=CEXM,xlab="Age",las=1,ylab="% of population",main="MPA-based management: Age structure")
  #hist(pmin(tH[-1]/Cmsy,2.5),breaks=30,xlim=c(0,2.5),xlab="Catch relative to MSY",main="Traditional management: Fishery catch")
  #hist(pmin(rH[-1]/Cmsy,2.5),breaks=30,xlim=c(0,2.5),col=4,xlab="Catch relative to MSY",main="MPA-based management: Fishery catch")
}



popdyn=function(n,e=0,k=1,r=1.5) n*exp(r*(1-n/k) + e)
fishplot=function(FmsyFact=1,esd=0,ssd=0,catchtype=2,tplot=30,trun=100*tplot,r=0.9,StockFreq=1){
  set.seed(1); errH=rnorm(trun,1,ssd); E=rnorm(trun,1,esd)-1; 
  nl=seq(0,1.1,len=1e3); Y=popdyn(nl,0,r=r)-nl; Bmsy=nl[which.max(Y)]; 
  Cmsy=max(Y); Fmsy=Cmsy/Bmsy; Nt=Bmsy+0.05; Ht=0; 
  StockRef=rep(seq(1,trun,by=StockFreq),each=StockFreq);
  if(catchtype==1) for(t in 2:trun){ H=errH[t]*FmsyFact*Cmsy; Ht=c(Ht,H); Nt=c(Nt,popdyn(Nt[t-1],E[t],r=r)-H); }
  # if(catchtype==2) for(t in 2:trun) Nt=c(Nt,popdyn(Nt[t-1],E[t])*errH[t]*(1-FmsyFact*Cmsy/(Cmsy+Bmsy)))
  # if(catchtype==2) for(t in 2:trun) Nt=c(Nt, pmax(popdyn(Nt[t-1],E[t]) - Nt[t-1]*errH[t]*FmsyFact*Fmsy),0.05)
  if(catchtype==2) for(t in 2:trun){ H=Nt[StockRef[t-1]]*errH[t]*FmsyFact*Fmsy;
  Ht=c(Ht,H); Nt=c(Nt, pmax((popdyn(Nt[t-1],r=r) - H)*exp(E[t]*Nt[t-1]),0.05)); }
  Nt=pmax(Nt,0); Nt[is.na(Nt)]=0; par(mfrow=c(1,1)); # layout(t(c(1,1,2)))
  plot(Nt[1:tplot],ylim=c(0,max(c(1,max(Nt[1:tplot])))),type="o",las=1,pch=16,ylab="Proportion of carrying capacity",main=paste0("Population biomass,  ",100*round(mean(Nt<Bmsy),2),"% years < Bmsy"),xlab="Years"); 
  polygon(c(0,2*tplot,2*tplot,0),c(Bmsy,Bmsy,0,0),col=rgb(1,0.3,0.4,alpha=0.2),border=0)
  text((9/30)*tplot,0.1,expression('B<B'['MSY']*', "Overfished"'),cex=1.35,col=2); # hist(Ht[-1]/Cmsy,xlim=c(0,2.5),xlab="",breaks=10,main="Fishery yield, relative to MSY"); abline(v=1,col=2,lwd=2);
}

