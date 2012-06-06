#-- source functions for Cumulative Exposure Model _ 11 May 2012 --#

#-- code speed-ups (June 03, 2012) --#
  #-- 1) replace all () with {} --#

DemogFun1<-function(k){
  out<-ifelse(k<=355,"Lamb","Ewe")
  return(out)
}

StatusFun1<-function(k){
  out<-ifelse(k=="C","I",k)
  return(out)
}

birth.mod.fun<-function(t){
  b<-ifelse(t %% 365==210,1,0)
  return(b)
}

aging.mod.fun<-function(t){
  b<-ifelse(t %% 365==243,1,0)
  return(b)
}

lambtransmission.mod.fun<-function(t){
  b<-ifelse(t %% 365 >= 210 & t %% 365 <=300, 1, 0)
  return(b)
}  

mother.fun<-function(k,NewMom){
  j<-ifelse(k %in% NewMom,1,0)
  return(j)
}  

lamb.check.fun<-function(k,CurrentMoms){
  j<-ifelse(k %in% CurrentMoms,1,0)
  return(j)
}

ewe.check.fun<-function(k,CurrentLambs){
  j<-ifelse(k %in% CurrentMoms,1,0)
  return(j)
}

grp.split.fun<-function(t){
  gsplit<-ifelse(t %% 365 ==60, 1, 0)
  return(gsplit)
}

grp.merge.fun<-function(t){
  gmerge<-ifelse(t %% 365 ==0,1,0)
  return(gmerge)
}


#transmission.fun<-function(season,temp.in,temp,contactnumber,Lambcontactnumber,LambOpen,LambTransmissionProb,chronicdose=chronicdose){
#  SpatGrp.in<-ifelse(season==0,temp.in$SpatGrp,temp.in$SpatGrpSeason)
#  if(temp.in$DemogGrp=="Ewe"){
#    if(season==0){
#      k<-subset(temp, DemogGrp=="Ewe")
#      g<-dim(k)[1]
#    } else {
#      k<-subset(temp, SpatGrpSeason=SpatGrp.in)
#      g<-dim(k)[1]
#    }
#    samp<-sample(1:g,contactnumber,replace=T)
#    contactset<-ifelse(g>=1,sum(k$SheddingRate[samp]),0)
#    contacts<-k$Status[samp]
#    contactIDs<-k$ID[samp]
#    groupsize<-g
#    lamb<-subset(temp,DemogGrp=="Lamb" & Mother==temp.in$ID)
#    contactL<-ifelse(dim(lamb)[1]==0,0,ifelse(lamb$Status=="I",1,0))  
#    #-- needs to be an indicator for whether her lamb has PN...?
#    
#    if(rbinom(1,1,prob=(1-exp(-(contactset))))==1){
#      DisStat.out<-"E"   
#      NewCt.out<-{temp.in$Count+1}
#      NewDAI.out<-exp(-contactset)
#      NewShedRate.out<-chronicdose
#    } else {
#      DisStat.out<-"S"
#      NewCt.out<-temp.in$Count
#      NewShedRate.out<-0
#      NewDAI.out<-0
#    }
#    
#  } else {    #-- for lambs
#    if(season==0){
#      l<-subset(temp, DemogGrp=="Lamb")
#      g<-dim(l)[1]
#    } else l<-subset(temp, SpatGrpSeason=SpatGrp.in)
#      g<-dim(l)[1]
#    #l<-subset(temp, SpatGrp==SpatGrp.in & DemogGrp=="Lamb")
#    samp<-sample(1:g,Lambcontactnumber,replace=T)
#    contactset<-ifelse(LambOpen==0,0,ifelse(g==0,
#                                            0,
#                                            sum(l$SheddingRate[samp])))
#    contacts<-l$Status[samp]
#    contactIDs<-l$ID[samp]
#    groupsize<-g
#    ewe<-subset(temp,DemogGrp=="Ewe" & ID==temp.in$Mother)
#    contactE<- ifelse(dim(ewe)[1]==0,0,ifelse(ewe$Status=="E"|ewe$Status=="I"|ewe$Status=="C",1,0))  
#    DisStat.out<-ifelse(contactE!=0,"E",
#                        ifelse(rbinom(1,1,prob=(1-exp(-(LambTransmissionProb*contactset))))==1,
#                               "E","S"))
#    NewCt.out<-ifelse(DisStat.out=="E",1,0)
#    NewShedRate.out<-ifelse(DisStat.out=="E",1,ifelse(DisStat.out=="S",0,chronicdose))
#    NewDAI.out<-1 
#  }
#  return(list(DiseaseStatus=DisStat.out,NewCount=NewCt.out,NewSheddingRate=NewShedRate.out))
#         #,NewDAI.out=NewDAI.out,contactset=contacts,contactIDs=contactIDs,groupsize=groupsize))
#}



transmission.fun<-function(season,temp.in,temp,contactnumber,Lambcontactnumber,LambOpen,LambTransmissionProb,chronicdose=chronicdose){
#  k<-switch(temp.in$Status,
#            S = transmission.fun(season,temp.in=temp.in,temp=temp,contactnumber=contactnumber,Lambcontactnumber=Lambcontactnumber,LambOpen=LambOpen,LambTransmissionProb=LambTransmissionProb,chronicdose=chronicdose),
#            E = acutechronic.fun(temp.in,temp,chronicdose=chronicdose,xi=xi,rho=rho),
#            I = acute.fun(Alpha,eta,temp.in),
#            C = chronic.fun(temp.in,temp,Gamma,chronicdecrease,tau,chronicdose),
#            R = recovered.fun(temp.in,Nu))
#  return(k) 
  
  SpatGrp.in<-ifelse(season==0,temp.in$SpatGrp,temp.in$SpatGrpSeason)
  k<-switch(temp.in$DemogGrp,
            Ewe = transmission.fun.adult(season,temp.in=temp.in,temp=temp,contactnumber=contactnumber,chronicdose=chronicdose),
            Lamb = transmission.fun.lamb(season,temp.in=temp.in,temp=temp,Lambcontactnumber=Lambcontactnumber,LambOpen=LambOpen,LambTransmissionProb=LambTransmissionProb)
            )
  return(k)
}

transmission.fun.adult<-function(season,temp.in=temp.in,temp=temp,contactnumber=contactnumber,chronicdose=chronicdose){
#  if(temp.in$DemogGrp=="Ewe"){
#    k<-ifelse(season==0,subset(temp,DemogGrp=="Ewe"),subset(temp,SpatGrpSeason=SpatGrp.in & DemogGrp=="Ewe"))
    if(season==0){
      k<-subset(temp, DemogGrp=="Ewe")
##      g<-dim(k)[1]
    } else {
#      k<-subset(temp, SpatGrpSeason=SpatGrp.in)
      k<-subset(temp, SpatGrpSeason=SpatGrp.in & DemogGrp=="Ewe")
    }
    g<-dim(k)[1]
    samp<-sample(1:g,contactnumber,replace=T)
    contactset<-ifelse(g>=1,sum(k$SheddingRate[samp]),0)
    contacts<-k$Status[samp]
    contactIDs<-k$ID[samp]
    groupsize<-g
    lamb<-subset(temp,DemogGrp=="Lamb" & Mother==temp.in$ID)
    contactL<-ifelse(dim(lamb)[1]==0,0,ifelse(lamb$Status=="I",1,0))  
    #-- needs to be an indicator for whether her lamb has PN...?
    
    if(rbinom(1,1,prob=(1-exp(-(contactset))))==1){
      DisStat.out<-"E"   
      NewCt.out<-{temp.in$Count+1}
      NewDAI.out<-exp(-contactset)
      NewShedRate.out<-chronicdose
    } else {
      DisStat.out<-"S"
      NewCt.out<-temp.in$Count
      NewShedRate.out<-0
      NewDAI.out<-0
    }
return(list(DiseaseStatus=DisStat.out,NewCount=NewCt.out,NewSheddingRate=NewShedRate.out))
  } 


transmission.fun.lamb<-function(season,temp.in=temp.in,temp=temp,Lambcontactnumber=Lambcontactnumber,LambOpen=LambOpen,LambTransmissionProb=LambTransmissionProb){
#else {    #-- for lambs
#  l<-ifelse(season==0,subset(temp,DemogGrp=="Lamb"),subset(temp,SpatGrpSeason=SpatGrp.in & DemogGrp=="Lamb"))
    if(season==0){
      l<-subset(temp, DemogGrp=="Lamb")
     } 
#    else l<-subset(temp, SpatGrpSeason=SpatGrp.in)
    else l<-subset(temp, SpatGrpSeason=SpatGrp.in & DemogGrp=="Lamb")
    
    g<-dim(l)[1]
    #l<-subset(temp, SpatGrp==SpatGrp.in & DemogGrp=="Lamb")
    samp<-sample(1:g,Lambcontactnumber,replace=T)
    contactset<-ifelse(LambOpen==0,0,ifelse(g==0,
                                            0,
                                            sum(l$SheddingRate[samp])))
    contacts<-l$Status[samp]
    contactIDs<-l$ID[samp]
    groupsize<-g
    ewe<-subset(temp,DemogGrp=="Ewe" & ID==temp.in$Mother)
    contactE<- ifelse(dim(ewe)[1]==0,0,ifelse(ewe$Status=="E"|ewe$Status=="I"|ewe$Status=="C",1,0))  
    DisStat.out<-ifelse(contactE!=0,"E",
                        ifelse(rbinom(1,1,prob=(1-exp(-(LambTransmissionProb*contactset))))==1,
                               "E","S"))
    NewCt.out<-ifelse(DisStat.out=="E",1,0)
    NewShedRate.out<-ifelse(DisStat.out=="E",1,0)
#    NewShedRate.out<-ifelse(DisStat.out=="E",1,ifelse(DisStat.out=="S",0,chronicdose))
    NewDAI.out<-1 
#  }
  return(list(DiseaseStatus=DisStat.out,NewCount=NewCt.out,NewSheddingRate=NewShedRate.out))
  #,NewDAI.out=NewDAI.out,contactset=contacts,contactIDs=contactIDs,groupsize=groupsize))
}

acutechronic.fun<-function(temp.in,temp,chronicdose,xi,rho,chronicdecrease){
  change<-ifelse(rbinom(1,1,prob=xi)==1,1,0)
  if(change==1){
    if(temp.in$Count==1 | is.na(temp.in$Count)=="TRUE"){
      k<-rbinom(1,1,prob=rho)
      DisStat.out<-ifelse(k==1,"I","C")
    }
    else {
      DisStat.out<-"C"
    } 
    NewCt.out<-temp.in$Count
    NewShedRate.out<-ifelse(DisStat.out=="I",1,chronicdose)
    NewDAI.out<-temp.in$DoseAtInfection
  } else {
    DisStat.out<-"E"
    NewCt.out<-temp.in$Count
    NewShedRate.out<-chronicdose
    NewDAI.out<-temp.in$DoseAtInfection
  }
  
  return(list(DiseaseStatus=DisStat.out,NewCount=NewCt.out,NewSheddingRate=NewShedRate.out))
    #,NewDAI.out=NewDAI.out))
}

acute.fun<-function(Alpha,eta,temp.in){
     DiseaseStatus<-ifelse(rbinom(1,1,prob=Alpha/(1-Alpha)*(eta))==1,"R","I")
     NewCount<-temp.in$Count
     NewSheddingRate<-ifelse(DiseaseStatus=="I",1,0)
     return(list(DiseaseStatus=DiseaseStatus,NewCount=NewCount,NewSheddingRate=NewSheddingRate))
}

chronic.fun<-function(temp.in,temp,Gamma,chronicdecrease,tau,chronicdose){
  currentgammaC<-chronicdecrease*temp.in$Count*Gamma 
  DiseaseStatus<-ifelse(rbinom(1,1,prob=tau)==1,"I",ifelse(rbinom(1,1,prob=Gamma)==1,"R","C"))
  NewCount<-temp.in$Count
  NewSheddingRate<-ifelse(DiseaseStatus=="C",chronicdose,ifelse(DiseaseStatus=="I",1,0))
  return(list(DiseaseStatus=DiseaseStatus,NewCount=NewCount,NewSheddingRate=NewSheddingRate))
  }

recovered.fun<-function(temp.in,Nu){
  DiseaseStatus<-ifelse(rbinom(1,1,prob=Nu)==1,"S","R")
  NewCount<-temp.in$Count
  NewSheddingRate<-0
  return(list(DiseaseStatus=DiseaseStatus,NewCount=NewCount,NewSheddingRate=NewSheddingRate))
}

survival.fun<-function(temp,individs,PNEweSurvProbs,ChronicEweSurvProbs,EweSurvProbs,Recr,GammaLamb,SLSdecrease){
  SurvivalStat.out<-rep(NA,individs)
  for(j in 1:individs){
    if(temp$DemogGrp[j]=="Lamb"){
      SurvivalStat.out[j]<-ifelse(rbinom(1,1,ifelse(temp$Status[j]=="I",((1-(1-Recr)/365)-(1-(1-SLSdecrease)/(1/GammaLamb))),(1-(1-Recr)/365)))==1,1,0)
    }
    else if(temp$Status[j]=="I"){
      SurvivalStat.out[j]<-ifelse(rbinom(1,1,PNEweSurvProbs[temp$Age[j]])==1,1,0)
    } else if(temp$Status[j]=="C"){
      SurvivalStat.out[j]<-ifelse(rbinom(1,1,ChronicEweSurvProbs[temp$Age[j]])==1,1,0)
    } else SurvivalStat.out[j]<-ifelse(rbinom(1,1,EweSurvProbs[temp$Age[j]])==1,1,0)
  }
  return(SurvivalStat.out)
}

cause.fun<-function(temp,individs,Alpha){
  CauseOut<-rep(NA,individs)
  for(j in 1:individs){
    temp.in<-temp[j,]  
    if(temp.in$DemogGrp=="Lamb"){
      CauseOut[j]<-ifelse(temp.in$Status=="I" & temp.in$StillAlive==0,
                          "PN",ifelse(temp.in$StillAlive==0,
                                      "Other","Alive"))
    }
    else CauseOut[j]<-ifelse(temp.in$StillAlive==1, "Alive",
                             ifelse(temp.in$Status=="I" & rbinom(1,1,prob=Alpha)==1,
                                    "PN","Other"))
  }
  return(CauseOut)
}

birth.fun<-function(BirthWindow,temp,BirthRate,SexRatio){
  NumPotentialMoms<-ifelse(BirthWindow==1,
                           dim(subset(temp,temp$SexRef==0 & temp$StillAlive==1 & temp$HasLamb==0))[1],0)
  PotentialMoms<-subset(temp,temp$SexRef==0 & temp$StillAlive==1 & temp$HasLamb==0)
  
  #-- pick set of potential moms to give birth. --#
  NumNewMoms<-BirthRate*NumPotentialMoms
  momsample<-sample(1:NumPotentialMoms,size=floor(NumNewMoms))
  NewMoms<-PotentialMoms[momsample,]$ID
  MomSpatGrpSeason<-PotentialMoms[momsample,]$SpatGrpSeason
  MomGroupsize<-PotentialMoms[momsample,]$Groupsize
  temp$HasLamb<-apply(as.matrix(temp$ID),1,function(x) mother.fun(x,NewMom=NewMoms))
  NewBirths<-ifelse(BirthWindow==1,length(NewMoms),0)
  
  NewRows<-data.frame(matrix(NA,nrow=NewBirths,ncol=dim(temp)[2]))
  names(NewRows)<-names(temp)
  if(NewBirths!=0){
    NewRows$SexRef<-rbinom(NewBirths,1,SexRatio)
    NewRows$Status<-rep("S",NewBirths)
    NewRows$StillAlive<-rep(1,NewBirths)
    NewRows$ID<-(max(temp$ID)+1):(max(temp$ID)+NewBirths)
    NewRows$Age<-rep(0,NewBirths)
    NewRows$DemogGrp<-rep("Lamb",NewBirths)
    NewRows$Mother<-NewMoms
    NewRows$HasLamb<-rep(0,NewBirths)
    NewRows$SheddingRate<-rep(0,NewBirths)
    NewRows$DoseAtInfection<-rep(NA,NewBirths)
    NewRows$SpatGrp<-1  #-- they need to be in the same spatial groups as their moms...
    NewRows$SpatGrpSeason<-MomSpatGrpSeason
    NewRows$Groupsize<-MomGroupsize
  }
  return(NewRows)
}

analysis.fun<-function(loop.output){
  persistence<-table(unlist(lapply(loop.output,is.null)))["FALSE"]
    #-- expand this function to include annual adult mortality rates and annual --#
    #-- recruitment rates --#
#  LambMort<-EweMort<-ceiling(persistence/365)
#    for(i in 1:length(LambMort)){
#      LambMort[i]<-sum(dim(loop.output[[i]]$DeathMat)[1]
#    }
  
#  N<-ChronicCount<-AcuteCount<-SCount<-RCount<-ECount<-mortality<-rep(NA,persistence)
#  contactsets<-groupsize<-vector("list",persistence)
  N<-ChronicCount<-AcuteCount<-SCount<-RCount<-ECount<-rep(NA,persistence)
  groupsize<-vector("list",persistence)
  
  for(j in 2:persistence){
#    mortality[j]<-dim(loop.output[[j]]$DeathMat)[1]
    N[j]<-dim(loop.output[[j]]$TimestepData)[1]
    ChronicCount[j]<-ifelse(is.na(table(loop.output[[j]]$TimestepData$Status)["C"])==TRUE,
                            0,table(loop.output[[j]]$TimestepData$Status)["C"])
    AcuteCount[j]<-ifelse(is.na(table(loop.output[[j]]$TimestepData$Status)["I"])==TRUE,
                          0,table(loop.output[[j]]$TimestepData$Status)["I"])
    SCount[j]<-ifelse(is.na(table(loop.output[[j]]$TimestepData$Status)["S"])==TRUE,
                      0,table(loop.output[[j]]$TimestepData$Status)["S"])
    RCount[j]<-ifelse(is.na(table(loop.output[[j]]$TimestepData$Status)["R"])==TRUE,
                      0,table(loop.output[[j]]$TimestepData$Status)["R"])
    ECount[j]<-ifelse(is.na(table(loop.output[[j]]$TimestepData$Status)["E"])==TRUE,
                      0,table(loop.output[[j]]$TimestepData$Status)["E"])
#    contactsets[[j]]<-loop.output[[j]][[3]]
    groupsize[[j]]<-table(loop.output[[j]]$TimestepData$Groupsize)
  }
  return(list(N=N,ECount=ECount,ChronicCount=ChronicCount,AcuteCount=AcuteCount,SCount=SCount,RCount=RCount,persistence=persistence,groupsize=groupsize))
#  return(list(N=N,ECount=ECount,ChronicCount=ChronicCount,AcuteCount=AcuteCount,SCount=SCount,RCount=RCount,mortality=mortality,persistence=persistence,contactsets=contactsets,groupsize=groupsize))
  #  return(N=N,ECount=ECount,ChronicCount=ChronicCount,AcuteCount=AcuteCount,SCount=SCount,RCount=RCount,mortality=mortality,persistence=persistence,contactsets=contactsets,groupsize=groupsize)
}



wrapper.fun<-function(ParamMat){
  model=ParamMat$model
  timesteps=ParamMat$timesteps
  BirthRate=ParamMat$BirthRate
  SexRatio=ParamMat$SexRatio
  tau = ParamMat$tau
  eta = ParamMat$eta
  rho = ParamMat$rho
  xi = ParamMat$xi
  n=ParamMat$n
  Recr=ParamMat$Recr
  SLSdecrease=ParamMat$SLSdecrease
  GammaLamb=ParamMat$GammaLamb
  Nu=ParamMat$Nu
  Gamma=ParamMat$Gamma
  Alpha=ParamMat$Alpha
  AlphaChronic=ParamMat$AlphaChronic
  contactnumber=ParamMat$contactnumber
  LambcontactnumberIn=ParamMat$Lambcontactnumber
  chronicdose=ParamMat$chronicdose
  LambTransmissionProb=ParamMat$LambTransmissionProb
  PropRecovered=ParamMat$PropRecovered
  chronicdecrease=ParamMat$chronicdecrease
  ngroups=ParamMat$ngroups
  
  SimOut<-IndividualSIR(model=model,timesteps=timesteps,
                        BirthRate=BirthRate,
                        SexRatio=SexRatio,
                        tau=tau,
                        eta=eta,
                        rho=rho,
                        xi=xi,
                        n=n,
                        Recr=Recr,
                        SLSdecrease=SLSdecrease,
                        GammaLamb=GammaLamb,
                        Nu=Nu,
                        Gamma=Gamma,
                        Alpha=Alpha,
                        AlphaChronic=AlphaChronic,
                        contactnumber=contactnumber,
                        LambcontactnumberIn=LambcontactnumberIn,
                        chronicdose=chronicdose,
                        LambTransmissionProb=LambTransmissionProb,
                        PropRecovered=PropRecovered,
                        chronicdecrease=chronicdecrease,
                        ngroups=ngroups)
  
  return(SimOut)
}


#-- disease status update function --#

disease.status<-function(temp.in,temp,individs,Alpha,eta,xi,rho,season,Lambcontactnumber,contactnumber,LambOpen,LambTransmissionProb,chronicdose,Nu,chronicdecrease,Gamma,tau){ 
  k<-switch(temp.in$Status,
   S = transmission.fun(season,temp.in=temp.in,temp=temp,contactnumber=contactnumber,Lambcontactnumber=Lambcontactnumber,LambOpen=LambOpen,LambTransmissionProb=LambTransmissionProb,chronicdose=chronicdose),
   E = acutechronic.fun(temp.in,temp,chronicdose=chronicdose,xi=xi,rho=rho),
   I = acute.fun(Alpha,eta,temp.in),
   C = chronic.fun(temp.in,temp,Gamma,chronicdecrease,tau,chronicdose),
   R = recovered.fun(temp.in,Nu))
  return(k)
}

require(stats)  

#-- should switch this to with. --#.
IndividualSIR<-function(model=model,timesteps=timesteps,
                        BirthRate=BirthRate,
                        eta=eta,
                        xi=xi,
                        rho=rho,
                        tau=tau,
                        SexRatio=SexRatio,
                        n=n,
                        Recr=Recr,
                        SLSdecrease=SLSdecrease,
                        GammaLamb=GammaLamb,
                        Nu=Nu,
                        Gamma=Gamma,
                        Alpha=Alpha,
                        AlphaChronic=AlphaChronic,
                        contactnumber=contactnumber,
                        LambcontactnumberIn=LambcontactnumberIn,
                        chronicdose=chronicdose,
                        LambTransmissionProb=LambTransmissionProb,
                        PropRecovered,
                        chronicdecrease,
                        ngroups=ngroups){
  
  StorageList<-vector("list",timesteps)
  AdultMort<-LambMort<-LambBirths<-rep(NA,timesteps)
  #  StorageList<-matrix(NA,ncol=17,nrow=timesteps*(N+100))
  ID<-1:n		
  OldEweSurvProbs<-c(rep(.833,3),rep(.945,4),rep(.850,13),0)	
  SexRef<-rep(0,n)	
  
  PNEweSurvProbs<-c(rep(1-(1-OldEweSurvProbs*(1-Alpha))/365,each=365),0)	
  ChronicEweSurvProbs<-c(rep(1-(1-OldEweSurvProbs*(1-AlphaChronic))/365,each=365),0)	
  EweSurvProbs<-c(rep(1-(1-OldEweSurvProbs)/365,each=365),0)	
  RefTS1<-as.data.frame(cbind(ID,SexRef))
  
  InitEweSurvProb<-c(.833,.833^2,.833^3,.833^3*.945,.833^3*.945^2,.833^3*.945^3,.833^3*.945^4,.833^3*.945*.850,.833^3*.945*.850^2,.833^3*.945*.850^3,.833^3*.945*.850^4,.833^3*.945*.850^5,.833^3*.945*.850^6,.833^3*.945*.850^7,.833^3*.945*.850^8,.833^3*.945*.850^9,.833^3*.945*.850^10,.833^3*.945*.850^11,.833^3*.945*.850^12,.833^3*.945*.850^13,0)
  NormEweSurvProb<-InitEweSurvProb/sum(InitEweSurvProb)
  
  RefTS1$OldAge<-sample(1:21,n,NormEweSurvProb,replace=T)-1	
  RefTS1$Age<-(RefTS1$OldAge)*365 + 182
  RefTS1$StillAlive<-rbinom(dim(RefTS1)[1],1,.9)	
  RefTS1$Status<-c(rep("I",1),rep("R",floor(PropRecovered*n)),rep("S",n-floor(PropRecovered*n)-1))	
  RefTS1$Cause<-rep(0,n)	
  RefTS1$SheddingRate<-c(rep(1,1),rep(0,n-1))   
  RefTS1$Count<-c(rep(1,1),rep(0,n-1))     
  RefTS1$Mother<-RefTS1$DemogGrp<-rep(NA,n)
  
  momoptions<-subset(RefTS1,Age>=365)$ID
  
  RefTS1$Mother<-sapply(RefTS1$Age,function(x) {ifelse(x>=365,NA,sample(momoptions,1))})
  RefTS1$DemogGrp<-sapply(RefTS1$Age,function(x) {ifelse(x>=365,"Ewe","Lamb")})
#  for(i in 1:n){
#    RefTS1$Mother[i]<-ifelse(RefTS1$Age[i]>=365,NA,sample(subset(RefTS1,Age>=365)$ID,1))
#    RefTS1$DemogGrp[i]<-ifelse(RefTS1$Age[i]>=365,"Ewe","Lamb")
#  }
  RefTS1$SpatGrp<-rep(1,n)
  RefTS1$SpatGrpSeason<-floor(runif(n,min=0,max=ngroups))
  RefTS1$DoseAtInfection<-rep(NA,n)
  RefTS1$Groupsize<-rep(n,n)
  StorageList[[1]][[1]]<-RefTS1
  
  LambWindow<-c(210,300)	
  
  #----------------------------------------------------#
  #------ Individual status update for-loop -----------#
  #----------------------------------------------------#
  flag=1
  for(i in 2:timesteps){
    #      for(i in 2:209){
    
    if(flag<1)
      break;
    
    #-- 1) Create a storage object for information from the previous timestep.
    temp<-StorageList[[i-1]][[1]]
    
    Age<-as.numeric(as.character(temp$Age))
    Sex<-temp$SexRef
    #-- Classify all individuals to a demographic group, Lamb, Ewe or Ram.
    DemogGrpNew<-sapply(Age,function(x) DemogFun1(x))
    
    temp$DemogGrp<-DemogGrpNew
    individs<-dim(temp)[1]
    
    #-- Loop through all individuals to update disease statuses.
    DiseaseStatus<-NewCount<-NewSheddingRate<-NewDoseAtInfection<-groupsizeNew<-rep(NA,individs)
    
    EweGroup<-subset(temp,DemogGrp=="Ewe")
    LambGroup<-subset(temp,DemogGrp=="Lamb")
    
    LambOpen<-lambtransmission.mod.fun(i)
    Lambcontactnumber<-min(dim(LambGroup)[1],LambcontactnumberIn)
    season<-ifelse(i %% 365 <=60, 0,1)
    
    #-- store contactsets --#
    contactsets<-vector("list",individs)
    
    ptm<-proc.time()
    for(j in 1:individs){
      temp.in<-temp[j,]
      k<-disease.status(temp.in,temp,individs,Alpha,eta,xi,rho,season,Lambcontactnumber,contactnumber,LambOpen,LambTransmissionProb,chronicdose,Nu,chronicdecrease,Gamma,tau)
      DiseaseStatus[j]<-k$DiseaseStatus
      NewCount[j]<-k$NewCount
      NewSheddingRate[j]<-k$NewSheddingRate
      
      if(season==0){
        groupsizeNew[j]<-individs
      } else groupsizeNew[j]<-dim(subset(temp,SpatGrpSeason==temp$SpatGrpSeason[j]))[1]
    }
    
    print(proc.time()-ptm)
    
    temp$Status<-DiseaseStatus
    temp$Count<-NewCount
    temp$SheddingRate<-NewSheddingRate
    temp$Groupsize<-groupsizeNew
        
    #-- 3) Update individual ages and "Alive" statuses
    ptm2<-proc.time()
    NewAge<-as.numeric(as.character(temp$Age))+1
    temp$Age<-NewAge
    temp$StillAlive<-survival.fun(temp,individs,PNEweSurvProbs,ChronicEweSurvProbs,EweSurvProbs,Recr,GammaLamb,SLSdecrease)
    temp$Cause<-cause.fun(temp,individs,Alpha)
    print(proc.time()-ptm2)
    #-- 4) Generate new individuals through birth process.
    ptm3<-proc.time()
    BirthWindow<-birth.mod.fun(i)
    NumPotentialMoms<-ifelse(BirthWindow==1,dim(subset(temp,temp$SexRef==0 & temp$StillAlive==1 & temp$HasLamb==0))[1],0)
    PotentialMoms<-subset(temp,temp$SexRef==0 & temp$StillAlive==1 & temp$HasLamb==0)
    
    #-- pick set of potential moms to give birth. --#
    NumNewMoms<-BirthRate*NumPotentialMoms
    momsample<-sample(1:NumPotentialMoms,size=floor(NumNewMoms))
    NewMoms<-PotentialMoms[momsample,]$ID
    MomSpatGrpSeason<-PotentialMoms[momsample,]$SpatGrpSeason
    MomGroupsize<-PotentialMoms[momsample,]$Groupsize
    temp$HasLamb<-apply(as.matrix(temp$ID),1,function(x) mother.fun(x,NewMom=NewMoms))
    NewBirths<-ifelse(BirthWindow==1,length(NewMoms),0)
    
    NewRows<-data.frame(matrix(NA,nrow=NewBirths,ncol=dim(temp)[2]))
    names(NewRows)<-names(temp)
    if(NewBirths!=0){
      NewRows$SexRef<-rbinom(NewBirths,1,SexRatio)
      NewRows$Status<-rep("S",NewBirths)
      NewRows$StillAlive<-rep(1,NewBirths)
      NewRows$ID<-(max(temp$ID)+1):(max(temp$ID)+NewBirths)
      NewRows$Age<-rep(0,NewBirths)
      NewRows$DemogGrp<-rep("Lamb",NewBirths)
      NewRows$Mother<-NewMoms
      NewRows$HasLamb<-rep(0,NewBirths)
      NewRows$SheddingRate<-rep(0,NewBirths)
      NewRows$DoseAtInfection<-rep(NA,NewBirths)
      NewRows$SpatGrp<-1  #-- they need to be in the same spatial groups as their moms...
      NewRows$SpatGrpSeason<-MomSpatGrpSeason
      NewRows$Groupsize<-MomGroupsize
    }
    
    #    time3[i]<-(proc.time()-ptm3)$elapsed
    print(proc.time()-ptm3)
    #-- 5) Make death mat of who died in this timestep, and how.
    deathnames<-c("StillAlive","DiseaseStatus","CauseOfDeath","DemogGrp","Age")
    DeathMata<-as.data.frame(cbind(temp$StillAlive, temp$Status, temp$Cause,temp$DemogGrp,temp$Age))
    names(DeathMata)<-deathnames
    DeathMat<-subset(DeathMata,StillAlive==0)
    TimestepData=data.frame(rbind(subset(temp,StillAlive==1),NewRows))
    
    names(TimestepData)<-names(temp)
    
    #-- 6) Store temp in StorageList.
#    StorageList[[i]]<-list(TimestepData=TimestepData,DeathMat=data.frame(DeathMat),contactsets=contactsets)
    StorageList[[i]]<-list(TimestepData=TimestepData)
    levels(StorageList[[i]][[1]]$Status)<-c("S","E","I","C","R")
    flag<-as.numeric(ifelse(is.na(table(StorageList[[i]][[1]]$Status)["I"])==TRUE,0,1))+as.numeric(ifelse(is.na(table(StorageList[[i]][[1]]$Status)["C"])==TRUE,0,1))+as.numeric(ifelse(is.na(table(StorageList[[i]][[1]]$Status)["E"])==TRUE,0,1))
    
    AdultMort[i]<-dim(subset(DeathMat,DemogGrp=="Ewe"))[1]
    LambMort[i]<-dim(subset(DeathMat,DemogGrp=="Lamb"))[1]
    LambBirths[i]<-dim(NewRows)[1]
    
    print(i)
  }
  
  #-- Analyze for-loop output --#
  ptm4<-proc.time()
  output.analysis<-analysis.fun(StorageList)
  persistence<-output.analysis$persistence; N<-output.analysis$N
  ChronicCount<-output.analysis$ChronicCount;AcuteCount<-output.analysis$AcuteCount;SCount<-output.analysis$SCount
  RCount<-output.analysis$RCount; mortality<-output.analysis$mortality; ECount<-output.analysis$ECount
#  contactsets<-output.analysis$contactsets; 
  groupsize<-output.analysis$groupsize
  print(proc.time()-ptm4)
  
  return(list(persistence=persistence,  #-- scalar
              ECount=ECount,#-- vector (variable length)
              ChronicCount=ChronicCount,#-- vector (variable length)
              AcuteCount=AcuteCount,
              SCount=SCount,
              RCount=RCount,
              N=N,
              AdultMort,
              LambMort,
              LambBirths,
#              mortality=mortality,
#              contactsets=contactsets,
              groupsize=groupsize
              #              time1=time1,
              #              time2=time2,
              #              time3=time3,
              # 		              StorageList=StorageList #-- 
  ))
}


#-- take a look at plyr. reshape package as well --#
#-- write a big csv for data storage; write list for everyone else --#
#-- help(package="plyr")
#-- 
#-- gc() #-- garbage collection
#-- Ideas for output clean-up: --#
#-- write each element of the list to its own storage object.  
#-- write each storage object once at the END of the simulation.
