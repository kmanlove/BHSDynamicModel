#-- source functions for Cumulative Exposure Model _ 11 May 2012 --#

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
 
      
transmission.fun<-function(season,temp.in,temp,contactnumber,Lambcontactnumber,LambOpen,LambTransmissionProb,chronicdose=chronicdose){
  SpatGrp.in<-ifelse(season==0,temp.in$SpatGrp,temp.in$SpatGrpSeason)
          if(temp.in$DemogGrp=="Ewe"){
            if(season==0){
              k<-subset(temp, DemogGrp=="Ewe")
            } else k<-subset(temp, SpatGrpSeason=SpatGrp.in)
                        
              samp<-sample(1:dim(k)[1],contactnumber,replace=T)
              contactset<-ifelse(dim(k)[1]>=1,sum(k$SheddingRate[samp]),0)
              contacts<-k$Status[samp]
              contactIDs<-k$ID[samp]
              groupsize<-dim(k)[1]
              lamb<-subset(temp,DemogGrp=="Lamb" & Mother==temp.in$ID)
              contactL<-ifelse(dim(lamb)[1]==0,0,ifelse(lamb$Status=="I",1,0))  
              #-- needs to be an indicator for whether her lamb has PN...?
 
          if(rbinom(1,1,prob=(1-exp(-(contactset))))==1){
            DisStat.out<-"E"   
            NewCt.out<-temp.in$Count+1
            NewDAI.out<-exp(-contactset)
            NewShedRate.out<-chronicdose
              } else {
                  DisStat.out<-"S"
                  NewCt.out<-temp.in$Count
                  NewShedRate.out<-0
                  NewDAI.out<-0
                }
              
          } else {    #-- for lambs
            if(season==0){
              l<-subset(temp, DemogGrp=="Lamb")
            } else l<-subset(temp, SpatGrpSeason=SpatGrp.in)
 
              #l<-subset(temp, SpatGrp==SpatGrp.in & DemogGrp=="Lamb")
              samp<-sample(1:dim(l)[1],Lambcontactnumber,replace=T)
              contactset<-ifelse(LambOpen==0,0,ifelse(dim(l)[1]==0,
                                 0,
                                 sum(l$SheddingRate[samp])))
              contacts<-l$Status[samp]
              contactIDs<-l$ID[samp]
              groupsize<-dim(l)[1]
              ewe<-subset(temp,DemogGrp=="Ewe" & ID==temp.in$Mother)
              contactE<- ifelse(dim(ewe)[1]==0,0,ifelse(ewe$Status=="E"|ewe$Status=="I"|ewe$Status=="C",1,0))  
              DisStat.out<-ifelse(contactE!=0,"E",
                                       ifelse(rbinom(1,1,prob=(1-exp(-(LambTransmissionProb*contactset))))==1,
                                              "E","S"))
              NewCt.out<-ifelse(DisStat.out=="E",1,0)
              NewShedRate.out<-ifelse(DisStat.out=="E",1,ifelse(DisStat.out=="S",0,chronicdose))
              NewDAI.out<-1 
              }
              return(list(DisStat.out=DisStat.out,NewCt.out=NewCt.out,NewShedRate.out=NewShedRate.out,NewDAI.out=NewDAI.out,contactset=contacts,contactIDs=contactIDs,groupsize=groupsize))
    		}

      
acutechronic.fun<-function(temp.in,temp,chronicdose,xi,rho){
    change<-ifelse(rbinom(1,1,prob=xi)==1,1,0)
           if(change==1){
            k<-rbinom(1,1,prob=rho)
            DisStat.out<-ifelse(k==1,"I","C")
            NewCt.out<-temp.in$Count
            NewShedRate.out<-ifelse(DisStat.out=="I",1,chronicdose)
            NewDAI.out<-temp.in$DoseAtInfection
          } else {
            DisStat.out<-"E"
            NewCt.out<-temp.in$Count
            NewShedRate.out<-chronicdose
            NewDAI.out<-temp.in$DoseAtInfection
          }
        
    return(list(DisStat.out=DisStat.out,NewCt.out=NewCt.out,NewShedRate.out=NewShedRate.out,NewDAI.out=NewDAI.out))
}
      
      
survival.fun<-function(temp,PNEweSurvProbs,ChronicEweSurvProbs,EweSurvProbs,Recr,GammaLamb,SLSdecrease){
  SurvivalStat.out<-rep(NA,dim(temp)[1])
  	for(j in 1:dim(temp)[1]){
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
            
cause.fun<-function(temp,Alpha){
  	CauseOut<-rep(NA,dim(temp)[1])
			for(j in 1:dim(temp)[1]){
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
  
  N<-ChronicCount<-AcuteCount<-SCount<-RCount<-ECount<-mortality<-rep(NA,persistence)
  contactsets<-groupsize<-vector("list",persistence)
  
	for(j in 2:persistence){
    mortality[j]<-dim(loop.output[[j]]$DeathMat)[1]
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
    contactsets[[j]]<-loop.output[[j]][[3]]
    groupsize[[j]]<-table(loop.output[[j]]$TimestepData$Groupsize)
	}
  return(list(N=N,ECount=ECount,ChronicCount=ChronicCount,AcuteCount=AcuteCount,SCount=SCount,RCount=RCount,mortality=mortality,persistence=persistence,contactsets=contactsets,groupsize=groupsize))
}
      
      

wrapper.fun<-function(ParamMat){
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

  SimOut<-IndividualSIR(timesteps=timesteps,
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
      
disease.status<-function(x,temp,temp.names,Alpha,eta,xi,rho,season,Lambcontactnumber,contactnumber,LambOpen,LambTransmissionProb,chronicdose,Nu,chronicdecrease,Gamma,tau){
      temp.in<-as.data.frame(t(x))
      colnames(temp.in)<-temp.names
      
        #-- Transmission
        if(temp.in$Status=="S"){
          transmit.out<-transmission.fun(season,temp.in=temp.in,temp=temp,contactnumber=contactnumber,Lambcontactnumber=Lambcontactnumber,LambOpen=LambOpen,LambTransmissionProb=LambTransmissionProb,chronicdose=chronicdose)
          DiseaseStatus<-transmit.out$DisStat.out
          NewCount<-transmit.out$NewCt.out
          NewDoseAtInfection<-transmit.out$NewDAI.out
          NewSheddingRate<-transmit.out$NewShedRate.out
          contactsets<-transmit.out$contactset
    		}
        
        #-- from incubatory to acute or chronic --#
        else if (temp.in$Status=="E"){
          acutechronic.out<-acutechronic.fun(temp.in,temp,chronicdose=chronicdose,xi=xi,rho=rho)
          DiseaseStatus<-acutechronic.out$DisStat.out
          NewCount<-acutechronic.out$NewCt.out
          NewSheddingRate<-acutechronic.out$NewShedRate.out
          NewDoseAtInfection<-acutechronic.out$NewDAI.out
        }
        
        #-- Recovery from Acute #-- maybe add death from acute as well. 
        else if(temp.in$Status=="I"){  #-- can either recover or stay in I. 
  					DiseaseStatus<-ifelse(rbinom(1,1,prob=Alpha/(1-Alpha)*(eta))==1,"R","I")
            NewCount<-temp.in$Count
            NewSheddingRate<-ifelse(DiseaseStatus=="I",1,0)
					}
        #-- need to add transition from acute to chronic. 
        
        #-- Recovery/retention in Chronic #-- maybe add death from chronic as well. 
        else if(temp.in$Status=="C"){  
            currentgammaC<-chronicdecrease*temp.in$Count*Gamma 
						DiseaseStatus<-ifelse(rbinom(1,1,prob=tau)==1,"A",ifelse(rbinom(1,1,prob=Gamma)==1,"R","C"))
            NewCount<-temp.in$Count
            NewSheddingRate<-ifelse(DiseaseStatus=="C",chronicdose,ifelse(DiseaseStatus=="A",1,0))
					}
        #-- Waning from recovered back to S --#
        else{
          DiseaseStatus<-ifelse(rbinom(1,1,prob=Nu)==1,"S","R")
          NewCount<-temp.in$Count
          NewSheddingRate<-0
        }  
        #-- record groupsize --#
        if(season==0){
          groupsizeNew<-dim(temp)[1]
        } else groupsizeNew<-dim(subset(temp,SpatGrpSeason==temp.in$SpatGrpSeason))[1]
      return(DiseaseStatus)
      }
      
  