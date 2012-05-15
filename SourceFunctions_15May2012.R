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
  b<-ifelse(t %% 365==209,1,0)
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
 
      
transmission.fun<-function(season,temp.in,temp,contactnumber,lambcontactnumber,LambOpen,LambTransmissionProb,chronicdose=chronicdose){
        SpatGrp.in<-ifelse(season==0,temp.in$SpatGrp,temp.in$SpatGrpSeason)
          if(temp.in$DemogGrp=="Ewe"){
              k<-subset(temp, SpatGrp==SpatGrp.in & DemogGrp=="Ewe")
              contactsetE<-sum(k$SheddingRate[sample(1:dim(k)[1],contactnumber,replace=T)])
              lamb<-subset(temp,DemogGrp=="Lamb" & Mother==temp.in$ID)
              contactL<-ifelse(dim(lamb)[1]==0,0,ifelse(lamb$Status=="I",1,0))  
              #-- needs to be an indicator for whether her lamb has PN...?
 
          if(rbinom(1,1,prob=(1-exp(-(contactsetE))))==1){
            # DiseaseStatus[j]<-ifelse(rbinom(1,1,prob=(1-exp(-(contactsetE))))==1,"I","C")   
            DisStat.out<-"E"   
            #-- need propchronic to reflect dosage
            NewCt.out<-temp.in$Count+1
            NewDAI.out<-exp(-contactsetE)
#            NewSheddingRate[j]<-ifelse(DiseaseStatus[j]=="I",1,ifelse(DiseaseStatus[j]=="S",0,chronicdose))
            NewShedRate.out<-chronicdose
              } else {
                  DisStat.out<-"S"
                  NewCt.out<-temp.in$Count
                  NewShedRate.out<-0
                  NewDAI.out<-0
                }
              
          } else {    #-- for lambs
              l<-subset(temp, SpatGrp==SpatGrp.in & DemogGrp=="Lamb")
              contactsetL<-ifelse(LambOpen==0,0,ifelse(dim(l)[1]==0,
                                 0,
                                 sum(l$SheddingRate[sample(1:dim(l)[1],Lambcontactnumber,replace=T)])))
              ewe<-subset(temp,DemogGrp=="Ewe" & ID==temp.in$Mother)
              contactE<- ifelse(dim(ewe)[1]==0,0,ifelse(ewe$Status=="E"|ewe$Status=="I"|ewe$Status=="C",ewe$SheddingRate,0))  
              DisStat.out<-ifelse(contactE!=0,"I",
                                       ifelse(rbinom(1,1,prob=(1-exp(-(LambTransmissionProb*contactsetL))))==1,
                                              "I","S"))
              NewCt.out<-ifelse(DisStat.out=="I",1,0)
              NewShedRate.out<-ifelse(DisStat.out=="I",1,ifelse(DisStat.out=="S",0,chronicdose))
              NewDAI.out<-1 
              }
              return(list(DisStat.out=DisStat.out,NewCt.out=NewCt.out,NewShedRate.out=NewShedRate.out,NewDAI.out=NewDAI.out))
    		}

      
acutechronic.fun<-function(temp.in,temp,chronicdose,xi){
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
      
      
survival.fun<-function(temp,PNLambSurvProb,PNEweSurvProbs,ChronicEweSurvProbs,EweSurvProbs,LambSurvProb){
    SurvivalStat.out<-rep(NA,dim(temp)[1])
  	for(j in 1:dim(temp)[1]){
#      temp.in<-temp[j,]
#       if(temp.in$DemogGrp=="Lamb"){
#				SurvivalStat.out[j]<-ifelse(rbinom(1,1,ifelse(temp.in$Status=="I",PNLambSurvProb,LambSurvProb))==1,1,0)
#			}
#			else{
#					SurvivalStat.out[j]<-ifelse(rbinom(1,1,prob=(ifelse(temp.in$Status=="I",
#                                   PNEweSurvProbs[temp.in$Age],
#                                   ifelse(temp.in$Status=="C",
#                                         ChronicEweSurvProbs[temp.in$Age],
#                                         EweSurvProbs[temp.in$Age]))))==1,1,0)
#			}
          if(temp$DemogGrp[j]=="Lamb"){
  			SurvivalStat.out[j]<-ifelse(rbinom(1,1,ifelse(temp$Status[j]=="I",PNLambSurvProb,LambSurvProb))==1,1,0)
			}
			else if(temp$Status[j]=="I"){
        SurvivalStat.out[j]<-ifelse(rbinom(1,1,PNEweSurvProbs[temp$Age[j]])==1,1,0)
			} else if(temp$Status[j]=="C"){
        SurvivalStat.out[j]<-ifelse(rbinom(1,1,ChronicEweSurvProbs[temp$Age[j]])==1,1,0)
			} else SurvivalStat.out[j]<-ifelse(rbinom(1,1,EweSurvProbs[temp$Age[j]])==1,1,0)
#					SurvivalStat.out[j]<-ifelse(rbinom(1,1,ifelse(temp$Status[j]=="I",PNEweSurvProbs[temp$Age[j]],ifelse(temp$Status[j]=="C", ChronicEweSurvProbs[temp$Age[j]],EweSurvProbs[temp$Age[j]])))==1,1,0)
#			}
  	}
    return(SurvivalStat.out)
}
      
#      survival.out<-survival.fun(temp=temp,PNLambSurvProb,PNEweSurvProbs,ChronicEweSurvProbs,EweSurvProbs,LambSurvProb)
      
cause.fun<-function(temp,InPNAdultSurvAdj){
  	CauseOut<-rep(NA,dim(temp)[1])
			for(j in 1:dim(temp)[1]){
        temp.in<-temp[j,]  
  			if(temp.in$DemogGrp=="Lamb"){
					CauseOut[j]<-ifelse(temp.in$Status=="I" & temp.in$StillAlive==0,
                           "PN",ifelse(temp.in$StillAlive==0,
                                       "Other","Alive"))
				}
				else CauseOut[j]<-ifelse(temp.in$StillAlive==1, "Alive",
                              ifelse(temp.in$Status=="I" & rbinom(1,1,prob=InPNAdultSurvAdj)==1,
                                     "PN","Other"))
			}
    return(CauseOut)
}
      
birth.fun<-function(i,temp,BirthRate){
    BirthWindow<-birth.mod.fun(i)
    NumPotentialMoms<-ifelse(BirthWindow==1,
                             dim(subset(temp,temp$SexRef==0 & temp$StillAlive==1 & temp$HasLamb==0))[1],0)
    PotentialMoms<-subset(temp,temp$SexRef==0 & temp$StillAlive==1 & temp$HasLamb==0)
    
    #-- pick set of potential moms to give birth. --#
    NumNewMoms<-BirthRate*NumPotentialMoms
    momsample<-sample(1:NumPotentialMoms,size=floor(NumNewMoms))
    NewMoms<-PotentialMoms[momsample,]$ID
    MomSpatGrpSeason<-PotentialMoms[momsample,]$SpatGrpSeason
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
		}
  return(NewRows)
}
    
analysis.fun<-function(loop.output){
  persistence<-table(unlist(lapply(loop.output,is.null)))["FALSE"]
  
  N<-ChronicCount<-AcuteCount<-SCount<-RCount<-mortality<-rep(NA,persistence)
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
	}
  return(list(N=N,ChronicCount=ChronicCount,AcuteCount=AcuteCount,SCount=SCount,RCount=RCount,mortality=mortality,persistence=persistence))
}