require(stats)	
source("work/StatProjects/Raina/sheep/Papers/DynamicModel/Code/IndividualTrackingModel/TwoSeasonModels_11May2012/TwoSeasonsGitRepo/SourceFunctions_15May2012")

#-- test changes for git. 
#-- test 2.

IndividualSIR<-function(timesteps=timesteps,
                        BirthRate=BirthRate,
                        SexRatio=SexRatio,
                        n=n,
                        K=K,
                        InPNAdultSurvAdj=InPNAdultSurvAdj,
                        InPNLambSurvProb=InPNLambSurvProb,
                        InLambSurvProb=InLambSurvProb,
                        Gamma=Gamma,Nu=Nu,
                        GammaChronic=GammaChronic,
                        Alpha=Alpha,
                        AlphaLamb=AlphaLamb,
                        AlphaChronic=AlphaChronic,
                        contactnumber=contactnumber,
                        LambForce=LambForce,
                        LambcontactnumberIn=LambcontactnumberIn,
                        chronicdose=chronicdose,
                        LambTransmissionProb=LambTransmissionProb,
                        LambEweTranmissionProb,
                        PropRecovered,
                        chronicdecrease){
  
	StorageList<-vector("list",timesteps)	
	ID<-1:n		
	OldEweSurvProbs<-c(rep(.833,3),rep(.945,4),rep(.850,13),0)	
  PNOldEweSurvProbs<-c(rep(.833,3),rep(.945,4),rep(.850,13),0)*InPNAdultSurvAdj	
  ChronicOldEweSurvProbs<-c(rep(.833,3),rep(.945,4),rep(.850,13),0)*AlphaChronic	
  OldLambSurvProb<-.6	
	SexRef<-rep(0,n)	
	Gamma=Gamma

	PNEweDeathProbs<-(1-PNOldEweSurvProbs)/365
	ChronicEweDeathProbs<-(1-ChronicOldEweSurvProbs)/365
	EweDeathProbs<-(1-OldEweSurvProbs)/365
	LambDeathProb<-(1-InLambSurvProb)/365
	PNLambDeathProb<-	((1-InPNLambSurvProb)/365)^(1/(1/Gamma))

	PNEweSurvProbs<-c(rep(1-PNEweDeathProbs,each=365),0)	
  ChronicEweSurvProbs<-c(rep(1-ChronicEweDeathProbs,each=365),0)	
  EweSurvProbs<-c(rep(1-EweDeathProbs,each=365),0)	
  LambSurvProb<-1-LambDeathProb	
  PNLambSurvProb<-1-PNLambDeathProb
	RefTS1<-as.data.frame(cbind(ID,SexRef))

  InitEweSurvProb<-c(.833,.833^2,.833^3,.833^3*.945,.833^3*.945^2,.833^3*.945^3,.833^3*.945^4,.833^3*.945*.850,.833^3*.945*.850^2,.833^3*.945*.850^3,.833^3*.945*.850^4,.833^3*.945*.850^5,.833^3*.945*.850^6,.833^3*.945*.850^7,.833^3*.945*.850^8,.833^3*.945*.850^9,.833^3*.945*.850^10,.833^3*.945*.850^11,.833^3*.945*.850^12,.833^3*.945*.850^13,0)
	NormEweSurvProb<-InitEweSurvProb/sum(InitEweSurvProb)

	RefTS1$OldAge<-sample(1:21,n,NormEweSurvProb,replace=T)-1	
	RefTS1$Age<-RefTS1$OldAge*365	
	RefTS1$StillAlive<-rbinom(dim(RefTS1)[1],1,.9)	
  RefTS1$Status<-c(rep("I",1),rep("R",floor(PropRecovered*n)),rep("S",n-floor(PropRecovered*n)-1))	
  RefTS1$Cause<-rep(0,n)	
  RefTS1$SheddingRate<-c(rep(1,1),rep(0,n-1))   
  RefTS1$Count<-c(rep(1,1),rep(0,n-1))   
  RefTS1$Mother<-rep(NA,length(RefTS1$Age)) 
  RefTS1$SpatGrp<-rep(1,length(RefTS1$Age))
  RefTS1$DoseAtInfection<-rep(NA,length(RefTS1$Age))
  
	StorageList[[1]][[1]]<-RefTS1

	LambWindow<-c(210,300)	
  
	#----------------------------------------------------#
	#------ Individual status update for-loop -----------#
	#----------------------------------------------------#
flag=1
	for(i in 2:timesteps){
    if(flag<1)
      break;
      
	#-- 1) Create a storage object for information from the previous timestep.
	temp<-StorageList[[i-1]][[1]]

		Age<-as.numeric(as.character(temp$Age))
		Sex<-temp$SexRef
			#-- Classify all individuals to a demographic group, Lamb, Ewe or Ram.
			temp$DemogGrp<-rep(NA,dim(temp)[1])
			temp$DemogGrp<-sapply(Age,DemogFun1)

    #-- Loop through all individuals to update disease statuses.
		DiseaseStatus<-rep(NA,dim(temp)[1])
      NewCount<-rep(NA,dim(temp)[1])
      NewSheddingRate<-rep(NA,dim(temp)[1])
      NewDoseAtInfection<-rep(NA,dim(temp)[1])

    EweGroup<-subset(temp,DemogGrp=="Ewe")
    LambGroup<-subset(temp,DemogGrp=="Lamb")
    
    SpatGrpNew<-rep(NA,dim(temp)[1])
    SpatGrpNew<-ifelse(t %% 365 <= 60, rep(1,dim(temp)[1]),ifelse(t %% 365 == 61, ceiling(runif(dim(temp)[1],min=0,max=ngroups)),temp$SpatGrp))
    temp$SpatGrp<-SpatGrpNew
    
    LambOpen<-lambtransmission.mod.fun(i)
    Lambcontactnumber<-min(dim(LambGroup)[1],LambcontactnumberIn)
    
    for(j in 1:dim(temp)[1]){
        temp.in<-temp[j,]
        #-- Transmission
        if(temp$Status[j]=="S"){
          transmit.out<-transmission.fun(temp.in,temp)
          DiseaseStatus[j]<-transmit.out$DisStat.out
          NewCount[j]<-transmit.out$NewCt.out
          NewDoseAtShedding[j]<-transmit.out$NewDAI.out
          NewSheddingRate[j]<-transmit.out$NewShedRate.out
          
#          if(temp$DemogGrp[j]=="Ewe"){
#              k<-subset(temp, SpatGrp==SpatGrp[j] & DemogGrp=="Ewe")
#              contactsetE<-sum(k$SheddingRate[sample(1:dim(k)[1],contactnumber)])
#              lamb<-subset(temp,DemogGrp=="Lamb" & Mother==temp$ID[j])
#              contactL<-ifelse(dim(lamb)[1]==0,0,ifelse(lamb$Status=="I",1,0))  
              #-- needs to be an indicator for whether her lamb has PN...?
# 
#          if(rbinom(1,1,prob=(1-exp(-(contactsetE))))==1){
            # DiseaseStatus[j]<-ifelse(rbinom(1,1,prob=(1-exp(-(contactsetE))))==1,"I","C")   
#            DiseaseStatus[j]<-"E"   
            #-- need propchronic to reflect dosage
#            NewCount[j]<-temp$Count[j]+1
#            NewDoseAtInfection[j]<-exp(-contactsetE)
#            NewSheddingRate[j]<-ifelse(DiseaseStatus[j]=="I",1,ifelse(DiseaseStatus[j]=="S",0,chronicdose))
#            NewSheddingRate[j]<-chronicdose
#              } else {
#                  DiseaseStatus[j]<-"S"
#                  NewCount[j]<-temp$Count[j]
#                  NewSheddingRate[j]<-0
#                  NewDoseAtInfection[j]<-0
#                }
              
#          } else {    #-- for lambs
#              l<-subset(temp, SpatGrp==SpatGrp[j] & DemogGrp=="Lamb")
#              contactsetL<-ifelse(LambOpen==0,0,ifelse(dim(l)[1]==0,
#                                 0,
#                                 sum(l$SheddingRate[sample(1:dim(l)[1],Lambcontactnumber)])))
#              ewe<-subset(temp,DemogGrp=="Ewe" & ID==temp$Mother[j])
#              contactE<- ifelse(dim(ewe)[1]==0,0,ifelse(ewe$Status=="E"|ewe$Status=="I"|ewe$Status=="C",ewe$SheddingRate,0))  
#              DiseaseStatus[j]<-ifelse(contactE!=0,"I",
#                                       ifelse(rbinom(1,1,prob=(1-exp(-(LambTransmissionProb*contactsetL))))==1,
#                                              "I","S"))
#              NewCount[j]<-ifelse(DiseaseStatus[j]=="I",1,0)
#              NewSheddingRate[j]<-ifelse(DiseaseStatus[j]=="I",1,ifelse(DiseaseStatus[j]=="S",0,chronicdose))
#              NewDoseAtInfection[j]<-1 
#              }
 	 		}
        
        #-- from incubatory to acute or chronic --#
        else if (temp$Status[j]=="E"){
          change<-ifelse(rbinom(1,1,prob=xi)==1,1,0)
          if(change==1){
            k<-rbinom(1,1,prob=rho)
            DiseaseStatus[j]<-ifelse(k==1,"I","C")
            NewCount[j]<-temp$Count[j]
            NewSheddingRate[j]<-ifelse(DiseaseStatus[j]=="I",1,chronicdose)
            NewDoseAtInfection[j]<-temp$DoseAtInfection[j]
          } else {
            DiseaseStatus[j]<-"E"
            NewCount[j]<-temp$Count[j]
            NewSheddingRate[j]<-chronicdose
            NewDoseAtInfection[j]<-temp$DoseAtInfection[j]
          }
        }
        
        #-- Recovery from Acute #-- maybe add death from acute as well. 
        else if(temp$Status[j]=="I"){  #-- can either recover or stay in I. 
						DiseaseStatus[j]<-ifelse(rbinom(1,1,prob=Gamma)==1,"R","I")
            NewCount[j]<-temp$Count[j]
            NewSheddingRate[j]<-ifelse(DiseaseStatus[j]=="I",1,0)
					}
        #-- need to add transition from acute to chronic. 
        
        #-- Recovery/retention in Chronic #-- maybe add death from chronic as well. 
        else if(temp$Status[j]=="C"){  
            currentgammaC<-chronicdecrease*temp$Count[j]*GammaChronic   
						DiseaseStatus[j]<-ifelse(rbinom(1,1,prob=1-currentgammaC)==1,"C","R")
            NewCount[j]<-temp$Count[j]
            NewSheddingRate[j]<-ifelse(DiseaseStatus[j]=="C",chronicdose,0)
					}
        #-- Waning from recovered back to S --#
        else{
          DiseaseStatus[j]<-ifelse(rbinom(1,1,prob=Nu)==1,"S","R")
          NewCount[j]<-temp$Count[j]
          NewSheddingRate[j]<-0
        }  
      }
    
		temp$Status<-DiseaseStatus
    temp$Count<-NewCount
    temp$SheddingRate<-NewSheddingRate
	#-- 3) Update individual ages and "Alive" statuses
		NewAge<-as.numeric(as.character(temp$Age))+1
		temp$Age<-NewAge

		SurvivalStatus<-rep(NA,dim(temp)[1])
		for(j in 1:dim(temp)[1]){
			if(temp$DemogGrp[j]=="Lamb"){
				SurvivalStatus[j]<-ifelse(rbinom(1,1,ifelse(temp$Status[j]=="I",PNLambSurvProb,LambSurvProb))==1,1,0)
			}
			else{
					SurvivalStatus[j]<-ifelse(rbinom(1,1,ifelse(temp$Status[j]=="I",
                                   PNEweSurvProbs[temp$Age[j]],
                                   ifelse(temp$Status[j]=="C",
                                         ChronicEweSurvProbs[temp$Age[j]],
                                         EweSurvProbs[temp$Age[j]])))==1,1,0)
			}
		}
		temp$StillAlive<-SurvivalStatus

		Cause<-rep(NA,dim(temp)[1])
			for(j in 1:dim(temp)[1]){
				if(temp$DemogGrp[j]=="Lamb"){
					Cause[j]<-ifelse(DiseaseStatus[j]=="I" & SurvivalStatus[j]==0,
                           "PN",ifelse(SurvivalStatus[j]==0,
                                       "Other","Alive"))
				}
				else Cause[j]<-ifelse(SurvivalStatus[j]==1, "Alive",
                              ifelse(DiseaseStatus[j]=="I" & rbinom(1,1,prob=InPNAdultSurvAdj)==1,
                                     "PN","Other"))
			}
		temp$Cause<-Cause

	#-- 5) Generate new individuals through birth process.
    BirthWindow<-birth.mod.fun(i)
		NumPotentialMoms<-ifelse(BirthWindow==1,
                             dim(subset(temp,temp$SexRef==0 & temp$StillAlive==1 & temp$HasLamb==0))[1],0)
    PotentialMoms<-subset(temp,temp$SexRef==0 & temp$StillAlive==1 & temp$HasLamb==0)
    
    #-- pick set of potential moms to give birth. --#
    NumNewMoms<-BirthRate*K/(K+dim(temp)[1])*NumPotentialMoms
    momsample<-sample(1:NumPotentialMoms,size=floor(NumNewMoms))
    NewMoms<-PotentialMoms[momsample,]$ID
    MomSpatGrps<-PotentialMoms[momsample,]$SpatGrp
    temp$HasLamb<-apply(as.matrix(temp$ID),1,function(x) mother.fun(x,NewMom=NewMoms))
    NewBirths<-ifelse(BirthWindow==1,length(NewMoms),0)

	#-- 6) Make death mat of who died in this timestep, and how.
  	deathnames<-c("StillAlive","DiseaseStatus","CauseOfDeath","DemogGrp","Age")
    DeathMata<-as.data.frame(cbind(temp$StillAlive, temp$Status, temp$Cause,temp$DemogGrp,temp$Age))
		names(DeathMata)<-deathnames
    DeathMat<-subset(DeathMata,StillAlive==0)

	#-- 7) Characterize new births
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
      NewRows$SpatGrp<-MomSpatGrps  #-- they need to be in the same spatial groups as their moms...

		}
	#-- 8) Store temp in StorageList.
	StorageList[[i]]<-list(TimestepData=data.frame(rbind(subset(temp,StillAlive==1),NewRows)),
                         DeathMat=data.frame(DeathMat))
  levels(StorageList[[i]][[1]]$Status)<-c("S","I","C","R")
    flag<-as.numeric(ifelse(is.na(table(StorageList[[i]][[1]]$Status)["I"])==TRUE,0,1))
    +as.numeric(ifelse(is.na(table(StorageList[[i]][[1]]$Status)["C"])==TRUE,0,1))
	}
	persistence<-table(unlist(lapply(StorageList,is.null)))["FALSE"]

  N<-ChronicCount<-AcuteCount<-SCount<-RCount<-mortality<-rep(NA,persistence)
	for(j in 2:persistence){
    mortality[j]<-dim(StorageList[[j]]$DeathMat)[1]
		N[j]<-dim(StorageList[[j]]$TimestepData)[1]
		ChronicCount[j]<-ifelse(is.na(table(StorageList[[j]]$TimestepData$Status)["C"])==TRUE,
                            0,table(StorageList[[j]]$TimestepData$Status)["C"])
		AcuteCount[j]<-ifelse(is.na(table(StorageList[[j]]$TimestepData$Status)["I"])==TRUE,
                          0,table(StorageList[[j]]$TimestepData$Status)["I"])
  	SCount[j]<-ifelse(is.na(table(StorageList[[j]]$TimestepData$Status)["S"])==TRUE,
                      0,table(StorageList[[j]]$TimestepData$Status)["S"])
    RCount[j]<-ifelse(is.na(table(StorageList[[j]]$TimestepData$Status)["R"])==TRUE,
                      0,table(StorageList[[j]]$TimestepData$Status)["R"])
	}

	return(list(persistence=persistence,
              ChronicCount=ChronicCount,
              AcuteCount=AcuteCount,
              SCount=SCount,
              RCount=RCount,
              N=N,
              mortality=mortality))
}
