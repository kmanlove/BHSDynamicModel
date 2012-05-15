require(stats)	

IndividualSIR<-function(timesteps=timesteps,
                        BirthRate=BirthRate,
                        SexRatio=SexRatio,
                        n=n,
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
 	 		}
        
        #-- from incubatory to acute or chronic --#
        else if (temp$Status[j]=="E"){
          acutechronic.out<-acutechronic.fun(temp.in,temp)
          DiseaseStatus[j]<-acutechronic.fun$DisStat.out
          NewCount[j]<-acutechronic.fun$NewCt.out
          NewSheddingRate[j]<-acutechronic.fun$NewShedRate.out
          NewDoseAtInfection[j]<-acutechronic.fun$NewDAI.out
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
    temp$StillAlive<-survival.fun(temp)
    temp$Cause<-cause.fun(temp)

	#-- 4) Generate new individuals through birth process.
    NewBirths<-birth.fun(i,temp)

	#-- 5) Make death mat of who died in this timestep, and how.
  	deathnames<-c("StillAlive","DiseaseStatus","CauseOfDeath","DemogGrp","Age")
    DeathMata<-as.data.frame(cbind(temp$StillAlive, temp$Status, temp$Cause,temp$DemogGrp,temp$Age))
		names(DeathMata)<-deathnames
    DeathMat<-subset(DeathMata,StillAlive==0)
    
	#-- 6) Store temp in StorageList.
  StorageList[[i]]<-list(TimestepData=data.frame(rbind(subset(temp,StillAlive==1),NewBirths)),
                         DeathMat=data.frame(DeathMat))
  levels(StorageList[[i]][[1]]$Status)<-c("S","E","I","C","R")
    flag<-as.numeric(ifelse(is.na(table(StorageList[[i]][[1]]$Status)["I"])==TRUE,0,1))
    +as.numeric(ifelse(is.na(table(StorageList[[i]][[1]]$Status)["C"])==TRUE,0,1))
	}
  
  #-- Analyze for-loop output --#
  output.analysis<-analysis.fun(StorageList)
  persistence<-output.analysis$persistence; N<-output.analysis$N
  ChronicCount<-output.analysis$ChronicCount;AcuteCount<-output.analysis$AcuteCount;SCount<-output.analysis$SCount
  RCount<-output.analysis$RCount; mortality<-output.analysis$mortality

	return(list(persistence=persistence,
              ChronicCount=ChronicCount,
              AcuteCount=AcuteCount,
              SCount=SCount,
              RCount=RCount,
              N=N,
              mortality=mortality))
}
