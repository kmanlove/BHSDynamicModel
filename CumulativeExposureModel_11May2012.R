require(stats)	

#-- should switch this to with. --#.
IndividualSIR<-function(timesteps=timesteps,
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
  RefTS1$Mother<-RefTS1$DemogGrp<-rep(NA,length(RefTS1$Age))
    for(i in 1:n){
      RefTS1$Mother[i]<-ifelse(RefTS1$Age[i]>=365,NA,sample(subset(RefTS1,Age>=365)$ID,1))
      RefTS1$DemogGrp[i]<-ifelse(RefTS1$Age[i]>=365,"Ewe","Lamb")
    }
  RefTS1$SpatGrp<-rep(1,length(RefTS1$Age))
  RefTS1$SpatGrpSeason<-floor(runif(length(RefTS1$Age),min=0,max=ngroups))
  RefTS1$DoseAtInfection<-rep(NA,length(RefTS1$Age))
  RefTS1$Groupsize<-rep(n,length(RefTS1$Age))
	StorageList[[1]][[1]]<-RefTS1

  time1<-rep(NA,timesteps)
  time2<-rep(NA,timesteps)
  time3<-rep(NA,timesteps)
  
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

    #-- Loop through all individuals to update disease statuses.
		DiseaseStatus<-rep(NA,dim(temp)[1])
    NewCount<-rep(NA,dim(temp)[1])
    NewSheddingRate<-rep(NA,dim(temp)[1])
    NewDoseAtInfection<-rep(NA,dim(temp)[1])
    groupsizeNew<-rep(NA,dim(temp)[1])

    EweGroup<-subset(temp,DemogGrp=="Ewe")
    LambGroup<-subset(temp,DemogGrp=="Lamb")
    
    LambOpen<-lambtransmission.mod.fun(i)
    Lambcontactnumber<-min(dim(LambGroup)[1],LambcontactnumberIn)
    season<-ifelse(i %% 365 <=60, 0,1)
    
    #-- store contactsets --#
    contactsets<-vector("list",dim(temp)[1])
    
    ptm<-proc.time()
    for(j in 1:dim(temp)[1]){
        temp.in<-temp[j,]
        #-- Transmission
        if(temp$Status[j]=="S"){
          transmit.out<-transmission.fun(season,temp.in,temp,contactnumber=contactnumber,Lambcontactnumber=Lambcontactnumber,LambOpen=LambOpen,LambTransmissionProb=LambTransmissionProb,chronicdose=chronicdose)
          DiseaseStatus[j]<-transmit.out$DisStat.out
          NewCount[j]<-transmit.out$NewCt.out
          NewDoseAtInfection[j]<-transmit.out$NewDAI.out
          NewSheddingRate[j]<-transmit.out$NewShedRate.out
          contactsets[[j]]<-transmit.out$contactset
 	 		}
        
        #-- from incubatory to acute or chronic --#
        else if (temp$Status[j]=="E"){
          if(temp$DemogGrp[j]=="Lamb"){
            DiseaseStatus[j]<-ifelse(rbinom(1,1,xi)==1,"I","E")
            NewCount[j]<-temp$Count[j]
            NewSheddingRate[j]<-1
            NewDoseAtInfection[j]<-temp$DoseAtInfection[j]
          } else{
          acutechronic.out<-acutechronic.fun(temp.in,temp,chronicdose=chronicdose,xi=xi,rho=rho)
          DiseaseStatus[j]<-acutechronic.out$DisStat.out
          NewCount[j]<-acutechronic.out$NewCt.out
          NewSheddingRate[j]<-acutechronic.out$NewShedRate.out
          NewDoseAtInfection[j]<-acutechronic.out$NewDAI.out
          }
        }
        
        #-- Recovery from Acute #-- maybe add death from acute as well. 
        else if(temp$Status[j]=="I"){  #-- can either recover or stay in I. 
  					DiseaseStatus[j]<-ifelse(rbinom(1,1,prob=Alpha/(1-Alpha)*(eta))==1,"R","I")
            NewCount[j]<-temp$Count[j]
            NewSheddingRate[j]<-ifelse(DiseaseStatus[j]=="I",1,0)
					}
        #-- need to add transition from acute to chronic. 
        
        #-- Recovery/retention in Chronic #-- maybe add death from chronic as well. 
        else if(temp$Status[j]=="C"){  
            currentgammaC<-chronicdecrease*temp$Count[j]*Gamma   
						DiseaseStatus[j]<-ifelse(rbinom(1,1,prob=tau)==1,"A",ifelse(rbinom(1,1,prob=Gamma)==1,"R","C"))
            NewCount[j]<-temp$Count[j]
            NewSheddingRate[j]<-ifelse(DiseaseStatus[j]=="C",chronicdose,ifelse(DiseaseStatus[j]=="A",1,0))
					}
        #-- Waning from recovered back to S --#
        else{
          DiseaseStatus[j]<-ifelse(rbinom(1,1,prob=Nu)==1,"S","R")
          NewCount[j]<-temp$Count[j]
          NewSheddingRate[j]<-0
        }  
        #-- record groupsize --#
        if(season==0){
          groupsizeNew[j]<-dim(temp)[1]
        } else groupsizeNew[j]<-dim(subset(temp,SpatGrpSeason==temp$SpatGrpSeason[j]))[1]
      }
    
     print(proc.time()-ptm)
    
		temp$Status<-DiseaseStatus
    temp$Count<-NewCount
    temp$SheddingRate<-NewSheddingRate
    temp$Groupsize<-groupsizeNew
    
#    table(temp$DemogGrp,temp$Status)
    
	#-- 3) Update individual ages and "Alive" statuses
    ptm2<-proc.time()
		NewAge<-as.numeric(as.character(temp$Age))+1
		temp$Age<-NewAge
    temp$StillAlive<-survival.fun(temp,PNEweSurvProbs,ChronicEweSurvProbs,EweSurvProbs,Recr,GammaLamb,SLSdecrease)
    temp$Cause<-cause.fun(temp,Alpha)
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
  StorageList[[i]]<-list(TimestepData=TimestepData,DeathMat=data.frame(DeathMat),contactsets=contactsets)
  levels(StorageList[[i]][[1]]$Status)<-c("S","E","I","C","R")
    flag<-as.numeric(ifelse(is.na(table(StorageList[[i]][[1]]$Status)["I"])==TRUE,0,1))+as.numeric(ifelse(is.na(table(StorageList[[i]][[1]]$Status)["C"])==TRUE,0,1))+as.numeric(ifelse(is.na(table(StorageList[[i]][[1]]$Status)["E"])==TRUE,0,1))

    print(i)
	}
  
  #-- Analyze for-loop output --#
  ptm4<-proc.time()
  output.analysis<-analysis.fun(StorageList)
  persistence<-output.analysis$persistence; N<-output.analysis$N
  ChronicCount<-output.analysis$ChronicCount;AcuteCount<-output.analysis$AcuteCount;SCount<-output.analysis$SCount
  RCount<-output.analysis$RCount; mortality<-output.analysis$mortality; ECount<-output.analysis$ECount
  contactsets<-output.analysis$contactsets; groupsize<-output.analysis$groupsize
  print(proc.time()-ptm4)
  
	return(list(persistence=persistence,  #-- scalar
              ECount=ECount,#-- vector (variable length)
              ChronicCount=ChronicCount,#-- vector (variable length)
              AcuteCount=AcuteCount,
              SCount=SCount,
              RCount=RCount,
              N=N,
              mortality=mortality,
              contactsets=contactsets,
              groupsize=groupsize,
              time1=time1,
              time2=time2,
              time3=time3,
              StorageList=StorageList #-- 
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
