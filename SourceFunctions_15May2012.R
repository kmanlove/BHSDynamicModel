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
 
      
transmission.fun<-function(temp.in,temp){
        SpatGrp.in<-temp.in$SpatGrp
          if(temp.in$DemogGrp=="Ewe"){
              k<-subset(temp, SpatGrp==SpatGrp.in & DemogGrp=="Ewe")
              contactsetE<-sum(k$SheddingRate[sample(1:dim(k)[1],contactnumber)])
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
                                 sum(l$SheddingRate[sample(1:dim(l)[1],Lambcontactnumber)])))
              ewe<-subset(temp,DemogGrp=="Ewe" & ID==temp.in$Mother)
              contactE<- ifelse(dim(ewe)[1]==0,0,ifelse(ewe$Status=="E"|ewe$Status=="I"|ewe$Status=="C",ewe$SheddingRate,0))  
              DisStat.out<-ifelse(contactE!=0,"I",
                                       ifelse(rbinom(1,1,prob=(1-exp(-(LambTransmissionProb*contactsetL))))==1,
                                              "I","S"))
              NewCt.out<-ifelse(DisStat.out=="I",1,0)
              NewShedRate.out<-ifelse(DisStat.out=="I",1,ifelse(DisStat.out=="S",0,chronicdose))
              NewDAI.out<-1 
              }
    		}
      return(list(DisStat.out=DisStat.out,NewCt.out=NewCt.out,NewShedRate.out=NewShedRate.out,NewDAI.out=NewDAI.out))
}
      
acutechronic.fun<-function(temp.in,temp){
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
        }
    return(list(DisStat.out=DisStat.out,NewCt.out=NewCt.out,NewShedRate.out=NewShedRate.out,NewDAI.out=NewDAI.out))
}