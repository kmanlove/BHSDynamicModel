#-- source functions for analysis of simulation output --#
#-- K. Manlove, 16 May 2012 --#

#-- duration.fun extracts the duration of every individual in every class they experienced over the course of the simulation --#

duration.fun<-function(loop.out){
  persist<-loop.out$persistence

small.StLi<-vector("list",persist-1)
for(i in 1:(persist-1)){
  small.StLi[[i]]<-loop.out$StorageList[[i+1]][[1]]
}

out<-do.call("rbind",small.StLi)

duration.list<-t<-vector("list",length(levels(as.factor(out$ID))))

for(i in 1:length(duration.list)){
  k<-subset(out,ID==i)
  duration.list[[i]]<-k
  duration.list[[i]]$switches<-duration.list[[i]]$timeinclass<-duration.list[[i]]$lasttime<-duration.list[[i]]$simend<-duration.list[[i]]$oldstatus<-rep(NA,dim(duration.list[[i]])[1])
  duration.list[[i]]$switches[1]<-duration.list[[i]]$timeinclass[1]<-0
  duration.list[[i]]$oldstatus[1]<-duration.list[[i]]$lasttime[1]<-duration.list[[i]]$simend[1]<-NA
  
  for(j in 2:dim(duration.list[[i]])[1]){
    duration.list[[i]]$switches[j]<-ifelse(j==dim(duration.list[[i]])[1],1,ifelse(duration.list[[i]]$Status[j]!=duration.list[[i]]$Status[j-1],1,0))
    duration.list[[i]]$timeinclass[j]<-ifelse(duration.list[[i]]$switches[j]==0,duration.list[[i]]$timeinclass[j-1]+1,0)
    duration.list[[i]]$oldstatus[j]<-ifelse(duration.list[[i]]$switches[j]==0,NA,as.character(duration.list[[i]]$Status[j-1]))
    duration.list[[i]]$lasttime[j]<-ifelse(duration.list[[i]]$switches[j]==0,NA,as.character(duration.list[[i]]$timeinclass[j-1]))
    duration.list[[i]]$simend[j]<-ifelse(j==dim(duration.list[[i]])[1],1,0)
  }
  
  t[[i]]<-subset(duration.list[[i]],switches==1)
}

switch.frame<-do.call("rbind",t)

R.dur<-as.numeric(as.character(subset(switch.frame,oldstatus=="R")$lasttime))
E.dur<-as.numeric(as.character(subset(switch.frame,oldstatus=="E")$lasttime))
I.dur<-as.numeric(as.character(subset(switch.frame,oldstatus=="I")$lasttime))
C.dur<-as.numeric(as.character(subset(switch.frame,oldstatus=="C")$lasttime))
S.dur<-as.numeric(as.character(subset(switch.frame,oldstatus=="S")$lasttime))

return(list(R.dur=R.dur,E.dur=E.dur,I.dur=I.dur,C.dur=C.dur,S.dur=S.dur))
}


#-- sir.plot builds the standard SIR plot: lines for each state with --#
  #-- y = count in state and x = time --#
sir.plot<-function(SimOut){
plot(SimOut$N~seq(1:length(SimOut$N)),type="l",ylim=c(0,max(na.omit(SimOut$N))))
lines(SimOut$AcuteCount~seq(1:length(SimOut$N)),type="l",col="red")
lines(SimOut$ChronicCount~seq(1:length(SimOut$N)),type="l",col="red",lty=2)
lines(SimOut$ECount~seq(1:length(SimOut$N)),type="l",col="purple",lty=1)
lines(SimOut$SCount~seq(1:length(SimOut$N)),type="l",col="blue",lty=1)
lines(SimOut$RCount~seq(1:length(SimOut$N)),type="l",col="green",lty=1)
}
