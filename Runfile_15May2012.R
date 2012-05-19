#-- runfile for CumulativeExposureModel_11May2012.R --#

source("work/StatProjects/Raina/sheep/Papers/DynamicModel/Code/IndividualTrackingModel/TwoSeasonModels_11May2012/TwoSeasonsGitRepo/ModelRunSource_15May2012.R")
source("work/StatProjects/Raina/sheep/Papers/DynamicModel/Code/IndividualTrackingModel/TwoSeasonModels_11May2012/TwoSeasonsGitRepo/CumulativeExposureModel_11May2012.R")

ParamMat<-read.table("work/StatProjects/Raina/sheep/Papers/DynamicModel/Code/IndividualTrackingModel/TwoSeasonModels_11May2012/TwoSeasonsGitRepo/ParamMat_17May2012.csv",header=T,sep="")
SimTest<-wrapper.fun(ParamMat)

#-- think about using lapply to read in output files. --#

source("work/StatProjects/Raina/sheep/Papers/DynamicModel/Code/IndividualTrackingModel/TwoSeasonModels_11May2012/TwoSeasonsGitRepo/OutputPlotSource_16May2012.R")

dur.out<-duration.fun(SimTest)
par(mfrow=c(1,1))
sir.plot(SimTest)

par(mfrow=c(2,3))
hist(dur.out$S.dur)
hist(R.dur)
hist(E.dur)
hist(I.dur)
hist(C.dur)



