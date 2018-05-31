library(changepoint)
library(imputeTS)

###Change the file location
testMWD<-read.csv2("C:/Users/veenas/Desktop/Latex-hmm/Codelibrary/MWD-1939-21.csv",sep=',')


testMWD<-apply(testMWD,2,as.numeric)
names(testMWD)<-c("Depth","PenRate","PercPress","FeedPress","FlushAirPress","RotPress","DampPress")
testMWD<-data.frame(testMWD)

################
###Deleting Rod change outliers, rod change will always happend after 3m depth
### mean-3Sd of feeder pressure removed
#######
MWD1<-subset(testMWD,Depth>3)
cutoff<-mean(MWD1$FeedPress)-3*sd(MWD1$FeedPress)
cleanmwd1<-subset(MWD1,FeedPress>cutoff)
MWD2<-subset(testMWD,Depth<3)
MWD3<-rbind(MWD2,cleanmwd1)

#########################
### Chancgepoint for pressure change, Changepoint package, 
####Change pont on percussion presure
#####
c1<-cpts(cpt.mean(MWD3$PercPress))
presschange<-MWD3$Depth[c1]
HardRock.MWD<-subset(MWD3,Depth>presschange)
MWD1<-HardRock.MWD

##############################
###Interpolation to make the data in even interval
##used package imputeTS in R.
###
a<-MWD1$Depth[1]
b<-MWD1$Depth[length(MWD1$Depth)]
Depth<-round(seq(from=a,to=b,by=0.01),2)
n<-length(Depth)
Missing.CleanMwd.1<-matrix(NA,n,7)
Missing.CleanMwd.1[,1]<-Depth

i=1
while (i <= n)
{
  for(d in 1:length(MWD1$Depth))
     {  if(Depth[i]==round(MWD1$Depth[d],2))
        {Missing.CleanMwd.1[i,]=as.matrix(round(MWD1[d,],2))}
     }
i=i+1
}
Missing.Mwd.1<-as.data.frame(Missing.CleanMwd.1)
names(Missing.Mwd.1)<-c("Depth","PenRate","PercPress","FeedPress","FlushAirPress","RotPress","DampPress")
Hmm.data<-apply(Missing.Mwd.1,2,na.interpolation)

###save data for HMM code, change the location 
write.csv(Hmm.data,"C:/Users/veenas/Desktop/Latex-hmm/Codelibrary/HMMdata.csv")
