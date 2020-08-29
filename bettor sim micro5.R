library(splines)
library(nnet)


#number of races in a simulation
rs=20000
#number of bettors on each race
bs=500
#number of races back that bettors remember
w=1
#standard deviation of fundamental value. This controls the signal to noise ratio in race outcomes. Higher values will make the odds better predictors of the outcomes.
sdStart=0
#number of horses in first race
hStart=10
#potential return rank multiplier. Turning this up makes bettors more strongly prefer high ranking bets.
rdm=1
#realized return extrapolation multiplier. Turning this up makes bettors more strongly prefer to bet on bets whose odds are close to the odds of bets that won recently.
em=0.005


#number of horses in a race
#hs=rep(hStart,rs)

hs<-rep(NA,rs)
hs[1]=hStart
for (i in 2:rs){
	hs[i]=ifelse(rbinom(1,1,.5)==0,hs[i-1],hs[i-1]+(rbinom(1,1,.5)*2-1))
	hs[i]<-ifelse(hs[i]<2,2,hs[i])
	hs[i]<-ifelse(hs[i]>18,18,hs[i])

}

hs<-rep(10,rs)

#this keeps track of races
race<-rep(1:rs,hs)

#set up the fundamentals for all bets in the world by drawing from a normal with mean 10 and sd set above before being truncated at 0.001.
rawProbs<-rnorm(sum(hs),0,1)^2
#normalize fundamentals in a race to sum to 1
denom<-rep(aggregate(rawProbs,list(race),sum)[,2],hs)
probs<-rawProbs/denom

#pick winners of races according to fundamental probabilities
winners<-aggregate(probs,list(race),function(x) min(which(cumsum(x)>runif(1))))
winner<-rep(0,sum(hs))
winner[c(0,cumsum(hs)[-rs])+winners[,2]]<-1

#set up holder for betting odds
odds<-matrix(NA,nrow=sum(hs),ncol=length(bs))
real<-matrix(NA,nrow=sum(hs),ncol=length(bs))

#loop over worlds
for (i in 1:length(bs)){

	#loop over races within world
	for (r in 1:rs){
		
		#select the current race.
		selT<-which(race==r)
		#check how many horses in this race
		h<-length(selT)
		#grab fundamental values.
		prob<-probs[selT]
		#bookmaker starts each race by putting 1 dollar on each horse. Effectively a negative vig. We can work in positive vig later.
		money<-rep(.00001,h)
		#pull out previous reference odds and previous realized returns
		#selP<-which(race<=(r-1)&race>=(r-runif(1,1,w)))
		#prevRef<-odds[selP]
		#prevReal<-prevRef*winner[selP]
				
		#loop over bettors in the market
		for (b in 1:bs[i]){

			#calculate current odds on current race
			odd<-sum(money)/money
			#build reference set that includes current race
			#ref<-c(prevRef,odd)
			#these are holders for potential return ranks and average realized returns from horses with similar odds for each horse in the current race.
			#ranks<-rep(NA,h)	
			#extraps<-rep(NA,h)	

			#loop over horses in race and calculate the potential return rank and realized return average for each horse in the race.
			#for (j in 1:h){
			#	ranks[j]<-mean(odd[j]>ref)
			#	#extraps[j]<-mean(prevReal[abs(log(odd[j]/prevRef))<=.3])		
			#}
			
			#fill in the average realized return for any bet whose odds are novel in the reference set
			#extraps[is.na(extraps)]<-mean(extraps,na.rm=T)
			#calculate the expected utility of betting on each horse
			#util<-prob*odd+rdm*ranks+em*extraps
			#util<-prob*odd+em*extraps
			#util<-prob*odd+rdm*ranks
			belief<-pnorm(rnorm(n=h,mean=qnorm(prob),sd=1))
			#util<-belief*odd+rdm*ranks
			util<-belief*odd
			#util<-prob*(odd^1.2)

			#put 1 dollar on the expected utility maximizing horse
			best<-which(util>1)
			money[best]<-money[best]+1
		}

		odds[selT,i]<-sum(money)/money
		print(r)

	}

	real[,i]<-winner*odds[,i]
}

#window width for analysis
wPrime=1
#number of windows
wPrimes=10

rankP<-lapply(1:length(bs),function(x) matrix(NA,nrow=sum(hs),ncol=wPrimes))
rankF<-lapply(1:length(bs),function(x) matrix(NA,nrow=sum(hs),ncol=wPrimes))
nRankP<-lapply(1:length(bs),function(x) matrix(NA,nrow=sum(hs),ncol=wPrimes))
nRankF<-lapply(1:length(bs),function(x) matrix(NA,nrow=sum(hs),ncol=wPrimes))
horses<-rep(hs, hs)

for (w in 1:length(bs)){
	odd<-odds[,w]
	rP<-rankP[[w]]
	rF<-rankF[[w]]
	nRP<-nRankP[[w]]
	nRF<-nRankF[[w]]

	for (i in wPrimes:(rs-wPrimes+1)){

		sel<-which(race==i)
		targ<-odd[sel]
		selR<-which(race>=(i-wPrimes)&race<=(i+wPrimes))
		refs<-odd[selR]
		races<-race[selR]
		nTarg<-length(targ)
		rawRankP1<-matrix(NA,nrow=nTarg,ncol=wPrimes)
		rawRankF1<-matrix(NA,nrow=nTarg,ncol=wPrimes)
		rawNRankP1<-matrix(NA,nrow=nTarg,ncol=wPrimes)
		rawNRankF1<-matrix(NA,nrow=nTarg,ncol=wPrimes)
		
		for (off in 1:wPrimes){

			refP<-refs[races<=(i-off*wPrime+wPrime)&races>(i-off*wPrime)]
			refF<-refs[races>=(i+off*wPrime-wPrime)&races<(i+off*wPrime)]
			nRefP<-length(refP)
			nRefF<-length(refF)			

			for (j in 1:nTarg){
				rawRankP1[j,off]<-(sum(targ[j]>refP)+.5*sum(targ[j]==refP))/nRefP
				rawRankF1[j,off]<-(sum(targ[j]>refF)+.5*sum(targ[j]==refF))/nRefF
			}

		}

		rP[sel,]<-rawRankP1
		rF[sel,]<-rawRankF1
		nRP[sel,]<-nRefP
		nRF[sel,]<-nRefF
		print(i)

	}
	rankP[[w]]<-rP
	rankF[[w]]<-rF

	nRankP[[w]]<-nRP
	nRankF[[w]]<-nRF

}

aveRank<-rowMeans(rankP[[1]])	
highest<-as.numeric(nRankP[[1]][,1]-rankP[[1]][,1]*nRankP[[1]][,1]<=.5)
summary(lm(real[,1]~I(1/odds[,1])))

summary(nRankP[[1]][,1])
summary(nRankF[[1]][,1])
summary(nRP[,1])


summary(lm(real[,1]~bs(1/odds[,1],100)+rankP[[1]]+rankF[[1]]))
summary(lm(real[,1]~bs(1/odds[,1],100)+aveRank))
summary(lm(real[,1]~bs(1/odds[,1],100)+rankP[[1]][,1]+highest))
summary(lm(real[,1]~bs(1/odds[,1],100)+rankP[[1]][,1]+aveRank))

summary(lm(real[,1]~bs(1/odds[,1],100)+factor(horses)+rankP[[1]]))
summary(lm(real[,1]~bs(1/odds[,1],20)*bs(horses,10)+aveRank))
summary(lm(real[,1]~bs(1/odds[,1],20)*bs(horses,10)+rankP[[1]]+rankF[[1]]))



length(rankP)

x<-lm(real[,1]~bs(rankP[[1]][,1],100))
x2<-lm(real[,1]~bs(rankP[[1]][,1],100)+bs(odds[,1],100))
anova(x,x2)

x<-lm(real[,1]~bs(aveRank,100))
x2<-lm(real[,1]~bs(aveRank,100)+bs(odds[,1],100))
anova(x,x2)

bins<-ceiling(rank(1/odds[,1])/length(winner)*100)
hold<-aggregate(cbind(real[,1],1/odds[,1]),list(bins),mean)
plot(y=hold[,2],x=hold[,3],col="red",xlim=c(0,1),bty="n",pch=19)


bins<-ceiling(rank(1/odds[,2])/length(winner)*100)
hold<-aggregate(cbind(real[,2],1/odds[,2]),list(bins),mean)
points(y=hold[,2],x=hold[,3],col="black",pch=19)
