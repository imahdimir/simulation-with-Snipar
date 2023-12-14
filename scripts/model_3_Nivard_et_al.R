wd = "/Users/mahdi/Dropbox/1-Git/simulation-with-Snipar/sim_no_indir"
setwd(wd)


d=read.table('direct_v1_obs.pgs.txt',header=T)
typeof(d)

dpar = d[sapply(d[,1],function(x) strsplit(x,'_')[[1]][1]=='22'),]

sib_mean <- aggregate(dpar$proband, by=list(dpar$FID), FUN=mean)

dpar$sib_mean = sim_mean[as.character(dpar$Group.1)]

pardiff = data.frame(IID=dpar[seq(1,dim(dpar)[1],2),'IID'],
                     sib_mean = (dpar[seq(1,dim(dpar)[1],2),'proband']+dpar[seq(2,dim(dpar)[1],2),'proband'])/2)

sib_mean = (dpar[seq(1,dim(dpar)[1],2),'proband'] + dpar[seq(2,dim(dpar)[1],2),'proband']) / 2
sim_mean = data.frame(FID=dpar[seq(1,dim(dpar)[1],2),'FID'],
                      sib_mean = (dpar[seq(1,dim(dpar)[1],2),'proband'] + dpar[seq(2,dim(dpar)[1],2),'proband']) / 2)
seq(1,dim(dpar)[1],2)

pardiff = rbind(pardiff,
                data.frame(IID=d[seq(2,dim(dpar)[1],2),'IID'],
                           sibdiff = d[seq(2,dim(dpar)[1],2),'proband']-d[seq(1,dim(dpar)[1],2),'proband'],
                           sibsum = d[seq(2,dim(dpar)[1],2),'proband']+d[seq(1,dim(dpar)[1],2),'proband']))

d=d[sapply(d[,1],function(x) strsplit(x,'_')[[1]][1]=='23'),]
d$patdiff = pardiff[match(d$FATHER_ID,pardiff$IID),'sibdiff']
d$patsum = pardiff[match(d$FATHER_ID,pardiff$IID),'sibsum']

y=read.table('phenotype.txt',header=F)

d$y=y[match(d[,2],y[,2]),3]

require(lme4)
rsib_bp = lmer(y~proband+patdiff+patsum+maternal+(1|FID),data=d)