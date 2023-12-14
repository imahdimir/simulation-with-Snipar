oo <- options(repos = "https://cran.r-project.org/")
install.packages("Matrix")
install.packages("lme4")
options(oo)



require(lme4)
require(glue)

wd_pref = "/Users/mahdi/Dropbox/1-Git/simulation-with-Snipar"


fit_model3 <- function(subdir, sum_or_direct, noise_level) {

  setwd(glue('{wd_pref}/{subdir}'))
  
  d=read.table(glue('{sum_or_direct}_{noise_level}obs.pgs.txt'), header=T)
  
  last_row = max(d[, 1])
  last_gen = strsplit(last_row, '_')[[1]][1]
  last_gen = as.integer(last_gen)
  
  dpar <- d[sapply(d[,1],function(x) strsplit(x,'_')[[1]][1]==toString(last_gen - 1)),]
  
  sib_mean <- aggregate(dpar$proband, by=list(dpar$FID), FUN=mean)
  sib_mean <- setNames(sib_mean, c('FID', 'sib_mean'))
  
  sib_mean <- merge(dpar, sib_mean)
  
  sib_mean <- sib_mean[, c('IID', 'sib_mean')]
  
  d <- d[sapply(d[,1],function(x) strsplit(x,'_')[[1]][1]==toString(last_gen)),]
  d <- d[, c('FID','IID', 'FATHER_ID', 'MOTHER_ID', 'proband', 'paternal', 'maternal')]
  
  d$pat_sib_mean <- sib_mean[match(d$FATHER_ID,sib_mean$IID),'sib_mean']
  d$mat_sib_mean <- sib_mean[match(d$MOTHER_ID,sib_mean$IID),'sib_mean']
  
  d$patdiff <- d$paternal - d$pat_sib_mean
  d$matdiff <- d$maternal - d$mat_sib_mean
  
  y <- read.table('phenotype.txt',header=F)
  
  d$y <- y[match(d$IID,y[,2]) , 3]
  
  pat_model = lmer(y ~ proband + patdiff + pat_sib_mean + maternal + (1 | FID), data=d)
  mat_model = lmer(y ~ proband + matdiff + mat_sib_mean + paternal + (1 | FID), data=d)
  
  write.csv(summary(pat_model)$coefficients , glue('model3_{noise_level}paternal.csv'))
  write.csv(summary(mat_model)$coefficients , glue('model3_{noise_level}maternal.csv'))
  
  return(list(pat_model, mat_model))

}


dirs = c('sim_no_indir', 'sim_no_indir_non_eq', 'sim_population', 'sim_population_non_eq')
noise_levels = c('', 'v1_', 'v10_')

cart_prod = expand.grid(dirs, noise_levels)

# drop rows of non_eq with no noise and v10 noise
cart_prod = cart_prod[-c(2,4,10,12), ]

cart_prod$s_or_d = ifelse(grepl('population', cart_prod$Var1), 'sum', 'direct')

ar = cart_prod


for (i in 1:nrow(ar)) {
  print(paste('Fitting model 3 for', ar[i,1], ar[i,2], ar[i,3]))
  x = fit_model3(ar[i,1], ar[i,3], ar[i,2])
  
  x1 = as.data.frame(summary(x[[1]])$coefficients)
  ar[i, 'patdiff'] = x1['patdiff', ]$Estimate
  ar[i, 'patdiff_se'] = x1['patdiff', ]['Std. Error']
  
  x2 = as.data.frame(summary(x[[2]])$coefficients)
  ar[i, 'matdiff'] = x2['matdiff', ]$Estimate
  ar[i, 'matdiff_se'] = x2['matdiff', ]['Std. Error']
}

library(openxlsx)
write.xlsx(ar, 'model3_results.xlsx')

x1 = as.data.frame(summary(x)$coefficients)
x1['patdiff', ]$Estimate
x1['patdiff', ]['Std. Error']


