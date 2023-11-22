library(glue)

wd = "/Users/mahdi/Dropbox/1-Git/snipar-simulate-save-last-3-gens/src/snipar-simulate-save-last-3-gens/sim"
setwd(wd)


add_k <- function(ro, st, t1){
  
  d = read.table(glue('direct_{st}obs.am_adj_pars.txt'), header = T)
  
  t1[ro, c('k', 'k_se')] <- d[d$parameter=='k', c('estimate', 'SE')]  
  
  return(t1)
  
}


add_delta <- function(ro, st, t1){
  
  d = read.table(glue('direct_{st}obs.2.effects.txt'))
  
  t1[ro, c('delta', 'delta_se')] <- d[d$V1 == 'proband', c('V2', 'V3')]
  
  return(t1)
  
}


add_alpha <- function(ro, st, t1) {
  
  A = as.matrix(c(0, 0, .5, .5))
  A = t(A)
  
  d = read.table(glue('direct_{st}obs.2.effects.txt'))
  
  al =  A %*% d$V2
  t1[ro, 'alpha'] <- al[1, 1]
  
  vc = read.table(glue('direct_{st}obs.2.vcov.txt'))
  vc = as.matrix(vc)
  
  se = A %*% vc %*% t(A)
  
  t1[ro, 'alpha_se'] <- sqrt(se[1, 1])
  
  return(t1)

}

add_diff_alpha <- function(ro, st, t1){
  
  A = as.matrix(c(0, 0, -1, 1))
  A = t(A)
  
  d = read.table(glue('direct_{st}obs.2.effects.txt'))
  
  al =  A %*% d$V2
  t1[ro, 'diff_alpha'] <- al[1, 1]
  
  vc = read.table(glue('direct_{st}obs.2.vcov.txt'))
  vc = as.matrix(vc)
  
  se = A %*% vc %*% t(A)
  
  t1[ro, 'diff_alpha_se'] <- sqrt(se[1, 1])
  
  return(t1)
  
  
}


add_eta_obs <- function(ro, st, t1){
  
  A = as.matrix(c(0, 0, .5, .5, 0, 0, 0, 0))
  A = t(A)
  
  d = read.table(glue('direct_{st}obs.3.effects.txt'))
  
  al =  A %*% d$V2
  t1[ro, 'eta'] <- al[1, 1]
  
  vc = read.table(glue('direct_{st}obs.3.vcov.txt'))
  vc = as.matrix(vc)
  
  se = A %*% vc %*% t(A)
  
  t1[ro, 'eta_se'] <- sqrt(se[1, 1])
  
  return(t1)
  
}


add_alpha_gp <- function(ro, st, t1){
  
  A = as.matrix(c(0, 0, 0, 0, .25, .25, .25, .25))
  A = t(A)
  
  d = read.table(glue('direct_{st}obs.3.effects.txt'))
  
  al =  A %*% d$V2
  t1[ro, 'alpha_gp'] <- al[1, 1]
  
  vc = read.table(glue('direct_{st}obs.3.vcov.txt'))
  vc = as.matrix(vc)
  
  se = A %*% vc %*% t(A)
  
  t1[ro, 'alpha_gp_se'] <- sqrt(se[1, 1])
  
  return(t1)

}


t1 <- data.frame(PGI = c("v0", "v1", "v10"))

mylist <- list('' ,'v1_', 'v10_')

for (i in seq_along(mylist)){
     
  print(paste(i, mylist[i]))
  
  t1 = add_k(i, mylist[i], t1)
  t1 = add_delta(i,mylist[i], t1)
  t1 = add_alpha(i, mylist[i], t1)
  t1 = add_diff_alpha(i, mylist[i], t1)
  t1 = add_eta_obs(i, mylist[i], t1)
  t1 = add_alpha_gp(i, mylist[i], t1)
  
}

write.table(t1, 'table1.txt', quote = FALSE, row.names = FALSE)







      