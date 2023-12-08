library(glue)
library(openxlsx)


add_k <- function(ro, st, t1, bc){
  
  d = read.table(glue('{bc}_{st}obs.am_adj_pars.txt'), header = T)
  
  t1[ro, c('k', 'k_se')] <- d[d$parameter=='k', c('estimate', 'SE')]  
  
  return(t1)
  
}


add_delta <- function(ro, st, t1, bc){
  
  d = read.table(glue('{bc}_{st}obs.2.effects.txt'))
  
  t1[ro, c('delta', 'delta_se')] <- d[d$V1 == 'proband', c('V2', 'V3')]
  
  return(t1)
  
}


add_alpha <- function(ro, st, t1, bc) {
  
  A = as.matrix(c(0, 0, .5, .5))
  A = t(A)
  
  d = read.table(glue('{bc}_{st}obs.2.effects.txt'))
  
  al =  A %*% d$V2
  t1[ro, 'alpha'] <- al[1, 1]
  
  vc = read.table(glue('{bc}_{st}obs.2.vcov.txt'))
  vc = as.matrix(vc)
  
  se = A %*% vc %*% t(A)
  
  t1[ro, 'alpha_se'] <- sqrt(se[1, 1])
  
  return(t1)

}

add_diff_alpha <- function(ro, st, t1, bc){
  
  A = as.matrix(c(0, 0, -1, 1))
  A = t(A)
  
  d = read.table(glue('{bc}_{st}obs.2.effects.txt'))
  
  al =  A %*% d$V2
  t1[ro, 'diff_alpha'] <- al[1, 1]
  
  vc = read.table(glue('{bc}_{st}obs.2.vcov.txt'))
  vc = as.matrix(vc)
  
  se = A %*% vc %*% t(A)
  
  t1[ro, 'diff_alpha_se'] <- sqrt(se[1, 1])
  
  return(t1)
  
  
}


add_eta <- function(ro, st, t1, bc){
  
  A = as.matrix(c(0, 0, .5, .5, 0, 0, 0, 0))
  A = t(A)
  
  d = read.table(glue('{bc}_{st}obs.3.effects.txt'))
  
  al =  A %*% d$V2
  t1[ro, 'eta'] <- al[1, 1]
  
  vc = read.table(glue('{bc}_{st}obs.3.vcov.txt'))
  vc = as.matrix(vc)
  
  se = A %*% vc %*% t(A)
  
  t1[ro, 'eta_se'] <- sqrt(se[1, 1])
  
  return(t1)
  
}

add_eta_p <- function(ro, st, t1, bc){
  
  A = as.matrix(c(0, 0, 1, 0, 0))
  A = t(A)
  
  d = read.table(glue('{bc}_{st}obs.3.paternal.effects.txt'))
  
  al =  A %*% d$V2
  t1[ro, 'eta_p'] <- al[1, 1]
  
  vc = read.table(glue('{bc}_{st}obs.3.paternal.vcov.txt'))
  vc = as.matrix(vc)
  
  se = A %*% vc %*% t(A)
  
  t1[ro, 'eta_p_se'] <- sqrt(se[1, 1])
  
  return(t1)
  
}

add_eta_m <- function(ro, st, t1, bc){
  
  A = as.matrix(c(0, 0, 1, 0, 0))
  A = t(A)
  
  d = read.table(glue('{bc}_{st}obs.3.maternal.effects.txt'))
  
  al =  A %*% d$V2
  t1[ro, 'eta_m'] <- al[1, 1]
  
  vc = read.table(glue('{bc}_{st}obs.3.maternal.vcov.txt'))
  vc = as.matrix(vc)
  
  se = A %*% vc %*% t(A)
  
  t1[ro, 'eta_m_se'] <- sqrt(se[1, 1])
  
  return(t1)
  
}


add_eta_imp <- function(ro, st, t1, bc){
  
  A = as.matrix(c(0, 0, .5, .5, 0, 0))
  A = t(A)
  
  d = read.table(glue('{bc}_{st}imp.3.effects.txt'))
  
  al =  A %*% d$V2
  t1[ro, 'eta_imp'] <- al[1, 1]
  
  vc = read.table(glue('{bc}_{st}imp.3.vcov.txt'))
  vc = as.matrix(vc)
  
  se = A %*% vc %*% t(A)
  
  t1[ro, 'eta_imp_se'] <- sqrt(se[1, 1])
  
  return(t1)
  
}


add_eta_p_imp <- function(ro, st, t1, bc){
  
  A = as.matrix(c(0, 0, 1, 0))
  A = t(A)
  
  d = read.table(glue('{bc}_{st}imp.3.paternal.effects.txt'))
  
  al =  A %*% d$V2
  t1[ro, 'eta_p_imp'] <- al[1, 1]
  
  vc = read.table(glue('{bc}_{st}imp.3.paternal.vcov.txt'))
  vc = as.matrix(vc)
  
  se = A %*% vc %*% t(A)
  
  t1[ro, 'eta_p_imp_se'] <- sqrt(se[1, 1])
  
  return(t1)
  
}

add_eta_m_imp <- function(ro, st, t1, bc){
  
  A = as.matrix(c(0, 0, 1, 0))
  A = t(A)
  
  d = read.table(glue('{bc}_{st}imp.3.maternal.effects.txt'))
  
  al =  A %*% d$V2
  t1[ro, 'eta_m_imp'] <- al[1, 1]
  
  vc = read.table(glue('{bc}_{st}imp.3.maternal.vcov.txt'))
  vc = as.matrix(vc)
  
  se = A %*% vc %*% t(A)
  
  t1[ro, 'eta_m_imp_se'] <- sqrt(se[1, 1])
  
  return(t1)
  
}

add_alpha_gp <- function(ro, st, t1, bc){
  
  A = as.matrix(c(0, 0, 0, 0, .25, .25, .25, .25))
  A = t(A)
  
  d = read.table(glue('{bc}_{st}obs.3.effects.txt'))
  
  al =  A %*% d$V2
  t1[ro, 'alpha_gp'] <- al[1, 1]
  
  vc = read.table(glue('{bc}_{st}obs.3.vcov.txt'))
  vc = as.matrix(vc)
  
  se = A %*% vc %*% t(A)
  
  t1[ro, 'alpha_gp_se'] <- sqrt(se[1, 1])
  
  return(t1)

}

add_alpha_gp_imp <- function(ro, st, t1, bc){
  
  A = as.matrix(c(0, 0, 0, 0, .5, .5))
  A = t(A)
  
  d = read.table(glue('{bc}_{st}imp.3.effects.txt'))
  
  al =  A %*% d$V2
  t1[ro, 'alpha_gp_imp'] <- al[1, 1]
  
  vc = read.table(glue('{bc}_{st}imp.3.vcov.txt'))
  vc = as.matrix(vc)
  
  se = A %*% vc %*% t(A)
  
  t1[ro, 'alpha_gp_imp_se'] <- sqrt(se[1, 1])
  
  return(t1)
  
}


beta_col = 'sum'

if (beta_col == 'direct') {
  wd = "/Users/mahdi/Dropbox/1-Git/snipar-simulate-save-last-3-gens/sim"
  setwd(wd)
}

if (beta_col == 'sum') {
  wd = "/Users/mahdi/Dropbox/1-Git/simulation-with-Snipar/sim_output_wt_indir_non_eq"
  setwd(wd)
}


t1 <- data.frame(PGI = c("v0", "v1", "v10"))

mylist <- list('' ,'v1_') #, 'v10_')

bc = beta_col

for (i in seq_along(mylist)){
     
  print(paste(i, mylist[i]))
  
  t1 = add_k(i, mylist[i], t1, bc)
  t1 = add_delta(i,mylist[i], t1, bc)
  t1 = add_alpha(i, mylist[i], t1, bc)
  t1 = add_diff_alpha(i, mylist[i], t1, bc)
  t1 = add_eta(i, mylist[i], t1, bc)
  t1 = add_eta_p(i, mylist[i], t1, bc)
  t1 = add_eta_m(i, mylist[i], t1, bc)
  t1 = add_eta_imp(i, mylist[i], t1, bc)
  t1 = add_eta_p_imp(i, mylist[i], t1, bc)
  t1 = add_eta_m_imp(i, mylist[i], t1, bc)
  t1 = add_alpha_gp(i, mylist[i], t1, bc)
  t1 = add_alpha_gp_imp(i, mylist[i], t1, bc)
  
}

write.xlsx(t1, 'table1.xlsx')



      