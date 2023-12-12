library(glue)


wd = "/Users/mahdi/Dropbox/1-Git/simulation-with-Snipar/sim_out_no_indir_non_eq"
setwd(wd)


keep_offspring_only <- function(prefx, st1, st2, last_gen_n){
  fn = glue("{prefx}_{st1}{st2}.pgs.txt")
  print(fn)
  
  d = read.table(fn , header = T)
  d = d[sapply(d$FID, function(x) strsplit(x, "_")[[1]][1]) == toString(last_gen_n),]
  
  fn = glue('{prefx}_{st1}{st2}_offspring_only.pgs.txt')
  print(fn)
  
  write.table(d, fn, quote = FALSE, row.names = FALSE)
}


isdir = TRUE
if (isdir){
  prefx = 'direct'
  
} else{
  prefx = 'sum'
}

noise <- list('v1_')   # '' ,'v1_', 'v10_')
obs_imp <- list('obs', 'imp')
cart_prod = expand.grid(noise, obs_imp)

last_n = 5


for (i in 1 : nrow(cart_prod)){
  
  keep_offspring_only(prefx, cart_prod[i, 1], cart_prod[i, 2], last_n)

}

