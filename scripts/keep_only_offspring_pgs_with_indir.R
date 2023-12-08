library(glue)

wd = "/Users/mahdi/Dropbox/1-Git/simulation-with-Snipar/sim_output_wt_indir"
setwd(wd)


keep_offspring_only <- function(st1, st2){
  fn = glue("sum_{st1}{st2}.pgs.txt")
  print(fn)
  
  d = read.table(fn , header = T)
  d = d[sapply(d$FID, function(x) strsplit(x, "_")[[1]][1]) == "23",]
  
  fn = glue('sum_{st1}{st2}_offspring_only.pgs.txt')
  print(fn)
  
  write.table(d, fn, quote = FALSE, row.names = FALSE)
}


noise <- list('' ,'v1_', 'v10_')
obs_imp <- list('obs', 'imp')

cart_prod = expand.grid(noise, obs_imp)


for (i in 1 : nrow(cart_prod)){
  
  keep_offspring_only(cart_prod[i, 1], cart_prod[i, 2])

}

