wd = "/Users/mahdi/Library/CloudStorage/Dropbox/1-Git/snipar-simulate-save-last-3-gens/src/snipar-simulate-save-last-3-gens/sim"
setwd(wd)

d = read.table('pedigree.txt', header = T)

dl = d[sapply(d$FID, function(x) strsplit(x, "_")[[1]][1]) == "23",]

cor(dl[, -c(1:8)])


dl = d[sapply(d$FID, function(x) strsplit(x, "_")[[1]][1]) == "22",]

cor(dl[, -c(1:8)])


dl = d[sapply(d$FID, function(x) strsplit(x, "_")[[1]][1]) == "21",]
cor(dl[, -c(1:8)])


d = read.table('direct_v1_obs.pgs.txt', header = T)



d = read.table('direct_obs.pgs.txt', header = T)
cor(d[,-c(1:4)])



