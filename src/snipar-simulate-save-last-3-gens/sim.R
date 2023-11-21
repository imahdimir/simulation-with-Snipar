wd = "/Users/mahdi/Dropbox/1-Git/snipar-simulate-save-last-3-gens/src/snipar-simulate-save-last-3-gens/sim"
setwd(wd)

d = read.table("direct_obs.pgs.txt", header = T)

d = d[sapply(d$FID, function(x) strsplit(x, "_")[[1]][1]) == "23",]

write.table(d, 'direct_obs_offspring.pgs.txt', quote = FALSE, row.names = FALSE)

cor(d[, -c(1:4)])


d = read.table("direct_v1_obs.pgs.txt", header = T)
d = d[sapply(d$FID, function(x) strsplit(x, "_")[[1]][1]) == "23",]

write.table(d, 'direct_v1_obs_offspring.pgs.txt', quote = FALSE, row.names = FALSE)

cor(d[, -c(1:4)])



d = read.table("direct_v10_obs.pgs.txt", header = T)

d = d[sapply(d$FID, function(x) strsplit(x, "_")[[1]][1]) == "23",]

write.table(d, 'direct_v10_obs_offspring.pgs.txt', quote = FALSE, row.names = FALSE)

cor(d[, -c(1:4)])



