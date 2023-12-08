
# no indir
wd = ""
setwd(wd)

# with indir
wd = "/Users/mahdi/Dropbox/1-Git/simulation-with-Snipar/sim_output_wt_indir_non_eq"
setwd(wd)

d = read.table('VCs.txt', header = T)

lstg = tail(d, 1)

eq_h = lstg['v_g'] / lstg['v_y']

print(paste('eq H. is: ', eq_h))

h2f = eq_h * (1- lstg['r_delta'])

print(paste('h2f: ', h2f))
