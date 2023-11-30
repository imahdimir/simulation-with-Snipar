wd = "/Users/mahdi/Dropbox/1-Git/snipar-simulate-save-last-3-gens/sim"
setwd(wd)


d = read.table('VCs.txt', header = T)

lstg = tail(d, 1)

eq_h = lstg['v_g'] / lstg['v_y']

h2f = eq_h * (1- lstg['r_delta'])

print(h2f)
