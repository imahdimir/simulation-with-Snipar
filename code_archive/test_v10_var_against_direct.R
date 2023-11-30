
# test for v10 variances being 11x of direct
d = read.table('causal_effects.txt', header = T)
v0 = var(d$direct)
print(v0)
v10 = var(d$direct_v10)
print(v10)
print(11 * v0 - v10)