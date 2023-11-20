# % true direct effect, direct effect + v1, direct effect + v10
pgs.py direct_obs --bed chr_@_combined --weights causal_effects.txt --beta_col direct --pedigree pedigree.txt --grandpar
pgs.py direct_v1_obs --bed chr_@_combined --weights causal_effects.txt --beta_col direct_v1 --pedigree pedigree.txt --grandpar
pgs.py direct_v10_obs --bed chr_@_combined --weights causal_effects.txt --beta_col direct_v10 --pedigree pedigree.txt --grandpar

# % I removed the parental gen from pgs files, keeping only offspring

# % 3 gen model
pgs.py direct_obs --pgs direct_obs_last_gen.pgs.txt --phenofile phenotype.txt --gen_models 1-3 --scale_phen --scale_pgs
pgs.py direct_v1_obs --pgs direct_v1_obs_last_gen.pgs.txt --phenofile phenotype.txt --gen_models 1-3 --scale_phen --scale_pgs
pgs.py direct_v10_obs --pgs direct_v10_obs_last_gen.pgs.txt --phenofile phenotype.txt --gen_models 1-3 --scale_phen --scale_pgs

