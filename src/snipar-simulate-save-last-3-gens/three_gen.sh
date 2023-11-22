# # combining gts
plink --bfile chr_1 --bmerge chr_1_par --out chr_1_off_and_par
plink --bfile chr_1_off_and_par --bmerge chr_1_gpar --out chr_1_combined

# # PGI
pgs.py direct_obs --bed chr_@_combined --weights causal_effects.txt --beta_col direct --pedigree pedigree.txt --grandpar
pgs.py direct_v1_obs --bed chr_@_combined --weights causal_effects.txt --beta_col direct_v1 --pedigree pedigree.txt --grandpar
pgs.py direct_v10_obs --bed chr_@_combined --weights causal_effects.txt --beta_col direct_v10 --pedigree pedigree.txt --grandpar

# # I removed the parental gen from pgs files, keeping only offspring

# # 3 gen model
pgs.py direct_obs --pgs direct_obs_offspring.pgs.txt --phenofile phenotype.txt --gen_models 1-3 --scale_phen --scale_pgs
pgs.py direct_v1_obs --pgs direct_v1_obs_offspring.pgs.txt --phenofile phenotype.txt --gen_models 1-3 --scale_phen --scale_pgs
pgs.py direct_v10_obs --pgs direct_v10_obs_offspring.pgs.txt --phenofile phenotype.txt --gen_models 1-3 --scale_phen --scale_pgs

# # k
pgs.py direct_obs --pgs direct_obs_offspring.pgs.txt --phenofile phenotype.txt --h2f 0.4123775,0
pgs.py direct_v1_obs --pgs direct_v1_obs_offspring.pgs.txt --phenofile phenotype.txt --h2f 0.4123775,0
pgs.py direct_v10_obs --pgs direct_v10_obs_offspring.pgs.txt --phenofile phenotype.txt --h2f 0.4123775,0

