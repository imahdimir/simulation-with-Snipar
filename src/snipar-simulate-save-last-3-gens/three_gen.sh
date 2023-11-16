# combining offspirng and parents genotypes
plink --bfile chr_1 --bmerge chr_1_par --out chr_1_off_and_par
plink --bfile chr_1_off_and_par --bmerge chr_1_gpar --out chr_1_combined

# % true direct effect
pgs.py direct_obs --bed chr_@_combined --weights causal_effects.txt --beta_col direct --pedigree pedigree.txt --grandpar
pgs.py direct_obs --pgs direct_obs.pgs.txt --phenofile phenotype.txt --gen_models 1-3

# % ture direct effect + noise(v1)
pgs.py direct_v1_obs --bed chr_@_combined --weights causal_effects.txt --beta_col direct_v1 --pedigree pedigree.txt --grandpar
pgs.py direct_v1_obs --pgs direct_v1_obs.pgs.txt --phenofile phenotype.txt --gen_models 1-3
