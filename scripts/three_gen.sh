# # I redone the simulation, adding whole last 3 gen pedigree to both imputed gts files.

cd sim

# # combining gts
plink --bfile chr_1 --bmerge chr_1_par --out chr_1_off_and_par
plink --bfile chr_1_off_and_par --bmerge chr_1_gpar --out chr_1_combined

# # create a new venv in order to install the grandpar branch
pyenv virtualenv 3.9.17 gpar
pyenv activate gpar

# # install the Snipar's grandpar branch, [this is a forked repo in which I have solved the fancy indexing error by np.nonzero()]
pip install https://github.com/imahdimir/snipar_h5py_indxing_err/archive/grandpar.zip

# # PGI using observed gts
pgs.py direct_obs --bed chr_@_combined --weights causal_effects.txt --beta_col direct --pedigree pedigree.txt --grandpar
pgs.py direct_v1_obs --bed chr_@_combined --weights causal_effects.txt --beta_col direct_v1 --pedigree pedigree.txt --grandpar
pgs.py direct_v10_obs --bed chr_@_combined --weights causal_effects.txt --beta_col direct_v10 --pedigree pedigree.txt --grandpar

# # PGI using imputed gpars gts
pgs.py direct_imp --bed chr_1_off_and_par --weights causal_effects.txt --beta_col direct --imp phased_impute_chr_1_par --grandpar
pgs.py direct_v1_imp --bed chr_1_off_and_par --weights causal_effects.txt --beta_col direct_v1 --imp phased_impute_chr_1_par --grandpar
pgs.py direct_v10_imp --bed chr_1_off_and_par --weights causal_effects.txt --beta_col direct_v10 --imp phased_impute_chr_1_par --grandpar

# # I kept only the offspring PGSs using 'keep_only_offspring_pgs.R'

# # compute h2f using 'compute_h2f.R'

# # k
pgs.py direct_obs --pgs direct_obs_offspring.pgs.txt --phenofile phenotype.txt --h2f 0.418828,0  --scale_phen --scale_pgs
pgs.py direct_v1_obs --pgs direct_v1_obs_offspring.pgs.txt --phenofile phenotype.txt --h2f 0.418828,0 --scale_phen --scale_pgs
pgs.py direct_v10_obs --pgs direct_v10_obs_offspring.pgs.txt --phenofile phenotype.txt --h2f 0.418828,0 --scale_phen --scale_pgs

# # 3 gen model using observed - both parents paternal and maternal 3 gen models
pgs.py direct_obs --pgs direct_obs_offspring.pgs.txt --phenofile phenotype.txt --gen_models 3 --scale_phen --scale_pgs
pgs.py direct_v1_obs --pgs direct_v1_obs_offspring.pgs.txt --phenofile phenotype.txt --gen_models 3 --scale_phen --scale_pgs
pgs.py direct_v10_obs --pgs direct_v10_obs_offspring.pgs.txt --phenofile phenotype.txt --gen_models 3 --scale_phen --scale_pgs

# # 3 gen model using imputed - both parents, paternal and maternal 3 gen models
pgs.py direct_imp --pgs direct_imp_offspring.pgs.txt --phenofile phenotype.txt --gen_models 3 --scale_phen --scale_pgs --gparsum
pgs.py direct_v1_imp --pgs direct_v1_imp_offspring.pgs.txt --phenofile phenotype.txt --gen_models 3 --scale_phen --scale_pgs --gparsum
pgs.py direct_v10_imp --pgs direct_v10_imp_offspring.pgs.txt --phenofile phenotype.txt --gen_models 3 --scale_phen --scale_pgs --gparsum

# # then I extract the table_1 using "table_1_simulation_results.R"
