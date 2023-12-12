""" 
    Simulation with Indirect Effects 

    """

pyenv virtualenv 3.9.17 sim_venv
pyenv activate sim_venv

python3.9 -m pip install --upgrade pip
pip install snipar==0.0.18

# cd to sim.py dir
# $ cd '/home/mahdi/Dropbox/1-Git/simulation-with-Snipar'

# run the simulation with indirect effects, I configured proper arguments in sim.py
# $ python sim.py

cd sim_out_wt_indir

# combining gts
plink --bfile chr_1 --bmerge chr_1_par --out chr_1_off_and_par
plink --bfile chr_1_off_and_par --bmerge chr_1_gpar --out chr_1_combined

# creating a new venv to install grandpar branch
pyenv virtualenv 3.9.17 gpar_venv
pyenv activate gpar_venv

# install the current Snipar's grandpar branch
python3.9 -m pip install --upgrade pip
pip install https://github.com/AlexTISYoung/snipar/archive/grandpar.zip

# I created a new column sum = direct + indirect in causal_effects.txt using archive_code/a_v10.py
# then I created sum_v1 and sum_v10 by adding noise to sum column


# PGI using observed gts
pgs.py sum_obs --bed chr_@_combined --weights causal_effects.txt --beta_col sum --pedigree pedigree.txt --grandpar
pgs.py sum_v1_obs --bed chr_@_combined --weights causal_effects.txt --beta_col sum_v1 --pedigree pedigree.txt --grandpar
pgs.py sum_v10_obs --bed chr_@_combined --weights causal_effects.txt --beta_col sum_v10 --pedigree pedigree.txt --grandpar

# PGI using imputed gpars gts
pgs.py sum_imp --bed chr_1_off_and_par --weights causal_effects.txt --beta_col sum --imp phased_impute_chr_1_par --grandpar
pgs.py sum_v1_imp --bed chr_1_off_and_par --weights causal_effects.txt --beta_col sum_v1 --imp phased_impute_chr_1_par --grandpar
pgs.py sum_v10_imp --bed chr_1_off_and_par --weights causal_effects.txt --beta_col sum_v10 --imp phased_impute_chr_1_par --grandpar

# I kept only the offspring PGSs using script/keep_only_offspring_pgs_with_indir.R

# compute h2f using 'compute_h2f.R'

# k
pgs.py sum_obs --pgs sum_obs_offspring_only.pgs.txt --phenofile phenotype.txt --h2f 0.221676198331546,0  --scale_phen --scale_pgs
pgs.py sum_v1_obs --pgs sum_v1_obs_offspring_only.pgs.txt --phenofile phenotype.txt --h2f 0.221676198331546,0 --scale_phen --scale_pgs
pgs.py sum_v10_obs --pgs sum_v10_obs_offspring_only.pgs.txt --phenofile phenotype.txt --h2f 0.221676198331546,0 --scale_phen --scale_pgs

# 3 gen model using observed - both parents paternal and maternal 3 gen models
pgs.py sum_obs --pgs sum_obs_offspring_only.pgs.txt --phenofile phenotype.txt --gen_models 3 --scale_phen --scale_pgs
pgs.py sum_v1_obs --pgs sum_v1_obs_offspring_only.pgs.txt --phenofile phenotype.txt --gen_models 3 --scale_phen --scale_pgs
pgs.py sum_v10_obs --pgs sum_v10_obs_offspring_only.pgs.txt --phenofile phenotype.txt --gen_models 3 --scale_phen --scale_pgs

# 3 gen model using imputed - both parents, paternal and maternal 3 gen models
pgs.py sum_imp --pgs sum_imp_offspring_only.pgs.txt --phenofile phenotype.txt --gen_models 3 --scale_phen --scale_pgs --gparsum
pgs.py sum_v1_imp --pgs sum_v1_imp_offspring_only.pgs.txt --phenofile phenotype.txt --gen_models 3 --scale_phen --scale_pgs --gparsum
pgs.py sum_v10_imp --pgs sum_v10_imp_offspring_only.pgs.txt --phenofile phenotype.txt --gen_models 3 --scale_phen --scale_pgs --gparsum

# then I extracted the table_1 result to fill up poulation effects rows using scripts/table_1_simulation_results.R