""" 
    Simulation with Indirect Effects 

    """

pyenv virtualenv 3.9.17 sim_venv
pyenv activate sim_venv

python3.9 -m pip install --upgrade pip
pip install snipar==0.0.18

# # Simulation

# I run the simulation with indirect effects by 4 random mating generations and 1 AM mating, 
# I configured proper arguments in sim.py

# $ cd '/home/mahdi/Dropbox/1-Git/simulation-with-Snipar'
# cd to sim.py dir

# $ python sim.py

# # Analysis

cd sim_out_no_indir_non_eq

# combining gts
plink --bfile chr_1 --bmerge chr_1_par --out chr_1_off_and_par
plink --bfile chr_1_off_and_par --bmerge chr_1_gpar --out chr_1_combined

# creating a new venv to install grandpar branch
pyenv virtualenv 3.9.17 gpar_venv
pyenv activate gpar_venv

# install the current Snipar's grandpar branch
python3.9 -m pip install --upgrade pip
pip install https://github.com/AlexTISYoung/snipar/archive/grandpar.zip


# PGI using observed gts
pgs.py direct_v1_obs --bed chr_@_combined --weights causal_effects.txt --beta_col direct_v1 --pedigree pedigree.txt --grandpar

# PGI using imputed gpars gts
pgs.py direct_v1_imp --bed chr_1_off_and_par --weights causal_effects.txt --beta_col direct_v1 --imp phased_impute_chr_1_par --grandpar

# kept only offspring PGSs, using script/keep_only_offspring_pgs_with_indir.R

# compute h2f using 'compute_h2f.R'

# k
pgs.py direct_v1_obs --pgs direct_v1_obs_offspring_only.pgs.txt --phenofile phenotype.txt --h2f 0.393986385925266,0 --scale_phen --scale_pgs

# 3 gen model using observed - both parents paternal and maternal 3 gen models
pgs.py direct_v1_obs --pgs direct_v1_obs_offspring_only.pgs.txt --phenofile phenotype.txt --gen_models 3 --scale_phen --scale_pgs

# 3 gen model using imputed - both parents, uni-parental
pgs.py direct_v1_imp --pgs direct_v1_imp_offspring_only.pgs.txt --phenofile phenotype.txt --gen_models 3 --scale_phen --scale_pgs --gparsum

# extracted the table 1 results using scripts/compute_sim_results_4_table1.R