"""

    """

from simulation_code.simulate import main

class WithIndirArgs :
    n_causal = 1000
    h2 = .5
    outprefix = 'sim_output_wt_indir/'
    nfam = 30 * 10 ** 3
    n_random = 0
    n_am = 22
    save_par_gts = True
    impute = True
    unphased_impute = False
    r_par = .5
    v_indir = .5
    r_dir_indir = .5
    beta_vert = 0

class NoIndirArgs :
    n_causal = 1000
    h2 = 0.5
    outprefix = 'sim_output_no_indir/'
    nfam = 30 * 10 ** 3
    n_random = 0
    n_am = 22
    save_par_gts = True
    impute = True
    unphased_impute = False
    r_par = 0.5
    v_indir = 0
    r_dir_indir = None
    beta_vert = 0

if __name__ == '__main__' :
    ar = WithIndirArgs()
    main(ar)
