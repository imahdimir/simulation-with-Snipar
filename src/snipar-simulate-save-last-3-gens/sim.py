"""

    """

import h5py
import numpy as np
from pysnptools.snpreader import Bed
from pysnptools.snpreader import SnpData
from snipar.ibd import write_segs_from_matrix
from snipar.simulate import am_indices
from snipar.simulate import compute_genetic_component
from snipar.simulate import compute_phenotype
from snipar.simulate import compute_phenotype_indirect
from snipar.simulate import compute_phenotype_vert
from snipar.simulate import compute_vcomps
from snipar.simulate import create_ped_output
from snipar.simulate import impute_all_fams_phased , impute_all_fams
from snipar.simulate import produce_next_gen
from snipar.simulate import produce_next_gen_unlinked
from snipar.simulate import random_mating_indices
from snipar.simulate import simulate_effects
from snipar.simulate import simulate_first_gen
from snipar.utilities import encode_str_array

class Args :
    n_causal = 1000
    h2 = 0.5
    outprefix = 'sim/'
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

class Population :
    unlinked = None
    haps = None
    new_haps = None
    ibd = None
    maps = None
    snp_ids = None
    alleles = None
    positions = None
    chroms = None
    f_inds = None
    m_inds = None
    y_ma = None
    y_fe = None
    y_m = None
    y_p = None
    total_matings = None
    pedcols = None
    ped = None
    vcols = None
    v = None
    a = None
    causal = None
    h2 = None
    h2_total = None
    delta_p = None
    delta_m = None
    g_ma = None
    g_fe = None
    eta_p = None
    eta_m = None
    eta_ma = None
    eta_fe = None
    par_ibd = None
    gpar_haps = None

def simulate_1st_gen(p , ar , maf , min_maf) :
    """ Simulate first generation """

    _f = simulate_first_gen
    _o = _f(ar.nfam , ar.n_causal , maf = maf , min_maf = min_maf)

    p.haps , p.maps , p.snp_ids , p.alleles , p.positions , p.chroms = _o

    return p

def gen_next_gen_haps(p: Population) :
    """ Generate haplotypes of new generation """
    p.new_haps = []
    p.ibd = []

    for cho in range(len(p.haps)) :
        in3 = p.haps[cho][: , 0 , : , :]
        in4 = p.haps[cho][: , 1 , : , :]

        if p.unlinked :
            _f = produce_next_gen_unlinked
            new_haps_chr , ibd_chr = _f(p.f_inds , p.m_inds , in3 , in4)
        else :
            in5 = p.maps[cho]
            _f = produce_next_gen
            new_haps_chr , ibd_chr = _f(p.f_inds , p.m_inds , in3 , in4 , in5)

        p.new_haps.append(new_haps_chr)
        p.ibd.append(ibd_chr)

    return p

def prepare_v_and_ped(p , ar) :
    pedcols = {
            'FID'           : None ,
            'IID'           : None ,
            'FATHER_ID'     : None ,
            'MOTHER_ID'     : None ,
            'SEX'           : None ,
            'PHENO'         : None ,
            'FATHER_PHENO'  : None ,
            'MOTHER_PHENO'  : None ,
            'DIRECT'        : None ,
            'FATHER_DIRECT' : None ,
            'MOTHER_DIRECT' : None ,
            }
    p.pedcols = list(pedcols.keys())

    if ar.v_indir == 0 :
        p.v = np.zeros((p.total_matings , 3))
        p.vcols = np.array(['v_g' , 'v_y' , 'r_delta']).reshape((1 , 3))

    else :
        p.v = np.zeros((p.total_matings , 8))
        p.vcols = np.array(['v_g' , 'v_y' , 'r_delta' , 'v_eg' , 'c_ge' ,
                            'r_eta' , 'r_delta_eta_c' ,
                            'r_delta_eta_tau']).reshape((1 , 8))
        p.pedcols += ['INDIRECT' , 'FATHER_INDIRECT' , 'MOTHER_INDIRECT']

    # variance components
    p.v[:] = np.nan

    # pedigree
    p.pedcols = np.array(p.pedcols)
    nro = p.total_matings * 2 * ar.nfam
    p.ped = np.zeros((nro , p.pedcols.shape[0]) , dtype = 'U30')

    return p

def sim_causal_effect_in_1st_gen_wt_no_indir(p , ar) :
    # Simulate direct effect component
    p.a , p.causal , p.h2 = simulate_effects(p.new_haps , ar.n_causal , ar.h2)

    # Compute parental phenotype
    _f = compute_phenotype
    _o = _f(p.haps , p.causal , p.a , 1 - p.h2)
    p.delta_p , p.delta_m , p.y_p , p.y_m = _o

    p.delta_p , p.delta_m = p.delta_p[p.f_inds] , p.delta_m[p.m_inds]
    p.y_p , p.y_m = p.y_p[p.f_inds] , p.y_m[p.m_inds]

    p = compute_phenotype_no_indir_no_beta_vert(p)

    p = rec_v_and_ped_wt_no_indir(p , ar , gen = 0)

    return p

def compute_phenotype_no_indir_no_beta_vert(p) :
    _o = compute_phenotype(p.new_haps , p.causal , p.a , 1 - p.h2)
    p.g_ma , p.g_fe , p.y_ma , p.y_fe = _o
    return p

def sim_causal_effect_in_1st_gen_wt_indir(p , ar) :
    # Compute indirect effect component
    _f = simulate_effects
    p.a , p.causal , p.h2_total = _f(p.new_haps ,
                                     ar.n_causal ,
                                     ar.h2 ,
                                     old_haps = p.haps ,
                                     v_indirect = ar.v_indir ,
                                     r_direct_indirect = ar.r_dir_indir ,
                                     father_indices = p.f_inds ,
                                     mother_indices = p.m_inds)

    # Compute parental phenotype
    _o = compute_phenotype(p.haps , p.causal , p.a[: , 0] , 1 - ar.h2)
    p.delta_p , p.delta_m , p.y_p , p.y_m = _o

    p.y_p , p.y_m = p.y_p[p.f_inds] , p.y_m[p.m_inds]

    p = compute_phenotype_wt_indir(p)

    p = rec_v_and_ped_wt_indir(p , ar , gen = 0)

    return p

def compute_phenotype_wt_indir(p) :
    print('Computing phenotypes')
    _o = compute_phenotype_indirect(p.new_haps ,
                                    p.haps ,
                                    p.f_inds ,
                                    p.m_inds ,
                                    p.causal ,
                                    p.a[: , 0] ,
                                    p.a[: , 1] ,
                                    1 - p.h2_total)
    p.g_ma , p.g_fe , p.y_ma , p.y_fe = _o[:4]
    p.eta_p , p.eta_m , p.delta_p , p.delta_m = _o[4 :]

    _f = compute_genetic_component
    p.eta_ma , p.eta_fe = _f(p.new_haps , p.causal , p.a[: , 1])

    return p

def sim_generations_wt_no_indir(p , ar) :
    for m in range(p.total_matings) :

        p = mate(p , ar , m)

        # generate haplotypes of new generation
        p.haps = p.new_haps
        p = gen_next_gen_haps(p)

        p.delta_p , p.delta_m = p.g_ma[p.f_inds] , p.g_fe[p.m_inds]
        p.y_p , p.y_m = p.y_ma[p.f_inds] , p.y_fe[p.m_inds]

        print('Computing phenotypes')
        if ar.beta_vert != 0 :
            _f = compute_phenotype_vert
            p.g_ma , p.g_fe , p.y_ma , p.y_fe = _f(p.new_haps ,
                                                   p.causal ,
                                                   p.a ,
                                                   1 - p.h2 ,
                                                   ar.beta_vert ,
                                                   p.y_p ,
                                                   p.y_m)
        else :
            p = compute_phenotype_no_indir_no_beta_vert(p)

        p = rec_v_and_ped_wt_no_indir(p , ar , m + 1)

        p = save_gpar_haps_and_par_ibd(p , m + 1)

    return p

def save_gpar_haps_and_par_ibd(p , gen) :
    if gen == p.total_matings - 2 :
        print('*** gpar gen ***')
        p.gpar_haps = p.new_haps

    if gen == p.total_matings - 1 :
        print('*** par gen ***')
        p.par_ibd = p.ibd

    return p

def sim_generations_wt_indir(p , ar) :
    for m in range(p.total_matings) :
        p = mate(p , ar , m)

        # generate haplotypes of new generation
        p.haps = p.new_haps
        p = gen_next_gen_haps(p)

        p = compute_phenotype_wt_indir(p)

        p = rec_v_and_ped_wt_indir(p , ar , m + 1)

        p = save_gpar_haps_and_par_ibd(p , m + 1)

    return p

def rec_v_and_ped_wt_no_indir(p , ar , gen) :
    p.v[gen , :] = compute_vcomps(p.g_ma ,
                                  p.g_fe ,
                                  p.y_ma ,
                                  p.y_fe ,
                                  delta_p = p.delta_p ,
                                  delta_m = p.delta_m)

    ro_strt = gen * 2 * ar.nfam
    ro_end = (gen + 1) * 2 * ar.nfam

    _f = create_ped_output
    p.ped[ro_strt :ro_end , :] = _f(p.g_ma ,
                                    p.g_fe ,
                                    p.y_ma ,
                                    p.y_fe ,
                                    p.y_p ,
                                    p.y_m ,
                                    p.delta_p ,
                                    p.delta_m ,
                                    father_indices = p.f_inds ,
                                    mother_indices = p.m_inds ,
                                    gen = gen + 1 ,
                                    header = False)

    return p

def rec_v_and_ped_wt_indir(p , ar , gen) :
    p.v[gen , :] = compute_vcomps(p.g_ma ,
                                  p.g_fe ,
                                  p.y_ma ,
                                  p.y_fe ,
                                  delta_p = p.delta_p ,
                                  delta_m = p.delta_m ,
                                  v_indirect = ar.v_indir ,
                                  eta_p = p.eta_p ,
                                  eta_m = p.eta_m)

    ro_strt = gen * 2 * ar.nfam
    ro_end = (gen + 1) * 2 * ar.nfam

    _f = create_ped_output
    p.ped[ro_strt :ro_end , :] = _f(p.g_ma ,
                                    p.g_fe ,
                                    p.y_ma ,
                                    p.y_fe ,
                                    p.y_p ,
                                    p.y_m ,
                                    p.delta_p ,
                                    p.delta_m ,
                                    eta_males = p.eta_ma ,
                                    eta_females = p.eta_fe ,
                                    eta_p = p.eta_p ,
                                    eta_m = p.eta_m ,
                                    father_indices = p.f_inds ,
                                    mother_indices = p.m_inds ,
                                    gen = gen + 1 ,
                                    header = False)

    return p

def mate_randomly(p , ar) :
    print('Mating Randomly')
    p.f_inds = random_mating_indices(ar.nfam)
    p.m_inds = random_mating_indices(ar.nfam)
    return p

def mate_assortatively(p , ar) :
    print('Matching assortatively')
    p.f_inds , p.m_inds = am_indices(p.y_ma , p.y_fe , ar.r_par)
    return p

def mate(p , ar , m) :
    print('Mating ' + str(m + 2))
    if m < ar.n_random :
        p = mate_randomly(p , ar)
    else :
        p = mate_assortatively(p , ar)
    return p

def add_v_and_ped_cols(p) :
    p.v = np.vstack((p.vcols , p.v))
    p.ped = np.vstack((p.pedcols , p.ped))
    return p

def forward_sim(p , ar) :
    """ Simulate population """

    # total matings after first gen of random mating
    p.total_matings = ar.n_random + ar.n_am

    p = prepare_v_and_ped(p , ar)

    # first generation by random mating
    print('Generating first generation by random-mating')
    p = mate_randomly(p , ar)

    # Generate haplotypes of new generation
    p = gen_next_gen_haps(p)

    if ar.v_indir == 0 :
        p = sim_causal_effect_in_1st_gen_wt_no_indir(p , ar)
        p = sim_generations_wt_no_indir(p , ar)

    else :
        p = sim_causal_effect_in_1st_gen_wt_indir(p , ar)
        p = sim_generations_wt_indir(p , ar)

    p = add_v_and_ped_cols(p)

    return p

def save_v_and_ped(p , ar) :
    print('Saving variance components')
    np.savetxt(ar.outprefix + 'VCs.txt' , p.v , fmt = '%s')

    print('Writing pedigree')
    np.savetxt(ar.outprefix + 'pedigree.txt' , p.ped , fmt = '%s')

def get_gen_inds(p , gen) :
    n_gen = str(gen)
    gen_inds = [x.split('_')[0] == n_gen for x in p.ped[: , 0]]
    return gen_inds

def save_phenotype_of_offspring_and_par(p , ar) :
    print('saving phenotypes of offspring and parents generations')

    _dc = {
            0 : (p.total_matings , '') ,
            1 : (p.total_matings - 1 , '_par') ,
            }

    for gen , gen_suf in _dc.values() :
        gen_inds = get_gen_inds(p , gen)
        phen_out = p.ped[gen_inds , :]
        _o = phen_out[: , [0 , 1 , 5]]
        _fn = 'phenotype' + gen_suf + '.txt'
        np.savetxt(ar.outprefix + _fn , _o , fmt = '%s')

def gen_bim_i(p , i) :
    bim_i = np.vstack((np.repeat(p.chroms[i] , p.snp_ids[i].shape[0]) ,
                       p.snp_ids[i] , p.maps[i] , p.positions[i] ,
                       p.alleles[i][: , 0] , p.alleles[i][: , 1])).T
    return bim_i

def save_bim_i(p , bim_i , i , gen_suf , ar) :
    np.savetxt(ar.outprefix + 'chr_' + str(p.chroms[i]) + gen_suf + '.bim' ,
               bim_i ,
               fmt = '%s' ,
               delimiter = '\t')

def gen_gts_chr(p , i , bim_i , gen_inds , gen_haps) :
    gts_chr = SnpData(iid = p.ped[gen_inds , 0 :2] ,
                      sid = p.snp_ids[i] ,
                      pos = bim_i[: , [0 , 2 , 3]] ,
                      val = np.sum(gen_haps[i] ,
                                   axis = 3 ,
                                   dtype = np.uint8).reshape((gen_haps[i].shape[
                                                                  0] * 2 ,
                                                              gen_haps[i].shape[
                                                                  2])))

    return gts_chr

def write_gts_as_bed(p , i , gts_chr , gen_suf , ar) :
    Bed.write(ar.outprefix + 'chr_' + str(p.chroms[i]) + gen_suf + '.bed' ,
              gts_chr ,
              count_A1 = True ,
              _require_float32_64 = False)

def save_gts_of_last_three_gens_then_impute_phased_parental_gts(p , ar) :
    print('Saving genotypes for last 3 generations')

    off_inds = get_gen_inds(p , p.total_matings)
    par_inds = get_gen_inds(p , p.total_matings - 1)
    gpar_inds = get_gen_inds(p , p.total_matings - 2)

    _dc = {
            0 : (off_inds , p.new_haps , '' , p.ibd) ,
            1 : (par_inds , p.haps , '_par' , p.par_ibd) ,
            2 : (gpar_inds , p.gpar_haps , '_gpar' , None) ,
            }

    for i in range(len(p.new_haps)) :
        print('Writing genotypes for chromosome ' + str(p.chroms[i]))

        bim_i = gen_bim_i(p , i)

        for g_inds , g_haps , g_suf , g_ibd in _dc.values() :
            save_bim_i(p , bim_i , i , g_suf , ar)

            gts_chr = gen_gts_chr(p , i , bim_i , g_inds , g_haps)

            write_gts_as_bed(p , i , gts_chr , g_suf , ar)

            if g_suf == '_gpar' :
                # skip imputaion of parents of grandparental generation
                continue

            print('Imputing parental genotypes and saving for gen' + g_suf)

            freqs = np.mean(gts_chr.val , axis = 0) / 2.0
            imp_ped = p.ped[g_inds , 0 :4]
            imp_ped = np.hstack((imp_ped , np.zeros((imp_ped.shape[0] , 2) ,
                                                    dtype = bool)))

            print('phased')
            hf = h5py.File(ar.outprefix + 'phased_impute_chr_' + str(
                    p.chroms[i]) + g_suf + '.hdf5' , 'w')

            phased_imp = impute_all_fams_phased(g_haps[i] , freqs , g_ibd[i])
            hf['imputed_par_gts'] = phased_imp
            del phased_imp
            hf['bim_values'] = encode_str_array(bim_i)
            hf['bim_columns'] = encode_str_array(np.array(['rsid' , 'map' ,
                                                           'position' ,
                                                           'allele1' ,
                                                           'allele2']))
            hf['pedigree'] = encode_str_array(imp_ped)
            hf['families'] = encode_str_array(imp_ped[0 : :2 , 0])
            hf.close()

            print('unphased')
            hf = h5py.File(ar.outprefix + 'unphased_impute_chr_' + str(
                    p.chroms[i]) + g_suf + '.hdf5' , 'w')

            uibd = g_ibd[i].copy()
            uibd[i] = np.sum(g_ibd[i] , axis = 2)

            imp = impute_all_fams(gts_chr , freqs , uibd[i])
            hf['imputed_par_gts'] = imp
            del imp
            hf['bim_values'] = encode_str_array(bim_i)
            hf['bim_columns'] = encode_str_array(np.array(['rsid' , 'map' ,
                                                           'position' ,
                                                           'allele1' ,
                                                           'allele2']))
            hf['pedigree'] = encode_str_array(imp_ped)
            hf['families'] = encode_str_array(imp_ped[0 : :2 , 0])
            hf.close()

def write_ibd_segs_of_offsrping_and_par(p , ar) :
    print('Write IBD segments of offspring and parents')

    off_inds = get_gen_inds(p , p.total_matings)
    par_inds = get_gen_inds(p , p.total_matings - 1)

    _dc = {
            0 : (off_inds , p.ibd , '') ,
            1 : (par_inds , p.par_ibd , '_par') ,
            }

    for g_inds , ibd , g_suf in _dc.values() :
        sibpairs = p.ped[g_inds , 1]
        sibpairs = sibpairs.reshape((int(sibpairs.shape[0] / 2) , 2))

        for i in range(len(p.haps)) :
            print('Writing IBD segments for chromosome ' + str(p.chroms[i]))

            ibd[i] = np.sum(ibd[i] , axis = 2)

            _fp = ar.outprefix + 'chr_'
            _fp += str(p.chroms[i]) + g_suf + '.segments.gz'

            _ = write_segs_from_matrix(ibd[i] ,
                                       sibpairs ,
                                       p.snp_ids[i] ,
                                       p.positions[i] ,
                                       p.maps[i] ,
                                       p.chroms[i] ,
                                       _fp)

def make_causal_out_and_save(p , ar) :
    causal_out = np.zeros((p.a.shape[0] , 5) , dtype = 'U30')
    snp_count = 0

    if ar.v_indir == 0 :
        causal_out = np.vstack((np.array(['SNP' , 'A1' , 'A2' , 'direct' ,
                                          'direct_v1']).reshape((1 , 5)) ,
                                causal_out))

        for i in range(len(p.haps)) :
            a_chr = p.a[snp_count :(snp_count + p.snp_ids[i].shape[0])]
            a_chr_v1 = a_chr + np.random.normal(0 , np.std(a_chr) , a_chr.shape)
            causal_out[snp_count :(snp_count + p.snp_ids[i].shape[0]) ,
            :] = np.vstack((p.snp_ids[i] , p.alleles[i][: , 0] ,
                            p.alleles[i][: , 1] , a_chr , a_chr_v1)).T

            snp_count += p.snp_ids[i].shape[0]


    else :
        causal_out = np.vstack((
                np.array(['SNP' , 'A1' , 'A2' , 'direct' , 'indirect']).reshape(
                        (1 , 5)) , causal_out))

        for i in range(len(p.haps)) :
            a_chr = p.a[snp_count :(snp_count + p.snp_ids[i].shape[0]) , :]
            causal_out[snp_count :(snp_count + p.snp_ids[i].shape[0]) ,
            :] = np.vstack((p.snp_ids[i] , p.alleles[i][: , 0] ,
                            p.alleles[i][: , 1] , a_chr[: , 0] ,
                            a_chr[: , 1])).T

            snp_count += p.snp_ids[i].shape[0]

    np.savetxt(ar.outprefix + 'causal_effects.txt' , causal_out , fmt = '%s')

def main(ar) :
    print('Simulating an initial generation by random-mating')
    print('Followed by ' + str(ar.n_random) + ' generations of random-mating')
    print('Followed by ' + str(ar.n_am) + ' generations of assortative mating')

    p = Population()
    p.unlinked = True  # since bgen arg is None

    p = simulate_1st_gen(p , ar , maf = None , min_maf = .05)

    # Perform simulation
    p = forward_sim(p , ar)

    save_v_and_ped(p , ar)

    save_phenotype_of_offspring_and_par(p , ar)

    save_gts_of_last_three_gens_then_impute_phased_parental_gts(p , ar)

    write_ibd_segs_of_offsrping_and_par(p , ar)

    make_causal_out_and_save(p , ar)

if __name__ == "__main__" :
    ar0 = Args()
    main(ar = ar0)
