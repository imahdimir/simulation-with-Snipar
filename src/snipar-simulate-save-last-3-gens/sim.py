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
from snipar.simulate import impute_all_fams
from snipar.simulate import impute_all_fams_phased
from snipar.simulate import produce_next_gen
from snipar.simulate import produce_next_gen_unlinked
from snipar.simulate import random_mating_indices
from snipar.simulate import simulate_effects
from snipar.simulate import simulate_first_gen
from snipar.utilities import encode_str_array

class Args :
    n_causal = 1000
    h2 = 0.5
    nfam = 30 * 10 ** 3
    n_random = 0
    n_am = 23
    save_par_gts = True
    impute = True
    r_par = 0.5
    v_indir = 0

    outprefix = 'sim/'

def gen_next_generation_haps(unlinked , haps , maps , f_inds , m_inds) :
    """ Generate haplotypes of new generation """
    new_haps = []
    ibd = []
    for chro in range(0 , len(haps)) :
        in3 = haps[chro][: , 0 , : , :]
        in4 = haps[chro][: , 1 , : , :]

        if unlinked :
            _fu = produce_next_gen_unlinked
            new_haps_chr , ibd_chr = _fu(f_inds , m_inds , in3 , in4)
        else :
            _fu = produce_next_gen
            in5 = maps[chro]
            new_haps_chr , ibd_chr = _fu(f_inds , m_inds , in3 , in4 , in5)

        new_haps.append(new_haps_chr)
        ibd.append(ibd_chr)

    return new_haps , ibd

def prepare_v_and_ped(v_indirect , total_matings , nfam) :
    # Pedigree columns
    pedcols = ['FID' , 'IID' , 'FATHER_ID' , 'MOTHER_ID' , 'SEX' , 'PHENO' ,
               'FATHER_PHENO' , 'MOTHER_PHENO' , 'DIRECT' , 'FATHER_DIRECT' ,
               'MOTHER_DIRECT']

    # Record variance component evolution and pedigree
    if v_indirect == 0 :
        V = np.zeros((total_matings + 1 , 3))
        V_header = np.array(['v_g' , 'v_y' , 'r_delta']).reshape((1 , 3))
    else :
        V = np.zeros((total_matings + 1 , 8))
        V_header = np.array(['v_g' , 'v_y' , 'r_delta' , 'v_eg' , 'c_ge' ,
                             'r_eta' , 'r_delta_eta_c' ,
                             'r_delta_eta_tau']).reshape((1 , 8))
        pedcols += ['INDIRECT' , 'FATHER_INDIRECT' , 'MOTHER_INDIRECT']

    V[:] = np.nan

    # Record full pedigree
    pedcols = np.array(pedcols)
    ped = np.zeros((nfam * 2 * (total_matings + 1) , pedcols.shape[0]) ,
                   dtype = 'U30')

    return pedcols , ped , V , V_header

def forward_sim(haps ,
                maps ,
                ngen_random ,
                ngen_am ,
                unlinked ,
                n_causal ,
                h2 ,
                v_indirect = 0 ,
                r_direct_indirect = 0 ,
                beta_vert = 0 ,
                r_y = None
                ) :
    """ Simulate population """

    nfam = haps[0].shape[0]

    print('Generating first generation by random-mating')
    # father and monther indices
    f_inds = random_mating_indices(nfam)
    m_inds = random_mating_indices(nfam)

    # Generate haplotypes of new generation
    _f = gen_next_generation_haps
    new_haps , ibd = _f(unlinked , haps , maps , f_inds , m_inds)

    # Matings
    total_matings = ngen_random + ngen_am

    _f = prepare_v_and_ped
    pedcols , ped , V , V_header = _f(v_indirect , total_matings , nfam)

    # Count generations
    a_count = 0

    # Simulate ngen_random generations of random-mating then ngen_am generations of assortative mating
    for gen in range(0 , total_matings) :
        if v_indirect == 0 :
            if a_count == 0 :
                # Simulate direct effect component
                a , causal , h2 = simulate_effects(new_haps , n_causal , h2)
                # Compute parental phenotype
                delta_p , delta_m , Y_p , Y_m = compute_phenotype(haps ,
                                                                  causal ,
                                                                  a ,
                                                                  1 - h2)
                delta_p , delta_m , Y_p , Y_m = delta_p[f_inds] , delta_m[
                    m_inds] , Y_p[f_inds] , Y_m[m_inds]
            else :
                delta_p , delta_m , Y_p , Y_m = G_males[f_inds] , G_females[
                    m_inds] , Y_males[f_inds] , Y_females[f_inds]
            # Compute phenotypes
            print('Computing phenotypes')
            if np.abs(beta_vert) > 0 and a_count > 0 :
                G_males , G_females , Y_males , Y_females = compute_phenotype_vert(
                        new_haps ,
                        causal ,
                        a ,
                        1 - h2 ,
                        beta_vert ,
                        Y_p ,
                        Y_m)
            else :
                G_males , G_females , Y_males , Y_females = compute_phenotype(
                        new_haps ,
                        causal ,
                        a ,
                        1 - h2)
        else :
            # Compute with indirect effects
            if a_count == 0 :
                # Compute indirect effect component
                a , causal , h2_total = simulate_effects(new_haps ,
                                                         n_causal ,
                                                         h2 ,
                                                         old_haps = haps ,
                                                         v_indirect = v_indirect ,
                                                         r_direct_indirect = r_direct_indirect ,
                                                         father_indices = f_inds ,
                                                         mother_indices = m_inds)
                # Compute parental phenotype
                delta_p , delta_m , Y_p , Y_m = compute_phenotype(haps ,
                                                                  causal ,
                                                                  a[: , 0] ,
                                                                  1 - h2)
                Y_p , Y_m = Y_p[f_inds] , Y_m[m_inds]
            # Compute phenotype
            print('Computing phenotypes')
            G_males , G_females , Y_males , Y_females , eta_p , eta_m , delta_p , delta_m = compute_phenotype_indirect(
                    new_haps ,
                    haps ,
                    f_inds ,
                    m_inds ,
                    causal ,
                    a[: , 0] ,
                    a[: , 1] ,
                    1 - h2_total)
            # Indirect effect component generation
            eta_males , eta_females = compute_genetic_component(new_haps ,
                                                                causal ,
                                                                a[: , 1])
        # Record variance components
        # Final offspring generation variance components
        if v_indirect == 0 :
            V[a_count , :] = compute_vcomps(G_males ,
                                            G_females ,
                                            Y_males ,
                                            Y_females ,
                                            delta_p = delta_p ,
                                            delta_m = delta_m)
        else :
            V[a_count , :] = compute_vcomps(G_males ,
                                            G_females ,
                                            Y_males ,
                                            Y_females ,
                                            delta_p = delta_p ,
                                            delta_m = delta_m ,
                                            v_indirect = v_indirect ,
                                            eta_p = eta_p ,
                                            eta_m = eta_m)
        # Record pedigree
        if v_indirect == 0 :
            ped[(a_count * 2 * nfam) :((a_count + 1) * 2 * nfam) ,
            :] = create_ped_output(G_males ,
                                   G_females ,
                                   Y_males ,
                                   Y_females ,
                                   Y_p ,
                                   Y_m ,
                                   delta_p ,
                                   delta_m ,
                                   father_indices = f_inds ,
                                   mother_indices = m_inds ,
                                   gen = a_count + 1 ,
                                   header = False)
        else :
            ped[(a_count * 2 * nfam) :((a_count + 1) * 2 * nfam) ,
            :] = create_ped_output(G_males ,
                                   G_females ,
                                   Y_males ,
                                   Y_females ,
                                   Y_p ,
                                   Y_m ,
                                   delta_p ,
                                   delta_m ,
                                   eta_males = eta_males ,
                                   eta_females = eta_females ,
                                   eta_p = eta_p ,
                                   eta_m = eta_m ,
                                   father_indices = f_inds ,
                                   mother_indices = m_inds ,
                                   gen = a_count + 1 ,
                                   header = False)
        a_count += 1

        # Match current generation into parent-pairs for next generation
        print('Mating ' + str(gen + 2))
        # Random mating
        if gen < ngen_random :
            print('Matching randomly')
            f_inds = random_mating_indices(nfam)
            m_inds = random_mating_indices(nfam)
        # Assortative mating
        if gen >= ngen_random :
            # Match assortatively
            print('Matching assortatively')
            f_inds , m_inds = am_indices(Y_males , Y_females , r_y)
            print('Parental phenotypic correlation: ' + str(round(
                    np.corrcoef(Y_males[f_inds] , Y_females[m_inds])[0 , 1] ,
                    4)))

        if gen == total_matings - 1 - 2 :
            print('*** gpar gen ***')
            haps_gpar = new_haps

        elif gen == total_matings - 1 - 1 :
            print('*** par gen ***')
            haps_par = new_haps
            ibd_par = ibd
            a_par = a

        elif gen == total_matings - 1 :
            print('*** offspring gen ***')
            haps_off = new_haps
            ibd_off = ibd
            a_off = a

            print('Sibling genotypic correlation: ' + str(round(
                    np.corrcoef(G_males , G_females)[0 , 1] , 4)))
            print('Sibling phenotypic correlation: ' + str(round(
                    np.corrcoef(Y_males , Y_females)[0 , 1] , 4)))

        else :
            # Generate haplotypes of new generation
            haps = new_haps
            _f = gen_next_generation_haps
            new_haps , ibd = _f(unlinked , haps , maps , f_inds , m_inds)

    V = np.vstack((V_header , V))
    ped = np.vstack((pedcols , ped))

    return haps_gpar , haps_par , ibd_par , a_par , haps_off , ibd_off , a_off , ped , V

def main(args) :
    pass

    ##
    print('Simulating an initial generation by random-mating')
    print('Followed by ' + str(args.n_random) + ' generations of random-mating')
    print('Followed by ' + str(args.n_am) + ' generations of assortative mating')

    unlinked = True  # since bgen arg is None

    ##
    # Simulate 1st gen
    haps , maps , snp_ids , alleles , positions , chroms = simulate_first_gen(
            args.nfam ,
            args.n_causal ,
            maf = None ,
            min_maf = 0.05)

    ##
    # Perform simulation
    haps_gpar , haps_par , ibd_par , a_par , haps_off , ibd_off , a_off , ped , V = forward_sim(
            haps ,
            maps ,
            args.n_random ,
            args.n_am ,
            unlinked ,
            args.n_causal ,
            args.h2 ,
            v_indirect = args.v_indir ,
            r_direct_indirect = args.r_dir_indir ,
            r_y = args.r_par ,
            beta_vert = args.beta_vert)

    ##
    # rm empty rows of ped
    print('ped: \n' , ped)

    ped = ped[~np.all(ped == '' , axis = 1)]

    print('ped: \n' , ped)

    ##
    print('Saving variance components')
    np.savetxt(args.outprefix + 'VCs.txt' , V , fmt = '%s')

    ##
    print('Writing pedigree')
    np.savetxt(args.outprefix + 'pedigree.txt' , ped , fmt = '%s')

    ##
    print('saving phenotypes of offspring and parents generations')

    # offspring
    n_last = ped[ped.shape[0] - 1 , 0].split('_')[0]
    print('n_last: ' , n_last)

    last_gen = [x.split('_')[0] == n_last for x in ped[: , 0]]
    phen_out = ped[last_gen , :]
    np.savetxt(args.outprefix + 'phenotype.txt' ,
               phen_out[: , [0 , 1 , 5]] ,
               fmt = '%s')

    # parents
    n_par = str(int(n_last) - 1)
    par_gen = [x.split('_')[0] == n_par for x in ped[: , 0]]
    phen_par = ped[par_gen , :]
    np.savetxt(args.outprefix + 'phenotype_par.txt' ,
               phen_par[: , [0 , 1 , 5]] ,
               fmt = '%s')

    # gpars
    n_gpar = str(int(n_last) - 2)
    gpar_gen = [x.split('_')[0] == n_gpar for x in ped[: , 0]]

    ##
    print('Saving genotypes for last 3 generations')

    _dc = {
            0 : (last_gen , haps_off , '' , ibd_off) ,
            1 : (par_gen , haps_par , '_par' , ibd_par) ,
            2 : (gpar_gen , haps_gpar , '_gpar' , None) ,
            }

    for i in range(len(haps_off)) :
        print('Writing genotypes for chromosome ' + str(chroms[i]))

        bim_i = np.vstack((np.repeat(chroms[i] , snp_ids[i].shape[0]) ,
                           snp_ids[i] , maps[i] , positions[i] ,
                           alleles[i][: , 0] , alleles[i][: , 1])).T

        for gen , gen_haps , gen_suf , ibd in _dc.values() :
            gts_chr = SnpData(iid = ped[gen , 0 :2] ,
                              sid = snp_ids[i] ,
                              pos = bim_i[: , [0 , 2 , 3]] ,
                              val = np.sum(gen_haps[i] ,
                                           axis = 3 ,
                                           dtype = np.uint8).reshape((gen_haps[
                                                                          i].shape[
                                                                          0] * 2 ,
                                                                      gen_haps[
                                                                          i].shape[
                                                                          2])))
            Bed.write(args.outprefix + 'chr_' + str(
                    chroms[i]) + gen_suf + '.bed' ,
                      gts_chr ,
                      count_A1 = True ,
                      _require_float32_64 = False)

            np.savetxt(args.outprefix + 'chr_' + str(
                    chroms[i]) + gen_suf + '.bim' ,
                       bim_i ,
                       fmt = '%s' ,
                       delimiter = '\t')

            if gen_suf == '_gpar' :
                continue  # skip imputation for grandparents

            print('Imputing parental genotypes and saving for gen' + gen_suf)
            freqs = np.mean(gts_chr.val , axis = 0) / 2.0
            imp_ped = ped[gen , 0 :4]
            imp_ped = np.hstack((imp_ped , np.zeros((imp_ped.shape[0] , 2) ,
                                                    dtype = bool)))

            print('phased')

            hf = h5py.File(args.outprefix + 'phased_impute_chr_' + str(
                    chroms[i]) + gen_suf + '.hdf5' , 'w')
            phased_imp = impute_all_fams_phased(gen_haps[i] , freqs , ibd[i])
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

    ##
    print('Write IBD segments')

    _dc = {
            0 : (last_gen , a_off , haps_off , ibd_off , '') ,
            1 : (par_gen , a_par , haps_par , ibd_par , '_par') ,
            }

    for gen , a , haps , ibd , gen_suf in _dc.values() :

        snp_count = 0
        sibpairs = ped[gen , 1]
        sibpairs = sibpairs.reshape((int(sibpairs.shape[0] / 2) , 2))

        if args.v_indir == 0 :
            causal_out = np.zeros((a.shape[0] , 5) , dtype = 'U30')
        else :
            causal_out = np.zeros((a.shape[0] , 5) , dtype = 'U30')

        for i in range(len(haps)) :
            print('gen: ' , gen_suf)
            print('Writing IBD segments for chromosome ' + str(chroms[i]))
            # Segments
            ibd1 = ibd.copy()

            if not args.unphased_impute :
                ibd1[i] = np.sum(ibd[i] , axis = 2)

            _ = write_segs_from_matrix(ibd1[i] ,
                                       sibpairs ,
                                       snp_ids[i] ,
                                       positions[i] ,
                                       maps[i] ,
                                       chroms[i] ,
                                       args.outprefix + 'chr_' + str(chroms[
                                                                         i]) + gen_suf + '.segments.gz')
            # Causal effects
            if args.v_indir == 0 :
                a_chr = a[snp_count :(snp_count + snp_ids[i].shape[0])]
                a_chr_v1 = a_chr + np.random.normal(0 ,
                                                    np.std(a_chr) ,
                                                    a_chr.shape)
                causal_out[snp_count :(snp_count + snp_ids[i].shape[0]) ,
                :] = np.vstack((snp_ids[i] , alleles[i][: , 0] ,
                                alleles[i][: , 1] , a_chr , a_chr_v1)).T
                if i == 0 :
                    causal_out = np.vstack((np.array(['SNP' , 'A1' , 'A2' ,
                                                      'direct' ,
                                                      'direct_v1']).reshape((1 ,
                                                                             5)) ,
                                            causal_out))
            else :
                a_chr = a[snp_count :(snp_count + snp_ids[i].shape[0]) , :]
                causal_out[snp_count :(snp_count + snp_ids[i].shape[0]) ,
                :] = np.vstack((snp_ids[i] , alleles[i][: , 0] ,
                                alleles[i][: , 1] , a_chr[: , 0] ,
                                a_chr[: , 1])).T
                if i == 0 :
                    causal_out = np.vstack((
                            np.array(['SNP' , 'A1' , 'A2' , 'direct' ,
                                      'indirect']).reshape((1 , 5)) ,
                            causal_out))
            snp_count += snp_ids[i].shape[0]

        np.savetxt(args.outprefix + 'causal_effects' + gen_suf + '.txt' ,
                   causal_out ,
                   fmt = '%s')

if __name__ == "__main__" :
    args1 = Args()
    main(args = args1)
