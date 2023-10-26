"""

    """

import h5py
import numpy as np
from pysnptools.snpreader import Bed
from pysnptools.snpreader import SnpData
from snipar.ibd import write_segs_from_matrix
from snipar.simulate import forward_sim
from snipar.simulate import impute_all_fams
from snipar.simulate import impute_all_fams_phased
from snipar.simulate import simulate_first_gen
from snipar.utilities import encode_str_array

class Args :
    n_causal = 1000
    h2 = 0.5
    outprefix = 'sim/'
    nfam = 30 * 10 ** 3
    impute = False
    n_random = 0
    n_am = 22
    save_par_gts = True
    unphased_impute = False
    r_par = 0.5
    v_indir = 0
    r_dir_indir = None
    beta_vert = 0

def save_phenotype_genotype(n_last ,
                            gen_out_suf ,
                            par_gen_out_suf ,
                            ped ,
                            chroms ,
                            snp_ids ,
                            maps ,
                            positions ,
                            alleles ,
                            new_haps ,
                            haps ,
                            ibd ,
                            a ,
                            args
                            ) :
    last_gen = [x.split('_')[0] == n_last for x in ped[: , 0]]
    phen_out = ped[last_gen , :]
    np.savetxt(args.outprefix + 'phenotype' + gen_out_suf + '.txt' ,
               phen_out[: , [0 , 1 , 5]] ,
               fmt = '%s')

    # Save to HDF5 file
    print('Saving genotypes')
    # save offspring genotypes
    for i in range(len(new_haps)) :
        print('Writing genotypes for chromosome ' + str(chroms[i]))
        bim_i = np.vstack((np.repeat(chroms[i] , snp_ids[i].shape[0]) ,
                           snp_ids[i] , maps[i] , positions[i] ,
                           alleles[i][: , 0] , alleles[i][: , 1])).T
        gts_chr = SnpData(iid = ped[last_gen , 0 :2] ,
                          sid = snp_ids[i] ,
                          pos = bim_i[: , [0 , 2 , 3]] ,
                          val = np.sum(new_haps[i] ,
                                       axis = 3 ,
                                       dtype = np.uint8).reshape((new_haps[
                                                                      i].shape[
                                                                      0] * 2 ,
                                                                  new_haps[
                                                                      i].shape[
                                                                      2])))
        Bed.write(args.outprefix + 'chr_' + str(
                chroms[i]) + gen_out_suf + '.bed' ,
                  gts_chr ,
                  count_A1 = True ,
                  _require_float32_64 = False)
        np.savetxt(args.outprefix + 'chr_' + str(
                chroms[i]) + gen_out_suf + '.bim' ,
                   bim_i ,
                   fmt = '%s' ,
                   delimiter = '\t')

        if args.save_par_gts :
            par_gen = [x.split('_')[0] == str(int(n_last) - 1) for x in
                       ped[: , 0]]
            par_gts_chr = SnpData(iid = ped[par_gen , 0 :2] ,
                                  sid = snp_ids[i] ,
                                  pos = bim_i[: , [0 , 2 , 3]] ,
                                  val = np.sum(haps[i] ,
                                               axis = 3 ,
                                               dtype = np.uint8).reshape((haps[
                                                                              i].shape[
                                                                              0] * 2 ,
                                                                          haps[
                                                                              i].shape[
                                                                              2])))
            Bed.write(args.outprefix + 'chr_' + str(
                    chroms[i]) + par_gen_out_suf + '.bed' ,
                      par_gts_chr ,
                      count_A1 = True ,
                      _require_float32_64 = False)
            del par_gts_chr

            np.savetxt(args.outprefix + 'chr_' + str(
                    chroms[i]) + par_gen_out_suf + '.bim' ,
                       bim_i ,
                       fmt = '%s' ,
                       delimiter = '\t')

        # Imputed parental genotypes
        if args.impute or args.unphased_impute :
            print('Imputing parental genotypes and saving')

            freqs = np.mean(gts_chr.val , axis = 0) / 2.0
            imp_ped = ped[last_gen , 0 :4]
            imp_ped = np.hstack((imp_ped , np.zeros((imp_ped.shape[0] , 2) ,
                                                    dtype = bool)))
            # Phased
            if args.impute :
                hf = h5py.File(args.outprefix + 'phased_impute_chr_' + str(
                        chroms[i]) + gen_out_suf + '.hdf5' , 'w')
                phased_imp = impute_all_fams_phased(new_haps[i] ,
                                                    freqs ,
                                                    ibd[i])
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

            # Unphased
            if args.unphased_impute :
                hf = h5py.File(args.outprefix + 'unphased_impute_chr_' + str(
                        chroms[i]) + gen_out_suf + '.hdf5' , 'w')
                ibd[i] = np.sum(ibd[i] , axis = 2)
                imp = impute_all_fams(gts_chr , freqs , ibd[i])
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
        del gts_chr

    # Write IBD segments
    snp_count = 0
    sibpairs = ped[last_gen , 1]
    sibpairs = sibpairs.reshape((int(sibpairs.shape[0] / 2) , 2))
    if args.v_indir == 0 :
        causal_out = np.zeros((a.shape[0] , 5) , dtype = 'U30')
    else :
        causal_out = np.zeros((a.shape[0] , 5) , dtype = 'U30')
    for i in range(len(haps)) :
        print('Writing IBD segments for chromosome ' + str(chroms[i]))
        # Segments
        if not args.unphased_impute :
            ibd[i] = np.sum(ibd[i] , axis = 2)
        _ = write_segs_from_matrix(ibd[i] ,
                                   sibpairs ,
                                   snp_ids[i] ,
                                   positions[i] ,
                                   maps[i] ,
                                   chroms[i] ,
                                   args.outprefix + 'chr_' + str(chroms[
                                                                     i]) + gen_out_suf + '.segments.gz')
        # Causal effects
        if args.v_indir == 0 :
            a_chr = a[snp_count :(snp_count + snp_ids[i].shape[0])]
            a_chr_v1 = a_chr + np.random.normal(0 , np.std(a_chr) , a_chr.shape)
            causal_out[snp_count :(snp_count + snp_ids[i].shape[0]) ,
            :] = np.vstack((snp_ids[i] , alleles[i][: , 0] , alleles[i][: , 1] ,
                            a_chr , a_chr_v1)).T
            if i == 0 :
                causal_out = np.vstack((np.array(['SNP' , 'A1' , 'A2' ,
                                                  'direct' ,
                                                  'direct_v1']).reshape((1 ,
                                                                         5)) ,
                                        causal_out))
        else :
            a_chr = a[snp_count :(snp_count + snp_ids[i].shape[0]) , :]
            causal_out[snp_count :(snp_count + snp_ids[i].shape[0]) ,
            :] = np.vstack((snp_ids[i] , alleles[i][: , 0] , alleles[i][: , 1] ,
                            a_chr[: , 0] , a_chr[: , 1])).T
            if i == 0 :
                causal_out = np.vstack((
                        np.array(['SNP' , 'A1' , 'A2' , 'direct' ,
                                  'indirect']).reshape((1 , 5)) , causal_out))
        snp_count += snp_ids[i].shape[0]  # count
    np.savetxt(args.outprefix + 'causal_effects' + gen_out_suf + '.txt' ,
               causal_out ,
               fmt = '%s')

##
def main(args) :
    pass

    ##

    print('Simulating an initial generation by random-mating')
    print('Followed by ' + str(args.n_random) + ' generations of random-mating')
    print('Followed by ' + str(args.n_am) + ' generations of assortative mating')

    unlinked = True  # because bgen arg is None

    # simulate base generation
    haps , maps , snp_ids , alleles , positions , chroms = simulate_first_gen(
            args.nfam ,
            args.n_causal ,
            maf = None ,
            min_maf = 0.05)

    # Perform simulation with forward_sim
    new_haps , haps , father_indices , mother_indices , ibd , ped , a , v = forward_sim(
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

    # Save variance components
    print('Saving variance components to ' + args.outprefix + 'VCs.txt')
    np.savetxt(args.outprefix + 'VCs.txt' , v , fmt = '%s')

    # Save pedigree and fam files
    print('Writing pedigree and fam files')
    np.savetxt(args.outprefix + 'pedigree.txt' , ped , fmt = '%s')

    # save phenotype genotype for the last generation
    n_last = ped[ped.shape[0] - 1 , 0].split('_')[0]
    gen_out_suf = ''
    par_gen_out_suf = '_par'

    save_phenotype_genotype(n_last ,
                            gen_out_suf ,
                            par_gen_out_suf ,
                            ped ,
                            chroms ,
                            snp_ids ,
                            maps ,
                            positions ,
                            alleles ,
                            new_haps ,
                            haps ,
                            ibd ,
                            a ,
                            args)

    # save phenotype genotype for the second last generation (parents)
    print('*** Parents Generation ***')

    n_last = str(int(n_last) - 1)
    gen_out_suf = '_par'
    par_gen_out_suf = '_gpar'

    save_phenotype_genotype(n_last ,
                            gen_out_suf ,
                            par_gen_out_suf ,
                            ped ,
                            chroms ,
                            snp_ids ,
                            maps ,
                            positions ,
                            alleles ,
                            new_haps ,
                            haps ,
                            ibd ,
                            a ,
                            args)

##
if __name__ == "__main__" :
    args1 = Args()
    main(args = args1)
