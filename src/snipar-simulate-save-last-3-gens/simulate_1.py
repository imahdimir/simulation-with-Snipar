"""

    """

import h5py
import numpy as np
from pysnptools.snpreader import Bed
from pysnptools.snpreader import SnpData
from snipar.ibd import write_segs_from_matrix
from snipar.simulate import impute_all_fams
from snipar.simulate import impute_all_fams_phased
from snipar.utilities import encode_str_array

from simulate_0 import Args

def main(args) :
    pass

    ##
    import numpy as np

    np_fn = {
            'maps'      : None ,
            'snp_ids'   : None ,
            'alleles'   : None ,
            'positions' : None ,
            'chroms'    : None ,
            'haps_gpar' : None ,
            'haps_par'  : None ,
            'haps_off'  : None ,
            'ibd_par'   : None ,
            'ibd_off'   : None ,
            'a_par'     : None ,
            'a_off'     : None ,
            'ped'       : None ,
            'V'         : None
            }

    for ky in np_fn.keys() :
        np_fn[ky] = np.load('sim_0/' + ky + '.npy')

    maps = np_fn['maps']
    snp_ids = np_fn['snp_ids']
    alleles = np_fn['alleles']
    positions = np_fn['positions']
    chroms = np_fn['chroms']
    haps_gpar = np_fn['haps_gpar']
    haps_par = np_fn['haps_par']
    haps_off = np_fn['haps_off']
    ibd_par = np_fn['ibd_par']
    ibd_off = np_fn['ibd_off']
    a_par = np_fn['a_par']
    a_off = np_fn['a_off']
    ped = np_fn['ped']
    V = np_fn['V']

    ##
    print('Saving variance components')
    np.savetxt(args.outprefix + 'VCs.txt' , V , fmt = '%s')

    ##
    print('Writing pedigree')
    # np.savetxt(args.outprefix + 'pedigree.txt' , ped , fmt = '%s')
    print(ped[-1 : -5 , :])

    ##
    print('saving phenotypes of offspring and parents generations')

    # offspring
    n_last = ped[ped.shape[0] - 1 , 0].split('_')[0]
    print(n_last)

    ##

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
            (last_gen , haps_off , '' , ibd_off)    : None ,
            (par_gen , haps_par , '_par' , ibd_par) : None ,
            (gpar_gen , haps_gpar , '_gpar' , None) : None ,
            }

    for i in range(len(haps_off)) :
        print('Writing genotypes for chromosome ' + str(chroms[i]))

        bim_i = np.vstack((np.repeat(chroms[i] , snp_ids[i].shape[0]) ,
                           snp_ids[i] , maps[i] , positions[i] ,
                           alleles[i][: , 0] , alleles[i][: , 1])).T

        for gen , gen_haps , gen_suf , ibd in _dc.keys() :
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

            print('unphased')

            hf = h5py.File(args.outprefix + 'unphased_impute_chr_' + str(
                    chroms[i]) + gen_suf + '.hdf5' , 'w')
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

    ##
    print('Write IBD segments')

    _dc = {
            (last_gen , a_off , haps_off , ibd_off , '')    : None ,
            (par_gen , a_par , haps_par , ibd_par , '_par') : None ,
            }

    for gen , a , haps , ibd , gen_suf in _dc.keys() :

        snp_count = 0
        sibpairs = ped[gen , 1]
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
