"""

    """

import numpy as np
import pandas as pd

def add_causal_effect_cols() :
    pass

    ##
    fp = '/Users/mahdi/Dropbox/1-Git/simulation-with-Snipar/sim_output_wt_indir/causal_effects.txt'
    df = pd.read_csv(fp , delimiter = '\s' , header = None)

    ##
    df.columns = ['SNP' , 'A1' , 'A2' , 'direct' , 'indirect']

    ##
    df.to_csv(fp , index = False , sep = ' ')

def add_sum_column_and_noisy_columns() :
    pass

    ##
    fp = '/Users/mahdi/Dropbox/1-Git/simulation-with-Snipar/sim_output_wt_indir/causal_effects.txt'
    df = pd.read_csv(fp , delimiter = '\s')

    ##
    df['sum'] = df['direct'] + df['indirect']

    ##
    df['sum_v1'] = df['sum'] + np.random.normal(0 ,
                                                1 * np.std(df['sum']) ,
                                                len(df))

    ##
    df['sum_v10'] = df['sum'] + np.random.normal(0 ,
                                                 np.sqrt(10) * np.std(
                                                         df['sum']) ,
                                                 len(df))

    ##
    df.to_csv(fp , index = False , sep = ' ')

    ##

def nex() :
    pass
    ##

    ##
    ce = df.to_numpy()

    ##
    _cols = ['SNP' , 'A1' , 'A2' , 'direct' , 'direct_v1' , 'direct_v10']
    ce = np.vstack([_cols , ce])

    ##
    np.savetxt('ce.txt' , ce , fmt = '%s')

    ##
    fp = '/Users/mahdi/Dropbox/1-Git/snipar-simulate-save-last-3-gens/src/snipar-simulate-save-last-3-gens/sim/table1.txt'

    df = pd.read_csv(fp , delim_whitespace = True , index_col = 0)

    ##
    df = df.round(4)

    ##
    df.to_csv(
            '/Users/mahdi/Dropbox/1-Git/snipar-simulate-save-last-3-gens/src/snipar-simulate-save-last-3-gens/sim/table1_round4.txt')

    ##
    import h5py

    fp = '/Users/mahdi/Dropbox/1-Git/snipar-simulate-save-last-3-gens/src/snipar-simulate-save-last-3-gens/sim1/phased_impute_chr_1.hdf5'

    hf = h5py.File(fp , 'r')

    ##
    hf.keys()

    ##
    import pandas as pd

    df = pd.DataFrame(hf['pedigree'])

    ##
    df = df.applymap(lambda x : x.decode())

    ##
    df1 = df.drop_duplicates(subset = [0 , 1])

    ##

    ##
