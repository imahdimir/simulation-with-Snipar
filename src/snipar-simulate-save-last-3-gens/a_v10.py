"""

    """

import numpy as np
import pandas as pd

def main() :
    pass

    ##
    fp = '/Users/mahdi/Dropbox/1-Git/snipar-simulate-save-last-3-gens/src/snipar-simulate-save-last-3-gens/sim/causal_effects.txt'
    df = pd.read_csv(fp , delimiter = '\s')

    ##
    df['direct_v10'] = df['direct'] + np.random.normal(0 ,
                                                       np.sqrt(10) * np.std(
                                                               df['direct']) ,
                                                       len(df))

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
