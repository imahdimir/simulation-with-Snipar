"""

    """

import numpy as np
import pandas as pd


def main():
    pass

    ##
    fp = '/Users/mahdi/Dropbox/1-Git/snipar-simulate-save-last-3-gens/src/snipar-simulate-save-last-3-gens/sim/causal_effects.txt'
    df = pd.read_csv(fp, delimiter='\s')

    ##
    df['direct_v10'] = df['direct'] + np.random.normal(0, np.sqrt(10) * np.std(df['direct']), len(df))

    ##
    ce = df.to_numpy()

    ##
    _cols = ['SNP' , 'A1' , 'A2' , 'direct' , 'direct_v1' , 'direct_v10']
    ce = np.vstack([_cols, ce])

    ##
    np.savetxt('ce.txt', ce, fmt = '%s')

    ##
