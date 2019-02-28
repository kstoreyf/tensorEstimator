import pandas as pd
import numpy as np
from astropy.cosmology import LambdaCDM


def main():
    print "Loading data..."
    #sample = 'Full'
    #sample = 'Bright-no'
    sample = 'Dim-no'
    #fn = '../data/DR7-{}.ascii'.format(sample)
    fn = '../data/random-DR7-{}.ascii'.format(sample)

    df = pd.read_csv(fn, index_col=0)
    #cosmo = LambdaCDM(H0=70, Om0=0.25,Ode0=0.75)
    cosmo = LambdaCDM(H0=70, Om0=0.30,Ode0=0.70)

    write_comoving_dist(df, fn, cosmo)

def write_comoving_dist(df, fn, cosmo):

    dcm_col = 'dcm_Om0-{:.2f}'.format(cosmo.Om0)
    if dcm_col not in df.columns:
        print "Adding colum {}...".format(dcm_col)
        df[dcm_col] = df['z'].apply(get_comoving_dist, args=(cosmo,))
        df.to_csv(fn)
        print "Added dcm column {}, saved new dataframe to {}".format(dcm_col, fn)
        return df
    else:
        print "Dcm column {} already in dataframe {}, returning original".format(dcm_col, fn)
        return df



def get_comoving_dist(z, cosmo):
    comov = cosmo.comoving_distance(z)
    return comov.value*cosmo.h




if __name__=="__main__":
    main()