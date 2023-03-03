import numpy as np
from itertools import combinations
import itertools


def flatten(list2d):
    return list(itertools.chain.from_iterable(list2d))

def jitter_point(mean,std=0.15):
    return np.random.normal(mean,std)

def inverse_variance_mean(means,variances,axis=1):

    # variances = standard_devs**2
    # variances = standard_devs

    weighted_mean = np.nansum(means/variances,axis=axis)/(np.nansum(1/variances,axis=axis))

    weighted_variance = (np.nansum(1/variances,axis=axis))**(-1)

    return weighted_mean, weighted_variance

def standard_error(data):

    return np.std(data)/np.sqrt(len(data))