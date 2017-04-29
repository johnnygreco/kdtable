from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import numpy as np

__all__ = ['ra_dec_to_xyz', 
           'angular_dist_to_euclidean_dist',
           'euclidean_dist_to_angular_dist']

def ra_dec_to_xyz(ra, dec):
    """Convert ra & dec to Euclidean points
    Parameters
    ----------
    ra, dec : ndarrays
    Returns
    x, y, z : ndarrays
    """
    sin_ra = np.sin(ra * np.pi / 180.)
    cos_ra = np.cos(ra * np.pi / 180.)

    sin_dec = np.sin(np.pi / 2 - dec * np.pi / 180.)
    cos_dec = np.cos(np.pi / 2 - dec * np.pi / 180.)

    return (cos_ra * sin_dec,
            sin_ra * sin_dec,
            cos_dec)


def angular_dist_to_euclidean_dist(D, r=1):
    """convert angular distances to euclidean distances"""
    return 2 * r * np.sin(0.5 * D * np.pi / 180.)


def euclidean_dist_to_angular_dist(D, r=1):
    """convert euclidean distances to angular distances"""
    return 2*np.arcsin(D/r/2.)*180./np.pi
