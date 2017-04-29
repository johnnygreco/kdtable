"""
Modified utilities from astroML: 
https://github.com/astroML/astroML
"""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import numpy as np

__all__ = ['ra_dec_to_xyz', 
           'angular_dist_to_euclidean_dist',
           'euclidean_dist_to_angular_dist']


def ra_dec_to_xyz(ra, dec):
    """
    Convert ra & dec to Euclidean points

    Parameters
    ----------
    ra, dec : ndarrays
        RA and DEC in degrees.

    Returns
    -------
    x, y, z : ndarrays
    """
    sin_ra = np.sin(ra*np.pi/180.)
    cos_ra = np.cos(ra*np.pi/180.)

    sin_dec = np.sin(np.pi/2-dec*np.pi/180.)
    cos_dec = np.cos(np.pi/2-dec*np.pi/180.)

    return (cos_ra*sin_dec, sin_ra*sin_dec, cos_dec)


def angular_dist_to_euclidean_dist(theta, r=1):
    """
    Convert angular distances to Euclidean distances.

    Parameters
    ----------
    theta : float or ndarray 
        Angular distance.
    r : float, optional
        Radius. Unit sphere assumed as default.

    Returns
    -------
    d : float or ndarray
        Euclidean distance
    """
    d = 2*r*np.sin(0.5*theta*np.pi/180.)
    return d


def euclidean_dist_to_angular_dist(d, r=1):
    """
    Convert Euclidean distances to angular distances.

    Parameters
    ----------
    d : float or ndarray
        Euclidean distance
    r : float, optional
        Radius. Unit sphere assumed as default.

    Returns
    -------
    theta : float or ndarray
        Angular distance.
    """
    theta = 2*np.arcsin(d/r/2.)*180./np.pi
    return theta
