from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import os, copy
import numpy as np
from astropy.table import Table
from astropy.coordinates import SkyCoord
from sklearn.neighbors import KDTree
from .utils import ra_dec_to_xyz, angular_dist_to_euclidean_dist
from .utils import euclidean_dist_to_angular_dist

__all__ = ['KDTable']

class KDTable(object):

    def __init__(self, fn=None, format='csv', ra_col='ra', 
                 dec_col='dec', data=None):
        self.data = data if data else Table.read(fn, format=format)
        self._kdt = None
        self._skycoord = None
        self.ra_col = ra_col
        self.dec_col = dec_col
    
    @property
    def kdt(self):
        if self._kdt is None:
            xyz = ra_dec_to_xyz(
                self.data[self.ra_col], self.data[self.dec_col])
            xyz = np.asarray(xyz).T
            self._kdt = KDTree(xyz)
        return self._kdt

    @property
    def colnames(self):
        return self.data.colnames

    @property
    def skycoord(self):
        if self._skycoord is None:
            self._skycoord = SkyCoord(
                self.data[self.ra_col], self.data[self.dec_col], unit='deg')
        return self._skycoord

    def query_radius(self, coords, r, count_only=False, return_seps=True):
        """
        Search for sources around coords within circle of 
        radius r in degrees. 
        """
        coords = np.asarray(coords)
        if coords.shape == (2,):
            ra, dec = coords
            xyz = np.array(ra_dec_to_xyz(ra, dec)).T.reshape(1, -1)
        else:
            ra, dec = coords[:, 0], coords[:, 1]
            xyz = np.array(ra_dec_to_xyz(ra, dec)).T
        if count_only:
            kws = dict(count_only=True)
            results = self.kdt.query_radius(
                xyz, angular_dist_to_euclidean_dist(r), **kws)
        else: 
            kws = dict(count_only=False, return_distance=True)
            idx, dist = self.kdt.query_radius(
                xyz, angular_dist_to_euclidean_dist(r), **kws)
            if len(idx)==1:
                idx = idx[0]
                dist = euclidean_dist_to_angular_dist(dist[0])
            else:
                if return_seps:
                    for i in range(len(dist)):
                        dist[i] = euclidean_dist_to_angular_dist(dist[i])
            results = (idx, dist) if return_seps else idx
        return results

    def make_cuts(self, cuts):
        self.data = self.data[cuts]
        self._kdt = None
        self._skycoord = None
    
    def copy(self):
        return copy.deepcopy(self)

    def add_col(self, name, arr):
        self.data[name] = arr

    def __getitem__(self, cols):
        return self.data[cols]

    def __len__(self):
        return len(self.data)
