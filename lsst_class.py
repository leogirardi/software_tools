# -*- coding: utf-8 -*-
"""
Created on Fri Oct 26 21:53:46 2018

@author: rstreet
"""

from astropy import units as u
from astropy.coordinates import SkyCoord

class LSSTFootprint:
    """Class describing the parameters of a single LSST field of view"""
    
    def __init__(self,ra_centre=None,dec_centre=None):
        self.ra = ra_centre
        self.dec = dec_centre
        self.radius = 3.5/2.0  # deg radius field of view
        self.coord = SkyCoord(ra_centre+' '+dec_centre, unit=(u.hourangle, u.deg))
        self.pixel_scale = 0.2 # arcsec
    
    def count_stars_in_footprint(self,coords):
        """Method to count the number of stars within this footprint, 
        given a set of SkyCoords"""
        
        sep = self.coord.separation(coords)  
        
        stars_in_field = np.where( sep.degree < self.radius )[0]
        
        print('Number of stars within the LSST footprint: '+\
                    str(len(stars_in_field)))
