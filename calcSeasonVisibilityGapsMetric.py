################################################################################################
# Metric to evaluate the calcSeasonVisibilityGapsMetric
#
# Author - Rachel Street: rstreet@lco.global
################################################################################################
import numpy as np
import lsst.sims.maf.metrics as metrics
import lsst.sims.maf.slicers as slicers
import lsst.sims.maf.metricBundles as metricBundles
from lsst.sims.maf.metrics import BaseMetric, seasonMetrics
import healpy as hp
import calcVisitIntervalMetric

class calcSeasonVisibilityGapsMetric(BaseMetric):

    def __init__(self, cols=['fieldRA','observationStartMJD',],
                       metricName='calcVisitIntervalMetric',
                       **kwargs):
        
        """tau_obs is an array of minimum-required observation intervals for four categories
        of time variability"""
        
        self.tau_obs = np.array([2.0, 20.0, 73.0, 365.0])
        self.ra_col = 'fieldRA'
        self.mjdCol = 'observationStartMJD'

        super(calcSeasonVisibilityGapsMetric,self).__init__(col=cols, metricName=metricName)
    
    def calcSeasonGaps(ra, time):
        """Given the RA of a field pointing, and time of observation, calculate the length of 
        the gaps between observing seasons.

        Parameters
        ----------
        ra : float
            The RA (in degrees) of the point on the sky
        time : np.ndarray
            The times of the observations, in MJD
        Returns
        -------
        np.ndarray
            Time gaps in days between sequential observing seasons
        """

        seasons = seasonMetrics.calcSeason(ra, time)
        firstOfSeason, lastOfSeason = seasonMetrics.indSeasonEdges(seasons)
        ngaps = len(firstOfSeason)-1
        season_gaps = data[self.mjdCol][lastOfSeason[0:ngaps-1]] - \
                        data[self.mjdCol[firstOfSeason[1:ngaps]]

        
    def run(self, dataSlice, slicePoint=None):
        season_gaps = calcSeasonGaps(dataSlice[self.ra_col], dataSlice[self.mjdCol])
    
        # To avoid the intensive calculation of the exact visibility of every pointing 
        # for 365d a year, we adopt the pre-calculated values for an example field in 
        # the Galactic Bulge, which receives good, but not circumpolar, annual visibility.  
        total_time_visible_days = 1975.1256 / 24.0
        expected_gap = 365.24 - total_time_visible_days
    
        metric_values = [0.0]*len(self.tau_obs)
        for i,tau in enumerate(self.tau_obs):
            if tau >= expected_gap:
                metric_values[i] = 0.0
                for t in season_gaps:
                    metric_values[i] += calcVisitIntervalMetric(t, tau)
                metric_values[i] /= 10.0

            else:
                metric_values[i] = 1.0

        return metric_values