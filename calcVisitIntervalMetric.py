################################################################################################
# Metric to evaluate the calcVisitIntervalMetric
#
# Author - Rachel Street: rstreet@lco.global
################################################################################################
import numpy as np
import lsst.sims.maf.metrics as metrics
import lsst.sims.maf.slicers as slicers
import lsst.sims.maf.metricBundles as metricBundles
from lsst.sims.maf.metrics import BaseMetric
import healpy as hp

class calcVisitIntervalMetric(BaseMetric):

    def __init__(self, cols=['observationStartMJD',],
                       metricName='calcVisitIntervalMetric',
                       **kwargs):
        """tau_obs is an array of minimum-required observation intervals for four categories
        of time variability"""
        
        self.mjdCol = 'observationStartMJD'
        self.tau_obs = np.array([2.0, 20.0, 73.0, 365.0])

        super(calcVisitIntervalMetric,self).__init__(col=cols, metricName=metricName)


    def run(self, dataSlice, slicePoint=None):
        
        # Calculate the median time interval from the observation 
        # sequence in the dataSlice
        delta_t = np.median(dataSlice[self.mjdCol])
        
        # Decay constant for metric value relationship as a function of 
        # observation interval
        K = 1.0/self.tau_obs

        metric_values = [0.0]*len(self.tau_obs)
        for i,tau in enumerate(self.tau_obs):
            idx = np.where(delta_t <= tau)[0]
            if delta_t <= tau:
                metric_values[i] = 1.0
            else:
                metric_values[i] = np.exp(-K*(delta_t - tau))
        
        return metric_values