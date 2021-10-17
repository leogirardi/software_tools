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

        self.mjdCol = 'observationStartMJD'

        super(calcVisitIntervalMetric,self).__init__(col=cols, metricName=metricName)


    def run(self, dataSlice, slicePoint=None):
        
        # Array of minimum-required observation intervals for four categories
        # of time variability
        tau_obs = np.array([2.0, 20.0, 73.0, 365.0])
        
        # Calculate the median time interval from the observation 
        # sequence in the dataSlice
        delta_t = np.median(dataSlice[self.mjdCol])
        
        # Decay constant for metric value relationship as a function of 
        # observation interval
        K = 1.0/tau_obs

        metric = 0.0
        idx = np.where(delta_t <= tau_obs)[0]
        if delta_t <= tau_obs:
            metric = 1.0
        else:
            metric = np.exp(-K*(delta_t - tau_obs))
        
        return metric