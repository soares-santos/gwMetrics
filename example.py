import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import unittest

import lsst.sims.maf
import lsst.sims.maf.db as db
import lsst.sims.maf.metrics as metrics
import lsst.sims.maf.slicers as slicers
import lsst.sims.maf.stackers as stackers
import lsst.sims.maf.plots as plots
import lsst.sims.maf.metricBundles as metricBundles
import lsst.sims.maf.utils as utils

import numpy as np
import healpy as hp
import defs as astr
import metrics

import make_maps
import os

dirdb = '/data/des40.a/data/marcelle/lsst-gw/OperationsSimulatorBenchmarkSurveys/'
file = 'kraken_1042_sqlite.db'
outdir = './example-outputs/'
if not os.path.exists(outdir):
    os.makedirs(outdir)

print "Generating gw event sample..."
# generate one simulated gw event
##mjds=np.array([59808.])
##ras=np.array([np.rad2deg(0.15)])
##decs=np.array([np.rad2deg(-0.15)])
##maps,mjds = make_maps.event_sample(ras=ras,decs=decs,mjds=mjds)
# generate simulated gw event sample over 1 year
maps,mjds = make_maps.event_sample(nevents=50)

#sum_maps = sum(maps[i] for i in range(len(maps)))
#plt.clf()
#hp.mollview(maps[0])
#plt.savefig(outdir+'LIGO_test_map.png')
#plt.clf()
#hp.mollview(sum_maps)
#plt.savefig(outdir+'LIGO_MAPS_sum.png')

# evaluate opsim realization of LSST survey
opsdb = db.OpsimDatabase(dirdb+file)
resultsDb = db.ResultsDb(outDir=outdir)
print "Starting to compute AreaProb..."
metric = metrics.AreaProb()
metric.setEvents(maps,mjds)
print "Calling slicer..."
##slicer = slicers.HealpixSlicer(nside=4)
slicer = slicers.UniSlicer()
sql0 = 'filter = "i" and night < 365'
cmap = plt.get_cmap('winter')
plotDict0={'title': 'u-band. error: 0.020, z: 2.1. 1 year. opsim: kraken_1042.',  'logScale': True, 'xlabel': 'Bayes Factor','cmap': cmap,'percentileClip': None, 'colorMin':1.0, 'colorMax':100.0, 'cbarFormat':'%.1g','cbar_edge': True, 'nTicks': 10, 'aspect': 'auto','xextent': None, 'origin': None}
HealpixSkyMap = plots.HealpixSkyMap()
HealpixHistogram = plots.HealpixHistogram()
mb0 = metricBundles.MetricBundle(metric, slicer, sql0,  plotDict=plotDict0, plotFuncs=[HealpixSkyMap, HealpixHistogram])
mbD = {0:mb0}
bgroup = metricBundles.MetricBundleGroup(mbD, opsdb, outDir=outdir, resultsDb=resultsDb)
bgroup.runAll()


