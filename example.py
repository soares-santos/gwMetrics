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


dirdb = '/data/des40.a/data/marcelle/lsst-gw/OperationsSimulatorBenchmarkSurveys/'
file = 'kraken_1042_sqlite.db'
outdir = './example-outputs/'

mjds=np.array([1.,8.,15.])
ras=np.random.rand(3)*180.
decs=np.random.rand(3)*(-30.)
print ras,decs
maps=make_maps.gaussian2d_from_sample_map(ras,decs,sigma_ra=5.,sigma_dec=5.)

sum_maps = sum(maps[i] for i in range(len(maps)))

plt.clf()
hp.mollview(maps[0])
plt.savefig(outdir+'LIGO_test_map.png')
plt.clf()
hp.mollview(sum_maps)
plt.savefig(outdir+'LIGO_MAPS_sum.png')
#hp.write_map(outdir+'LIGO_test_map.fits',maps[0])


opsdb = db.OpsimDatabase(dirdb+file)
resultsDb = db.ResultsDb(outDir=outdir)
print "Starting to compute AreaProb..."
metric = metrics.AreaProb()
metric.setEvents(maps,mjds)
print "Calling slicer..."
slicer = slicers.HealpixSlicer(nside=4)
sql0 = 'filter = "i" and night < 50'
cmap = plt.get_cmap('winter')
plotDict0={'title': 'u-band. error: 0.020, z: 2.1. 1 year. opsim: kraken_1042.',  'logScale': True, 'xlabel': 'Bayes Factor','cmap': cmap,'percentileClip': None, 'colorMin':1.0, 'colorMax':100.0, 'cbarFormat':'%.1g','cbar_edge': True, 'nTicks': 10, 'aspect': 'auto','xextent': None, 'origin': None}
HealpixSkyMap = plots.HealpixSkyMap()
HealpixHistogram = plots.HealpixHistogram()
mb0 = metricBundles.MetricBundle(metric, slicer, sql0,  plotDict=plotDict0, plotFuncs=[HealpixSkyMap, HealpixHistogram])
mbD = {0:mb0}
bgroup = metricBundles.MetricBundleGroup(mbD, opsdb, outDir=outdir, resultsDb=resultsDb)
bgroup.runAll()


