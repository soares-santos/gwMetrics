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
dirdb = '/data/des40.a/data/marcelle/lsst-gw/OperationsSimulatorBenchmarkSurveys/'
file = 'kraken_1042_sqlite.db'
opsdb = db.OpsimDatabase(dirdb+file)
outdir = './example-outputs/'
resultsDb = db.ResultsDb(outDir=outdir)
metric = metrics.AstrometryMetric()
slicer = slicers.HealpixSlicer(nside=4)
sql0 = 'filter = "u" and night < 100'
cmap = plt.get_cmap('winter')
plotDict0={'title': 'u-band. error: 0.020, z: 2.1. 1 year. opsim: kraken_1042.',  'logScale': True, 'xlabel': 'Bayes Factor','cmap': cmap,'percentileClip': None, 'colorMin':1.0, 'colorMax':100.0, 'cbarFormat':'%.1g','cbar_edge': True, 'nTicks': 10, 'aspect': 'auto','xextent': None, 'origin': None}
HealpixSkyMap = plots.HealpixSkyMap()
HealpixHistogram = plots.HealpixHistogram()
mb0 = metricBundles.MetricBundle(metric, slicer, sql0,  plotDict=plotDict0, plotFuncs=[HealpixSkyMap, HealpixHistogram])
mbD = {0:mb0}
bgroup = metricBundles.MetricBundleGroup(mbD, opsdb, outDir=outdir, resultsDb=resultsDb)
bgroup.runAll()


