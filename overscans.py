import os
import sys

import csv
import argparse

if  os.environ.get('EOTEST_DIR') is None:# replace with your eotest directory
    os.environ['EOTEST_DIR'] = '/Users/duncan/lib/eotest/python' 
from lsst.afw.image import ImageF
import lsst.eotest.image_utils as imutils
from lsst.eotest.sensor.AmplifierGeometry import makeAmplifierGeometry

import matplotlib.pyplot as plt
from matplotlib import gridspec
import numpy as np
from scipy.optimize import curve_fit
from astropy.io.fits import getheader

def exp_fit(x, c, A, B):
    return A * np.exp(-x/c) + B 

def mean_overscans(imagefile, nexp=7, plots=None, verbose=False):
    result = [] # decay constant, amplitude, flux of image, residue of first overscan, sensor number, segment, run number
    for segment in range(1,16+1):
        header = getheader(imagefile)
        info = [header['LSST_NUM'],'Segment {}'.format(segment),header['RUNNUM'].split()[0]]
        image_untrimmed = ImageF(imagefile,segment)
        amp = makeAmplifierGeometry(imagefile)
        
        flat_mean = np.mean(imutils.trim(image_untrimmed, imaging=amp.imaging).array[:,-1]) 
        # mean value of the last column of the image
        
        image = imutils.trim(image_untrimmed, imaging=amp.serial_overscan)
        overscans_mean = [np.mean([image.array[i,j] for i in np.arange(len(image.array))]) for j in range(len(image.array[0]))]
        bias = np.mean(overscans_mean[5:-2])
        over_mean_subfit = overscans_mean[:nexp]
        params, cov = curve_fit(exp_fit, np.arange(nexp), over_mean_subfit, p0=(10, 10, 20000),bounds=([.1,0,0],[20,300,50000]))
        residue = params[1]/(flat_mean-bias)
        
        result.append([params[0],params[1],flat_mean,residue, *info])
        
        if verbose:
            print(params)
        if plots is None:
            continue
        fig = plt.figure(figsize=(10,10))
        plt.plot(over_mean_subfit, ls='none', marker='.')
        xfit = np.linspace(0,nexp-1, 50)
        plt.plot(xfit,[exp_fit(x, *params) for x in xfit])
        plt.title('Superflat Mean Serial Overscan in {0} {1} Run {2}'.format(*info))
        plt.figtext(0.5,0.5,('Decay constant: {p[0]:.03g} pixels \nAmplitude: {p[1]:.03g}'\
                    ' ADU\nImage flux: {0:.00f} ADU\nResidue in first overscan pixel: {1:.03%}').format(flat_mean,residue,p=params))
        plt.ylabel('ADU')
        plt.xlabel('Pixel')
        plt.legend(['Data','Fit'])
        fig.patch.set_facecolor('white')
        plt.savefig(('{0}/{1}_{2}_run{3}.png'.format(plots, *info)).replace(" ",""))
        plt.close(fig)

    return result

def main():

    parser = argparse.ArgumentParser(description='Retrieve overscan information from flat test images.')
    parser.add_argument('indir', type=str, nargs=1, help='Directory with input files')
    parser.add_argument('out', type=str, nargs=1, help='Output file')
    parser.add_argument('--npixels', type=int, nargs=1, help='Number of overscan pixels to fit with an exponential.')
    parser.add_argument('--plots', type=str, nargs='?',const='plots', help='Save plots of overscan fits.')

    args = parser.parse_args()

    indir = args.indir[0]
    outfile = args.out[0]
    npixels = args.npixels[0]
    makeplots = args.plots is not None
    if makeplots:
        plotsdir = args.plots
        os.makedirs(plotsdir, exist_ok=True)
    else:
        plotsdir = None
    flats_directory = os.fsencode(indir)
    data = [['decay constant (pixels)', 'amplitude (ADU)', 'flux of image (ADU)', 'residue of first overscan (ratio)', 'sensor number', 'segment', 'run number']]

    out = csv.writer(open(outfile,'w'))

    for root, dir, f in os.walk(flats_directory):
        for file in f:
            filename = os.path.join(root, file).decode('UTF-8');
            if filename.endswith('.fits'):
                print('Checking file: ' + filename)
                data = data + mean_overscans(filename, plots=plotsdir, nexp=npixels)
    out.writerows(data)
    return

if __name__ == "__main__": main()

