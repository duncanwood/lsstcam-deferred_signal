import argparse
import numpy as np
import os
if  os.environ.get('EOTEST_DIR') is None:# replace with your eotest directory
    os.environ['EOTEST_DIR'] = '/Users/duncan/lib/eotest/python'
from lsst.eotest.sensor.AmplifierGeometry import makeAmplifierGeometry
from lsst.afw.image import ImageF
import lsst.eotest.image_utils as imutils
from lsst.eotest.sensor.eperTask import *
import matplotlib.pyplot as plt
from astropy.io.fits import getheader

def flats(imagefile, verbose=False):
    result = {}
    ccd = makeAmplifierGeometry(imagefile)
    for segment in range(1,16+1):
        image_untrimmed = ImageF(imagefile,segment)#imutils.dm_hdu(segment))
        amp = ccd[segment]
        flatmean = imutils.mean(imutils.unbias_and_trim \
                (image_untrimmed, overscan=ccd.serial_overscan, imaging=ccd.imaging))
        result[segment] = flatmean
    return result

def sensor_CTI(flatname, verbose):
    result = []
    amps = imutils.allAmps(flatname)
    s_task = EPERTask()
    s_task.config.direction = 's'
    s_task.config.verbose = verbose
    s_task.config.cti = True
    scti, bias_ests = s_task.run(flatname, 2, amps, 2)

    return scti

def cti_vs_flux(indir, verbose=False):
    sctis = []
    fluxes = []
    flatdir = os.fsencode(indir)
    sensorname = None
    for root, d, f in os.walk(flatdir):
        for file in f:
            filename = os.path.join(root, file).decode('UTF-8');
            if filename.endswith('.fits'):
                if verbose: 
                    print('Checking file: ' + filename)
                header = getheader(filename)
                if not (header['IMGTYPE'].replace(" ","") == 'FLAT'): continue
                if sensorname is None:
                    sensorname = header['LSST_NUM']
                else:
                    if not (sensorname == header['LSST_NUM']):
                        continue
                sctis.append([cti.value for cti in sensor_CTI(filename, verbose).values()])
                fluxes.append(list(flats(filename).values()))
    return np.array(sctis), np.array(fluxes), sensorname

def plot_cti_vs_flux(sctis, fluxes, plotdir,sensor=""):
    import itertools
    order = np.argsort(fluxes, axis=0)
#    data = np.array([fluxes, sctis])
#    sdata = np.sort(data,axis=1)
    fig = plt.figure(figsize=(10,10))
    fig.patch.set_facecolor('white')

    x = np.take_along_axis(fluxes, order, axis=0)
    y = np.take_along_axis(sctis, order, axis=0)
    plt.plot(x,y,marker='.')
    plt.title('sCTI vs Flux in Flats {}'.format(sensor))
    plt.xlabel('Flux (ADU)')
    plt.ylabel('Serial CTI (ratio)')
    plt.legend(['Segment {}'.format(i) for i in range(1,16+1)])
    plt.savefig(plotdir)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description='Retrieve overscan information from flat test images.')
    parser.add_argument('indir', type=str, nargs=1, help='Directory with input files')
    #parser.add_argument('out', type=str, nargs=1, help='Output file')
    parser.add_argument('--plots', type=str, nargs='?',const='plots', help='Save plots of overscan fits.')
    parser.add_argument('--verbose',action="store_true", default=False, help='Print runtime information to stdout.')

    args = parser.parse_args()

    indir = args.indir[0]
    #outfile = args.out[0]
    verbose = args.verbose
    makeplots = args.plots is not None
    if makeplots:
        plotsdir = args.plots
    #    os.makedirs(plotsdir, exist_ok=True)
    else:
        plotsdir = None

    sctis, fluxes,sensorname = cti_vs_flux(indir, verbose)

    if plotsdir is not None:
        plot_cti_vs_flux(sctis, fluxes, plotsdir, sensorname)
    return 0 

if __name__ == "__main__": main()
