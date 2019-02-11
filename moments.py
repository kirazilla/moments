#!/usr/local/bin/python3
"""
Compute moments
Written by Monica Rainer
Python3
2019-02-11

Usage:
    moments.py -h
    moments.py
    moments.py <list> 
    moments.py <list> [--range=<values>]
    moments.py <list> <outname>
    moments.py <list> <outname> [--range=<values>]

Options:
    -h,--help           : show this screen
    list                : list of data files (ccf either FITS or txt files)
    outname             : root name of the output file, '_m1.txt', '_m2.txt', 
                          and '_m3.txt' will be added [Default: moments]
    --range=<values>    : borders of the line in km/s, written as --range=x1,x2
"""


#from __future__ import (absolute_import, division, print_function, unicode_literals)

from matplotlib import pyplot as plt
import sys
from docopt import docopt
import numpy as np
from astropy.io import fits


docopt_args = docopt(__doc__)

if not docopt_args['<list>']:
    sys.exit('No input files given')

lista = docopt_args['<list>']

out = docopt_args['<outname>']
if not out:
    out = 'moments'

name0 = '_'.join((out,'m0.txt'))
name1 = '_'.join((out,'m1.txt'))
name2 = '_'.join((out,'m2.txt'))
name3 = '_'.join((out,'m3.txt'))

#names = [line.rstrip('\n') for line in open(lista)]
names = [line.strip() for line in open(lista)]


rvs = []
ccfs = []
jds = []
def read_fits(fname):
    global rvs
    global ccfs
    global jds

    data, hea = fits.getdata(fname,0,header=True)
    try:
        ccfs.append(data[69])
    except:
        print('%s skipped!' % (fname))
        return

    jds.append(hea['HIERARCH TNG DRS BJD'])

    length = hea['NAXIS1']
    start = hea['CRVAL1']
    step = hea['CDELT1']
    stop = start + step*length
    wave = np.arange(start,stop,step)
    #print(len(wave))
    rvs.append(wave)


for name in names:
    try:
        read_fits(name)
    except:
        rv, flux = np.loadtxt(name, usecols=(0,1), unpack=True)
        rvs.append(rv)
        ccfs.append(flux)
        jds.append(0)

rvs = np.asarray(rvs)
ccfs = np.asarray(ccfs)
jds = np.asarray(jds)


def show_ccf(rvs,ccfs):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    for i in range(len(ccfs)):
        ax.plot(rvs[i],ccfs[i]/np.mean(ccfs[i]))

    coords = []

    def onclick(event):
        ix = event.xdata
        print('x = %f'%(ix))
        ax.axvline(x=ix)
        fig.canvas.draw()

        #global coords
        coords.append(ix)

        if len(coords) == 2:
            fig.canvas.mpl_disconnect(cid)
            #plt.close()

            #return coords

    cid = fig.canvas.mpl_connect('button_press_event', onclick)

    print("Click on the figure to define the line region.")

    plt.show(block=False)

    #while len(coords) != 2:
    #    pass


    return coords


if not docopt_args['--range']:
    coords = show_ccf(rvs,ccfs)


    while True:

        if len(coords) == 2:
            print(coords)
        flag = input('Do you accept these line limits (y/n, default y)?\n')
        if flag == 'n':
            plt.clf()
            coords = show_ccf(rvs,ccfs)
        else:
            break


else:
    coords = docopt_args['--range'].split(',')
    print(coords)

# compute moments!
mom0 = []
mom1 = []
mom2 = []
mom3 = []
errmom0 = []
errmom1 = []
errmom2 = []
errmom3 = []

for idx, ccf in enumerate(ccfs):

    # 1. normalize line
    x1 = np.searchsorted(rvs[idx],float(coords[0]))
    x2 = np.searchsorted(rvs[idx],float(coords[1]))

    idxarray = np.zeros(ccf.shape, dtype=bool)
    idxarray[0:x1-1] = True
    idxarray[x2+1:-1] = True

    #linfit = np.poly1d(np.polyfit(rvs[idx][0:x1-1,x2+1:-1],ccf[0:x1-1,x2+1:-1],1))
    linfit = np.poly1d(np.polyfit(rvs[idx][idxarray],ccf[idxarray],1))
    fitvalues = linfit(rvs[idx])
    norccf = 1.0 - (ccf/fitvalues)

    # 2. compute moments (https://www.aanda.org/articles/aa/full/2003/05/aa3122/node2.html)
    #m1 = np.true_divide(np.sum(np.multiply(rvs[idx][x1:x2],norccf[x1:x2])),np.sum(norccf[x1:x2]))
    #m2 = np.true_divide(np.sum(np.multiply(np.power(rvs[idx][x1:x2],2),norccf[x1:x2])),np.sum(norccf[x1:x2]))
    #m3 = np.true_divide(np.sum(np.multiply(np.power(rvs[idx][x1:x2],3),norccf[x1:x2])),np.sum(norccf[x1:x2]))

    ew = np.trapz(norccf[x1:x2]) # equivalent width
    m1 = np.true_divide(np.trapz(np.multiply(rvs[idx][x1:x2],norccf[x1:x2])), ew)
    m2 = np.true_divide( np.trapz( np.multiply( np.power(rvs[idx][x1:x2],2 ), norccf[x1:x2] ) ), ew)
    m3 = np.true_divide(np.trapz(np.multiply(np.power(rvs[idx][x1:x2],3),norccf[x1:x2])), ew)


    # 3. compute statistical uncertainties on moments
    # (http://www.ster.kuleuven.be/~zima/famias/famias_manual/Data_Manager.html#page:sigmammom)

    SNR = np.true_divide(np.mean(norccf[idxarray]),np.std(norccf[idxarray]))
    #print(SNR)

    with np.errstate(divide='ignore', invalid='ignore'):
        sigmaI = np.true_divide(SNR , np.sqrt(norccf[x1:x2]))
    sigmaI = np.nan_to_num(sigmaI)
    #print(sigmaI)


    #deltav02 =  np.sum(np.power(sigmaI,2))
    #deltav12 = np.sum(np.power(np.multiply(rvs[idx][x1:x2],sigmaI),2))
    #deltav22 = np.sum(np.power(np.multiply(np.power(rvs[idx][x1:x2],2),sigmaI),2))
    #deltav32 = np.sum(np.power(np.multiply(np.power(rvs[idx][x1:x2],3),sigmaI),2))

    deltav02 = np.trapz(np.power(sigmaI,2))
    deltav12 = np.trapz(np.power(np.multiply(rvs[idx][x1:x2],sigmaI),2))
    deltav22 = np.trapz(np.power(np.multiply(np.power(rvs[idx][x1:x2],2),sigmaI),2))
    deltav32 = np.trapz(np.power(np.multiply(np.power(rvs[idx][x1:x2],3),sigmaI),2))

    #errm1 = np.sqrt( np.true_divide( deltav12 , np.power( np.sum(norccf[x1:x2]) ,2) ) ) +
    #      \ np.multiply( deltav02 , np.power( np.true_divide 
    #      \ (np.sum( np.multiply( rvs[idx][x1:x2] , norccf[x1:x2] ) ),
    #      \ np.power( np.sum(norccf[x1:x2]) ,2 )  ) ,2 ) )

    errm0 = np.sqrt( np.true_divide( deltav02 , np.power( ew ,2) ) ) + \
            np.power ( deltav02 , 2 )

    errm1 = np.sqrt( np.true_divide( deltav12 , np.power( ew ,2) ) ) + \
            np.multiply( deltav02 , np.power( np.true_divide \
            (np.trapz( np.multiply( rvs[idx][x1:x2] , norccf[x1:x2] ) ), \
            np.power( ew ,2 )  ) ,2 ) )

    errm2 = np.sqrt( np.true_divide( deltav22 , np.power( ew ,2) ) ) + \
            np.multiply( deltav02 , np.power( np.true_divide \
            (np.trapz( np.multiply( np.power( rvs[idx][x1:x2] , 2) , norccf[x1:x2] ) ), \
            np.power( ew ,2 )  ) ,2 ) )

    errm3 = np.sqrt( np.true_divide( deltav32 , np.power( ew ,2) ) ) + \
            np.multiply( deltav02 , np.power( np.true_divide \
            (np.trapz( np.multiply( np.power( rvs[idx][x1:x2] , 3) , norccf[x1:x2] ) ), \
            np.power( ew ,2 )  ) ,2 ) )

    # 4. save moments and errors

    mom0.append(ew)
    mom1.append(m1)
    mom2.append(m2)
    mom3.append(m3)

    errmom0.append(errm0)
    errmom1.append(errm1)
    errmom2.append(errm2)
    errmom3.append(errm3)

mom0 = np.asarray(mom0)
mom1 = np.asarray(mom1)
mom2 = np.asarray(mom2)
mom3 = np.asarray(mom3)

errmom0 = np.asarray(errmom0)
errmom1 = np.asarray(errmom1)
errmom2 = np.asarray(errmom2)
errmom3 = np.asarray(errmom3)

out0 = np.vstack( ( np.flipud(jds), np.flipud(mom0), np.flipud(errmom0) ) )
out1 = np.vstack( ( np.flipud(jds), np.flipud(mom1), np.flipud(errmom1) ) )
out2 = np.vstack( ( np.flipud(jds), np.flipud(mom2), np.flipud(errmom2) ) )
out3 = np.vstack( ( np.flipud(jds), np.flipud(mom3), np.flipud(errmom3) ) )

np.savetxt(name0, np.transpose(out0), fmt='%.8f', header='# 1. BJD - 2. zero moment (EW) - 3. statistical uncertanties\n')
np.savetxt(name1, np.transpose(out1), fmt='%.8f', header='# 1. BJD - 2. first moment (RV) - 3. statistical uncertanties\n')
np.savetxt(name2, np.transpose(out2), fmt='%.8f', header='# 1. BJD - 2. second moment (variance) - 3. statistical uncertanties\n')
np.savetxt(name3, np.transpose(out3), fmt='%.8f', header='# 1. BJD - 2. third moment (skewness) - 3. statistical uncertanties\n')



