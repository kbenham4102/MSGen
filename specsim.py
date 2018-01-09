import numpy as np
import scipy.stats as stats

def gauss(x,mz,sigma,A):
    return A*np.exp(-(x-mz)**2/2/sigma**2)


def sig(respwr, mz):
        return mz/(2.355*respwr)


def peaklist(filename):
    arr = []
    pknums = []
    f = open(filename)
    for line in f:
        if 'Num peaks:' in line:
            x = line.split(' ')
            pknums.append(int(x.pop(-1)))
        if ('Name:' not in line and 'MW:' not in line and 'Comment:' not in line
         and 'Num' not in line and line != '\n'):
            arr1 = line.split('\t')
            arr2 = arr1.pop(-1)
            arr.append(arr1)
    arr = np.array(arr).astype(np.float)
    f.close()
    return arr

# specsim.genMS(mz, win, pkvar = 0.2, noiselevel = 5e-3, bsline = False,
# baselevel = 0.05, rp = 4000, points = 30000)
# Generate mass spectrum with dict input mz
# Inputs:
# mz - type: dict, format {'mz1':[m/z, int],.....} where m/z is the peak center
#                      and int is the signal height on a scale [0,1]
# win - spectral window on which to create spectrum, m/z 100 - m/z 1000 recomm.
#
# pkvar - scale of random samples used to scale positive peak heights, 0.2
# typical
#
# noiselevel - height of electronics noise
#
# bsline - if True, raised exponential baseline is added at beginning of
#           spectrum, default is False.
# baselevel - coefficient to determine prominence of raised baseline
#             realistic values ~ 0.05 - 0.2
# misalign - option to add slight random shift to m/z access, simulates ppm
#            error for peak center values
# rp - resolving power of the spectrum, recommended value = 4,000
#
# points - number of data points to use, recommended = 30,000
# RETURNS:
# sigsim - vector of simulated intensity
# x - vector of m/z values over which sigsim is calculated


def genMS(mz, win, pkvar = 0.2, noiselevel = 5e-3, bsline = False,
baselevel = 0.05, misalign = False, rp = 4000, points = 30000):


    mz0 = win[0]
    mzf = win[1]

    x = np.linspace(mz0,mzf, points)
    sigsim = np.zeros(points)

    # Add random shift to m/z axis
    if misalign == True:
        # Generate spectrum with Gaussian peaks and small random peak shift
        for key in mz:
            sigsim += gauss(x, mz[key][0], sig(rp,mz[key][0]+
                        np.random.normal(0,0.005)), mz[key][1])
    else:
        # Generate spectrum with Gaussian peaks
        for key in mz:
            sigsim += gauss(x, mz[key][0], sig(rp,mz[key][0]), mz[key][1])

    # Add baseline rise at low mass end
    if bsline == True:
        baseline = baselevel*np.exp(-1*np.linspace(0,5,points))
        sigsim += baseline

    # Add peak variance
    for i in range(points):
        if sigsim[i] > 0:
            sigsim[i] *= abs(1 + np.random.normal(0,pkvar))

    # Add electronic baseline noise
    noise = np.random.normal(0,1,points)*noiselevel
    sigsim += noise

    return x, sigsim

def file_sz(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    k = l.split(' ')
    k.pop(-1)
    return i + 1, len(k)

def importgendat(fileloc):
    points, n = file_sz(fileloc)
    arr = np.empty((points, n))
    f = open(fileloc)
    i = 0
    for line in f:
        arr1 = line.split(' ')
        b = arr1.pop(-1)
        arr2 = np.array(arr1).astype(np.float)
        arr[i,:] = arr2
        i += 1
    mz = arr[:,0]
    signal = arr[:,1:]
    return mz, signal
