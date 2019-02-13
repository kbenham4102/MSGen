import numpy as np
import scipy.stats as stats
import matplotlib.pyplot as plt
import specsim as s

# Create array of mz's and intensities in data file
fname = 'yeast1.txt'
arr = s.peaklist(fname)

# Determine boundaries of x axis
x0 = np.amin(arr[:,0])
xn = np.amax(arr[:,0])
win = [x0,xn]
# Note: max and min mz values must be within core set of protiens between
# populations.

# Normalize database intensities
arr[:,1] /= np.amax(arr[:,1])

pknames = np.arange(len(arr[:,0]))
mz = dict(zip(pknames, arr))


respwr = 50000
numsamp = 40
Enoise = 5e-4
pkv = 0.10
bsl = False
misal = False
pts = 30000
blev = 0.05

datamat = np.empty((30000, numsamp + 1))
for i in range(numsamp):
    x, datamat[:,i+1] = s.genMS(mz, win, rp = respwr, noiselevel = Enoise, pkvar = pkv,
                bsline = bsl, misalign = misal, points = pts, baselevel = blev)


datamat[:,0] = x



directory = 'C:/Users/Kevin/scripts/genyeastdat/'
filename = fname.split('.')[0] + '_' + str(int(pkv*100)) + 'p_' + str(numsamp) + 'samp.txt'

print("Saving...")
file = open(directory + 'HEADER_' + filename , 'w')
file.write("M/Z Window = %d m/z - %d m/z \n" % (x0, xn))
file.write("Number of data points = %d\n" %pts)
file.write("Number of spectra generated = %d\n" % numsamp)
file.write("Resolving Power = %d\n" % respwr)
file.write("Noiselevel = %1.3e\n" % Enoise)
file.write("Peak Variance Scale = %d\n" % pkv)
file.write("Added baseline asymmetry: %d\n" % bsl)
if bsl == True:
    file.write("Prefactor for baseline asymmetry: %d\n" % blev)
file.write("M/Z misalignemnt variation: %d\n" % misal)
file.close()


file = open(directory + filename, 'w')
#file.write(filename.partition('.')[0] + "\n")

for i in range(len(datamat[:,0])):
    file.writelines('%+1.5e ' % k for k in datamat[i,:])
    file.write('\n')
file.close()
