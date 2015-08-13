################################################################################
#									       #
#  Program to calculate the Conditional Entropy for a single pulsating star    # 
#	This program calculates and saves a file with periods and entropies    #
#	To run: type in terminal -> python3 CE.py			       #
#									       #
#  To  change the file: change the .dat file in the "#Load the data" section   #
#  To change the period range: change the values in the                        #
#                              "#Creates the period array" section	       # 
#  The standard precision is 0.0001 but for some stars we need a finner        #
#   precision, 0.00001 is usualy enough.  *precision is the period step!       #
#									       #
################################################################################

import numpy
from periodogram import find_period, rephase, get_phase
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def hc(bins, row, col, size):
    return numpy.sum((lambda p: p*numpy.log(numpy.sum(bins[i,:])/size/p) \
                                if p > 0 else 0)(bins[i][j] / size)
                     for i in row for j in col) if size > 0 else numpy.PINF

#Load the data
data = numpy.ma.array(data=numpy.loadtxt('/path/to/star/file'),
                      mask=None, dtype=float)

#Creates the period array
periods = numpy.arange(0.1, 1.00001, 0.00001) #period range (p_min, p_max, step)

#Set the number of rows and columns to calculate entropy
xbins = 10  	#keep this fix
ybins = 5	#keep this fix
row = numpy.linspace(0, xbins-1, xbins)
col = numpy.linspace(0, ybins-1, ybins)

#For loop to calculate entropy
entropy = []
for p in periods:
    r = rephase(data, p, 0)
    r.T[1] = (r.T[1]-numpy.min(r.T[1]))/(numpy.max(r.T[1])-numpy.min(r.T[1]))
    bins, binsX, binsY = numpy.histogram2d(r.T[0], r.T[1], [xbins, ybins],
                                           [[0,1], [0,1]])
    ent = hc(bins, row, col, len(r.T[1]))
    #print(p, ent)
    entropy.append(ent)

#Save the period and entropy into a file
#numpy.savetxt('period_entropy.dat',
#              numpy.dstack((periods, entropy))[0],
#              fmt='%s')

#Find the position of the minimum entropy to get the correspondent period
minimum = numpy.argmin(entropy)
e_min = entropy[minimum]
p_min = periods[minimum]
print(p_min)
#Print the minimum entropy and the correspondent period
#print('\n', p_min, e_min)

'''
#plot the entropy against periods
fig = plt.figure()
plt.plot(periods, entropy, 'r')
plt.plot(periods, entropy, 'k+', markersize=12)
fig.suptitle('OGLE-LMC-CEP-0010')
plt.xlabel('Periods')
plt.ylabel('Entropy')
#fig.savefig('0010_test11.png')
'''

