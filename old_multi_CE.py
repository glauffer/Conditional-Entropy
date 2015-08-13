################################################################################
#									       #
#  Program to calculate the Conditional Entropy for all the files in a 
#    directory  	                                                       # 
#	This program calculates and saves a file with file name and its period #
#	To run: type in terminal -> python3 multi_CE.py			       #
#									       #
#  To  change the files directory: change the path in the 
#                                  "#Folder Location" section                  #
#  To  change the results directory: change the path in the 
#                                  "#Results Location" section                 #
#  To change the period range: change the values in the                        #
#                              "#Creates the period array" section	       # 
#  The standard precision is 0.0001 but for some stars we need a finner        #
#   precision, 0.00001 is usualy enough.  *precision is the period step!       #
#									       #
################################################################################

import os
import numpy
from periodogram import rephase, get_phase
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def hc(bins, row, col, size):	#changed here, original: bins[:,j]
    return numpy.sum((lambda p: p*numpy.log(numpy.sum(bins[i,:])/size/p) \
                                if p > 0 else 0)(bins[i][j] / size)
                     for i in row for j in col) if size > 0 else numpy.PINF

#Folder Location
path = '/path/to/star/data'

os.chdir(path)

#Results Location
out = '/output/path/'
files = os.listdir(path)
files.sort()

ent_data = []
ce_period = []

#Loop to calculate entropy for all the files in the directory
for star in files:
	data = numpy.ma.array(data=numpy.loadtxt(star), mask=None, dtype=float)
	periods = numpy.arange(0.1, 32.0001, 0.0001)
	xbins = 10
	ybins = 5
	row = numpy.linspace(0, xbins-1, xbins)
	col = numpy.linspace(0, ybins-1, ybins)

	entropy = []
	
	for p in periods:
		r = rephase(data, p, 0)
		r.T[1] = (r.T[1]-numpy.min(r.T[1]))/(numpy.max(r.T[1])-numpy.min(r.T[1]))
		bins, binsX, binsY = numpy.histogram2d(r.T[0], r.T[1], [xbins, ybins],[[0,1], [0,1]])
		ent = hc(bins, row, col, len(r.T[1]))
		#print(star, p, ent)
		entropy.append(ent)

	name_of_star = star + ".txt"
	#name_of_plot = star + ".png"

	numpy.savetxt(os.path.join(out, name_of_star),
			numpy.dstack((periods, entropy))[0],
			fmt='%s')

#Find the position of the minimum entropy to get the correspondent period
	minimum = numpy.argmin(entropy)
	e_min = entropy[minimum]
	p_min = periods[minimum]
	line = star, p_min, e_min
	ce_per = star, p_min
	ent_data.append(line)
	ce_period.append(ce_per)
	#numpy.savetxt(os.path.join(out, 'CE_RRAB_finner_period_remaining_stars.dat'), ce_period, fmt='%s')

#Print the minimum entropy and the correspondent period
	#print('\n', star, p_min, e_min)

#plot the entropy against periods
	#fig = plt.figure()
	#plt.plot(periods, entropy, 'r')
	#plt.plot(periods, entropy, 'k+', markersize=12)
	#fig.suptitle(star)
	#plt.xlabel('Periods')
	#plt.ylabel('Entropy')
	#fig.savefig(os.path.join(out,name_of_plot))

numpy.savetxt(os.path.join(out, 'Allfiles_result'), ent_data, fmt='%s')

