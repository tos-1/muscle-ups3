'''
@PHDTHESIS{2015PhDT.......115S,
       author = {{Sissom}, Daniel J.},
        title = "{Early Growth in a Perturbed Universe: Exploring Dark Matter Halo Populations in 2LPT and ZA Simulations}",
     keywords = {Astrophysics;Astronomy, Astrophysics - Cosmology and Nongalactic Astrophysics, Astrophysics - Astrophysics of Galaxies},
       school = {Vanderbilt University},
         year = 2015,
        month = jan,
       adsurl = {https://ui.adsabs.harvard.edu/abs/2015PhDT.......115S},
      adsnote = {Provided by the SAO/NASA Astrophysics Data System}
}
'''

import sys
import struct
import os
import numpy as N

def read_bgc2 ( filename ) :
	offset = 4
	groupoffset = 8
	particleoffset = 8
	headersize = 1024
	groupsize = 4*8 + 10*4
	particlesize = 1*8 + 6*4
	headerformat = '=Q 16q 19d'
	groupformat = '=2q 2Q 10f'
	particleformat = '=q 6f '
	#print " Reading " + filename + " ... "
	fd = open ( filename , 'rb')
	bin_string = fd.read()
	fd.close()
	#print " Finished reading file "
	bin_string = bin_string[offset:]

	# Header stuff
	header_bin = bin_string[:headersize]
	header_pad = headersize - 36*8
	header = list( struct.unpack( headerformat , header_bin[:-header_pad]) )

	# Group stuff
	ngroups = header[8]
	#print ' ngroups  = ' , ngroups
	groupstart = headersize + groupoffset
	groupend = groupstart + ngroups * groupsize
	group_bin = bin_string[groupstart:groupend]
	group = []
	for i in range ( ngroups ):
		group.append( list( struct.unpack( groupformat , group_bin[ i * groupsize :( i +1) * groupsize ]) ) )	

	# Particle stuff
	particlestart = headersize + groupoffset + ngroups * groupsize + particleoffset
	particle_bin = bin_string[particlestart:]
	particle = []
	p_start = 0
	for i in range ( ngroups ):
		npart = group[i][2]
		particle.append([])
		for j in range ( npart ):
			particle[ i ].append( list( struct.unpack( particleformat , particle_bin[ p_start : p_start + particlesize ]) ) )
			p_start += particlesize
		p_start += particleoffset

	#print " Finished parsing bgc2 file "
	return header , group , particle
