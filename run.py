'''
This is the file to call from a bash script to submit a job. 
It makes sense to leave the flag pos == False (default) to save the positions on a binary.

E.g.
>> python run.py -scheme=muscleups -box=100. -ng=128 -z_pk=50 -redshift=0.
'''

from muscleups import muscleups
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-cosmology',default = 'ehu', help = 'How to generate the linear power spectrum. The options are Eisenstein and Hu or CLASS')
parser.add_argument('-hubble', default=0.7, help = 'normalized hubble rate',type = float)
parser.add_argument('-omega_b',default = 0.0225, help = 'physical density of baryons',type = float)
parser.add_argument('-Omega_cdm',default = 0.25, help = 'cosmological cdm density',type = float)
parser.add_argument('-ns',default = 0.96, help = 'spectral index',type = float)
parser.add_argument('-s8',default = 0.8, help = 'power spectrum normalization',type = float) 
parser.add_argument('-z_pk',default = 50., help = 'redshift at which the initial power spectrum is computed',type = float)
parser.add_argument('-redshift',default = 0., help = 'redshift at which the output is computed',type = float)
parser.add_argument('-ng',default = 64, help = 'number of particles per side',type = int)
parser.add_argument('-box',default = 64., help = 'size of the box in Mpc/h',type = float)
parser.add_argument('-sigmaalpt',default = 4., help = 'interpolating scale between 2lpt and sc, in Mpc/h',type = float)
parser.add_argument('-scheme',default = '2lpt', help = 'scheme among which to choose the evolution. The options are zeld,2lpt,sc,hmlpt,muscle')
parser.add_argument('-smallscheme',default = None, help = 'selecting this activates alpt. It works only with sc and muscle, while 2lpt on large scales is automatically set')
parser.add_argument('-threads',default = 1., help = 'number of threads used by pyfftw',type = int)
parser.add_argument('-extra',default ='', help = 'initial stringany for the output fileroot',type = str)
parser.add_argument('-seed',default =1, help = 'seed of the random number generator of initial conditions',type = int)
parser.add_argument('-exact_pk',default = True, help = 'flag to not fix linear mean value of fourier amplitudes of density',action="store_false")
parser.add_argument('-makeic',default = False, help = 'flag to write the parameter file for Gadget2. If z_pk!=redshift an error is raised',action="store_true")
parser.add_argument('-return_pos',default = False, help = 'flag to returns the array of pos, instead of writing a binary',action="store_true")


args = vars(parser.parse_args())

print('the scheme you chose was', args['scheme'])

MSC = muscleups( cosmology = args['cosmology'], h = args['hubble'] , omega_b = args['omega_b'] , Omega_cdm = args['Omega_cdm'], ns = args['ns'] , sigma8 =args['s8'] , z_pk = args['z_pk'], redshift = args['redshift'], sigmaalpt = args['sigmaalpt'], scheme =args['scheme'], ng = args['ng'] ,boxsize = args['box'] ,smallscheme = args['smallscheme'] ,
makeic = args['makeic'], return_pos=args['return_pos'],threads = args['threads'] , extra_info = args['extra'] ,seed=args['seed'], exact_pk=args['exact_pk'] )

MSC.generate()
