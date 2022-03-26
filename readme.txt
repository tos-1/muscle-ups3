This repo contains analytical Lagrangian algorithms to perform approximate gravitational evolution of cold dark matter. It performs Zel'dovich, 2lpt, alpt, MUSCLE  and MUSCLEUPS. For more details, see https://ui.adsabs.harvard.edu/abs/2021MNRAS.505.2999T/abstract.

The cosmological module cosmo.py generates the linear power spectrum, which is produced with Eisentein and Hu fitting functions by default, if CLASS is not installed. Alternatively, one can use CLASS (https://lesgourg.github.io/class_public/class.html version 2.9.3), 

readgadget.py, readsnap.py are taken from Pylians (https://pylians3.readthedocs.io/en/master/gadget.html), they are not necessary to run the code, only to read the Gadget binaries in which the particles are saved, and compute the basic statistics in make_statistics.py.

See the quickstart.ipynb notebook for an example.

It is easy to run the code in a virtual enviroment

>> virtualenv -p /usr/local/bin/python3.6 venv #(pyfftw creates problems with python3.8)
>> source venv/bin/activate
>> pip install numpy scipy colossus pyfftw h5py multiprocess 
#optionals
#>> pip install matplotlib ipython
#>> pip install jedi==0.17.2 (downgrade jedi because it clashes with ipython)

to install
>> make

# if you want to use a notebook
>> pip install ipykernel
>> python -m ipykernel install --user --name=muscleups

from terminal you can use the run.py file

If you find this code useful, please cite https://ui.adsabs.harvard.edu/abs/2021MNRAS.505.2999T/abstract
