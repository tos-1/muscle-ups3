# MUSCLEUPS
This repo contains analytical Lagrangian algorithms to perform approximate gravitational evolution of cold dark matter. It performs Zel'dovich, 2lpt, alpt, MUSCLE  and MUSCLEUPS. For more details, see the [paper](https://ui.adsabs.harvard.edu/abs/2021MNRAS.505.2999T/abstract).

The cosmological module _cosmo/cosmo.py_ generates the linear power spectrum, which is produced with Eisentein and Hu fitting functions by default, if CLASS is not installed. Alternatively, it will call [CLASS](https://lesgourg.github.io/class_public/class.html). We tested it with version 2.9.3 

_utils/readgadget.py_, _utils/readsnap.py_ are taken from [Pylians](https://pylians3.readthedocs.io/en/master/gadget.html), they are not necessary to run the code, only to read the Gadget binaries in which the particles are saved, e.g if you want to use your own initial conditions or to compute basic statistics in _utils/make_statistics.py_.

## Installation

It is easy to run the code in a virtual enviroment. You can install it from your command line and activate it
```bash
virtualenv -p /usr/local/bin/python3.6 venv #(pyfftw creates problems with python3.8)
source venv/bin/activate
pip install numpy scipy colossus pyfftw h5py multiprocess 
pip install matplotlib ipython
pip install jedi==0.17.2 #(downgrade jedi because it clashes with ipython)
```
to compile, in setup.py change the path of _include_dirs/_ and _library_dirs/_ to point to your gsl libraries, then
```bash
make
```
To use the code in a notebook, you must create a kernel
```bash
pip install ipykernel
python -m ipykernel install --user --name=muscleups
```

## Usage
from terminal you can use the _run.py_ file
```python
python run.py -paramfile='params.py'
```

The _paramfile.py_ contains the parameters of the simulations. More information are in the documentation of _muscleups.py_. You can see all the option arguments and their explanation

```python
python run.py -h
```

and in the _run.py_. You can look at _quickstart.ipynb_ notebook for an interactive example.

MUSCLEUPS forces the collapse of halo particles with respect to traditional Lagrangian algorithms of analytical dark matter displacement
![This is an image](https://github.com/tos-1/muscle-ups3/blob/master/images/pos.png)

This results in an improved matter power spectrum.
![This is an image](https://github.com/tos-1/muscle-ups3/blob/master/images/pk.png)

One can show that also the cross correlation with N-body increases. For more details see the [paper](https://ui.adsabs.harvard.edu/abs/2021MNRAS.505.2999T/abstract).

## License
[MIT](https://choosealicense.com/licenses/mit/)
