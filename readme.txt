>> virtualenv -p /usr/local/bin/python3.6 venv #(pyfftw creates problems with python3.8)
>> source venv/bin/activate
>> pip install numpy scipy colossus pyfftw h5py 
#optionals
#>> pip install matplotlib ipython
#>> pip install jedi==0.17.2 (downgrade jedi because it clashes with ipython)

to install
>> make

# if you want to use a notebook
>> pip install ipykernel
>> python -m ipykernel install --user --name=muscleups

to run it interatively look at the notebook

from terminal you can use the run.py file
