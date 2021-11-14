from distutils.core import setup, Extension
import numpy.distutils.misc_util
import os

bins = os.getcwd() + '/ext/'
if not os.path.exists(bins):
    os.mkdir(bins)
    print('creating the folder for binaries', bins)

extra_compile_args = ['-std=c99', '-Wall']

c_ext1 = Extension("ext/_halomodel",
		  sources=["src/_halomodel.c","src/halomodel.c","src/pnum.c","src/bresenham.c","src/search.c","src/sort.c"],
                  include_dirs=['/home/federico/gsl/include'],
                  library_dirs=['/home/federico/gsl/lib'],
                  extra_compile_args=extra_compile_args,
                  extra_link_args=['-lgsl', '-lm', '-lgslcblas'])
c_ext2 = Extension("ext/_mas", ["src/_mas.c","src/mas.c"], extra_compile_args=extra_compile_args)
c_ext3 = Extension("ext/_Paint", ["src/_Paint.c","src/Paint.c"], extra_compile_args=extra_compile_args)
c_ext4 = Extension("ext/_masHalos", ["src/_masHalos.c","src/masHalos.c"], extra_compile_args=extra_compile_args)

setup(
    ext_modules=[c_ext1, c_ext2, c_ext3, c_ext4 ],
    include_dirs=numpy.distutils.misc_util.get_numpy_include_dirs()
)
