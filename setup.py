import os, sys, numpy
from distutils.core import setup, Extension

IRMSD = Extension('minimsmbuilder/rmsdcalc',
                  sources = ["IRMSD/theobald_rmsd.c",
                             "IRMSD/rmsd_numpy_array.c"],
                  extra_compile_args=["-std=c99","-O2","-shared","-msse2",
                                      "-msse3","-fopenmp"],
                  extra_link_args=['-lgomp'],
                  include_dirs = [numpy.get_include(),
                                  os.path.join(numpy.get_include(), 'numpy')]
                 )

setup(name='minimsmbuilder',
      author = 'Robert McGibbon',
      packages = ['minimsmbuilder'],
      package_dir = {'minimsmbuilder': 'lib'},
      ext_modules = [IRMSD],
      )