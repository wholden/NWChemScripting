from setuptools import setup

setup(name='NWChemScripting',
      version='0.1',
      description='Code for scripting calculations on NWChem.',
      author='',
      packages=['NWChemScripting'],
      requires=[],
      scripts=[
              'bin/center_xyz.py',
              'bin/do_vtc_xes.py',
              'bin/make_bq_charges.py',
              'bin/process_dat_to_spectrum.py',
              'bin/nw_spectrum_vtc_wespecmod.py'
      ]
      )
