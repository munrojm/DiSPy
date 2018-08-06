from numpy.distutils.core import setup, Extension


extension_irr = Extension(
    name = 'DiSPy.irreps',
    extra_f77_compile_args=['-std=legacy'],
    sources=['DiSPy/irreps.pyf','DiSPy/pir_dat_module.f']
)

setup(
    name='DiSPy',
    version='0.1.0',
    author='Jason M. Munro',
    author_email='munrojm@psu.edu',
    packages=['DiSPy'],
    scripts=['bin/dispy'],
    license='LICENSE.txt',
    description='Utility to apply the distortion symmetry method.',
    long_description=open('README.rst').read(),
    provides=['dispy'],
    requires=[
        "numpy",
        "ase",
        "spglib"],
    ext_modules=[extension_irr],
)
