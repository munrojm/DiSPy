from setuptools import find_packages
from numpy.distutils.core import setup, Extension


extension_irr = Extension(
    name="DiSPy.core.stokes",
    extra_f90_compile_args=["-std=legacy"],
    sources=["DiSPy/core/stokes.pyf", "DiSPy/core/pir_dat_module.f90"],
)

setup(
    name="DiSPy",
    version="0.3",
    author="Jason M. Munro",
    author_email="jmunro@lbl.gov",
    packages=find_packages(exclude=["tests", "*.tests", "*.tests.*", "tests.*"]),
    package_dir={"DiSPy": "DiSPy"},
    package_data={"DiSPy": ["SPG_dict.txt", "PIR_data.txt"]},
    scripts=["scripts/dispy"],
    license="LICENSE.txt",
    description="Utility to apply the distortion symmetry method.",
    long_description=open("README.rst").read(),
    provides=["dispy"],
    install_requires=["pymatgen==2022.0.8"],
    python_requires=">=3.8",
    ext_modules=[extension_irr],
)
