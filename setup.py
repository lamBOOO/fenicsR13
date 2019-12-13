"Setup file"

from setuptools import setup

def readme(): # pylint: disable=missing-function-docstring
    with open('README.rst') as file:
        return file.read()

setup(
    name="fenicsR13",
    version="1.0",
    description="FEniCS (FEM) Solver for Non.-Eq.-Gases Based on R13 Equations",
    long_description=readme(),
    url="https://git.rwth-aachen.de/lamBOO/fenicsR13",
    author="Lambert Theisen, Manuel Torrilhon",
    author_email="lambert.theisen@rwth-aachen.de, mt@mathcces.rwth-aachen.de",
    license="None",
    packages=["fenicsR13"],
    zip_safe=False,
    entry_points={
        "console_scripts": [
            "fenicsR13=fenicsR13.fenicsR13:main",
            "geoToH5=fenicsR13.geoToH5:geo_to_h5"
        ],
    }
)
