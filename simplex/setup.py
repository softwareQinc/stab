from skbuild import setup
from pathlib import Path

curdir = Path(__file__).parent
long_description = (curdir / "pysimplex.md").read_text()

setup(
    name="pysimplex",
    version="1.4.0",
    description="A fast Clifford circuit simulator",
    long_description=long_description,
    long_description_content_type="text/markdown",
    author="Alec Edgington",
    author_email="alec.edgington@cambridgequantum.com",
    license="Apache 2",
    url="https://github.com/CQCL/simplex",
    python_requires=">=3.8",
    classifiers=[
        "Environment :: Console",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "License :: OSI Approved :: Apache Software License",
        "Operating System :: POSIX :: Linux",
        "Operating System :: MacOS :: MacOS X",
        "Operating System :: Microsoft :: Windows",
        "Intended Audience :: Developers",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering",
    ],
    packages=["pysimplex"],
    cmake_install_dir="pysimplex",
)
