from setuptools import setup, find_packages

setup(
    name="gut-ibm-tools",
    version="0.1.0",
    description="Analysis toolkit for GutIBM Enterobacteriaceae simulation",
    packages=find_packages(),
    python_requires=">=3.9",
    install_requires=[
        "numpy>=1.21",
        "h5py>=3.0",
        "scipy>=1.7",
    ],
    extras_require={
        "viz": ["matplotlib>=3.5"],
        "dev": ["pytest>=7.0", "ruff>=0.1"],
    },
)
