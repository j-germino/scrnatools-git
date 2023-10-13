import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="scrnatools",
    version="1.11",
    author="Joe Germino",
    author_email="joe.germino@ucsf.edu",
    description="Tools for single cell RNA sequencing pipelines",
    long_description="scrnatools provides helper plotting and analysis functions built on top of the scanpy and scVI packages to perform common tasks in scRNA-seq analysis pipelines. Extensive description of general best practices and detailed default pipeline examples can be found in the scanpy and scvi-tools documentation. For an example walkthrough of a typical scRNA-seq pipeline and the added functionality provided by scrnatools see the jupyter notebooks provided at https://github.com/j-germino/scrnatools-git/tree/main/examples. These can be uploaded to Google drive and run within a Google colab environemnt, using a GPU session. API documentation is available at: https://scrnatools-git.readthedocs.io/en/latest/.",
    url="https://github.com/j-germino/sc-rna-tools",
    packages=setuptools.find_packages(),
    install_requires=[
        "scanpy",
        "scrublet",
        "scvi-tools",
        "matplotlib<3.6",
        "pandas<2.0.0",
        "numpy",
        "seaborn",
        "scikit-misc",
        "leidenalg",
        "anndata",
        "scipy",
        "numba<0.57",
        "numpy<1.24",
        "starlette<0.28.0"
    ],
    python_requires='>=3, <3.11',
)
