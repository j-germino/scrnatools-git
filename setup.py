import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="scrnatools",
    version="0.4",
    author="Joe Germino",
    author_email="joe.germino@ucsf.edu",
    description="Tools for single cell RNA sequencing pipelines",
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
        "datetime",
        "anndata",
        "scipy",
        "numba<0.57",
        "numpy<1.24",
        "starlette<0.28.0"
    ],
    python_requires='>=3, <3.11',
)
