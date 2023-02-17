import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="scrnatools",
    version="0.3.3",
    author="Joe Germino",
    author_email="joe.germino@ucsf.edu",
    description="Tools for single cell RNA sequencing pipelines",
    url="https://github.com/j-germino/sc-rna-tools",
    packages=setuptools.find_packages(),
    install_requires=[
        "scanpy",
        "scrublet",
        "scvi-tools",
        "matplotlib",
        "pandas",
        "numpy",
        "seaborn",
        "memory-profiler",
        "leidenalg",
        "datetime",
        "anndata",
        "scipy",
    ]
)
