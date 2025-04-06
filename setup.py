from setuptools import setup, find_packages

setup(
    name="patient_deduplicator",
    version="0.1.0",
    packages=find_packages(),
    install_requires=[
        "pandas>=1.3.0",
        "numpy>=1.20.0",
        "matplotlib>=3.4.0",
        "pyarrow>=6.0.0",  # For parquet support
    ],
    author="Almog Kedem",
    author_email="almog.kedem@mail.huji.ac.il",
    description="A package for identifying potential duplicate patients across different studies",
    keywords="bioinformatics, oncology, patient deduplication",
    python_requires=">=3.8",
)