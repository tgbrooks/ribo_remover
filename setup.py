from setuptools import setup, find_packages

setup(
    name="ribo-remover",
    version="0.1.0",
    #packages=find_packages(where="src"),
    packages = ["ribo_remover"],
    package_dir={"": "src"},
    package_data={
        "ribo_remover": ["data/*", "data/blast_db/*"],
    },
    python_requires=">=3.10",
)
