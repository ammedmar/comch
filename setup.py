import setuptools
from pathlib import Path

with open("README.md", "r") as fh:
    long_description = fh.read()

version_namespace = {}
exec(Path("comch/_version.py").read_text(), version_namespace)

notebook_requirements = [
    "ipykernel",
    "jupyterlab",
    "nbformat",
]

docs_requirements = [
    "sphinx",
    "sphinx-rtd-theme",
]

setuptools.setup(
    name="comch",
    version=version_namespace["__version__"],
    author="Anibal M. Medina-Mardones",
    author_email="ammedmar@gmail.com",
    description="A specialized computer algebra system for the study of commutativity up to coherent homotopies",
    license="MIT",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/ammedmar/comch",
    project_urls={
        "Documentation": "https://comch.readthedocs.io/en/latest/",
        "Source": "https://github.com/ammedmar/comch",
        "Tracker": "https://github.com/ammedmar/comch/issues",
    },
    packages=setuptools.find_packages(exclude=["*.tests", "*.tests.*"]),
    install_requires=[],
    extras_require={
        "docs": docs_requirements,
        "notebooks": notebook_requirements,
        "dev": docs_requirements + notebook_requirements + [
            "build",
            "twine",
        ],
    },
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Science/Research",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3 :: Only",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Programming Language :: Python :: 3.12",
        "Programming Language :: Python :: 3.13",
        "Operating System :: OS Independent",
        "Topic :: Scientific/Engineering :: Mathematics",
    ],
    python_requires=">=3.10",
)
