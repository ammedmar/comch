import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="chain_level_steenrod_operations",
    version="0.0.1",
    author="Anibal M. Medina-Mardones",
    author_email="anibalm3@gmail.com",
    description="A package to study Steenrod"
               +"operations at the chain level",
    long_description=long_description,
    long_description_content_type="text/markdown",
    #url="https://github.com/pypa/sampleproject",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)