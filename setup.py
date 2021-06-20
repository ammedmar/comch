import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="comch",
    version="0.1.2",
    author="Anibal M. Medina-Mardones",
    author_email="ammedmar@gmail.com",
    description="A package to study commutativity"
               +"up to coherent homotopies",
    long_description=long_description,
    long_description_content_type="text/markdown",
    #url="https://github.com/ammedmar/comch",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3',
)
