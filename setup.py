import setuptools

with open("docs/README.md", "r") as fh:
   long_description = fh.read()

setuptools.setup(
    name="pyScatSpheres",
    version="1.0.2rc2",
    author="Tarik Ronan Drevon",
    author_email="tarik.drevon@stfc.ac.uk",
    description="Scattering of array of spheres, scalar theory",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://pypi.org/project/pyScatSpheres",
    project_urls={
        'Documentation': 'http://88.123.115.137:8000/projects/scattering/',
    },
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License ",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
    install_requires=['numpy','scipy','pandas','colorama','py3nj','matplotlib','easygui',''],
)
