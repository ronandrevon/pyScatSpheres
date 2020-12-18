import setuptools

with open("README.md", "r") as fh:
   long_description = fh.read()

setuptools.setup(
    name="pyScatSpheres",
    version="1.0.3.dev3",
    author="Tarik Ronan Drevon",
    author_email="tarik.drevon@stfc.ac.uk",
    description="Scattering of array of spheres, scalar theory",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://pypi.org/project/pyScatSpheres",
    project_urls={
        'Documentation': 'https://pyscatspheres.readthedocs.io/en/latest/',
        'Source':'https://github.com/ronandrevon/pyScatSpheres',
    },
    keywords=["diffraction","spherical","scalar","wave","array"],
    python_requires='>=3.6',
    packages=['pyScatSpheres'],
    include_package_data=True,
    package_data={'':['data*.pkl']},
    install_requires=['numpy','scipy','pandas','colorama','py3nj','matplotlib','easygui'],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License ",
        "Operating System :: OS Independent",
    ],
)
