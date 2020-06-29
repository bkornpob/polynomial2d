import setuptools

setuptools.setup(
    name="polynomial2d",
    version="1.2.0",
    author="Kornpob Bhirombhakdi",
    author_email="kbhirombhakdi@stsci.edu",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    url="https://github.com/bkornpob/polynomial2d",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3"
        ,"License :: OSI Approved :: MIT License"
        ,"Operating System :: OS Independent"
    ],
    python_requires='>=3'
)
