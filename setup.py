# filepath: GeneFinder1/setup.py
from setuptools import setup, find_packages

setup(
    name='GenePointer',
    version='1.0',
    packages=find_packages(),
    install_requires=[
        # List your dependencies here
    ],
    author='Aleksander Roosimaa',
    author_email='aleksander.roosimaa@ut.ee',
    description='A pipeline for finding/identifying resistance associated genes in bacterial genomes',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    url='',  # Replace with your GitHub URL
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
    ],
    python_requires='>=3.8',
)