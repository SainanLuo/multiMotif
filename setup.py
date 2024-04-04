from setuptools import setup, find_packages

# Read the contents of your README file
with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    name='multiMotif',
    version='1.2.0',
    packages=find_packages(),
    install_requires=[
        'biopython>=1.78',
        'matplotlib>=3.3',
        'numpy>=1.19',
        'pandas>=1.1',
        'plotly>=5.18.0'
    ],
    author='Sainan Luo',
    description='A Computational Tool for Variable Motif scanning and Sequence-based Relative Position Visualization of Search Results in Sequences.',
    long_description=long_description,  # The new line
    long_description_content_type='text/markdown',  # The new line
    url='https://github.com/SainanLuo/multiMotif',
    entry_points={
        'console_scripts': [
            'multiMotif=multiMotif.main:main',
        ],
    },
)
