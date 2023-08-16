from setuptools import setup, find_packages

setup(
    name='variamotif',
    version='0.1',
    packages=find_packages(),
    install_requires=[
        'biopython>=1.78',
        'matplotlib>=3.3',
        'numpy>=1.19',
        'pandas>=1.1'
    ],
    author='Sainan Luo',
    description='A Computational Tool for Variable Motif scanning and Sequence-based Relative Position Visualization of Search Results in Sequences.',
    url='https://github.com/Luo-Sainan/VariaMotif',
    entry_points={
        'console_scripts': [
            'variamotif=variamotif.main:main',
        ],
    },
)
