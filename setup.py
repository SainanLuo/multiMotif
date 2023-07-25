from setuptools import setup, find_packages

setup(
    name='variamotif',
    version='0.1',
    packages=find_packages(),
    install_requires=[
        #
    ],
    author='Sainan Luo',
    description='A Computational Tool for Variable Motif scanning and Sequence-based Relative Position Visualization of Search Results in Sequences.',
    url='https://github.com/Luo-Sainan/VariaMotif',
    entry_points={
        'console_scripts': [
            'variamotif=variamotif_cli:main',
        ],
    },
)
