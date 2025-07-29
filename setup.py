from setuptools import setup, find_packages

setup(
    name='pathwas',
    version='0.1.0',
    description='Toolbox for pathway-level association studies',
    packages=find_packages(),
    install_requires=[
        'numpy',
        'pandas',
        'scipy',
        'scikit-allel',
        'mygene',
        'gseapy',
    ],
    entry_points={
        'console_scripts': [
            'pathwas=pathwas.cli:main'
        ]
    },
)