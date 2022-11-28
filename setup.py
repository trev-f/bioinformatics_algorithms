from setuptools import find_packages, setup

setup(
    name='bioinformatics_textbook',
    packages=find_packages(),
    version='0.1.0',
    description='Notes and solutions for Bioinformatics Algorithms 3rd ed. by Compeau and Pevzner.',
    author='Trevor F. Freeman',
    license='MIT',
    install_requires=[
        'Click',
    ],
    entry_points={
        'console_scripts': [
            'bioinformatics-textbook = bioinformatics_textbook.cli:cli',
        ],
    },
)
