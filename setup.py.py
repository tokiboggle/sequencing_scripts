from setuptools import setup, find_packages

setup(
    name='p450-mutagenesis-analysis',
    version='0.1.0',
    packages=find_packages(),
    install_requires=[
        'biopython>=1.81',
        'numpy>=1.22.0'
    ],
    entry_points={
        'console_scripts': [
            'p450-analysis=src.main:main'
        ]
    },
    author='Your Name',
    description='P450 Enzyme Saturation Mutagenesis Analysis Pipeline',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
)
