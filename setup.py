from setuptools import setup, find_packages

setup(
    name='cemba-data',
    use_scm_version=True,
    setup_requires=['setuptools_scm'],
    author='Hanqing Liu',
    author_email='hanliu@salk.edu',
    description='Pipelines for single nucleus methylome and multi-omic dataset.',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    url='https://github.com/lhqing/cemba_data',
    license='MIT',
    classifiers=[
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.7",
    ],
    packages=find_packages(exclude=('doc',)),
    include_package_data=True,
    package_data={
        '': ['*.txt', '*.tsv', '*.csv', '*.fa', '*Snakefile', '*ipynb']
    },
    install_requires=['pandas>=1.0',
                      'numpy',
                      'seaborn==0.10',
                      'matplotlib',
                      'papermill',
                      'dnaio',
                      'pysam'],
    entry_points={
        'console_scripts': ['yap=cemba_data.__main__:main',
                            'yap-internal=cemba_data._yap_internal_cli_:internal_main'],
    }
)
