from setuptools import setup, find_packages

yap_version = '1.0.0'

setup(
    name='yap',
    version=yap_version,
    author='Hanqing Liu',
    author_email='hanliu@salk.edu',
    packages=find_packages(),
    description='Pipelines for single nucleus methylome and multi-omic dataset.',
    long_description=open('README.md').read(),
    include_package_data=True,
    install_requires=['pandas>=1.0', 'numpy', 'seaborn', 'matplotlib', 'papermill', 'dnaio', 'pysam'],
    package_data={
        '': ['*.txt', '*.tsv', '*.csv', '*.fa', '*Snakefile', '*ipynb']
    },
    entry_points={
        'console_scripts': ['yap=cemba_data.__main__:main',
                            'yap-internal=cemba_data._yap_internal_cli_:internal_main'],
    }
)

if __name__ == '__main__':
    f = open("cemba_data/__init__.py", 'w')
    f.write(f"__version__ = '{yap_version}'\n")
    f.close()
