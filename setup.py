from setuptools import setup, find_packages


setup(name='BACILLUSme',
      version='0.0.9',
      description='ME-model data files specific to Bacillus subtilis str. 168',
      author='Juan D. Tibocha-Bonilla',
      url='https://github.com/jdtibochab/bacillusme',
      install_requires=['Biopython', 'setuptools', 'cobra<=0.5.11', 'xlrd',
                        'pandas', 'six', 'scipy', 'numpy'],
      packages=find_packages(),
      package_data={'': ['building_data/*', 'me_models/*']},
      license='MIT')
