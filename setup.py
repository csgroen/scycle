from setuptools import setup, find_packages
from scycle import __version__

setup(name='scycle',
      version= __version__,
      description='Cell cycle pseudotime in single-cell RNA-seq',
      url='http://github.com/csgroen/scycle',
      author=['Clarice Groeneveld', 'Andrei Zinovyev', 'Aziz Fouche'],
      author_email='clarice.groeneveld@ligue-cancer.net',
      license='MIT',
      packages=find_packages(),
      install_requires =[
          'scanpy',
          'anndata',
          'pandas',
          'numpy',
          'plotnine',
          'plotly',
          'typing',
          'POT',
          'scrublet',
          'sica @ git+https://github.com/ncaptier/Stabilized_ICA',
          'transmorph',
      ],
      dependency_links=[],
      zip_safe=False)
