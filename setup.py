from setuptools import setup, find_packages

setup(name='scycle',
      version='0.1.10',
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
      ],
      zip_safe=False)