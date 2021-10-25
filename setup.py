from setuptools import setup, find_packages

setup(name='scycle',
      version= '0.1.17',
      description='Cell cycle pseudotime in single-cell RNA-seq',
      url='http://github.com/csgroen/scycle',
      author=['Clarice Groeneveld', 'Andrei Zinovyev', 'Aziz Fouche'],
      author_email='clarice.groeneveld@curie.fr',
      license='MIT',
      packages=find_packages(),
      install_requires =[
          'scanpy',
          'elpigraph-python',
          'anndata',
          'pandas',
          'numpy',
          'matplotlib',
          'plotnine',
          'plotly',
          'typing',
          'scrublet',
          'requests',
          'stabilized-ica',
          'transmorph',
          'scikit-misc'
      ],
      dependency_links=[],
      zip_safe=False)
