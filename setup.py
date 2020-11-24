from setuptools import setup, find_packages

setup(name='scycle',
      version='0.1.5',
      description='Cell cycle pseudotime in single-cell RNA-seq',
      url='http://github.com/csgroen/scycle',
      author=['Clarice Groeneveld', 'Andrei Zinovyev'],
      author_email='clarice.groeneveld@ligue-cancer.net',
      license='MIT',
      packages=find_packages(),
      install_requires =[
          'scanpy',
          'anndata',
          'pandas',
          'numpy',
          'plotnine',
          'typing',
      ],
      include_package_data = True,
      package_data ={'': ['data/*.pkl']},
      zip_safe=False)