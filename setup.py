try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

    config = {
              'description': 'A Box model for Glacier Bay based on salinity conservation',
              'author': 'Chuning Wang',
              'url': '.',
              'download_url': '.',
              'author_email': 'chuning@esm.rutgers.edu',
              'version': '0.1',
              'install_requires': ['nose'],
              'packages': ['box_gb'],
              'scripts': [],
              'name': 'projectname'
             }

    setup(**config)
