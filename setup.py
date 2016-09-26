try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

try:
	from distutils.command.build_py import build_py_2to3 \
		as build_py
except ImportError:
	from distutils.command.build_py import build_py

config = {
    'description': 'hydradna',
    'author': 'Justin Vrana',
    'url': '',
    'download_url': '',
    'author_email': 'justin.vrana@gmail.com',
    'version': '0.1.0',
    'install_requires': [],
    'packages': ['hydradna'],
    'scripts': [],
    'name': 'hydradna',
    'license': 'Copyright University of Washington'
}

setup(**config)
