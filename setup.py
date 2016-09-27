try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

try:
	from distutils.command.build_py import build_py_2to3 \
		as build_py
except ImportError:
	from distutils.command.build_py import build_py


def readme():
    with open('README.md') as f:
        return f.read()


setup(
    name='jdna',
    version='0.1',
    packages=['jdna'],
    url='https://github.com/jvrana/jdna',
    license='MIT',
    author='Justin Vrana',
    author_email='justin.vrana@gmail.com',
    description='A(nother) pythonic DNA manipulator',
    long_description=readme(),
    keywords='dna biology cloning',
    tests_require=['pytest'],
)
