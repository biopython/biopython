======================
Biopython Installation
======================

This document describes how to install Biopython in your system. To make things as simple as possible, 
it basically assumes you have nothing related to Python or Biopython on your computer and want a working 
installation of Biopython. Biopython supports Mac, UNIX/Linux and Windows machines, but the 
steps to install biopython vary across the platforms.


C Compiler
==========

1) Unix/Linux
-------------

GCC is recommende as the C compiler, which is usually available as part of the standard set of packages on any Unix or Linux system.

2) Mac OS X
------------

Install Apple’s XCode suite from the App Store, and then from the Xcode options also install the optional command line utilities.

3) Windows
-----------

We recommend you install Biopython and its dependencies using the provided pre-compiled biopython wheels. i.e. You don’t need a C compiler. 
See subsequent section for more details.

Installing Python
=================

Python is a interpreted, interactive object-oriented programming language and the home for all things python is `https://www.python.org. <https://www.python.org/>`_

To install Python in your system, follow the instructions given in the "Properly Install Python" section of `The Hitchhiker's Guide to Python <http://docs.python-guide.org/en/latest/>`_ based on the platform you are working on. 

Biopython is designed to work with Python 2.7, Python 3.4, Python 3.5 and 3.6. As of Biopython 1.68, support for Python 3.3 has been removed.

Installing Biopython in your system
===================================

UNIX/Linux
----------

For installing Biopython, you will require a software package called ``pip`` to do so. If you are using the latest 3.x version of Python, you will have the package preinstalled in your UNIX/Linux system. You can check the version of the software by doing::

  $ pip --version
  pip 9.0.1 from /usr/lib/python3.6/site-packages (python 3.6)
  
If you don't see an output like the one above, you do not have pip installed in the system. To install pip in your UNIX/Linux system, you can do::

  $ python -m ensurepip
  
Once ``pip`` is installed, Go to the terminal and type::

  $ pip install biopython
  
This command will fetch the latest version of biopython from Python Package Index and start the installation procedure. The command will also install any of the required dependencies of Biopython which are not available in your system, like `Numpy <http://www.numpy.org/>`_.

Once the installation is done, you can check if the installation is working by firing up the python interpreter in the terminal and doing the following::
  
  $ python
  Python 3.6.1+ (default, Mar 30 2016, 22:46:26) 
  [GCC 5.3.1 20160330] on linux
  Type "help", "copyright", "credits" or "license" for more information.
  >>> import Bio
  >>> Bio.__version__
  1.69
  
If you see an output like the one above, you have successfully installed Biopython in your system. Otherwise you will notice an error.

If you wish to install a specific version of Biopython for the system, you can type in the terminal::

  $ pip install biopython==1.68
  
This command will install 1.68 version of Biopython in your system. Once you have biopython installed, you can proceed ahead.

Windows
-------

For installing Biopython, you will require a software package called ``pip`` to do so. If you are using the latest 3.x version of Python, you will have the package preinstalled in your UNIX/Linux system. You can check the version of the software by doing::

  $ pip --version
  pip 9.0.1 from /usr/lib/python3.6/site-packages (python 3.6)
  
If you don't see an output like the one above, you do not have pip installed in the system. To install pip in your UNIX/Linux system, you can do::

  $ python -m ensurepip
  
Once ``pip`` is installed, Go to the terminal and type::

  $ pip install biopython
  
This command will fetch latest version of biopython from PyPI and start the installation. The command will also install any of the required dependencies of Biopython which are not available in your system, like `Numpy <http://www.numpy.org/>`_.

Once the installation is done, you can check if the installation is working by firing up the python interpreter in the terminal and doing the following::
  
  $ python
  Python 3.6.1+ (default, Mar 30 2016, 22:46:26) 
  [GCC 5.3.1 20160330] on linux
  Type "help", "copyright", "credits" or "license" for more information.
  >>> import Bio
  >>> Bio.__version__
  1.69
  
If you see an output like the one above, you have successfully installed Biopython in your system. Otherwise you will notice an error.

If you wish to install a specific version of Biopython for the system, you can type in the terminal::

  $ pip install biopython==1.68
  
This command will install 1.68 version of Biopython in your system. Once you have biopython installed, you can proceed ahead.


Mac OS X
--------

To install Biopython in Mac OS X, you will need to have the package ``pip`` installed in your system.

You can install `pip` in Mac OS X using the following command::

  $ python -m ensurepip

This command will install the software package in your system. To proceed further, do::

  $ pip install biopython
  
This will install the biopython package in your system. For checking the installation or installing a specific version of biopython in your system, follow the steps mentioned in the UNIX/Linux section for installing biopython in your system

Bonus: Installing Biopython using Anaconda
==========================================

For Anaconda users, the Biopython source files are available in the `conda-forge <https://conda-forge.org/>`_ channel. To install Biopython, just do::

  $ conda install -c conda-forge biopython
  
You will see an output something like this::

  Fetching package metadata ...........
  Solving package specifications: ..........
  
  Package plan for installation in environment Users/username/Miniconda2:
  
  The following packages will be downloaded:
  
  package                    |            build
  ---------------------------|-----------------
  conda-env-2.6.0            |                0          498 B
  mkl-2017.0.3               |                0       126.3 MB
  numpy-1.13.1               |           py27_0         3.3 MB
  biopython-1.70             |      np113py27_0         2.1 MB
  conda-4.3.22               |           py27_0         520 KB
  ------------------------------------------------------------
                                         Total:       132.2 MB
                                         
  The following NEW packages will be INSTALLED:

    biopython: 1.68-np113py36_0
    mkl:       2017.0.3-0
    numpy:     1.13.1-py36_0
  
  Proceed ([y]/n)?
  
Type "y" and the installation will start. Once installation is finished, you can check if the installation worked properly by doing the steps mentioned in the UNIX/Linux section above.

Making use of Anaconda platform will ensure that you get all the dependencies required to run biopython in your system successfully. To get Anaconda on your system, you can refer to the instructions given here- `https://www.continuum.io/downloads <https://www.continuum.io/downloads>`_
