======================
Biopython Installation
======================

This document describes installing Biopython on your computer. To make things as simple as possible, 
it basically assumes you have nothing related to Python or Biopython on your computer and want to 
end up with a working installation of Biopython. Biopython supports Mac, UNIX/Linux and Windows machines, but the 
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

We recommend you install Biopython and its dependencies using the provided pre-compiled Windows Installers. i.e. You don’t need a C compiler. 
See subsequent section for more details.

Installing Python
=================

Python is a interpreted, interactive object-oriented programming language and the home for all things python is `https://www.python.org. <https://www.python.org/>`_

Biopython is designed to work with Python 2.7, Python 3.4, Python 3.5 and 3.6. As of Biopython 1.68, support for Python 3.3 has been removed.

1) Installing Python in UNIX/Linux
----------------------------------

Python 2.7 and Python 3.5 usually comes pre-installed with most modern Linux systems like Ubuntu, Fedora. If you would like to check 
the version present in your system, you can type the command::

  $ python --version
  
This command will show the version of python2 present in your system. To check the exact version of python3 interpreter, type::

  $ python3 --version
    
This command will show the version of python3 present in your system. If you aren't able to find the Python software available 
in your Linux system, you can download Python from the site here- `https://www.python.org/downloads/source/ <https://www.python.org/downloads/source/>`_ and select the source gzipped tarball for the version you would like to install.
For example, If you want to install Python 3.6 in your system, Look for the latest stable Python 3.6 release and Click on the "Download Gzipped source tarball".

Once the download is done, navigate to the folder where the source tarball is stored. For example, if the tarball is stored in your Downloads folder, navigate to Downloads folder and do::

  $ tar -zxvf Python-3.6.1.tgz
  
Then enter into the created directory::

  $ cd Python-3.6.1
  
Now, start the build process by configuring everything into your system::

  $ ./configure
  
Build all of the files with::
  
  $ make

Finally, you’ll need to have root permissions on the system and then install everything::
  
  $ make install
  
If there were no errors and everything worked correctly, you should now be able to type ``python`` at a command prompt and enter into the python interpreter::

  $ python
  Python 2.7.5 (...)
  ...
  Type "help", "copyright", "credits" or "license" for more information.
  >>>
  
(The precise version text and details will depend on the version you installed and your operating system.)



2) Python installation on Windows
---------------------------------

Installation on Windows is most easily done using handy windows installers. As described above in the UNIX section, you should go to the webpage for the current stable version of Python to download this installer. At the current time, you’d go to `https://www.python.org/downloads/windows/. <https://www.python.org/downloads/windows/>`_ and depending on the type of system and version required, download Python-x.x.x.msi (where x.x.x is the exact version number, like 3.6.1).

The installer is an executable program, so you only need to double click it to run it. Then just follow the friendly instructions. On all newer Windows machines you’ll need to have Administrator privileges to do this installation.

3) Python installation on Mac OS X
----------------------------------

Apple includes python on Mac OS X, and while you can use this many people have preferred to install the latest version of python as well in parallel. We refer you to `https://www.python.org <https://www.python.org>`_ for more details, although otherwise the UNIX instructions apply.


Installing Biopython in your system
===================================

UNIX/Linux
----------

For installing Biopython, you will require a software package called ``pip`` to do so. If you are using the latest 3.x version of Python, you will have the package preinstalled in your UNIX/Linux system. You can check the version of the software by doing::

  $ pip --version
  pip 9.0.1 from /usr/lib/python3.6/site-packages (python 3.6)
  
If you don't see an output like the one above, you do not have pip installed in the system. To install pip in your UNIX/Linux system, follow the instructions given in the link here-`https://packaging.python.org/guides/installing-using-linux-tools/ <https://packaging.python.org/guides/installing-using-linux-tools/>`_


Once ``pip`` is installed, Go to your terminal and type::

  $ pip install biopython
  
This command will fetch the latest version of biopython from Python Package Index and start the installation procedure. The command will also install any of the required dependencies of Biopython which are not available in your system, like `Numpy <http://www.numpy.org/>`_.

Once the installation is done, you can check if the installation is working by firing up the python interpreter in the terminal and doing the following::
  
  $ python
  Python 3.5.1+ (default, Mar 30 2016, 22:46:26) 
  [GCC 5.3.1 20160330] on linux
  Type "help", "copyright", "credits" or "license" for more information.
  >>> import Bio
  >>> Bio.__version__
  1.69
  
If you see an output like the one above, you have successfully installed Biopython in your system. Otherwise you will notice an error.

If you wish to install a specific version of Biopython for the system, you can type in the terminal::

  $ pip install biopython==1.68
  
This command will install 1.68 version of Biopython in your system.

Mac OS X
--------

To install Biopython in Mac OS X, you will need to first go to `http://biopython.org/wiki/Download <http://biopython.org/wiki/Download>`_ and download the latest gzipped tarball from the link. Once it is done, Extract the tarball and then go to the terminal and type::

  $ python setup.py build 
  $ python setup.py test
  $ sudo python setup.py install

This will install the biopython package in your system. To check the installation, follow the steps mentioned in the UNIX/Linux section. To install a specific version of Biopython, search for the tarball in the biopython download link and follow the steps as it is.

