.. _`chapter:contributing`:

Where to go from here – contributing to Biopython
=================================================

Bug Reports + Feature Requests
------------------------------

Getting feedback on the Biopython modules is very important to us.
Open-source projects like this benefit greatly from feedback,
bug-reports (and patches!) from a wide variety of contributors.

The main forums for discussing feature requests and potential bugs are
the `Biopython mailing list <http://biopython.org/wiki/Mailing_lists>`__
and issues or pull requests on GitHub.

Additionally, if you think you’ve found a new bug, you can submit it to
our issue tracker at https://github.com/biopython/biopython/issues (this
replaced the older Open Bioinformatics Foundation hosted RedMine
tracker). This way, it won’t get buried in anyone’s Inbox and forgotten
about.

Mailing lists and helping newcomers
-----------------------------------

We encourage all our uses to sign up to the main Biopython mailing list.
Once you’ve got the hang of an area of Biopython, we’d encourage you to
help answer questions from beginners. After all, you were a beginner
once.

Contributing Documentation
--------------------------

We’re happy to take feedback or contributions - either via a bug-report
or on the Mailing List. While reading this tutorial, perhaps you noticed
some topics you were interested in which were missing, or not clearly
explained. There is also Biopython’s built-in documentation (the
docstrings, these are also
`online <http://biopython.org/DIST/docs/api>`__), where again, you may
be able to help fill in any blanks.

Contributing cookbook examples
------------------------------

As explained in Chapter :ref:`chapter:cookbook`,
Biopython now has a wiki collection of user contributed “cookbook”
examples, http://biopython.org/wiki/Category:Cookbook – maybe you can
add to this?

.. _`sec:maintain_dist`:

Maintaining a distribution for a platform
-----------------------------------------

We currently provide source code archives (suitable for any OS, if you
have the right build tools installed), and pre-compiled wheels via
https://github.com/biopython/biopython-wheels to cover the major
operating systems.

Most major Linux distributions have volunteers who take these source
code releases, and compile them into packages for Linux users to easily
install (taking care of dependencies etc). This is really great and we
are of course very grateful. If you would like to contribute to this
work, please find out more about how your Linux distribution handles
this. There is a similar process for conda packages via
https://github.com/conda-forge/biopython-feedstock thanks to the
conda-forge team.

Below are some tips for certain platforms to maybe get people started
with helping out:

Windows
   – You must first make sure you have a C compiler on your Windows
   computer, and that you can compile and install things (this is the
   hard bit - see the Biopython installation instructions for info on
   how to do this).

RPMs
   – RPMs are pretty popular package systems on some Linux platforms.
   There is lots of documentation on RPMs available at
   http://www.rpm.org to help you get started with them. To create an
   RPM for your platform is really easy. You just need to be able to
   build the package from source (having a C compiler that works is thus
   essential) – see the Biopython installation instructions for more
   info on this.

   To make the RPM, you just need to do:

   .. code:: console

      $ python -m build
      $ # Use the generated .tar.gz to create RPM with rpmbuild and a SPEC file

Macintosh
   – Since Apple moved to Mac OS X, things have become much easier on
   the Mac. We generally treat it as just another Unix variant, and
   installing Biopython from source is just as easy as on Linux. The
   easiest way to get all the GCC compilers etc installed is to install
   Apple’s X-Code. We might be able to provide click and run installers
   for Mac OS X, but to date there hasn’t been any demand.

Once you’ve got a package, please test it on your system to make sure it
installs everything in a good way and seems to work properly. Once you
feel good about it, make a pull request on GitHub and write to our
`Biopython mailing list <http://biopython.org/wiki/Mailing_lists>`__.
You’ve done it. Thanks!

Contributing Unit Tests
-----------------------

Even if you don’t have any new functionality to add to Biopython, but
you want to write some code, please consider extending our unit test
coverage. We’ve devoted all of
Chapter :ref:`chapter:testing` to this topic.

Contributing Code
-----------------

There are no barriers to joining Biopython code development other than
an interest in creating biology-related code in Python. The best place
to express an interest is on the Biopython mailing lists – just let us
know you are interested in coding and what kind of stuff you want to
work on. Normally, we try to have some discussion on modules before
coding them, since that helps generate good ideas – then just feel free
to jump right in and start coding!

The main Biopython release tries to be fairly uniform and interworkable,
to make it easier for users. You can read about some of (fairly
informal) coding style guidelines we try to use in Biopython in the
contributing documentation at http://biopython.org/wiki/Contributing. We
also try to add code to the distribution along with tests (see
Chapter :ref:`chapter:testing` for more info on the
regression testing framework) and documentation, so that everything can
stay as workable and well documented as possible (including docstrings).
This is, of course, the most ideal situation, under many situations
you’ll be able to find other people on the list who will be willing to
help add documentation or more tests for your code once you make it
available. So, to end this paragraph like the last, feel free to start
working!

Please note that to make a code contribution you must have the legal
right to contribute it and license it under the Biopython license. If
you wrote it all yourself, and it is not based on any other code, this
shouldn’t be a problem. However, there are issues if you want to
contribute a derivative work - for example something based on GPL or
LPGL licensed code would not be compatible with our license. If you have
any queries on this, please discuss the issue on the mailing list or
GitHub.

Another point of concern for any additions to Biopython regards any
build time or run time dependencies. Generally speaking, writing code to
interact with a standalone tool (like BLAST, EMBOSS or ClustalW) doesn’t
present a big problem. However, any dependency on another library - even
a Python library (especially one needed in order to compile and install
Biopython like NumPy) would need further discussion.

Additionally, if you have code that you don’t think fits in the
distribution, but that you want to make available, we maintain Script
Central (http://biopython.org/wiki/Scriptcentral) which has pointers to
freely available code in Python for bioinformatics.

Hopefully this documentation has got you excited enough about Biopython
to try it out (and most importantly, contribute!). Thanks for reading
all the way through!
