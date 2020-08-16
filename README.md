This is a fork of original SYMORO repository. As the oiginal repository does not seems to be maintained anymore, this one offers minimal maintenance and will try to deal with bug fixes and pull requests.

SYMORO
======

SYOMRO is a software package for SYmbolic MOdeling of RObots.

This software package is developed as part of the OpenSYMORO project by
the robotics team at [IRCCyN][lk:irccyn] under the supervision of Wisama
Khalil.

For details on the algorithms used, please see [the paper][lk:hal]
published in the AIM 2014 conference.

Requirements
------------

NOT TESTED WITH ALL VERSIONS OF DEPENDENCIES. You should stick to the lower required version to avoid issues.

+ python (>= 2.7) For python 3 version, choose the branch python3 forked from https://github.com/cmigliorini/symoro
+ sympy (>= 0.7.6)
+ numpy (>= 1.6.1)
+ wxPython (>= 2.8.12) NOT COMPATIBLE WITH >=4.0.0 
+ PyOpenGL (>= 3.0.1b2)


Getting Started
---------------
+ Install on windows: have a look at [Setup][lk:setup] but does not seem to work first try.

+ Install on ubuntu (tested 14.04 and 20.04):

  `sudo apt install python-setuptools`
  `sudo apt install python-pip`
  `sudo pip install sympy==0.7.6`
  `sudo appt install python-dev`
  `sudo pip install numpy==1.6.1`
  `sudo apt install python-wxgtk2.8 #(pythron-wxgtk3.0 also tested)`
  `sudo pip install PyOpenGL==3.0.1b2straight away`

  


Licence
-------
See [LICENCE][lk:licence].

# Known issues

- [Setup][lk:setup] from setup.py fails on ubuntu and windows. Major issue apparently comes from wxPython install

- Jacobian and Jacobian determinant computation fails if intermediate frame == 0. 
- IGM - Pieper method fails

