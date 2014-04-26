#Datacard Viewer


Datacard Viewer is Graphical User Interface tool for exploration of Higgs measurements in datacards.

The Higgs searches and measurements are encoded in complex likelihood functions that comprise
thousands parameters. For convenience of extracting results, this information is stored in ROOT
(RooFit) format. On top of that there is a layer of ASCII configuration files describing the contents of the
ROOT files, which is parsed using Python. Diagnostic information about these inputs is also produced in
form of ROOT or text files. The amount of information has become too large to inspect manually with a
text editor or the root browser. 


##Instalation

###Setting up a Web server

First of all you, you need to install Flask(A Python Web Microframework). 
More information: [Flask installation](http://flask.pocoo.org/docs/installation/#installation)

**Requirements:** Python 2.6 or higher.

It can be done by:

1) simply typing in terminal:

   ~~~ sh
   $ pip install Flask
   ~~~
2) building and installing from the source:

to build:
   ~~~ sh
   $ python libs/flask-0.11/setup.py build
   ~~~

to install
   ~~~ sh
   $ sudo python libs/flask-0.11/setup.py install
   ~~~

###Configuring the server

To configure your web server go to [run.py](run.py) and change IP address. 

###Launching the server

To configure your web server type:

   ~~~ sh
   $ python run.py
   ~~~
##Documentation

+ [Documentation](DOCUMENTATION.md)
+ [Tasks](TASKS.md)

