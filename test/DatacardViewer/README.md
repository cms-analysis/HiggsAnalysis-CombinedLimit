#Datacard Viewer


Datacard Viewer is Graphical User Interface tool for exploration of Higgs measurements in datacards.

The Higgs searches and measurements are encoded in complex likelihood functions that comprise
thousands parameters. For convenience of extracting results, this information is stored in ROOT
(RooFit) format. On top of that there is a layer of ASCII configuration files describing the contents of the
ROOT files, which is parsed using Python. Diagnostic information about these inputs is also produced in
form of ROOT or text files. The amount of information has become too large to inspect manually with a
text editor or the root browser. 


##Installation

###General Requirements

+ Python 2.6 or higher.
+ Flask 0.11 or higher.

###Flask

Flask can be installed in two different ways:

1) using [Flask Local and System-Wide installation guides](http://flask.pocoo.org/docs/installation/#installation).

2) or building and installing Flask System-Wide from the source:

   ~~~ sh
   $ cd libs/flask-0.11
   $ python setup.py build
   $ sudo python setup.py install
   ~~~

###Configuring the server

The server automatically configures it self to: "yourExternalIPAddress":5000.
To configure your web server manually go to [run.py](run.py) and change IP address and/or port. 

###Launching the server

To launch your web server go back to the DatacardViewer folder and type:

   ~~~ sh
   $ python run.py
   ~~~
If it doesn't work, go to [run.py](run.py) and configure it manually.

##Documentation

+ [Documentation](DOCUMENTATION.md)
+ [Tasks](TASKS.md)
