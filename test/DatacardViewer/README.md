#Datacard Viewer


Datacard Viewer is a Graphical User Interface tool for exploration of Higgs measurements in datacards.

##How to use it

The tool works like this:

+ Install Flask server on network.
+ Type server's address in a browser.
+ Upload datacards and their ROOT files you want view.
+ View them.

##Installation

###General Requirements

Server side (libraries):
+ Python 2.6 or higher.
+ Flask 0.11 or higher.

Client side (browser):
+ Safari, Chrome, Firefox, Opera.
+ IE8, IE9 or newer.

###Flask

Flask can be installed in two different ways:

1) using [Flask Local and System-Wide installation guides](http://flask.pocoo.org/docs/installation/#installation).

2) or building and installing System-Wide Flask from the source:

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
