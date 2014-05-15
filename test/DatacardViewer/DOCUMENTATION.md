#DatacardViewer technical documentation

The technical documentation reflects the DatacardViewer technology stack, architecture, structure and documentation.

# Table of Contents

* [Overview](#overview)
  * [About](#about)
  * [Technology stack](#technology-stack)
  * [Architecture and structure](#architecture-and-structure)
* [Documentation](#documentation)


#Overview

##About

The Higgs searches and measurements are encoded in complex likelihood functions that comprise
thousands parameters. For convenience of extracting results, this information is stored in ROOT
(RooFit) format. On top of that there is a layer of ASCII configuration files describing the contents of the
ROOT files, which is parsed using Python. Diagnostic information about these inputs is also produced in
form of ROOT or text files. The amount of information has become too large to inspect manually with a
text editor or the root browser. 

Datacard Viewer is a Graphical User Interface tool for exploration of Higgs measurements in datacards.
[More information](presentation/).

##Technology stack

Libraries used to create this tool and links.

**Server side**:
+ [Flask](http://flask.pocoo.org/docs/)

**Client side**:
+ [AngularJS](https://angularjs.org/)
+ [Bootstrap](http://getbootstrap.com/)
+ [Handsontable](http://handsontable.com/)
+ [jQuery](http://jquery.com/)
+ [JSRootIO](http://root.cern.ch/gitweb?p=rootjs.git;a=summary)

##Architecture and structure

**DatacardViewer Architecture model:**
![Alt text](presentation/DatacardViewer_schema.png?raw=true "DatacardViewer Architecture")

**DatacardViewer structure:**
```
app
|-- app.py
|-- static
    |-- css
    |-- img
    |-- js
        |-- app.js, controllers.js, etc.
    |-- lib
        |-- angular
            |-- angular.js, etc.
    |-- partials
|-- templates
    |-- index.html
```
This structure is a mix of Flask and AngularJS structures.

#Documentation

+ [Main javascript documentation for future developers](documentation/)
