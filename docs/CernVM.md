## Standalone use inside CernVM

Standalone by adding the `CMS` group to the CVMFS Configuration.
A minimal `CernVM` working context setup can be found in the CernVM Marketplace
under `Experimental/HiggsCombine` or at
[https://cernvm-online.cern.ch/context/view/9ee5960ce4b143f5829e72bbbb26d382](https://cernvm-online.cern.ch/context/view/9ee5960ce4b143f5829e72bbbb26d382).
At least 2GB of disk space should be reserved on the virtual machine for
Combine to work properly.

# Available machines for standalone combine

The standalone version can be easily compiled via _CVMFS_ as it relies on dependencies which are already installed at _/cvmfs/cms.cern.ch/_. Access to _/cvmfs/cms.cern.ch/_ can be obtained from lxplus machines or via [`CernVM`](https://cernvm-online.cern.ch/). The only requirement will be to add the _CMS_ group to the CVMFS configuration as shown in the picture

![](cvmsf_config.png)

At least 2GB of disk space should be reserved on the virtual machine for combine to work properly. A minimal CernVM working context setup can be found in the CernVM Marketplace under [`Experimental/HiggsCombine`](https://cernvm-online.cern.ch/context/view/9ee5960ce4b143f5829e72bbbb26d382). 

To use this predefined context, first locally launch the CernVM (eg you can use the .ova with VirtualBox, by downloading from [here](http://cernvm.cern.ch/releases/production/cernvm4-micro-2020.07-1.ova) and launching the downloaded file. You can click on "pair an instance of CernVM" from the cernvm-online dashboard, which displays a PIN. In the VirtualBox terminal, pair the virtual machine with this PIN code (enter in the terminal using #PIN eg `#123456`. After this, you will be asked again for username (use `user`) and then a password (use `hcomb`).

In case you do not want to use the cvmfs area, you will need to adapt the location of the dependencies listed in both the Makefile and env_standalone.sh files.
