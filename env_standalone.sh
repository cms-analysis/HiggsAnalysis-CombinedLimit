 source /cvmfs/cms.cern.ch/slc7_amd64_gcc820/external/gcc/8.2.0-pafccj/etc/profile.d/init.sh
 source /cvmfs/cms.cern.ch/slc7_amd64_gcc820/external/gsl/2.2.1-pafccj/etc/profile.d/init.sh
 source /cvmfs/cms.cern.ch/slc7_amd64_gcc820/external/boost/1.67.0-pafccj/etc/profile.d/init.sh
 source /cvmfs/cms.cern.ch/slc7_amd64_gcc820/external/xz/5.2.2-pafccj/etc/profile.d/init.sh
 source /cvmfs/cms.cern.ch/slc7_amd64_gcc820/external/tbb/2019_U3-pafccj/etc/profile.d/init.sh
 source /cvmfs/cms.cern.ch/slc7_amd64_gcc820/lcg/root/6.14.09-nmpfii2/bin/thisroot.sh

 export LD_LIBRARY_PATH=/cvmfs/cms.cern.ch/slc7_amd64_gcc820/cms/vdt/0.4.0-pafccj/lib:/cvmfs/cms.cern.ch/slc7_amd64_gcc820/external/pcre/8.37/lib/:${LD_LIBRARY_PATH}
 export PATH=${PATH}:${PWD}/exe:${PWD}/scripts
 export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${PWD}/lib
 export PYTHONPATH=${PYTHONPATH}:${PWD}/lib/python:${PWD}/lib

