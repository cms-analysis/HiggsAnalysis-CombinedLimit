# equivalent to CMSSW 14_1_X
. /cvmfs/cms.cern.ch/el9_amd64_gcc12/external/gcc/12.3.1-40d504be6370b5a30e3947a6e575ca28/etc/profile.d/init.sh
. /cvmfs/cms.cern.ch/el9_amd64_gcc12/lcg/root/6.30.07-8b1a11e1ef0e074fdfd44e162b27e71c/etc/profile.d/init.sh
# this also provides py3-six
. /cvmfs/cms.cern.ch/el9_amd64_gcc12/external/py3-pandas/2.0.1-7e0980e2f417a1fc2106b953723ea098/etc/profile.d/init.sh
. /cvmfs/cms.cern.ch/el9_amd64_gcc12/external/gsl/2.6-5e2ce72ea2977ff21a2344bbb52daf5c/etc/profile.d/init.sh
. /cvmfs/cms.cern.ch/el9_amd64_gcc12/external/tbb/v2021.9.0-7e31cf78e4a7495ff0337c859d1dc314/etc/profile.d/init.sh
. /cvmfs/cms.cern.ch/el9_amd64_gcc12/cms/vdt/0.4.3-793cee1e1edef0e54b2bd5cb1f69aec9/etc/profile.d/init.sh
. /cvmfs/cms.cern.ch/el9_amd64_gcc12/external/boost/1.80.0-87b5de10acd2f2c8a325345ad058b814/etc/profile.d/init.sh
. /cvmfs/cms.cern.ch/el9_amd64_gcc12/external/pcre/8.43-e34796d17981e9b6d174328c69446455/etc/profile.d/init.sh
. /cvmfs/cms.cern.ch/el9_amd64_gcc12/external/eigen/3bb6a48d8c171cf20b5f8e48bfb4e424fbd4f79e-3ca740c03e68b1a067f3ed0679234a78/etc/profile.d/init.sh
 export PATH=${PATH}:${PWD}/build/bin
 export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${PWD}/build/lib
 export PYTHONPATH=${PYTHONPATH}:${PWD}/build/lib/python
