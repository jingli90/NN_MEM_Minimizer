# NN_MEM_Minimizer
# branch anneal_perm

Download and setup MEM according to https://twiki.cern.ch/twiki/bin/viewauth/CMS/IPHCMEMCPP

Download the code
git clone -b anneal_perm https://github.com/jingli90/NN_MEM_Minimizer.git

Change the include path in Makefile and test.cpp to the MEM you installed
Change the path in test.cpp to the path of your path of config.cfg, weightInit.txt, and root files 
Change the path in test/test_Batch/RunBatch_test according to your working directory

Change the options in test.cpp if you want

make

cd test

./test

If you want to launch the jobs, go to test/test_Batch/
bsub -q 1nw -N RunBatch_test
