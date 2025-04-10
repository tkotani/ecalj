python --version
python ./auto/jobsubmit.py --inpath INPUT/testSGA --epath ${HOME}/bin/ --nqsub 2 --ncore 16 --bnd4all True --niter 2 --kkmesh 6 6 6 3 3 3 
# python ./auto/jobsubmit.py --inpath INPUT/testSGA --epath /home/k0413/k041303/bin/ --nqsub 1 --ncore 128 --bnd4all True --niter 2
#        auto/jobsubmit.py --inpath INPUT/testSGA --epath /home/stakano/ecalj/bin/ --nsub 2 --niter 2 --ncore 32 --bnd4all True --gw80 False --koption 6 8 10
