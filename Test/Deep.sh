#!/bin/sh
/bin/mkdir -p latestjobs histo
for ch in  elnu
do
  	for ds in tzq
        do
          	chds=${ch}_${ds}
                /bin/cat <<EOF > latestjobs/${chds}.sh
#!/bin/bash
export X509_USER_PROXY=x509up_u527
source ~/.bashrc
source /cvmfs/cms.cern.ch/cmsset_default.sh
source /cvmfs/grid.cern.ch/etc/profile.d/setup-cvmfs-ui.sh
./build/dptest calchisto ${ch} ${ds}
EOF
                /bin/cat <<EOF > latestjobs/${chds}.job
Executable    = /bin/bash
universe      = vanilla
Output        = latestjobs/dt_${chds}.out
Error         = latestjobs/dt_${chds}.err
Log           = latestjobs/dt_${chds}.log
Arguments     = latestjobs/${chds}.sh
Queue
EOF
                condor_submit latestjobs/${chds}.job
        done
done

