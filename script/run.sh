mkdir -p workdir/output
./bin/parallelRun -list script/run_step1.sh  2>&1|tee log.1
./bin/parallelRun -list script/run_step2.sh  2>&1|tee log.2
./bin/parallelRun -list script/run_step3.sh  2>&1|tee log.3
./bin/parallelRun -list script/run_step4.sh  2>&1|tee log.4
