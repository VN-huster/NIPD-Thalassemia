# NIPT-Thalassemia
NIPT-Thalassemia

# pipeline of parper

## Requirement
* perl
* R
* bgzip
* tabix
* vcf-merge (from [VCFtools](https://vcftools.github.io/perl_module.html))
* beagle

## Step0: prepare data
which not in git repo:
* `copy workdir/input/F*.vcf.gz`
* `copy workdir/input/F*.vcf.gz.tbi`
* `copy workdir/input/F*.cnv.list`
* `copy bin/`
* `copy db/`

create output dir:  
```
mkdir -p workdir/output
```

## All Step: run all steps in one shell
```
sh script/run.sh
```

## Step1: phasing parents
```
sh script/run_step1.sh
```
or 
```
./bin/parallelRun -list script/run_step1.sh
```
for parallel run

## Step2: NIPD
```
sh script/run_step2.sh
```
or 
```
./bin/parallelRun -list script/run_step2.sh
```
for parallel run

## Step3
```
sh script/run_step3.sh
```
or 
```
./bin/parallelRun -list script/run_step3.sh
```
for parallel run

## Step4
```
sh script/run_step4.sh
```
or 
```
./bin/parallelRun -list script/run_step4.sh
```
for parallel run

## Create table1
```
sh script/create.table1.sh
```

## Create table3
```
sh script/create.table3.sh
```

## Plot figure S2
after *step3*, we can plot figure S2 use:
```
sh script/plot.figureS2.sh
```
or
```
Rscript script/fetus_hap_consistency.alpha.R workdir/input/alpha.list workdir/output
Rscript script/fetus_hap_consistency.beta.R workdir/input/beta.list workdir/output
```

## Plot figure S3
```
sh script/plot.figureS3.sh
```
