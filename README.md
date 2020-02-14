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
get ref panel and copy to `db/` .   
get samples' vcfs and copy to `workdir/input` .  
create output dir:  
```
mkdir -p workdir/output
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

## Step2
```
./bin/parallelRun -list script/run_step2.sh
```

## Step3
```
./bin/parallelRun -list script/run_step3.sh
```

## Step4
```
./bin/parallelRun -list script/run_step4.sh
```

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
