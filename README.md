# NIPT-Thalassemia
NIPT-Thalassemia

# pipeline of parper

## step0
get ref panel and copy to `db/` .
get samples' vcf and copy to `workdir/input` .

## step1
```
./bin/parallelRun -list script/run_step1.sh
```

## step2
```
./bin/parallelRun -list script/run_step2.sh
```

## step4
```
./bin/parallelRun -list script/run_step3.sh
```

## step4
```
./bin/parallelRun -list script/run_step4.sh
```

## create table1
```
sh script/create.table1.sh
```

## create table2
```
sh script/create.table2.sh
```

## plot figure S2
```
sh script/plot.figureS2.sh
```

## plot figure S3
```
sh script/plot.figureS3.sh
```
