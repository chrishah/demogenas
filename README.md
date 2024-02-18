# demogenas
democratizing genome assembly

Usage:

You are going to need at least 2 file:
 - data file - a tab delimited file specifying the location and type of data for your samples (example: `data/testdata/test.data.tsv`)
 - config file - specifying parameters for the pipeline (example: `data/testdata/test.config.yaml`)

The example data file (`data/testdata/test.data.tsv`) specifies 5 hypothetical samples comprising different kinds of input data:
 - fastq_only - a sample for which multiple libraries were sequenced and all data come in fastq format
 - bam_only - a sample for which multiple libraries were sequenced and all data come in bam format
 - fastq_bam - a sample for which multiple libraries were sequenced and data come as fastq and bam format
 - fast5_only - a sample for which multiple libraries were sequenced and data comes as fast5
 - all_types - a sample for which multiple libraries were sequenced and data come in fastq, bam and fast5 format

Per default, `demogenas` will process (trim, errorcorrect, merge) and assemble all samples and datatypes in the datafile automatically, with the particular steps (trimmers, correctors, assemblers) as specified in the config file.

```bash
./demogenas -t local -m assemble --configfile=data/testdata/test.config.yaml --dry
``` 

If you want to process only selected samples in the data file you can specify the sample(s) of interest with the `--select=` option, e.g.:
```bash
# assemble one sample 'fastq_only'
./demogenas -t local -m assemble --configfile=data/testdata/test.config.yaml --dry --select=fastq_only

# assemble two samples 'fastq_only' and 'bam_only'
./demogenas -t local -m assemble --configfile=data/testdata/test.config.yaml --dry --select=fastq_only,bam_only
```

You can check the different behaviours of the pipeline for different sample/data types - note that `demogenas` will do different things depending on the data types provided, e.g.:
```bash
# only fastq
./demogenas -t local -m assemble --configfile=data/testdata/test.config.yaml --dry --select=fastq_only
# only bam
./demogenas -t local -m assemble --configfile=data/testdata/test.config.yaml --dry --select=bam_only
# combination of fastq and bam
./demogenas -t local -m assemble --configfile=data/testdata/test.config.yaml --dry --select=fastq_bam
# fast5 only
./demogenas -t local -m assemble --configfile=data/testdata/test.config.yaml --dry --select=fast5_only
# combination of fastq, bam and fast5
./demogenas -t local -m assemble --configfile=data/testdata/test.config.yaml --dry --select=all_types
```

Via the config file I control which steps are performed, e.g.:
```bash
# process illumina data and assemble with all relevant assemblers
./demogenas -t local -m assemble --configfile=data/testdata/test.config.yaml --dry --select=fastq_only
# process illumina data and assemble only with spades
./demogenas -t local -m assemble --configfile=data/testdata/test.config.spadesonly.yaml --dry --select=fastq_only
```

If you don't want to go all the way and assemble, there are other modes, such as:
 - `-m trim_illumina`
 - `-m correct_illumina`
 - `-m merge_illumina`
 - `-m eval_illumina`

```bash
# trim illumina data
./demogenas -t local -m trim_illumina --configfile=data/testdata/test.config.yaml --dry --select=fastq_only
# trim illumina data and evaluate (fastqc, kmer spectrum)
./demogenas -t local -m eval_illumina --configfile=data/testdata/test.config.yaml --dry --select=bam_only
# trim and correct illumina data
./demogenas -t local -m correct_illumina --configfile=data/testdata/test.config.yaml --dry --select=fastq_bam
# trim, correct and merge illumina dat
./demogenas -t local -m merge_illumina --configfile=data/testdata/test.config.yaml --dry --select=fastq_only
```

Note that for a sample comprising only fast5 data none of the illumina specific steps will be done. You can try.
```bash
./demogenas -t local -m trim_illumina --configfile=data/testdata/test.config.yaml --dry --select=fast5_only
./demogenas -t local -m correct_illumina --configfile=data/testdata/test.config.yaml --dry --select=fast5_only
./demogenas -t local -m merge_illumina --configfile=data/testdata/test.config.yaml --dry --select=fast5_only
```

We have extra modes to evaluate your assemblies - this can be done at any time, even if not all assemblies are finished yet.

First, run `prepare_assmblies` mode - this will gather all assemblies that are finished at this stage.
```bash
./demogenas -m prepare_assemblies -t local --configfile=data/testdata/test.config.yaml --select="fastq_bam"
```

Now, to evaluate with the methods as specified in the config file, run `evaluate_assemblies` like e.g. so:
```bash
./demogenas -m evaluate_assemblies -t local --configfile=data/testdata/test.config.yaml --select="fastq_bam"
```

