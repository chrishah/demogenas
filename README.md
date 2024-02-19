# demogenas
***democratizing genome assembly***

## Introduction

`demogenas` is a pipeline for denovo genome assembly and genome evaluation. It's written in Snakemake. 

In terms of dependencies you'll need:
 - `Snakemake` - best get via conda - note that we have been testing extensively only with Snakemake version 5.9.1 and we expect that some issues will arise with newer Snakemake versions.
 - `Singularity` - globally installed - 3.11.4 and newer should work. 

Futher you'll need to to clone this repository to get the context for the worklow and a number of scripts that will be used by the workflow.
```bash
git clone --recursive https://github.com/chrishah/demogenas.git
cd demogenas
```

The master wrappers script can be executed as follows. You're going to have to be in the repo to execute the script.
```bash
./demogenas
```

The wrapper script will call `snakemake` and trigger certain parts of the pipeline indicated via the `-m` option.

To run the pipeline successfully you are going to need at least 2 files:
 - data file - a tab delimited file specifying the location and type of data for your samples (example: `data/testdata/test.data.tsv`)
 - config file - specifying parameters for the pipeline (example: `data/testdata/test.config.yaml`)

The principal input data for now are:
 - Illumina paired end reads in fastq or bam format (specified in data file in the columns `f_fastq` and `r_fastq` or `bam`)
 - ONT reads in fast5 format (data file column 'fast5_dir')
 - long reads in fastq format (data file column 'long')

A given sample (data file column 'sample') can combine multiple datatypes specified in multiple lines with unique library names (data file column 'lib').

The example data file (`data/testdata/test.data.tsv`) specifies 5 hypothetical samples comprising different kinds of input data:
 - fastq_only - a sample for which multiple illumina libraries were sequenced and all data come in fastq format
 - bam_only - a sample for which multiple illumina libraries were sequenced and all data come in bam format
 - fastq_bam - a sample for which multiple illumina libraries were sequenced and data come as fastq and bam format
 - fast5_only - a sample for which multiple ONT libraries were sequenced and data comes as fast5
 - all_types - a sample for which multiple illumina and ONT libraries were sequenced and data come in fastq, bam and fast5 format

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

Here's the rulegraph of the workflow that would be exectued by the last command above.

<img src="https://github.com/chrishah/demogenas/blob/demo/rulegraph.all_types.pdf" height="500">


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
 - `-m eval_kmer_plot`
 - `-m kmer_filter`
 - `-m call_ont`

```bash
# trim illumina data
./demogenas -t local -m trim_illumina --configfile=data/testdata/test.config.yaml --dry --select=fastq_only
# trim and correct illumina data
./demogenas -t local -m correct_illumina --configfile=data/testdata/test.config.yaml --dry --select=fastq_bam
# trim, correct and merge illumina dat
./demogenas -t local -m merge_illumina --configfile=data/testdata/test.config.yaml --dry --select=fastq_only
# trim illumina data and evaluate (fastqc, kmer spectrum)
./demogenas -t local -m eval_kmer_plot --configfile=data/testdata/test.config.yaml --dry --select=fastq_bam
```

The latter command would just trim Illumina reads and calculate and plot kmer frequencies.

<img src="https://github.com/chrishah/demogenas/blob/demo/rulegraph.eval_kmer_plot.pdf" height="500">

Note that for a sample comprising only fast5 data none of the illumina specific steps will be done. You can try.
```bash
./demogenas -t local -m trim_illumina --configfile=data/testdata/test.config.yaml --dry --select=fast5_only
./demogenas -t local -m correct_illumina --configfile=data/testdata/test.config.yaml --dry --select=fast5_only
./demogenas -t local -m merge_illumina --configfile=data/testdata/test.config.yaml --dry --select=fast5_only
```

We have extra modes to evaluate your assemblies - this can be done at any time, even if not all assemblies are finished yet.

First, run `prepare_assemblies` mode - this will gather all assemblies that are finished at this stage (potentially restricted only to a single sample id via `--select=` option as below).
```bash
./demogenas -m prepare_assemblies -t local --configfile=data/testdata/test.config.yaml --select="fastq_bam"
```
Since this is just a demo and no assembly has actually been run yet the above command will not actually trigger any jobs. Snakemake will tell you that there's nothing to be done. However, if any assemblies for this particular sample had been completed the above command would have gathered them in a particular place now. Check out the content of this target directory:
```bash
ls -1 results/fastq_bam/assembly_evaluation/assemblies/
```
You'll see a list of files. These are just empty files now that ship with the repo for the purpose of this demo. Filenames as given by demogenas should be indicative of the origin of each file. The principal naming scheme is `<assembler-trimmer-correction-merger>.min<length>.fasta`. The file `platanus-trimgalore-bless-usearch-auto.min1000.fasta` for example was produced via platanus, based on reads trimmed with trimgalore, corrected with bless and merged with usesarch. The assembly has been filtered to retain only scaffolds of a minimum length of 1000bp. If one sticks to the principal naming scheme one can also put in external assemblies for subsequent evaluation, such as the file `something.min1000.fasta`. 


Now, to evaluate all assemblies finished at this moment with the methods as specified in the config file, run `evaluate_assemblies` like e.g. so:
```bash
./demogenas -m evaluate_assemblies -t local --configfile=data/testdata/test.config.yaml --select="fastq_bam" --dry
```


