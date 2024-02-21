# STORMqc

## Pipeline Overview

The STORMqc pipeline aims to provide a way to verify the quality of a STORM-seq library preparation and sequencing run.
It currently works for human data with the intention of adding mouse capabilities in the future. Broadly speaking, the
pipeline runs:

  - FastQC on both the read 1 and read 2 FASTQ files to get basic read quality metrics
  - Aligns the FASTQ files with STAR ("regular" STAR, not STARsolo) for getting counts to undesirable mapping locations
  (rRNA, mtDNA, ERCC spike-ins, and non-annotated genomic space)
  - `samtools stats` and `samtools flagstat` for alignment metrics
  - MultiQC is run over the processed files to collate everything
  - Specific metrics (including, but not limited to, average base quality, number of sequences, and unique mapping
  percentage) are plotted in a 384-well plate layout for visualization of spatially correlated failures

The general steps for running the pipeline are:

  1. Download repository and create resources (only needs to be done once)
  2. After resources have been created, create files to run pipeline
  3. Run pipeline

More details are provided below on each of these steps.

## Download Repository and Create Resources

This process only needs to be done once to get setup for all experiments.

### Download Repository

To download the STORMqc git repository from [GitHub](https://github.com/huishenlab/STORMqc), run
```
cd /path/to/where/you/want/the/qc/directory/to/go
git clone git@github.com:huishenlab/STORMqc.git <name of what you want your directory to be>
```

### Create Resources

**Note, the resources must be created before you can start creating files to run the pipeline.**

#### Requirements for creating resources

The following are the tools you will need access to in order to create the resources. Many of them are likely to be
available on your typical unix/linux distribution, but some may need to be installed or acquired yourself (e.g.,
samtools, bedtools, and STAR).

  - wget
  - python
  - samtools
  - GNU awk
  - bedtools
  - STAR (loaded via linux module)

#### Submitting the resource creator script

Say you've named your directory `storm_fun`, your first step in `storm_fun` is to create the reference resources needed
for performing the scRNA-seq alignments:
```
cd storm_fun/resources
sbatch gather_rna_resources.slurm
```
This process will take an hour or two and will prepare both human (hg38) and mouse (mm10) resources. Right now, the
pipeline is not setup to actually process mouse data, but the resources are there for easing the transition process once
mouse capabilities are added in the future.

#### Resources created / downloaded

  - FASTA and GTF for ERCC spike-in sequences
  - Human files
    - GENCODE primary assembly FASTA and GTF (GRCh38)
    - FANTOM enhancer location BED (hg38)
    - RepeatMasker repeat sequence location BED (hg38)
  - Mouse files
    - GENCODE primary assembly FASTA and GTF (GRCm39)
    - FANTOM enhancer location BED (mm10)
    - RepeatMasker repeat sequence location BED (mm10)
  - Merged files
    - Human + ERCC FASTA and GTF
    - Human annotated space BED (Human GTF annotation locations + ERCC GTF annotation locations + enhancers + repeats)
    - Mouse + ERCC FASTA and GTF
    - Mouse annotated space BED (Mouse GTF annotation locations + ERCC GTF annotation locations + enhancers + repeats)
  - STAR indexes
    - Merged human + ERCC with merged human GTF
    - Merged mouse + ERCC with merged mouse GTF
  - Gene ID lists
    - rRNA gene IDs (mouse and human)
    - mtDNA gene IDs (mouse and human)
    - ERCC gene IDs (mouse and human)

## Create Pipeline Input Files

Once you've the resource files, you'll be able to start creating input files for the pipeline. If you haven't created
the resources and try running the Python script to create the input files, it will produce an error message reminding
you to create said resources.

Whenever you get a new dataset, you'll want to first store those in a directory that is accessible to the STORMqc
pipeline (i.e., you won't be able to store the new data on a remote server and run the pipeline on your local machine).
Once you've downloaded the data, you'll want to create (or modify) the input file for the Python script.

The input file is a TSV with two columns, an experiment name column and the full path to the FASTQs you just downloaded.
An example of the contents of an input file is:
```
test_data	/path/to/my/STORM/fastqs
test_data2	/path/to/my/newer/STORM/fastqs
```
Notice that there are two entries shown. You can have 1 or more entries in your input file, where each entry is a unique
experiment ID and directory path. This gives you flexibility to either process multiple datasets in one go or to have
one input file that you add new rows to each time you get more data.

For each row in a file, the Python script will check to make sure a directory with the experiment ID doesn't already
exist. If it does, the script will default to asking you to confirm you want to overwrite the pipeline files in that
directory. If you confirm that you want them overwritten, only the pipeline files will be overwritten, not any processed
files generated by the Snakemake pipeline. You can also specify from the command line to always overwrite or never
overwrite any reoccurring experimental IDs. One final thing to note, only the experiment IDs are required to be unique.
Therefore, if you accidentally repeat the same FASTQ directory with different experiment IDs, the script won't catch it.

To create the pipeline files, run these commands (starting from the top directory of your pipeline):
```
cd bin
python make_pipeline_files.py \
    [-m|--min-read-count 100] \
    [-O|--overwrite ask/always/never] \
    [-V|--visualize-config /path/to/viz_config.yaml] \
    input_file.tsv
```
The `-m|--min-read-count` option allows you to adjust the minimum number of reads needed in a FASTQ to be included when
running the pipeline. The default is 100 reads. Note, the code will assume that a file larger than 10 kB is large
enough to be included in the pipeline. This is a reasonable size if the minimum read count is 100, but if you increase
that number, you may have additional wells processed that you didn't want. The `-O|--overwrite` option provides the
ability to toggle how to handle repeated experiment IDs. The `-V|--visualize-config` option allows you to provide your
own YAML file for controlling the y-axis limits of the quality control (QC) plots. The default YAML file can be found in
`bin/default_visualization_limits.yaml`. If you plan to change the QC plot limits, it is highly suggested that you make
a copy of `bin/default_visualization_limits.yaml` and use this new file as your input to `make_pipeline_files.py`. This
will allow you to maintain the default values while still having a file that you can adjust as needed.

The generator script will create all directories and files needed for running the pipeline. These will be placed in the
`results` directory:
```
results
|- experiment_id_1
|  |- analysis
|  |  |- L001
|  |  |- L002
|  |  |- ...
|  |- benchmark
|  |  |- L001
|  |  |- L002
|  |  |- ...
|  |- log
|  |  |- L001
|  |  |- L002
|  |  |- ...
|  |- pipeline_files
|  |  |- L001_config.yaml
|  |  |- L001_rerun.slurm
|  |  |- L001_samples.tsv
|  |  |- L001_submit.slurm
|  |  |- L001_too_few_reads.txt
|  |  |- ...
|- experiment_id_2
|  |- ...
```

Regarding the contents of each experiment directory, the `analysis` directory is where all processed outputs from the
pipeline will be stored. The `benchmark` directory includes time and memory statistics from the different jobs that are
run in the pipeline, while the `log` directory is for any error output when running the pipeline. The script will create
subdirectories in each of the output directories and run files for each lane found in the data directory (given by
`L001`/`L002`/etc. in the read name). The run files are as follows: `*_config.yaml` is the configuration file for
running the pipeline, `*_samples.tsv` is the samplesheet based on the FASTQ file name, `*_submit.slurm` is the script
that is submitted to the SLURM job scheduler for running the whole pipeline, `*_rerun.slurm` is the script that is
submitted if you want to only rerun making the QC plots, and `*_too_few_reads.txt` is a list of wells that weren't
included in the samplesheet because they didn't meet the minimum read count requirement in at least one of the FASTQ
files for that well.

## Running the Pipeline

Now that you've created the pipeline files, it's simply a matter of submitting the slurm scripts to the scheduler and
sitting back while the jobs run. You can do that by running (from the top directory of the pipeline):
```
sbatch results/experiment_id_1/pipeline_files/L001_submit.slurm
```

If you have multiple lanes and/or multiple experiments, you can submit multiple scripts at one time, but note that you
will be competing with yourself for resources, so I don't suggest submitting more than a couple at a time.
```
sbatch results/experiment_id_1/pipeline_files/L001_submit.slurm
sbatch results/experiment_id_1/pipeline_files/L002_submit.slurm
```

### Rerunning Only the Quality Control Plot Rule

There may be times when you want to regenerate the QC plots with different y-axis limits. For these instances, there is
a SLURM script that is provided that will only rerun the rule for generating the QC plots. After adjusting your limits
YAML config file, run (or similar for other lanes):
```
sbatch results/experiment_id_1/pipeline_files/L001_rerun.slurm
```
Assuming you haven't deleted any files that were created by the pipeline, this will only regenerate the PDF of the QC
plots.

## Processed Output from Pipeline

As mentioned previously all processed outputs form the pipeline can be found in the `results/experiment_id/analysis`
directory. Inside this directory, you will see a directory for each lane of data (`L001`/`L002`/etc.) where each of
these directories contain the following (`*` represents that every sample will have a corresponding file):
```
L001
|- fastqc
|  |- *_R1_fastqc.html
|  |- *_R1_fastqc.zip
|  |- *_R2_fastqc.html
|  |- *_R2_fastqc.zip
|- multiqc
|  |- multiqc_report_data/
|  |- multiqc_report.html
|- plots
|  |- platetools_data.tsv
|  |- quality_control_plots.pdf
|- renamed_fastqs
|  |- *_R1.fastq.gz
|  |- *_R2.fastq.gz
|- samtools
|  |- *.stats
|  |- *.flagstat
|- star
|  |- *.non_annotated.tsv
|  |- *Log.final.out
|  |- *Log.out
|  |- *Log.progress.out
|  |- *ReadsPerGene.out.tab
|  |- *SJ.out.tab
```

### Processed Output Descriptions

  - `fastqc`: FastQC output files for reads 1 and 2 for each sample ID processed
  - `multiqc`: output from MultiQC
  - `plots`: additional figures for quick checks of STORM-seq run
  - `renamed_fastqs`: symlinks to raw data - provides consistent naming scheme for Snakemake
  - `samtools`: stats and flagstat output from `samtools`
  - `star`: log and counts files from STAR, also includes read counts aligned to non-annotated space (non-annotated
  space considered as genomic space that does *not* fall in the annotation GTF (both ERCC and your species of interest),
  RepeatMasker defined space, and FANTOM5 enhancer defined space.
