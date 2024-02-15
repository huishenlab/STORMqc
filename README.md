# STORMqc
QC pipeline for STORM data

# Download Repository and Create Resources

This process only needs to be done once to get setup for all experiments.

## Download Repository

To download the STORMqc git repository from [GitHub](https://github.com/huishenlab/STORMqc), run
```
cd /path/to/where/you/want/the/qc/directory/to/go
git clone git@github.com:huishenlab/STORMqc.git <name of what you want your directory to be>
```

## Create Resources

### Requirements for creating resources

The following are the tools you will need access to in order to create the resources. Many of them are likely to be
available on your typical unix/linux distribution, but some may need to be installed or acquired yourself (e.g.,
samtools, bedtools, and STAR).

  - wget
  - python
  - samtools
  - GNU awk
  - bedtools
  - STAR (loaded via linux module)

### Submitting the resource creator script

Say you've named your directory `storm_fun`, your first step in `storm_fun` is to create the reference resources needed
for performing the scRNA-seq alignments:
```
cd storm_fun/resources
sbatch gather_rna_resources.slurm
```
This process will take an hour or two and will prepare both human (hg38) and mouse (mm10) resources. Right now, the
pipeline is not setup to actually process mouse data, but the resources are there for easing the transition process once
mouse capabilities are added in the future.

**Note, the resources must be created before you can proceed to the next step.**

# Pipeline Overview
