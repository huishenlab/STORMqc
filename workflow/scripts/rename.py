import pandas as pd
import sys
import os

def create_symlink(out_dir, fastq_dir, fq, samp, n_read):
    if len(fq.split(',')) > 1:
        print('WARNING: More than one FASTQ file given. Taking only the first file!')
        fq = fq.split(',')[0]

    fq_path = os.path.join(fastq_dir, fq)
    if os.path.exists(fq_path):
        new_name = os.path.join(out_dir, f'{samp}_R{n_read}.fastq.gz')
        os.symlink(fq_path, new_name)
        print(f'{fq_path} successfully symlinked to {new_name}')
    else:
        print(f'ERROR: {fq_path} not found in {fastq_dir}')
        sys.exit(1)

def rename_fastqs(samplesheet, fastq_dir, out_dir):
    samples = pd.read_csv(samplesheet, sep='\t')

    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    for idx, row in samples.iterrows():
        samp = row['sample']
        fq1 = row['fq1']
        fq2 = row['fq2']

        print(f'Renaming paired-end reads for: {samp}')

        create_symlink(out_dir, fastq_dir, fq1, samp, 1) # read 1
        create_symlink(out_dir, fastq_dir, fq2, samp, 2) # read 2

with open(snakemake.log[0], 'w') as fh:
    sys.stderr = sys.stdout = fh
    rename_fastqs(
        snakemake.params['samplesheet'],
        snakemake.params['fastq_dir'],
        snakemake.output['sym_dir'],
    )
