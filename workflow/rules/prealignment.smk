#----------------------------------------------------------------------------------------------------------------------
# Notes on where variables are defined, if not defined in prealignment.smk
#
# LOG       - workflow/Snakefile
# ANALYSIS  - workflow/Snakefile
# BENCHMARK - workflow/Snakefile
# config    - workflow/Snakefile
#
#----------------------------------------------------------------------------------------------------------------------

checkpoint rename_fastq_files:
    output:
        sym_dir = directory(f'{ANALYSIS}/renamed_fastqs'),
    params:
        samplesheet = config['samplesheet'],
        fastq_dir = config['fastq_dir'],
    log:
        f'{LOG}/rename/rename.log',
    threads: 1,
    resources:
        mem_gb = config['hpc_parameters']['memory']['small'],
        time = config['hpc_parameters']['runtime']['short'],
    conda:
        '../envs/python.yaml'
    script:
        '../scripts/rename.py'

def get_renamed_fastqs(wildcards):
    cp_output = checkpoints.rename_fastq_files.get().output.sym_dir

    return list(expand(cp_output + '/' + wildcards.sample + '_R{read}.fastq.gz', read = [1, 2]))

rule fastqc:
    input:
        get_renamed_fastqs
    output:
        f'{ANALYSIS}/fastqc/{{sample}}_R1_fastqc.html',
        f'{ANALYSIS}/fastqc/{{sample}}_R1_fastqc.zip',
        f'{ANALYSIS}/fastqc/{{sample}}_R2_fastqc.html',
        f'{ANALYSIS}/fastqc/{{sample}}_R2_fastqc.zip',
    params:
        dir = f'{ANALYSIS}/fastqc',
    log:
        stdout = f'{LOG}/fastqc/{{sample}}.out',
        stderr = f'{LOG}/fastqc/{{sample}}.err',
    benchmark:
        f'{BENCHMARK}/fastqc/{{sample}}.txt'
    threads: config['hpc_parameters']['threads']['fastqc'],
    resources:
        mem_gb = config['hpc_parameters']['memory']['small'],
        time = config['hpc_parameters']['runtime']['short'],
    conda:
        '../envs/babraham.yaml',
    shell:
        '''
        mkdir -p {params.dir}
        fastqc --outdir {params.dir} --threads {threads} {input} 2> {log.stderr} 1> {log.stdout}
        '''

rule trim_reads:
    input:
        get_renamed_fastqs
    output:
        f'{ANALYSIS}/trim_reads/{{sample}}_R1_val_1.fq.gz',
        f'{ANALYSIS}/trim_reads/{{sample}}_R2_val_2.fq.gz',
    params:
        outdir = f'{ANALYSIS}/trim_reads',
    log:
        stdout = f'{LOG}/trim_reads/{{sample}}.out',
        stderr = f'{LOG}/trim_reads/{{sample}}.err',
    benchmark:
        f'{BENCHMARK}/trim_reads/{{sample}}.txt',
    conda:
        '../envs/babraham.yaml'
    threads: config['hpc_parameters']['threads']['trim']
    resources:
        mem_mb = config['hpc_parameters']['memory']['small'],
        time = config['hpc_parameters']['runtime']['short'],
    shell:
        """
        trim_galore \
            --output_dir {params.outdir} \
            --cores {threads} \
            --illumina \
            --length 36 \
            --paired \
            --fastqc \
            {input} \
            2> {log.stderr} 1> {log.stdout}
        """
