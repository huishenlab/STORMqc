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

# TODO: Refactor to minimize repetition of code
def get_renamed_fastqs(wildcards):
    cp_output = checkpoints.rename_fastq_files.get().output.sym_dir

    return list(expand(cp_output + '/' + wildcards.sample + '_R{read}.fastq.gz', read = [1, 2]))

def get_renamed_read1(wildcards):
    cp_output = checkpoints.rename_fastq_files.get().output.sym_dir

    return cp_output + '/' + wildcards.sample + '_R1.fastq.gz'

def get_renamed_read2(wildcards):
    cp_output = checkpoints.rename_fastq_files.get().output.sym_dir

    return cp_output + '/' + wildcards.sample + '_R2.fastq.gz'

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
