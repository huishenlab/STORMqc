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

checkpoint sortmerna_assets:
    output:
        db_dir = directory(f'{ANALYSIS}/sortmerna_database')
    log:
        f'{ANALYSIS}/sortmerna_database/sortmerna_database.log'
    threads: 1
    resources:
        mem_gb = config['hpc_parameters']['memory']['small'],
        time = config['hpc_parameters']['runtime']['short'],
    conda:
        '../envs/sortmerna.yaml'
    shell:
        '''
        mkdir -p {output.db_dir}
        wget \
            --directory-prefix {output.db_dir} \
            https://github.com/biocore/sortmerna/releases/download/v4.3.4/database.tar.gz
        tar -xvf {output.db_dir}/database.tar.gz -C {output.db_dir}
        rm -f {output.db_dir}/database.tar.gz

        sortmerna \
            --index 1 \
            --workdir {output.db_dir} \
            --threads {threads} \
            -ref {output.db_dir}/smr_v4.3_fast_db.fasta 
        '''

def get_sortmerna_db(wildcards):
    dir_loc = checkpoints.sortmerna_assets.get().output.db_dir
    return f'{dir_loc}/smr_v4.3_fast_db.fasta'

def get_sortmerna_db_idx(wildcards):
    dir_loc = checkpoints.sortmerna_assets.get().output.db_dir
    return f'{dir_loc}/idx'

rule sortmerna_run:
    input:
        fq1 = get_renamed_read1,
        fq2 = get_renamed_read2,
        db = get_sortmerna_db,
        db_idx = get_sortmerna_db_idx,
    output:
        f'{ANALYSIS}/sortmerna/{{sample}}/{{sample}}.log',
        f'{ANALYSIS}/sortmerna/{{sample}}/{{sample}}.sam.gz',
        work_dir = directory(f'{ANALYSIS}/sortmerna/{{sample}}'),
    params:
        sample = '{sample}',
    log:
        stdout = f'{LOG}/sortmerna/{{sample}}.out',
        stderr = f'{LOG}/sortmerna/{{sample}}.err',
    threads: config['hpc_parameters']['threads']['sortmerna'],
    resources:
        mem_gb = config['hpc_parameters']['memory']['small'],
        time = config['hpc_parameters']['runtime']['short'],
    conda:
        '../envs/sortmerna.yaml'
    shell:
        '''
        sortmerna \
            --workdir {output.work_dir} \
            --aligned {output.work_dir}/{params.sample} \
            --sam \
            --threads {threads} \
            --idx-dir {input.db_idx} \
            -ref {input.db} \
            --reads {input.fq1} \
            --reads {input.fq2} 2> {log.stderr} 1> {log.stdout}
        '''
