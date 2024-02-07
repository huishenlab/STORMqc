#----------------------------------------------------------------------------------------------------------------------
# Notes on where variables are defined, if not defined in prealignment.smk
#
# LOG       - workflow/Snakefile
# ANALYSIS  - workflow/Snakefile
# BENCHMARK - workflow/Snakefile
# config    - workflow/Snakefile
# get_renamed_fastqs - workflow/rules/prealignment.smk
#
#----------------------------------------------------------------------------------------------------------------------

rule star_align:
    input:
        get_renamed_fastqs
    output:
        f'{ANALYSIS}/star/{{sample}}Log.final.out',
        f'{ANALYSIS}/star/{{sample}}Log.out',
        f'{ANALYSIS}/star/{{sample}}Log.progress.out',
        f'{ANALYSIS}/star/{{sample}}ReadsPerGene.out.tab',
        f'{ANALYSIS}/star/{{sample}}SJ.out.tab',
        bam = f'{ANALYSIS}/star/{{sample}}Aligned.sortedByCoord.out.bam',
        bai = f'{ANALYSIS}/star/{{sample}}Aligned.sortedByCoord.out.bam.bai',
    params:
        star_idx = config['star']['index'],
        prefix = f'{ANALYSIS}/star/{{sample}}',
    log:
        f'{LOG}/star/{{sample}}.log',
    benchmark:
        f'{BENCHMARK}/star/{{sample}}.txt',
    threads: config['hpc_parameters']['threads']['star'],
    resources:
        mem_gb = config['hpc_parameters']['memory']['large'],
        time = config['hpc_parameters']['runtime']['long'],
    conda:
        '../envs/star.yaml',
    shell:
        '''
        STAR \
            --runThreadN {threads} \
            --genomeDir {params.star_idx} \
            --readFilesIn {input} \
            --readFilesCommand zcat \
            --outFileNamePrefix {params.prefix} \
            --outSAMtype BAM SortedByCoordinate \
            --clip5pNbases 0 14 \
            --quantMode GeneCounts > {log}

        samtools index -@ {threads} {output.bam} >> {log}
        '''

rule samtools_stat:
    input:
        f'{ANALYSIS}/star/{{sample}}Aligned.sortedByCoord.out.bam',
    output:
        f'{ANALYSIS}/samtools/{{sample}}.stats',
    log:
        f'{LOG}/samtools/{{sample}}.stats.log',
    benchmark:
        f'{BENCHMARK}/samtools/{{sample}}.stats.txt',
    threads: config['hpc_parameters']['threads']['samtools'],
    resources:
        mem_gb = config['hpc_parameters']['memory']['small'],
        time = config['hpc_parameters']['runtime']['short'],
    conda:
        '../envs/star.yaml',
    shell:
        '''
        samtools stats -@ {threads} {input} 1> {output} 2> {log}
        '''

rule samtools_flagstat:
    input:
        f'{ANALYSIS}/star/{{sample}}Aligned.sortedByCoord.out.bam',
    output:
        f'{ANALYSIS}/samtools/{{sample}}.flagstat',
    log:
        f'{LOG}/samtools/{{sample}}.flagstat.log',
    benchmark:
        f'{BENCHMARK}/samtools/{{sample}}.flagstat.txt',
    threads: config['hpc_parameters']['threads']['samtools'],
    resources:
        mem_gb = config['hpc_parameters']['memory']['small'],
        time = config['hpc_parameters']['runtime']['short'],
    conda:
        '../envs/star.yaml',
    shell:
        '''
        samtools flagstat -@ {threads} {input} 1> {output} 2> {log}
        '''
