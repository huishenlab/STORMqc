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
        bam = temp(f'{ANALYSIS}/star/{{sample}}Aligned.sortedByCoord.out.bam'),
        bai = temp(f'{ANALYSIS}/star/{{sample}}Aligned.sortedByCoord.out.bam.bai'),
    params:
        star_idx = config['star']['index'],
        prefix = f'{ANALYSIS}/star/{{sample}}',
    log:
        f'{LOG}/star/{{sample}}.star_align.log',
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

rule non_annotated_space:
    input:
        f'{ANALYSIS}/star/{{sample}}Aligned.sortedByCoord.out.bam',
    output:
        n_reads = f'{ANALYSIS}/star/{{sample}}.non_annotated.tsv',
        subset = temp(f'{ANALYSIS}/star/{{sample}}Aligned.sortedByCoord.subsetted.bam'),
    params:
        annot = config['resources']['annotated_space'],
    log:
        f'{LOG}/star/{{sample}}.non_annotated_space.log',
    benchmark:
        f'{BENCHMARK}/star/{{sample}}.non_annotated_space.txt',
    threads: 1
    resources:
        mem_gb = config['hpc_parameters']['memory']['small'],
        time = config['hpc_parameters']['runtime']['short'],
    conda:
        '../envs/star.yaml',
    shell:
        '''
        # Pull out proper-pair mapped reads that are not flagged as not primary, duplicate, QC-fail or supplementary
        samtools view -f 3 -F 3840 -o {output.subset} {input} 2> {log}

        # total number of reads from subsetted BAM
        echo -e "total_reads\t`samtools view -c {output.subset}`" 1> {output.n_reads} 2>> {log}

        # number of reads that have no overlap with annotated space
        echo -e "non_annotated\t`bedtools intersect -v -a {output.subset} -b {params.annot} | samtools view | wc -l`" 1>> {output.n_reads} 2>> {log}
        '''
