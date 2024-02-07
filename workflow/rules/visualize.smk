#----------------------------------------------------------------------------------------------------------------------
# Notes on where variables are defined, if not defined in prealignment.smk
#
# LOG       - workflow/Snakefile
# SAMPLES   - workflow/Snakefile
# ANALYSIS  - workflow/Snakefile
# BENCHMARK - workflow/Snakefile
# config    - workflow/Snakefile
#
#----------------------------------------------------------------------------------------------------------------------

def get_multiqc_params(wildcards):
    indirs = f'{ANALYSIS}/fastqc {ANALYSIS}/star {ANALYSIS}/samtools'

    return indirs

rule multiqc:
    input:
        # FASTqc
        expand(f'{ANALYSIS}/fastqc/{{samples.sample}}_R{{read}}_fastqc.zip', read=[1, 2], samples=SAMPLES.itertuples()),
        # STAR
        expand(f'{ANALYSIS}/star/{{samples.sample}}Log.final.out', samples=SAMPLES.itertuples()),
        expand(f'{ANALYSIS}/star/{{samples.sample}}ReadsPerGene.out.tab', samples=SAMPLES.itertuples()),
    output:
        directory(f'{ANALYSIS}/multiqc/multiqc_report_data'),
        f'{ANALYSIS}/multiqc/multiqc_report_data/multiqc_data.json',
        f'{ANALYSIS}/multiqc/multiqc_report.html',
    params:
        dirs = get_multiqc_params,
        out_dir = f'{ANALYSIS}/multiqc',
    log:
        f'{LOG}/multiqc/multiqc.log'
    benchmark:
        f'{BENCHMARK}/multiqc/multiqc.txt',
    threads: 1
    resources:
        mem_gb = config['hpc_parameters']['memory']['small'],
        time = config['hpc_parameters']['runtime']['short']
    conda:
        '../envs/python.yaml'
    shell:
        '''
        multiqc -f -o {params.out_dir} -n multiqc_report.html {params.dirs} 2> {log}
        '''

rule platetools_data:
    input:
        multiqc_data = f'{ANALYSIS}/multiqc/multiqc_report_data/multiqc_data.json',
        read_counts = expand(f'{ANALYSIS}/star/{{samples.sample}}ReadsPerGene.out.tab', samples=SAMPLES.itertuples()),
    output:
        platetools_data = f'{ANALYSIS}/plots/platetools_data.tsv',
    params:
        rrna = config['gene_ids']['rrna'],
        mito = config['gene_ids']['mito'],
        ercc = config['gene_ids']['ercc'],
    log:
        f'{LOG}/plots/platetools_data.log',
    benchmark:
        f'{BENCHMARK}/plots/platetools_data.txt',
    threads: 1
    resources:
        mem_gb = config['hpc_parameters']['memory']['small'],
        time = config['hpc_parameters']['runtime']['short'],
    conda:
        '../envs/python.yaml'
    script:
        '../scripts/collect_platetools_data.py'
