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
        non_annotated = expand(f'{ANALYSIS}/star/{{samples.sample}}.non_annotated.tsv', samples=SAMPLES.itertuples()),
    output:
        platetools_data = f'{ANALYSIS}/plots/platetools_data.tsv',
    params:
        rrna = config['resources']['rrna_gene_ids'],
        mito = config['resources']['mito_gene_ids'],
        ercc = config['resources']['ercc_gene_ids'],
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

rule make_plots:
    input:
        pt_data = f'{ANALYSIS}/plots/platetools_data.tsv',
    output:
        pdf = f'{ANALYSIS}/plots/quality_control_plots.pdf',
    log:
        log_file = f'{LOG}/plots/make_plots.log',
    benchmark:
        f'{BENCHMARK}/plots/make_plots.txt',
    threads: 1
    resources:
        mem_gb = config['hpc_parameters']['memory']['small'],
        time = config['hpc_parameters']['runtime']['short'],
    conda:
        '../envs/r.yaml'
    script:
        '../scripts/create_qc_plots.R'
