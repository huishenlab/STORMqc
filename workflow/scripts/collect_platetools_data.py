from functools import reduce
import pandas as pd
import numpy as np
import json
import sys
import re
import os

def read_file(fname):
    with open(fname, 'r') as fh:
        data = json.load(fh)
    return data

def parse_name(s):
    pieces = s.split('_')

    m1 = re.match('^[A-P][0-2][0-9]$', pieces[0])
    if m1 is None:
        print(f'ERROR: Unknown identifier: {s}. Must be in the form: (well)_(cell id)_(read #)')
        sys.exit(1)
    well = pieces[0]

    read = None
    id = '_'.join(pieces[1:])
    if 'R1' in pieces or 'R2' in pieces:
        read = pieces[-1]
        id = '_'.join(pieces[1:-1])

    return (well, id, read)

def parse_qual_data(l):
    total = 0
    numerator = 0
    q30 = 0

    for pair in l:
        if pair[0] >= 30:
            q30 += pair[1]

        total += pair[1]
        numerator += pair[0] * pair[1]

    return (round(numerator/total, 2), round(100* q30 / total, 2))

def parse_per_sequence_qual_scores(data):
    score_list = data['report_plot_data']['fastqc_per_sequence_quality_scores_plot']['datasets'][0]

    indexes = {}
    for idx, d in enumerate(score_list):
        base = d['name'].replace('_R1', '').replace('_R2', '')
        read = 'R1' if 'R1' in d['name'] else 'R2'

        try:
            indexes[base][read] = idx
        except:
            indexes[base] = {read: idx}

    qual_data = {'well': [], 'cell_id': [], 'qual_avg_R1': [], 'qual_q30_R1': [], 'qual_avg_R2': [], 'qual_q30_R2': []}
    for base, dic in indexes.items():
        well, id, _ = parse_name(base)

        avg_qual_R1, q30_qual_R1 = parse_qual_data(score_list[dic['R1']]['data'])
        avg_qual_R2, q30_qual_R2 = parse_qual_data(score_list[dic['R2']]['data'])

        qual_data['well'].append(well)
        qual_data['cell_id'].append(id)
        qual_data['qual_avg_R1'].append(avg_qual_R1)
        qual_data['qual_q30_R1'].append(q30_qual_R1)
        qual_data['qual_avg_R2'].append(avg_qual_R2)
        qual_data['qual_q30_R2'].append(q30_qual_R2)

    return pd.DataFrame(qual_data)

def parse_general_stats(data):
    sub = data['report_saved_raw_data']['multiqc_general_stats']

    uniq = set([s.replace('_R1', '').replace('_R2', '') for s in sub.keys()])

    read_data = {'well': [], 'cell_id': [], 'n_reads_mil_R1': [], 'n_reads_mil_R2': [], 'percent_uniq_map': []}
    for name in uniq:
        well, id, _ = parse_name(name)
        read_data['well'].append(well)
        read_data['cell_id'].append(id)

        read_data['n_reads_mil_R1'].append(
            round(sub[f'{name}_R1']['FastQC_mqc-generalstats-fastqc-total_sequences'] / 1000000.0, 2)
        )
        read_data['n_reads_mil_R2'].append(
            round(sub[f'{name}_R2']['FastQC_mqc-generalstats-fastqc-total_sequences'] / 1000000.0, 2)
        )
        read_data['percent_uniq_map'].append(
            sub[name]['STAR_mqc-generalstats-star-uniquely_mapped_percent']
        )

    return pd.DataFrame(read_data)

def parse_read_counts(count_files, params):
    rrna = pd.read_csv(params['rrna'], sep=' ', names=['id', 'name'])
    mito = pd.read_csv(params['mito'], sep=' ', names=['id', 'name'])
    ercc = pd.read_csv(params['ercc'], sep=' ', names=['id', 'name'])

    count_data = {'well': [], 'cell_id': [], 'rrna': [], 'mito': [], 'ercc': []}
    for fname in count_files:
        _, name = os.path.split(fname)
        well, id, _ = parse_name(name.replace('ReadsPerGene.out.tab', ''))
        df = pd.read_csv(fname, sep='\t', names=['id', 'unstranded', 'forward', 'reverse'], header=None, skiprows=4)
        total_reads = df['reverse'].sum()

        df_rrna = pd.merge(rrna, df, on=['id'])
        df_mito = pd.merge(mito, df, on=['id'])
        df_ercc = pd.merge(ercc, df, on=['id'])

        count_data['well'].append(well)
        count_data['cell_id'].append(id)
        count_data['rrna'].append(round(100 * df_rrna['reverse'].sum() / total_reads, 4))
        count_data['mito'].append(round(100 * df_mito['reverse'].sum() / total_reads, 4))
        count_data['ercc'].append(round(100 * df_ercc['reverse'].sum() / total_reads, 4))

    return pd.DataFrame(count_data)

def main(multiqc_data, count_files, params, oname):
    data = read_file(multiqc_data)

    datasets = []
    datasets.append(parse_per_sequence_qual_scores(data))
    datasets.append(parse_general_stats(data))
    datasets.append(parse_read_counts(count_files, params))

    df = reduce(lambda x, y: pd.merge(x, y, on=['well', 'cell_id']), datasets)
    df.to_csv(oname, sep='\t', na_rep='NA', index=False)

with open(snakemake.log[0], 'w') as fh:
    sys.stderr = sys.stdout = fh
    main(
        snakemake.input['multiqc_data'],
        snakemake.input['read_counts'],
        snakemake.params,
        snakemake.output['platetools_data'],
    )
