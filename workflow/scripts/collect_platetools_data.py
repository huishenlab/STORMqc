import pandas as pd
import numpy as np
import json
import sys

def read_file(fname):
    with open(fname, 'r') as fh:
        data = json.load(fh)
    return data

def parse_name(s):
    pieces = s.split('_')

    if len(pieces) == 3:
        well = pieces[0]
        id = pieces[1]
        read = pieces[2]

        return (well, id, read)

    print(f'ERROR: Unknown identifier: {s}. Must be in the form: (well)_(cell id)_(read #)')
    sys.exit(1)

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

    qual_data = {'well': [], 'cell_id': [], 'read': [], 'qual_avg': [], 'qual_q30': []}
    for i in range(len(score_list)):
        well, id, read = parse_name(score_list[i]['name'])

        avg_qual, q30_qual = parse_qual_data(score_list[i]['data'])

        qual_data['well'].append(well)
        qual_data['cell_id'].append(id)
        qual_data['read'].append(read)
        qual_data['qual_avg'].append(avg_qual)
        qual_data['qual_q30'].append(q30_qual)

    return pd.DataFrame(qual_data)

def parse_general_stats(data):
    sub = data['report_saved_raw_data']['multiqc_general_stats']

    read_data = {'well': [], 'cell_id': [], 'read': [], 'n_reads_mil': []}
    for name, dic in sub.items():
        well, id, read = parse_name(name)
        n_reads = dic['FastQC_mqc-generalstats-fastqc-total_sequences']

        read_data['well'].append(well)
        read_data['cell_id'].append(id)
        read_data['read'].append(read)
        read_data['n_reads_mil'].append(round(n_reads / 1000000.0, 2))

    return pd.DataFrame(read_data)

def main(fname, oname):
    data = read_file(fname)

    qual_data = parse_per_sequence_qual_scores(data)
    gen_data  = parse_general_stats(data)

    df = qual_data.merge(gen_data, on=['well', 'cell_id', 'read'])
    df.to_csv(oname, sep='\t', na_rep='NA', index=False)

with open(snakemake.log[0], 'w') as fh:
    sys.stderr = sys.stdout = fh
    main(
        snakemake.input['multiqc_data'],
        snakemake.output['platetools_data'],
    )
