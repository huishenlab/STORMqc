from functools import reduce
import pandas as pd
import numpy as np
import json
import sys
import re
import os

def read_json_file(fname):
    """Read JSON file into a dict.

    Inputs -
        fname - str, filename
    Returns -
        dict
    """
    with open(fname, 'r') as fh:
        data = json.load(fh)
    return data

def parse_name(s):
    """Parse name string into its well, id, and read number.

    Inputs -
        s - str, name string for a well
    Returns -
        tuple - (str, str, str)
    """
    pieces = s.split('_')

    # Make sure the name string has the well in the right location
    m1 = re.match('^[A-P][0-2][0-9]$', pieces[0])
    if m1 is None:
        print(f'ERROR: Unknown identifier: {s}. Must be in the form: (well)_(cell id)_(read #)')
        sys.exit(1)

    well = pieces[0]

    # Default to case where R1 and R2 are not in the name string
    # Can then check and change if needed
    read = None
    id = '_'.join(pieces[1:])
    if 'R1' in pieces or 'R2' in pieces:
        read = pieces[-1]
        id = '_'.join(pieces[1:-1])

    return (well, id, read)

def parse_qual_data(l):
    """Parse list of quality score data.

    Inputs -
        l - list, elements are 2-element lists with qual score and corresponding count
    Returns -
        tuple, (float, float)
    """
    total     = 0 # total number of reads
    numerator = 0 # sum of score*count for weighted averaging
    q30       = 0 # number of reads with score >= 30

    for pair in l:
        if pair[0] >= 30:
            q30 += pair[1]

        total += pair[1]
        numerator += pair[0] * pair[1]

    return (round(numerator/total, 2), round(100 * q30 / total, 2))

def parse_per_sequence_qual_scores_dict(score_list):
    """Parse dict with name and data entries.

    Inputs -
        score_list - dict, follows format of [{'name': '(well)_(id)_(read)', 'data': [[score, count], ...]
    Returns -
        pd.DataFrame
    """
    # Reads 1 and 2 may not always be next to one another, so scan through list and pull out indexes for each read
    indexes = {}
    for idx, d in enumerate(score_list):
        base = d['name'].replace('_R1', '').replace('_R2', '')
        read = 'R1' if 'R1' in d['name'] else 'R2'

        try:
            indexes[base][read] = idx
        except:
            indexes[base] = {read: idx}

    # Calculate stats and setup for creating DataFrame
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

def parse_per_sequence_qual_scores_json(data):
    """Wrapper around parse_per_sequence_qual_scores_dict for JSON dict input.

    Inputs -
        data - dict, from multiqc_data.json
    Returns -
        pd.DataFrame
    """
    score_list = data['report_plot_data']['fastqc_per_sequence_quality_scores_plot']['datasets'][0]

    return parse_per_sequence_qual_scores_dict(score_list)

def parse_per_sequence_qual_scores_file(fname):
    """Wrapper around parse_per_sequence_qual_scores_dict for text file input.

    Inputs -
        fname - str, file name
    Returns -
        pd.DataFrame
    """
    # If reading from a text file, we need to put the file contents into the correct data structure for
    # parse_per_sequence_qual_scores_dict
    with open(fname, 'r') as fh:
        file_contents = fh.read()

    # File should have an even number of rows (a score and a count line for each sample)
    # Verify that's the case before continuing
    file_contents = file_contents.strip('\n').split('\n')
    if len(file_contents) % 2 != 0:
        print('ERROR: malformed file: {fname}')
        sys.exit(1)

    # Parse file contents into expected data structure
    score_list = []
    for i in range(0, len(file_contents)-1, 2):
        scores = file_contents[i].split('\t')
        counts = file_contents[i+1].split('\t')
        if len(scores) != len(counts):
            print('ERROR: different number of entries in score and count lines: {fname}')

        name = counts[0]
        data = [[float(scores[j]), float(counts[j])] for j in range(1, len(scores))]
        score_list.append({'name': name, 'data': data})

    return parse_per_sequence_qual_scores_dict(score_list)

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
    """Parse counts matrix file to find counts to undesired locations.

    Inputs -
        count_files - list, file names for input
        params - dict, contains file names that include undesired locations
    Returns -
        pd.DataFrame
    """
    # rRNA, mitochondrial DNA, and ERCC gene names
    rrna = pd.read_csv(params['rrna'], sep=' ', names=['id', 'name'])
    mito = pd.read_csv(params['mito'], sep=' ', names=['id', 'name'])
    ercc = pd.read_csv(params['ercc'], sep=' ', names=['id', 'name'])

    # Parse input files
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

def parse_non_annotated(infiles):
    """Parse text files with info on reads aligning to non-annotated space.

    Inputs -
        infiles - list, file names for input
    Returns -
        pd.DataFrame
    """
    out_data = {'well': [], 'cell_id': [], 'non_annot': []}
    for fname in infiles:
        _, name = os.path.split(fname)
        well, id, _ = parse_name(name.replace('.non_annotated.tsv', ''))

        with open(fname, 'r') as fh:
            data = fh.read()

        m = re.search(r'total_reads\s+(\d+)\s+non_annotated\s+(\d+)', data, re.MULTILINE)
        if m is None:
            print(f'ERROR: Could not find necessary data from {fname}')
            sys.exit(1)

        percent = round(100 * int(m.group(2)) / int(m.group(1)), 2)

        out_data['well'].append(well)
        out_data['cell_id'].append(id)
        out_data['non_annot'].append(percent)

    return pd.DataFrame(out_data)

def main(input, params, oname):
    # Read MultiQC data JSON file
    data = read_json_file(input['multiqc_data'])

    # Figure out how we can get the quality score data
    qual_scores = None
    if 'fastqc_per_sequence_quality_scores_plot' in data['report_plot_data'].keys():
        qual_scores = parse_per_sequence_qual_scores_json(data)
    else:
        qual_scores = parse_per_sequence_qual_scores_file(
            input['multiqc_data'].replace(
                'multiqc_data.json',
                'mqc_fastqc_per_sequence_quality_scores_plot_1.txt'
            )
        )

    # Collect datasets for merging
    datasets = []
    datasets.append(qual_scores)
    datasets.append(parse_general_stats(data))
    datasets.append(parse_read_counts(input['read_counts'], params))
    datasets.append(parse_non_annotated(input['non_annotated']))

    df = reduce(lambda x, y: pd.merge(x, y, on=['well', 'cell_id']), datasets)
    df.to_csv(oname, sep='\t', na_rep='NA', index=False)

# Run script
with open(snakemake.log[0], 'w') as fh:
    sys.stderr = sys.stdout = fh
    main(
        snakemake.input,
        snakemake.params,
        snakemake.output['platetools_data'],
    )
