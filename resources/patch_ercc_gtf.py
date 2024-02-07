#!/usr/bin/env python

with open('ERCC92.gtf') as exon_gtf, open('ERCC92.genes.patched.gtf', 'w') as gene_gtf:
    for line in exon_gtf:
        f = line.strip().split('\t')
        #f[0] = f[0].replace('-','_')  # required for RNA-SeQC/GATK (no '-' in contig name)

        attr = f[8]
        if attr[-1]==';':
            attr = attr[:-1]
        attr = dict([i.split(' ') for i in attr.replace('"','').split('; ')])
        # add gene_name, gene_type
        attr['gene_name'] = attr['gene_id']
        attr['gene_type'] = 'ercc_control'
        attr['gene_status'] = 'KNOWN'
        attr['level'] = 2
        for k in ['id', 'type', 'name', 'status']:
            attr['transcript_'+k] = attr['gene_'+k]

        attr_str = []
        for k in ['gene_id', 'transcript_id', 'gene_type', 'gene_status', 'gene_name',
            'transcript_type', 'transcript_status', 'transcript_name']:
            attr_str.append('{0:s} "{1:s}";'.format(k, attr[k]))
        attr_str.append('{0:s} {1:d};'.format('level', attr['level']))
        f[8] = ' '.join(attr_str)

        # write gene, transcript, exon
        gene_gtf.write('\t'.join(f[:2]+['gene']+f[3:])+'\n')
        gene_gtf.write('\t'.join(f[:2]+['transcript']+f[3:])+'\n')
        f[8] = ' '.join(attr_str[:2])
        gene_gtf.write('\t'.join(f[:2]+['exon']+f[3:])+'\n')

