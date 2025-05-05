#!usr/local/bin/python3

import sys
import re
import pandas as pd

def parse_ipr(ipr_file):
    ipr = pd.read_table(ipr_file, header=None)
    ipr = ipr.iloc[:, [0, 4, 11, 13]]
    ipr.columns = ['query', 'pfam_id', 'ipr_id', 'Go_term']
    ipr[['pfam_id', 'ipr_id', 'Go_term']] = ipr[['pfam_id', 'ipr_id', 'Go_term']].apply(lambda col: col.str.replace('-', '', regex=False))
    return ipr

def parse_blastp(blastp_results):
    blastp = pd.read_table(blastp_results, header=None)
    blastp = blastp.iloc[:, [0, 8]]
    blastp.columns = ['query', 'attributes']
    blastp[['gene_name', 'species_name', 'protein_description']] = blastp['attributes'].apply(extract_info)


    return blastp

def extract_info(attr):
    gene, species, protein_desc = '', '', ''
    note_match = re.search(r'Note=(.*)', attr)
    if note_match:
        note = note_match.group(1)
        gene_match = re.search(r'GN=(.*?)\s+PE=', note)
        if gene_match:
            gene = gene_match.group(1).strip()
        species_match = re.search(r'OS=(.*?)\s+OX=', note)
        if species_match:
            species = species_match.group(1).strip()
        parts = note.split('|')
        if len(parts) >= 5:
            protein_words = re.search(r'(.*?)\s+OS=', parts[4]).group(1).strip().split()
            if len(protein_words) > 1:
                protein_desc = ' '.join(protein_words[1:])

    return pd.Series([gene, species, protein_desc])

def build_new_attr(row):
        attr = row['attribute']
        if pd.notna(row.get('ipr_id')) and row['ipr_id']:
            attr += f';Dbxref=InterPro:{row["ipr_id"]}'
        if pd.notna(row.get('pfam_id')) and row['pfam_id']:
            attr += f';Pfam:{row["pfam_id"]}'
        if pd.notna(row.get('Go_term')) and row['Go_term']:
            attr += f';Ontology_term={row["Go_term"]}'
        if pd.notna(row.get('gene_name')) and row['gene_name']:
            note_str = f'Similar to {row["gene_name"]}: {row["protein_description"]} ({row["species_name"]})'
            attr += f';Note={note_str}'
        return attr


def main():
    gff_file, blastp_results, ipr_file = sys.argv[1:]
    gff = pd.read_table(gff_file, names=['seqname', 'source', 'feature',
                                         'start', 'end', 'score',
                                         'strand', 'frame', 'attribute'])
    blastp = parse_blastp(blastp_results)
    ipr = parse_ipr(ipr_file)
    gff_mrna = gff[gff['feature'] == 'mRNA'].copy()
    gff_mrna['ID'] = gff_mrna['attribute'].str.extract(r'ID=([^;]+)')
    merged = gff_mrna.merge(ipr.groupby('query').agg({
        'ipr_id': lambda x: ','.join(sorted(set(x.dropna()))),
        'pfam_id': lambda x: ','.join(sorted(set(x.dropna()))),
        'Go_term': lambda x: ','.join(sorted(set(x.dropna())))
    }).reset_index(), left_on='ID', right_on='query', how='left')
    merged = merged.merge(blastp, left_on='ID', right_on='query', how='left')
    merged['attribute'] = merged.apply(build_new_attr, axis=1)
    merged.to_csv("merged", sep = '\t', index= False)

    gff.loc[gff['feature'] == 'mRNA', 'attribute'] = merged['attribute'].values
    gff.to_csv('annotated.gff', sep='\t', index=False, header=False)

if __name__ == '__main__':
