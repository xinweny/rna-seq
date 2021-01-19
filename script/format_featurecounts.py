import pandas as pd
import argparse
import re

def main():
    # Set up command line parser
    parser = argparse.ArgumentParser(description='Simplify featureCounts output')
    parser.add_argument('-i', required=True, help='Path to featureCounts table.')
    parser.add_argument('-o', required=True, help='Output file name.')

    args = parser.parse_args()

    fc_df = pd.read_csv(args.i, header=0, comment='#', sep='\t')

    # Drop columns
    fc_df = fc_df.drop(columns=['Chr', 'Start', 'End', 'Strand'])

    # Remove version from Ensembl gene id
    fc_df['Geneid'] = fc_df['Geneid'].apply(lambda x: re.sub(r"\.[0-9]+", '', x))

    # Reformat column name
    new_cols = ['Geneid', 'Length'] + [col.split('/')[-1].replace('.bam', '') for col in fc_df.columns[2:]]
    fc_df.columns = new_cols

    # Save to output
    fc_df.to_csv(args.o, sep='\t', index=False)

# Execute main()
if __name__ == '__main__':
    main()
