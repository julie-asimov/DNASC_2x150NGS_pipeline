import argparse
import os
import pandas as pd


def get_mode(row):
    if row['status'].lower() == 'fail':
        if row['Depth'] <= 40:
            return 'Seq Fail'
        elif row['%Mapped'] < 90:
            return '%Mapped Fail'
        elif 90 <= row['%Mapped'] < 98:
            return 'Ambiguous'
    


def create_echo_work_list(data, output_file):
    # Group by 'STOCK_ID'
    grouped_data = data.groupby('STOCK_ID')

    echo_wl = pd.DataFrame()  # Initialize the echo_wl DataFrame

    for stock_id, group in grouped_data:
        if 'pass' not in group['status'].values and 'Seq Fail' in group['Mode'].values:
            for _, row in group[group['Mode'] == 'Seq Fail'].iterrows():
                for _ in range(3):
                    source_well_coord = row['384MP_WELL_COORD']
                    source_well_coord = source_well_coord if not pd.isna(source_well_coord) else ""
                    echo_row = {
                        'Source Plate Barcode': row['384MP_PLATE_ID'],
                        'Source Plate Type': '384PP_AQ_BP',
                        'Source Well': "{" + str(source_well_coord) + "}",
                        'Destination Plate Barcode': '0',
                        'Destination Well': '',
                        'Transfer Volume': '1500',
                        '': row['SAMPLE_NAME']
                    }
                    echo_wl = echo_wl.append(echo_row, ignore_index=True)

    if not echo_wl.empty:
        # Export the DataFrame to CSV
        echo_wl.to_csv(output_file, index=False)


def summary_file(metadata_csv: str, output_dir: str) -> None:
    """
    Generate a summary file and a list of passed wells from sequencing data.
    Args:
        metadata_csv (str): Path to metadata CSV file.
        output_dir (str): Path to output directory.
    Returns:
        None
    """

    # Set paths to files
    qc_file = os.path.join(output_dir, "qc.csv")
    metrics_file = os.path.join(output_dir, "metrics.csv")

    # Load data into dataframes
    df_metrics = pd.read_csv(metrics_file)
    df_metadata = pd.read_csv(metadata_csv)
    df_qc = pd.read_csv(qc_file)
    
    # Merge dataframes
    df = pd.merge(df_metrics, df_metadata, on='SEQ_NAME')
#    print(df)
    df = pd.merge(df, df_qc, left_on='INDEX_ID', right_on='filename')
    print(df)

    # Filter by criteria
    df['status'] = (
        (df['Coverage'] == 100) &
        df['SNPs'].isna() &
        df['Missing_Positions'].isna() &
        (df['%Mapped'] > 90) &
        df['Homopolymer_SNPs'].isna() &
        df['Softclip_pileup'].isna() &
        df['Backbone_Softclip_pileup'].isna() &
        df['Backbone_SNPs'].isna() &
        df['Backbone_homopolymer_SNPs'].isna() &
        (df['total_filtered_reads'] >= 1000)
    ).map({True:'pass', False:'fail'})

    
# Apply the function to each row of the dataframe for rows with 'fail' in the 'Status' column
    df['Mode'] = df[df['status'].str.lower() == 'fail'].apply(get_mode, axis=1)


# Get passed samples and sample with the highest coverage
    df_pass = df[df['status'] == 'pass']
    #make the sample with the highest depth of coverage available
    df['max'] = df_pass.groupby('STOCK_ID')['CONCENTRATION_NGUL'].transform('max')
    # Create a mask for selecting samples with highest coverage
    mask = (df['status'] == 'pass') & (df['CONCENTRATION_NGUL'] == df['max'])

    # Label selected samples
    df.loc[mask, 'selected'] = 'yes'

    #Add columns
    df['Comments'] = pd.Series(dtype=str)
    df['Scientist'] = pd.Series(dtype=str)

    # Format data
    df = df.drop(['max'], axis=1)

    df = df.reindex(columns=[
       'Comment', 'Scientist', 'Mode', 'status', 'selected', 'ASSOCIATED_WELLS','SEQ_NAME', 'Not_Enough_Reads','Plasmid_ID', 'SNPs',
       'Missing_Positions', 'Softclip_pileup', 'Homopolymer_SNPs', 'Coverage',
       'Depth', '%Mapped', 'Low_coverage_regions','Lacz_coverage', 'Backbone_Found', 'Backbone_SNPs',
       'Backbone_homopolymer_SNPs', 'Backbone_Softclip_pileup',
       'Known_BB_SNPs', 'REQUEST_ID','INDEX_ID', 'HISTORICAL_AVAILABLE_GLY','WELL_ID', 'PLATE_ID',
       'STOCK_ID', 'PROCESS_ID', 'CONCENTRATION_NGUL', 'WELL_NUMBER', 'SOURCE','DNA_TO_ASSEMBLE', 'DNA_WELLS_TO_ASSEMBLE',
       'PLASMID_ALIAS', 'AVAILABLE', 'SEQ_CONFIRMED', 'SAMPLE_NAME',
       'Antibiotic', 'ASSEMBLY_PLATE_ID', 'ASSEMBLY_WELL_ID',
       'ASSEMBLY_WELL_COORD', 'ASSEMBLY_TYPE', 'AGAR_PLATE_ID', 'AGAR_WELL_ID',
       'AGAR_WELL_COORD', 'COLOR', 'BG_COLOR', 'COMP_CELL',
       'GLYCEROL_PLATE_ID', 'GLYCEROL_WELL_ID', 'GLYCEROL_WELL_COORD',
       'GLYCEROL_AVAILABILITY', 'GLYCEROL_SEQ_CONFIRMED', '96MP_PLATE_ID',
       '96MP_WELL_ID', '96MP_AVAILABLE', '96MP_SEQ_CONFIRMED',
       '96MP_WELL_COORD', '384MP_PLATE_ID', '384MP_WELL_ID',
       '384MP_WELL_COORD', '384MP_AVAILABLE', '384MP_SEQ_CONFIRMED','filename', 'total_raw_Reads', 'total_raw_bases',
       'total_filtered_reads', 'total_filtered_bases', 'percent_filtered'])


    # Write output files
    metadata_csv_name = os.path.splitext(os.path.basename(metadata_csv))[0]
    new_file_name = metadata_csv_name + '_results.csv'
    df.to_csv(os.path.join(output_dir, new_file_name), sep=',', index=False, header=True)

    #df_take = df[df['selected'] == 'yes']
    df_take = df[(df['selected'] == 'yes') & (df['HISTORICAL_AVAILABLE_GLY'] == 'TRUE')]
    listfil = df_take['SEQ_NAME'].str.rsplit('_', 1).str.get(0).str.replace('_', '-')
    take = listfil.tolist()

    df_wells = df_metadata[df_metadata['SAMPLE_NAME'].isin(take)]
    df_wells = df_wells[df_wells['ASSOCIATED_WELLS'].str.contains('well').fillna(False)]
    well_ids = ','.join(df_wells['ASSOCIATED_WELLS'])

    new_file_name_well = metadata_csv_name + '_AvailableWells.txt'
    with open(os.path.join(output_dir, new_file_name_well), 'w') as f:
        f.write(well_ids) 

    #get list of seq confirmed
    listfil = df_pass['SEQ_NAME'].str.rsplit('_', 1).str.get(0).str.replace('_', '-')
    take_sc = listfil.tolist()

    df_wells = df_metadata[df_metadata['SAMPLE_NAME'].isin(take_sc)]
    df_wells = df_wells[df_wells['ASSOCIATED_WELLS'].str.contains('well').fillna(False)]
    well_ids = ','.join(df_wells['ASSOCIATED_WELLS'])

    new_file_name_well = metadata_csv_name + '_SeqConfwells.txt'
    with open(os.path.join(output_dir, new_file_name_well), 'w') as f:
        f.write(well_ids)



# Assuming you already have the 'data' DataFrame with required columns
    output_file_path = os.path.join(output_dir, metadata_csv_name + '_SeqRepeat_WL.csv')
    create_echo_work_list(df, output_file_path)


def main():
    parser = argparse.ArgumentParser(description='Generate summary file and list of passed wells from sequencing data.')
    parser.add_argument('metadata_csv', type=str, help='Path to metadata CSV file.')
    parser.add_argument('output_dir', type=str, help='Path to output directory.')
    args = parser.parse_args()

    summary_file(args.metadata_csv, args.output_dir)


if __name__ == '__main__':
    main()
