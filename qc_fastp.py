import os
import json
import pandas as pd
import re


def qc_fastp(output_dir):
    # get names of all json files found in folder into pandas

    fastp_dir=os.path.join(output_dir, "fastp")

    json_files = [pos_json for pos_json in os.listdir(fastp_dir) if pos_json.endswith('.json')]



    # here I define my pandas Dataframe with the columns I want to get from the json
    jsons_data = pd.DataFrame(columns=['filename', 'total_raw_Reads', 'total_raw_bases', 'total_filtered_reads', 'total_filtered_bases','percent_filtered'])

    # we need both the json and an index number so use enumerate(). Read in json file one by one
    for index, js in enumerate(json_files):
        with open(os.path.join(fastp_dir, js)) as json_file:
            json_text = json.load(json_file)
            #grab stats from json files
            filename = re.sub('(_L001).*',"", js)
            total_raw_reads = json_text['summary']['before_filtering']['total_reads']
            total_raw_bases = json_text['summary']['before_filtering']['total_bases']
            total_filtered_reads = json_text['summary']['after_filtering']['total_reads']
            total_filtered_bases = json_text['summary']['after_filtering']['total_bases']
            try:
                percent_filtered = (1-(total_filtered_reads/total_raw_reads))*100 #if there are zero reads this will scream
            except:
                percent_filtered = 0

            # here I push a list of data into a pandas DataFrame at row given by 'index'
            jsons_data.loc[index] = [filename,total_raw_reads,total_raw_bases,total_filtered_reads,total_filtered_bases,percent_filtered]

    #output
    jsons_data.to_csv(os.path.join(output_dir, 'qc.csv'), sep=',', index=False, header=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser("requires: (1) output folder")
    parser.add_argument('output_dir')
    args = parser.parse_args()

    qc_fastp(args.output_dir)


#    try:
 #       file_path = os.path.join(args.fastp_dir, qc.csv)
  #      check_file_size(file_path)
   # except Exception as e:
     #   print(f"Error checking file: {filename}")
    #    print(e)

