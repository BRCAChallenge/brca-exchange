#!/usr/bin/env python
"""
This script downloads, merges and postprocesses the preprocessed ENIGMA data
from Synapse, and updates the ENIGMA_combined_preprocessed_hg38 file
"""
import argparse
import os
import re
import shutil
import subprocess
import synapseclient




def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-u", "--username", help="Synapse username")
    parser.add_argument("-p", "--password", help="Synapse password")
    parser.add_argument("-s", "--synapse_entity", default="syn7188267",
                        help="Synapse ID for the folder of preprocessed data")
    parser.add_argument("-w", "--work_dir", 
                        default="../brca/pipeline-data/data",
                        help="Parent directory for the working directory")
    parser.add_argument('-a', "--artifacts_dir",
                        help='Artifacts directory with pipeline artifact files.')
    parser.add_argument("-k", "--keep_tmp_dir", default  = False,
                        help="Keep the temporary directory (do not delete)")
    args = parser.parse_args()

    #
    # Set up a temp workspace and log in to Synapse
    temp_dir =  args.work_dir + "/enigma/"
    os.makedirs(temp_dir)
    syn = synapseclient.Synapse()
    syn.login(args.username, args.password)
    
    #
    # Download all the ENIGMA_last_updated*hg38.tsv files, the preprocessed files
    query = 'SELECT id, name FROM entity WHERE parentId=="%s"' \
        % (args.synapse_entity)
    results = syn.query(query)
    for item in results['results']:
        if re.search("ENIGMA_last_updated(.)*hg38.tsv", item['entity.name']):
            print "fetching", item['entity.name']
            syn.get(item['entity.id'], downloadLocation=temp_dir)
    
    #
    # Merge the preprocessed files, postprocess the merged file, and upload
    # the resulting file to Synapse
    merge_cmd = "python enigma-merge_hg38.py -i %s -o %s" % (temp_dir, temp_dir)
    subprocess.check_call(merge_cmd, shell=True)
    postprocess_cmd = "enigma_postprocess.py -i %s -o %s -a %s" \
        % (temp_dir + "/ENIGMA_combined_hg38.tsv", 
           temp_dir + "/ENIGMA_combined_postprocessed_hg38.tsv", 
           args.artifacts_dir)
    subprocess.check_call(postprocess_cmd, shell=True)
    output_entity = synapseclient.File(temp_dir + "/ENIGMA_combined_postprocessed_hg38.tsv",
                                       parent=args.synapse_entity)
    output_entity = syn.store(output_entity)

    # 
    # Cleanup...
    if not args.keep_tmp_dir:
        shutil.rmtree(temp_dir)
               
    

if __name__ == "__main__":
    main()
