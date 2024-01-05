# parse Kraken2 results
# get all hamming_distace and LCA pairs
# get all LCA taxids
# get all LCA taxonomies

import os
import argparse
import concurrent.futures
import pandas as pd
from collections import Counter
import csv
import subprocess
from shlex import quote
import time
import re

parser = argparse.ArgumentParser(description="Parse Kraken2 results\nget all hamming_distace and LCA pairs\nget all LCA taxids\nget all LCA taxonomies")
parser.add_argument('--kraken2_res_folder', type=str, required=True,
                    help='The path to the folder containing kragen2 results (output of run_kraken2.py)',
                    default = "/hpc/projects/data.science/duo.peng/Cas13_design_for_virus/cas13a_guide_design_py/top_guides_out/1_neighbors_kraken2.out/")
parser.add_argument('-n', '--num_workers', type=int, default=32, help="number of workers for parallel processing of kraken2 result parsing")
args = parser.parse_args()

def chunk_list(lst, chunk_size):
    """
    Breaks a list into chunks of the specified size.
    The last chunk may have fewer items if the list size is not divisible by the chunk size.
    """
    return [lst[i:i+chunk_size] for i in range(0, len(lst), chunk_size)]

def taxid2taxonomy(taxid="10239"):
    """get the taxonomy of a taxid
    input: taxid, output: taxonomy text
    """
    api_call = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=taxonomy&id=" + taxid + "&api_key=ae0efc082ed8a9e97c68724ad915598d7308&rettype=native&retmode=text"
    command = ["curl", api_call]
    command = ' '.join(quote(s) for s in command)
    completed_process = subprocess.run(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    output = completed_process.stdout.decode('utf-8')
    return output.replace("\n", "")

def taxidlist2taxonomy(taxid=["10239","9606","11118"]):
    """get the taxonomy of a taxid
    input: taxid, output: dictionary of taxid:taxonomy
    """
    taxid_query = ",".join(taxid)
    api_call = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=taxonomy&id=" + taxid_query + "&api_key=ae0efc082ed8a9e97c68724ad915598d7308&rettype=native&retmode=text"
    command = ["curl", api_call]
    command = ' '.join(quote(s) for s in command)
    completed_process = subprocess.run(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    output = completed_process.stdout.decode('utf-8')
    #split the output into a list of taxonomies
    fields = re.split('\n\d+\.', output)
    if len(fields) != len(taxid):
        print("warning: number of taxonomies does not match the number of taxids")
        print(f"taxid: {taxid}")
        return None
    return dict(zip(taxid, fields))

def worker_taxonomy(taxid):
    start_time = time.time()
    taxonomy = taxid2taxonomy(taxid = taxid)
    end_time = time.time()
    elapsed_time = end_time - start_time
    if elapsed_time < 1.2:
        time.sleep(1 - elapsed_time + 0.2)
    return [taxid,taxonomy]

def worker_readcsv(filepath):
    df = pd.read_csv(filepath, sep='\t', header=None)
    df.columns = ["classification", "sequence_ID", "taxonomy_ID", "sequence_length", "LCA"]
    df["hd"] = df["sequence_ID"].apply(lambda x: x.split("_")[1])
    _hd_LCA = df.apply(lambda x: (x["hd"] + "_" + x["LCA"]), axis=1)
    return _hd_LCA

def read_kraken2_res(kraken2_res_folder):
    filelist = os.listdir(kraken2_res_folder)
    filepathlist = [os.path.join(kraken2_res_folder, filename) for filename in filelist]
    with concurrent.futures.ProcessPoolExecutor(max_workers=args.num_workers) as executor:
        results = executor.map(worker_readcsv, filepathlist)
    return results

def main():
    #############################################################################################
    # read the kraken2 results and get a sense of the distribution of LCA over hamming distance##
    #############################################################################################
    #read hd_LCA results
    print("Parsing all Kraken2 results and extracting hd LCA columns...")
    hd_LCA = read_kraken2_res(args.kraken2_res_folder)
    concatenated_list = []
    for lst in hd_LCA:
        concatenated_list.extend(lst)
    
    # get unique taxids
    print("Getting unique taxids...")
    counted_elements = Counter(concatenated_list)
    unique_taxids = list(set([value.split(':')[0].split('_')[1] for value in counted_elements.keys()]))
    unique_taxids.remove("0") # remove the taxid 0
    #translate unique_taxids to taxonomy
    print("Translating unique taxids to taxonomy...")
    chunked_taxids = chunk_list(unique_taxids, 100)
    taxid2taxonomy_dict = {}
    failed_taxids = []
    for l in chunked_taxids:
        start_time = time.time()
        mydict = taxidlist2taxonomy(taxid=l)
        if mydict is None:
            mydict = taxidlist2taxonomy(taxid=l) # try again
        if mydict is not None:
            taxid2taxonomy_dict.update(mydict)
        else:
            failed_taxids.extend(l)
            print(f"failed to get taxonomy for taxids: {l}")
        end_time = time.time()
        elapsed_time = end_time - start_time
        if elapsed_time < 0.3:
            time.sleep(0.3)
    taxid2taxonomy_dict["0"]="NA"
    # write the taxid2taxonomy_dict to a file
    print("Writing the taxid2taxonomy_dict to a file...")
    outfile="taxid2taxonomy_dict.txt"
    with open(outfile, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile, delimiter='\t')
        writer.writerow(['taxid', 'taxonomy'])
        for key, value in taxid2taxonomy_dict.items():
            writer.writerow([key, value.replace("\n","")])
    print(f"Done! wrote to {outfile}")

    # Create a list of unique values and their counts + taxonomy
    unique_counts = [(value, count, taxid2taxonomy_dict[value.split(':')[0].split('_')[1]].replace("\n","")) for value, count in counted_elements.items()]
    # Write the list to the file, add taxnomy 
    print("Writing the hammingDistance_LCA_taxonomy_counts to the file...")
    outfile="hd_LCA_taxonomy_counts.txt"
    with open(outfile, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile, delimiter='\t')
        writer.writerow(['hammingDistance_taxid', 'Count', "taxonomy"])  # Write the header
        writer.writerows(unique_counts)  # Write the rows
    print(f"Done! wrote to {outfile}")

if __name__ == "__main__":
    main()
