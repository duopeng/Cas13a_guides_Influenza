import os
import argparse
import concurrent.futures
import pandas as pd
import traceback

parser = argparse.ArgumentParser(description='process the tiled windows and count the coverage of segments and subtypes for each guide')
parser.add_argument('--folder_path_prefix', type=str, required=True,
                    help='The prefix of the folder path, for example Cas13a_SARS_CoV2_chunk_ will match both Cas13a_SARS_CoV2_chunk_1 and Cas13a_SARS_CoV2_chunk_2')
parser.add_argument('--genome_df', type=str, required=True, help='The path to the genome dataframe')
parser.add_argument('-n', '--num_workers', type=int, default=10)
parser.add_argument('-o', '--outputdir', type=str, default=".")
args = parser.parse_args()

#load genome dataframe and generate genomeID-to-segment lookup dictionary
df = pd.read_csv(args.genome_df, sep="\t")
genomeID2segment = dict(zip(df["sequence_ID"], df["segment"]))
genomeID2spp = dict(zip(df["sequence_ID"], df["spp"]))
genomeID2subtype = dict(zip(df["sequence_ID"], df["subtype"]))
genomeID2strain = dict(zip(df["sequence_ID"], df["strain"]))

def isNaN(num):
    return num!= num

for key in genomeID2subtype.keys():
    if isNaN((genomeID2subtype[key])):
        genomeID2subtype[key] = "NA"



def main():
    print("processing chunk folders...")
    folder_list = get_folder_list(args.folder_path_prefix)
    #process each folder in parallel
    with concurrent.futures.ProcessPoolExecutor(max_workers=args.num_workers) as executor:
        results = executor.map(get_guide_dict, folder_list)

    for result in results:
        print(result)

    # merge all guide_dict files
    print("Merging guide frequency...")
    # the guides are arranged by segment -> guide -> count
    guide_dict = {} # guide_dict[segment][guide] = count
    guide_info = {}
    for folder_path in folder_list:
        parent_dir = os.path.dirname(folder_path)
        folder_name = os.path.basename(folder_path)
        with open(os.path.join(parent_dir,f"{folder_name}.guide_table.tsv"), 'r') as file:
            for line in file:
                line = line.strip()
                fields = line.split("\t")
                field0, guide, seg_cov_str = fields[0], fields[1], fields[2]         
                #update guide_dict and guide_info
                update_dict2(guide_dict, field0, guide, seg_cov_str)
                guide_info[guide] = fields[3:]

    #write guide_dict to file, order by dict value in decending order
    with open(f"{args.outputdir}/guides_table.tsv", 'w') as file:
        for key1 in guide_dict.keys():
            for key2, value in sorted(guide_dict[key1].items(), key=lambda item: item[1], reverse=True):
                info = "\t".join(guide_info[key2])
                file.write(f"{key1}\t{key2}\t{value}\t{info}\n")
    print("Finished merging guide frequency")


def update_dict(mydict, key1, key2, value):
    '''
    update a dictionary with two lvl keys, increment the value by value
    '''
    if key1 in mydict.keys():
        if key2 in mydict[key1].keys():
            mydict[key1][key2] += value
        else:
            mydict[key1][key2] = value
    else:
        mydict[key1] = {key2: value}
    return mydict

def update_dict2(mydict, key1, key2, value):
    '''
    update a dictionary with two lvl keys, increment the value by value (value is in the format of id1:count1|id2:count2|)
    '''
    if key1 in mydict.keys():
        if key2 in mydict[key1].keys():
            mydict[key1][key2] = merge_id_count_array(mydict[key1][key2], value)
        else:
            mydict[key1][key2] = value
    else:
        mydict[key1] = {key2: value}
    return mydict

def merge_id_count_array(arr1,arr2):
    '''
    merge two id:count arrays, each id:count array is in the format of id1:count1|id2:count2|...
    '''
    # store arr1 in a dictionary
    dict1 = {}
    for item in arr1.split("|"):
        fields = item.split(":")
        if len(fields) >= 2:
            dict1[fields[0]] = int(fields[1])
    # add arr2 to dict1
    for item in arr2.split("|"):
        fields = item.split(":")
        if len(fields) >= 2:
            if fields[0] in dict1.keys():
                dict1[fields[0]] += int(fields[1])
            else:
                dict1[fields[0]] = int(fields[1])
    # convert dict1 back to array
    arr1 = "|".join([f"{key}:{dict1[key]}" for key in dict1.keys()])
    return arr1

def get_guide_dict(folder_path):
    '''
    process a chunk-folder, scan all subdirectories, and reads targets.txt and windows.txt
    '''
    print(f"Counting guide frequency in {folder_path}")
    try:
        guide_dict = {} # guide_dict[seq_id][guide] = count
        segment_dict = {} # segment_dict[segment][guide] = count
        subtype_dict = {} # subtype_dict[guide][subtype] = count
        strain_dict = {} # strain_dict[guide][strain] = count
        guide_windows = {}
        for lv2_dir in os.listdir(folder_path):

            # get guide information
            seq_id = lv2_dir
            segment = genomeID2segment.get(lv2_dir, "NA")
            subtype = genomeID2subtype.get(lv2_dir, "NA")
            strain = genomeID2strain.get(lv2_dir, "NA")

            # read the target file
            target_file = os.path.join(folder_path, lv2_dir, "targets.txt")
            if os.path.isfile(target_file): #check if target.txt exists
                with open(target_file, 'r') as file: # read target.txt, update guide_dict
                    #skip the first line
                    next(file)
                    for line in file:
                        line = line.strip()
                        guide_dict = update_dict(guide_dict, segment, line, 1)  # NOTE: this dict determines how guides are stored (e.g. by segment, seq_id or by strain)
                        segment_dict = update_dict(segment_dict, segment, line, 1) # update segment : guide count
                        subtype_dict = update_dict(subtype_dict, line, subtype, 1) # update guide : subtype count
                        strain_dict = update_dict(strain_dict, line, strain, 1) # update guide : strain count

            # read the windows file
            windows_file = os.path.join(folder_path, lv2_dir, "windows.txt")
            if os.path.isfile(windows_file): #check if windows.txt exists
                with open(windows_file, 'r') as file: # read windows.txt, update guide_windows
                    #skip the first line
                    next(file)
                    for line in file:
                        line = line.strip()
                        fields = line.split("\t")
                        guide_windows[fields[2]] = line
        #write guide_dict to file, order by dict value in decending order
        parent_dir = os.path.dirname(folder_path)
        folder_name = os.path.basename(folder_path)
        with open(os.path.join(parent_dir,f"{folder_name}.guide_table.tsv"), 'w') as file, open(os.path.join(parent_dir,f"{folder_name}.guide_table_log.txt"), 'w') as log_file:
            i = 0
            for key1 in guide_dict.keys(): # key1 = strain
                for key2 in guide_dict[key1]: # key2 = guide

                    #calculate subtype coverage
                    subtype_coverage = list(set(subtype_dict[key2]))
                    if "NA" in subtype_coverage:
                        subtype_coverage.remove("NA")
                    subtype_coverage = str(len(subtype_coverage))

                    #get segment coverage, format: seg1:count | seg2:count | ...
                    guide = key2 # guide
                    valid_keys = []
                    for key in segment_dict.keys(): # Filter out keys that cannot be converted to integers
                        try:
                            # Attempt to convert the key to an integer
                            int_key = int(key)
                            valid_keys.append(int_key)
                        except (ValueError, TypeError):
                            # Ignore keys that cannot be converted to integers
                            pass
                    sorted_keys = sorted(valid_keys, key=lambda x: int(x))
                    segment_coverage = '|'.join([f'{seg}:{segment_dict[seg][guide]}' for seg in sorted_keys if guide in segment_dict[seg]])

                    #calculate strain coverage
                    strain_coverage = list(set(strain_dict[key2]))
                    if "NA" in strain_coverage:
                        strain_coverage.remove("NA")
                    strain_coverage = str(len(strain_coverage))

                    i += 1
                    if i % 10000 == 0:
                        log_file.write(f"processed {i} guides\n")

                    # write guide to file
                    file.write(f"{key1}\t{key2}\t{segment_coverage}\t{subtype_coverage}\t{strain_coverage}\t{guide_windows[key2]}\n") if key2 in guide_windows.keys() else file.write(f"{key1}\t{key2}\t{guide_dict[key1][key2]}\t\t\t\t\t\t\t\t\n")
        return f"Finished processing {folder_path}"
    except Exception as e:
        # Handle the exception
        # You can choose to raise a custom exception, log the error, etc.
        traceback_msg = traceback.print_exc()
        return f"Error occurred: {str(e)}\n traceback: {traceback_msg}"

def get_folder_list(folder_path_prefix):
    """get a list of folders that contains the prefix
    Args:
        folder_path_prefix (_type_): path to the parent folder to be scanned
    Returns:
        list: a list of folder paths
    """
    folder_list = []
    parent_dir = os.path.dirname(folder_path_prefix)
    prefix = os.path.basename(folder_path_prefix)
    for lv1_dir in os.listdir(parent_dir):
        if prefix in lv1_dir and os.path.isdir(os.path.join(parent_dir,lv1_dir)): # go through all folders that contains prefix
            folder_list.append(os.path.join(parent_dir, lv1_dir))
    return folder_list


if __name__ == "__main__":
    main()