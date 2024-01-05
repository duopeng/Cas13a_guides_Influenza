import re

#define functions to parse gene, segment names, subtypes types
def parse_subtypes(text):
    """parse the subtypes
    """
    match1 = re.search(r"(H\dN\d)", text, re.IGNORECASE)
    if match1:
        return match1[1].upper()
    else:
        return None
def parse_gene_name(text):
    """Parse gene name from the description"""
    gene_name = None
    
    match1 = re.search(r"(hamagglutinin|heamagglutinin|hemegglutinin|hemagglutin|hemaggulutinin|haemagglutinin|hemagglutinin|neuraminidase|matrix protein|nucleoprotein|neuraminidase|matrix protein|polymerase protein 2|polymerase protein 3|membrane protein 1|membrane protein 2|nonstructural protein)", text, re.IGNORECASE)
    if match1:
        gene_name = match1[1]
        return gene_name
    
    match2 = re.search(r"virus(.+) (?:gene|mrna|rna|pseudogene)", text, re.IGNORECASE)
    match3 = re.search(r"(?:gene|mrna|rna|pseudogene) for ([\d\w ]+)[,\.\(]", text, re.IGNORECASE)


    to_search = ""
    if match2 and match3:
        to_search = match2[1] + match3[1]
    elif match2:
        to_search = match2[1]
    elif match3:
        to_search = match3[1]
    
    match_gene_symbol = re.search(r"(PB2|PB 2|Pb2| pb2 |PB1-F2|Pb1|PB1|PB 1| pb1 |PA| pa |HA|NP|NA|M1|M2|NS1|NEP|NS2|NS|P3| ns2 | ns1 | pb2 | p3 | ns )", to_search)
    
    if match_gene_symbol:
        gene_name = match_gene_symbol[1]
        return gene_name
    else:
        return None

def parse_segment_name(text):
    segment_name = None
    match1 = re.search(r"(?:segment|seg|segment:) ([\d])", text, re.IGNORECASE)
    
    if match1:
        segment_name = match1[1]
    return segment_name

def update_dict(mydict, key, value):
    if key in mydict:
        mydict[key].append(value)
    else:
        mydict[key] = [value]
    return mydict

def parse_spp(text):
    match = re.search(r'(Influenza\s.\svirus)', text, re.IGNORECASE)
    if match:
        spp_name = match[1].lower()
        return spp_name
    else:
        return None
    
# define a dictionary to collapse gene names
gene_name_collapse = {
 ' ns ': '8',
 ' pa ': '3',
 ' pb1 ': '2',
 ' pb2 ': '1',
 'HA': '4',
 'HEMAGGLUTIN': '4',
 'Haemagglutinin': '4',
 'Hemagglutin': '4',
 'M1': '7',
 'M2': '7',
 'Matrix Protein': '7',
 'Matrix protein': '7',
 'NA': '6',
 'NEP': '8',
 'NP': '5',
 'NS': '8',
 'NS1': '8',
 'NS2': '8',
 'Neuraminidase': '6',
 'Nonstructural Protein': '8',
 'Nonstructural protein': '8',
 'Nucleoprotein': '5',
 'P3': '1',
 'PA': '3',
 'PB 1': '2',
 'PB1': '2',
 'PB1-F2': '2',
 'PB2': '1',
 'haemagglutinin': '4',
 'hamagglutinin': '4',
 'heamagglutinin': '4',
 'hemagglutin': '4',
 'hemaggulutinin': '4',
 'hemegglutinin': '4',
 'matrix protein': '7',
 'membrane protein 1': '7',
 'neuraminidase': '6',
 'nonstructural protein': '8',
 'nucleoprotein': '5',
 'polymerase protein 2': '1',
 'polymerase protein 3': '1'
}

genome_pattern = re.compile(r'\((?P<strain_name>[\w\d\s\-\/]+)\((?P<subtype>H\dN\d)\)\)')
genome_pattern2 = re.compile(r'\((?P<strain_name>[\w\d\s\-\/]+)\)')

def extract_genome_strain(input_str):
    # returns: strain_name, subtype, gene_name, seg_name 
    match = re.search(genome_pattern, input_str)
    match2 = re.search(genome_pattern2, input_str)
    if match:
        #print("Strain name:", match.group("strain_name")) 
        #print("Subtype:", match.group("subtype"))  
        gene_name = parse_gene_name(input_str)
        seg_name = parse_segment_name(input_str)
        return match.group("strain_name"), match.group("subtype"), gene_name, seg_name
    elif match2:
        #print("Strain name:", match2.group("strain_name")) 
        gene_name = parse_gene_name(input_str)
        seg_name = parse_segment_name(input_str)
        return match2.group("strain_name"), None, gene_name, seg_name
    else:
        #print("No match!!")
        return None, None, None, None
    
def dict_increment(mydict, key):
    if key in mydict:
        mydict[key] += 1
    else:
        mydict[key] = 1
    return mydict

def dict_update(mydict, k1, k2, v):
    if k1 in mydict:
        mydict[k1][k2] = v
    else:
        mydict[k1] = {k2:v}
    return mydict

def desc2seg(row):
    desc = row["description"]
    segment_name = parse_segment_name(desc)
    if segment_name:
        return segment_name
    else:
        gene_name = parse_gene_name(desc)
        if gene_name:
            segment_name = gene_name_collapse[gene_name]
            return segment_name
        else:
            return None
        
def desc2subtype(row):
    desc = row["description"]
    subtype_name = parse_subtypes(desc)
    if subtype_name:
        return subtype_name
    else:
        return None
    
def desc2strain(row):
    desc = row["description"]
    strain_name, subtype_name, gene_name, segment_name = extract_genome_strain(desc)
    if strain_name and len(strain_name) > 7 and " form" not in strain_name and "form " not in strain_name:
        return strain_name
    else:
        return None