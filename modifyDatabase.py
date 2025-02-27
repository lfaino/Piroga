from Bio import SeqIO
from ete4 import NCBITaxa
from Bio import Entrez
from tqdm import tqdm
from collections.abc import Iterable




Entrez.email = 'luigi.faino@uniroma1.it'
Entrez.api_key= "50c28bfb47188bbc324aeb93d6e6e3a1a909"

ncbi = NCBITaxa()
#ncbi.update_taxonomy_database()
#file = '/data3/paperCrosara/database/rrnDB-5.8_16S_rRNA/rrnDB-5.8_16S_rRNA.fasta'
#file_out = '/data3/paperCrosara/database/rrnDB-5.8_16S_rRNA/rrnDB-5.8_16S_rRNA.monica.fasta'
#file = '/data3/paperCrosara/database/sh_general_release_dynamic_s_all_04.04.2024_dev.fasta'
#file_out = '/data3/paperCrosara/database/sh_general_release_dynamic_s_all_04.04.2024_dev.monica.fasta'
# file = '/data3/paperCrosara/database/silva/species_taxid.fasta'
# file_out = '/data3/paperCrosara/database/silva/species_taxid.monica.fasta'

file_out = "/data3/paperCrosara/database/Refseq_NCBI/taxonomy/nucl_gb.accession2taxid"

def get_taxid(species):
    species = species.replace(" ", "+").strip()
    search = Entrez.esearch(term = species, db = "taxonomy", retmode = "xml")
    record = Entrez.read(search)
    return record['IdList'][0]

def ncbi_fetch(name):
    taxid = get_taxid(name)
    return(taxid)

def is_part_and_not_last(strings_key, target_string_dict): #(key, dict_rank)
    taxid_to_remove = []
    for key in target_string_dict:
        if strings_key in target_string_dict[key] and strings_key != target_string_dict[key][-1]:
            taxid_to_remove.append(strings_key)
    if len(taxid_to_remove):
        return (taxid_to_remove)
    else:
        return ""

def get_desired_ranks(taxid):

    value_taxid = taxid
    taxid_dict = {}
    try:
        lineage = ncbi.get_lineage(value_taxid)
        #lineage_string = "_".join(map(str, lineage))
        lineage2ranks = ncbi.get_rank(lineage)
        ranks2lineage = dict((rank, value_taxid) for (value_taxid, rank) in lineage2ranks.items())
        #ranks2lineage_select = {'{}'.format(rank): ranks2lineage.get(rank, '<not present>') for rank in desired_ranks}
        for key in ranks2lineage:
            if ranks2lineage[key] != "<not present>":
                ranks2lineage[key] = list(ncbi.get_taxid_translator([ranks2lineage[key]]).values())[0]
        taxid_dict[value_taxid] = ranks2lineage

        return(taxid_dict,lineage)
    except:
        taxid_dict = ""
        lineage_string = ""
        print(str(taxid))
        return(taxid_dict,lineage_string)

def remove_contained_strings(strings):
    filtered_strings = []

    for s in strings:
        # Check if the string is not contained in any other string
        if not any(s != other and s in other for other in strings):
            filtered_strings.append(s)

    return filtered_strings

# list_to_use = []
# #desired_ranks = ["species","genus","family","order","class","phylum","clade","superkingdom","subspecies","species subgroup","species group","strain"]
# list_taxid = []
# with open(file_out, 'r') as fh:
#     for line in fh:
#         id, acc, taxid, gi = line.rstrip().split("\t")
#         if not "taxid" in taxid:
#             list_taxid.append(taxid)


# list_uniq = list(set(list_taxid))
list_uniq= [189395]#,346,611301]

dict_rank = {}
section_all= []
dict_lineage= {}



for taxid in tqdm(list_uniq):
    rank, lineage2ranks= get_desired_ranks(taxid)#, desired_ranks)
    if rank != "":
        for data in lineage2ranks:
            #if data == 4751:
            dict_rank[taxid] = rank
            dict_lineage[taxid] = lineage2ranks
            break

# for key in dict_lineage:




for key in dict_rank:
    if isinstance(dict_rank[key], Iterable):
        for section in dict_rank[key]:
            for group in dict_rank[key][section]:
                section_all.append(group)
section_all_unique = list(set(section_all))
#section_all_unique.remove("no rank")

with open("/tmp/testTaxonomy.all.csv", "w") as f:
    f.write("taxid\t" + '\t'.join(section_all_unique) + "\n")
    for taxid_key in dict_rank:
        if isinstance(dict_rank[taxid_key], Iterable):
            for section in dict_rank[taxid_key]:
                ranks2lineage_select = {'{}'.format(rank): dict_rank[taxid_key][section].get(rank, '<not present>') for rank in section_all_unique}
                ordered_list = [ranks2lineage_select[key] for key in section_all_unique]
                ordered_list_taxid = [str(taxid_key)] + ordered_list
                final = '\t'.join(ordered_list_taxid)
                f.write(final + "\n")
list_key =[]
for key in dict_lineage:
    list_key.append(key)

for key in list_key:
    result = is_part_and_not_last(key, dict_lineage)
    if result != "":
        for remove_taxid in result:
            del dict_lineage[remove_taxid]

with open("/tmp/testTaxonomy.uniq.csv", "w") as fu:
    for key in dict_lineage:
        fu.write(str(key) + "\n")

print("done")
