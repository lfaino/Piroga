#! /usr/bin/python
import sys
import uuid
import random
import string
from matplotlib import pyplot as plt
import matplotlib
import subprocess as sb
import re
import Bio
from collections import Counter
import gzip
from Bio.Align.Applications import MafftCommandline
from itertools import chain
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
from Bio.Align.AlignInfo import SummaryInfo
import os
from Bio import Entrez
from Bio import SeqIO
from Bio.Blast.Applications import NcbimakeblastdbCommandline
from multiprocessing import Pool
from io import StringIO
from difflib import SequenceMatcher
from tqdm import tqdm
import tempfile
import argparse
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap
import csv
from ete3 import NCBITaxa
import requests
import xml.etree.ElementTree as ET

import warnings
warnings.filterwarnings("ignore")


matplotlib.style.use('ggplot')
KTIMPORTTAX = "ktImportTaxonomy -tax %s -o %s %s"
MINIMAP = "minimap2 -x ava-ont -t %s %s %s | gzip -1 > %s"
MINIMAP_S = "minimap2 -a -x map-ont --secondary=no -t %s %s %s "
MINIMAP_SY = "minimap2 -a -x map-ont -t %s %s %s "
SAMSORT = "samtools sort -o %s -O BAM"
SAMINDEX = "samtools index %s"
RACON = "racon -t %s -f %s %s %s"
PORECHOP = "porechop -i %s -t %s -o %s"
JELLYFISH_COUNT= "jellyfish count -m 100 -s 100M -t 10 -o /dev/stdout -C %s"
JELLYFISH_DUMP= "jellyfish dump -o %s /dev/stdin"
SPADES = "spades -s %s -o %s --only-assembler"
NANOPOLISHV = "nanopolish variants --consensus -o %s -w %s -r %s -b %s -g %s -t %s --min-candidate-frequency 0.1 -p 1 " \
             "--fix-homopolymers"
NANOPOLISHI = "nanopolish index -d %s %s "
NANOPOLISHVA = "nanopolish vcf2fasta -g %s %s"
BWAI = "bwa index %s"
BWA = "bwa mem -t %s %s %s"
BCFTOOLS_MP = "bcftools mpileup -Ou -f %s %s"
BCFTOOLS_CALL = "bcftools call -mv -o %s --ploidy 1"
BCFTOOLS_IN = "bcftools index %s"
BCFTOOLS_CO = "bcftools consensus -f %s -o %s %s"
FREEBAYES = "freebayes -f %s -p 1 %s"
KROCUS_D = "krocus_database_downloader --species %s"
BGZIP = "bgzip %s"
TABIX = "tabix -p vcf %s"

mlst_database = "https://pubmlst.org/data/dbases.xml"

ncbi = NCBITaxa()

def setting():
    parser = argparse.ArgumentParser(prog='quantify and detect pathogens', usage='%(prog)s [options]')
    parser.add_argument("-o","--output", nargs="?", default="output", help="output name")
    parser.add_argument("-r", "--retmax", nargs="?", type = int, default="100000", help="number of id to retrieve per itaration")
    parser.add_argument("--barcode", help="barcode to analyse", required=True)
    parser.add_argument("-f5","--folder_fast5", help="folder containing fast5 files",default="/data/sharedata/")
    parser.add_argument("-t", "--threads", nargs="?", type = int, default="10", help="number of threads")
    parser.add_argument("-d", "--database", nargs="?", choices=['Plant_Bacteria_Funghi.nal', 'nt.nal', 'ITS_16S_18S_28S_LSU_SSU.nal',
                                                                'Metazoa_Organism_OR_metazoa_All_Fields_.nal','Xanthomonas_genomes.nal','curatedXylellaDatabase.nal',
                                                                'customITSreduced.nal','customITSreducedReduced.nal','BacteriaFungiVitis.nal',
                                                                '50_4000RefSeq_biomol_genomic_ITS_16S_18S_28S_LSU_SSU.nal','bacterialGenome2024_01_24.fna.nal'],
                        default="",
                        help="database name; can be a fasta file or a subset of ncbi database; Default is Plant_Bacteria_Funghi")
    parser.add_argument("-de", "--database_external", nargs="?")
    parser.add_argument("-tmp", "--temp", nargs="?", default="/tmp/")
    parser.add_argument("--folder_fastq", help='folder where fastq file are located', required=True,
                        default="/data/sharedata/")#, widget='DirChooser')
    parser.add_argument("-e", "--email", nargs="?", default="", help="email for blast search", required=True)
    parser.add_argument("-s", "--search", nargs="?", default="", help="ncbi search to retrieve GIs")
    parser.add_argument("-sp", "--species", nargs="?", default="", help="PUBMLST species to use")
    parser.add_argument("-m", "--min",  default="10", help="minimum number of reads to plot a species as "
                                                           "fraction of total mappd reads [0-100]")
    parser.add_argument("-p", "--percentage", type = float, default="90", help="minimum identity to consider "
                                                                               "a valid alignment")
    parser.add_argument("-mr", "--min_reads", type = float, default="100", help="minimum number of reads to "
                                                                                "sequence to consider a sample in the "
                                                                                "analysis")
    parser.add_argument("-int", "--interval", default="100-20000", help="minimum and max read length"
                                                                                " to use in the correction step; default [100-20000]")
    parser.add_argument("-su", "--subset", default="0", type = int, help="number of reads to use for each barcode [0]")
    parser.add_argument("-mrc", "--num_reads_cor", type=int, default="0",
                        help="minimum number of reads to use in overlapping for reads correction [0]. 0 = no filtering")

    parser.add_argument("-a", "--assemble", action='store_true', help="assembled-reads", default=False)
    parser.add_argument("-ua", "--use_assembled", action='store_true', help="use assembled reads for "
                                                                            "discovery", default=False)
    parser.add_argument("-dh", "--database_history", action='store_true', help="", default=False)
    parser.add_argument("-c", "--correct", action='store_true', help="correct or not reads via racon",
                        default=False)
    parser.add_argument("-mafft", "--mafft", action='store_true', help="correct or not reads via racon",
                        default=False)
    parser.add_argument("-u", "--update", action='store_true', help="update database", default=False)#nargs="?", default="", required=True)
    parser.add_argument("-v",'--verbose', help='be verbose', dest='verbose', action='store_true', default=False)


    args = parser.parse_args()
    return args



def get_desired_ranks(taxid, desired_ranks):
    lineage = ncbi.get_lineage(taxid)
    lineage2ranks = ncbi.get_rank(lineage)
    ranks2lineage = dict((rank, taxid) for (taxid, rank) in lineage2ranks.items())
    return {'{}_id'.format(rank): ranks2lineage.get(rank, '<not present>') for rank in desired_ranks}

# def main(taxids, desired_ranks, path):
#     with open(path, 'w') as csvfile:
#         fieldnames = ['{}_id'.format(rank) for rank in desired_ranks]
#         writer = csv.DictWriter(csvfile, delimiter='\t', fieldnames=fieldnames)
#         writer.writeheader()
#         for taxid in taxids:
#             writer.writerow(get_desired_ranks(taxid, desired_ranks))


def best_align_match(samfile, iteration):
    if iteration > 0:
        samDict = {}
        count = 0
        filtered = samfile + ".reduced.paf"
        with gzip.open(samfile, "r") as sam, open(filtered, "w") as reduced:
            for line in sam:
                line = line.decode()
                paf = line.split("\t")
                if paf[0] in samDict:
                    if count < iteration:
                        reduced.write(line)
                        count += 1
                else:
                    samDict[paf[0]] = "0"
                    reduced.write(line)
                    count = 0
    else:
        filtered = samfile

    return filtered

def get_consensus_seq(filename):
    common_alignment = MultipleSeqAlignment(
        chain(*AlignIO.parse(filename, "fasta"))
    )
    summary = SummaryInfo(common_alignment)
    consensus = summary.dumb_consensus(0.7, "N")
    return consensus


def consensus_seq(file, fastq, set_count):
    print(fastq.name)
    record_dict = SeqIO.to_dict(SeqIO.parse(fastq.name, "fasta"))
    samDict = {}
    count = 0
    with gzip.open(file, "r") as sam:
        for line in sam:
            line = line.decode()
            paf = line.split("\t")
            if paf[0] in samDict:
                if count < set_count:
                    samDict[paf[0]].append((paf[5], paf[4]))
                    count += 1
            else:

                samDict[paf[0]] = [(paf[0], "+"), (paf[5], paf[4])]
                count = 0
    files_to_align = []
    for key in samDict:
        fasta = tempfile.NamedTemporaryFile(suffix=".fasta", delete=False, mode="w")
        list_record = samDict[key]
        for pair in list_record:
            id, strand = pair
            if id in record_dict:
                if strand == "+":
                    SeqIO.write(record_dict[id], fasta, "fasta")
                else:
                    SeqIO.write(record_dict[id].reverse_complement(), fasta, "fasta")
        files_to_align.append(fasta.name)
    result_list = []
    with Pool(processes=50) as pool:
        for result in tqdm(pool.imap(func=mafft, iterable=files_to_align), total=len(files_to_align)):
            result_list.append(result)
    outfile = fastq.name + ".fasta"
    with open(outfile, "w") as fh:
        for line in result_list:
            fh.write(line)

    print(outfile)
    return outfile


def mafft(fasta):
    aln = tempfile.NamedTemporaryFile(suffix=".fasta", delete=False, mode="w")
    mafft_cline = MafftCommandline(input=fasta)
    stdout, stderr = mafft_cline()
    aln.write(stdout)
    cons = get_consensus_seq(aln.name)
    with open(fasta) as f:
        first_line = f.readline()
    consensus = first_line.rstrip() + "\n" + str(cons) + "\n"
    return (consensus)


def convert_color(test):
    values = dict()
    for index, letter in enumerate(string.ascii_letters):
        values[letter] = index + 1
    complete_colors = []
    for name_sp in test:
        if " " in name_sp:
            name_a = name_sp.split(" ")
            genera = name_a[0]
            del name_a[0]
            specie = "".join(name_a)
            name = name_sp
        else:
            genera = name_sp
            specie = "unknown"
            name  = name_sp + "unknown"
        genera = re.sub("[^a-z]+", "", genera.lower())
        specie = re.sub("[^a-z]+", "", specie.lower())
        name = re.sub("[^a-z]+", "", name.lower())
        genera = list(genera)
        specie = list(specie)
        name = list(name)
        rgb = [genera, name, specie]
        number_rgb = []
        for i in rgb:
            number  = sum([values[l] for l in i])
            random.seed(number)
            number_rgb.append(random.randint(0,255)/255)
        complete_colors.append([name_sp ,tuple(number_rgb)])
    return(complete_colors)

def blast(elm):
    blastn_cline = Bio.Blast.Applications.NcbiblastnCommandline(db=elm[1], query=elm[0], evalue=0.001, out=elm[2] ,outfmt = "6 qseqid sseqid bitscore sscinames pident evalue staxids qlen" )
    try:
        blastn_cline()
        if not verbose:
            os.remove(elm[0])
    except ValueError:
        print(blastn_cline)
    return(elm[2])

def common_start(*s):
   l = len(s)
   if l == 0:
       return None
   elif l == 1:
       return s[0]
   start = s[0]
   while start != "" and not all(ii.startswith(start) for ii in s[1:]):
       start = start[:-1]

   return start

def download_fasta(elm):
    search, retmax, retstart, email = elm
    Bio.Entrez.email = email
    handle = Bio.Entrez.esearch(db="nucleotide", term=search, retmax=retmax, retstart=retstart)
    record = Bio.Entrez.read(handle)
    records = str("\n".join(record["IdList"]))
    return (records)

def counpute_count(dict_match):
    genera_dict = {}
    species_dict = {}
    id, match_blast = dict_match
    if len(match_blast) > 1:
        substring_counts = {}
        length_array = int(len(match_blast))
        for i in range(0, length_array):
            for j in range(i + 1, length_array):
                string1 = match_blast[i]
                string2 = match_blast[j]
                match = SequenceMatcher(None, string1, string2).find_longest_match(0, len(string1), 0, len(string2))
                matching_substring = string1[match.a:match.a + match.size]
                if (matching_substring not in substring_counts):
                    substring_counts[matching_substring] = 1
                else:
                    substring_counts[matching_substring] += 1
        species = sorted(substring_counts.items(), key=lambda x: (x[1], x[0]), reverse=True)
        specie_id = True
        genera_id = True
        for organism in species:
            if " " in organism[0] and organism[0][0].isupper() and organism[0] in match_blast:
                fields = organism[0].split(" ")
                genera = fields[0]
                specie = organism[0]
                if genera != "" and specie != "" and specie_id:
                    if organism[0] in species_dict and specie_id:
                        specie_id = False
                        species_dict[organism[0]].append(id)
                        genera_dict[genera].append(id)
                        genera_id = False
                    else:
                        specie_id = False
                        species_dict[organism[0]] = [id]
                        genera_id = False
                        genera_dict[genera] = [id]
                elif genera != "" and specie == "" and genera_id:
                    if genera in genera_dict:
                        genera_dict[genera].append(id)
                        genera_id = False
                    else:
                        genera_id = False
                        genera_dict[genera] = [id]
    else:
        fields = match_blast[0].split(" ")
        genera = fields[0]
        species = match_blast[0]
        if genera in genera_dict:
            genera_dict[genera].append(id)
        else:
            genera_dict[genera] = [id]
        if species in species_dict:
            species_dict[species].append(id)
        else:
            species_dict[species] = [id]
    return([species_dict, genera_dict])

def racon(fastq, threads, cwd, num_reads_cor, mafft):
    sam = tempfile.NamedTemporaryFile(dir=dirpathTEMP,suffix=".paf.gz", delete=False)
    reads = tempfile.NamedTemporaryFile(dir=dirpathTEMP,suffix=".fasta", delete=False)
    print("RUNNING MINIMAP")
    m = MINIMAP % (threads, fastq, fastq, sam.name)
    print(m)
    minimap = sb.Popen(m, shell=True, cwd=cwd, stderr=sb.PIPE)
    minimap.communicate()
    filtered = best_align_match(sam.name, num_reads_cor)
    print("RUNNING RACON")
    r = RACON % (threads, fastq, filtered , fastq)
    print(r)
    racon_cmd = sb.Popen(r, shell=True, cwd=cwd, stdout=reads, stderr=sb.PIPE)
    racon_cmd.communicate()
    if not verbose:
        os.remove(sam.name)
    if os.path.exists(filtered):
        if not verbose:
            os.remove(filtered)
    sam = tempfile.NamedTemporaryFile(dir=dirpathTEMP,suffix=".paf.gz",delete=False)
    print("RUNNING MINIMAP")
    m = MINIMAP % (threads, reads.name, fastq, sam.name)
    print(m)
    minimap = sb.Popen(m, shell=True, cwd=cwd, stderr=sb.PIPE)
    minimap.communicate()
    filtered = best_align_match(sam.name, num_reads_cor)
    print("RUNNING RACON")
    r = RACON % (threads, fastq, filtered, reads.name)
    print(r)
    output = tempfile.NamedTemporaryFile(dir=dirpathTEMP,suffix=".fasta", delete=False)
    racon_cmd = sb.Popen(r, shell=True, cwd=cwd, stdout=output, stderr=sb.PIPE)
    racon_cmd.communicate()
    sam = tempfile.NamedTemporaryFile(dir=dirpathTEMP,suffix=".paf.gz", delete=False)
    print("RUNNING MINIMAP")
    m = MINIMAP % (threads, output.name, output.name, sam.name)
    print(m)
    minimap = sb.Popen(m, shell=True, cwd=cwd, stderr=sb.PIPE)
    minimap.communicate()
    if num_reads_cor == 0:
        num_reads_cor = 100
    outputCons = consensus_seq(sam.name, output, num_reads_cor)
    if not verbose:
        os.remove(sam.name)
    if os.path.exists(filtered):
        if not verbose:
            os.remove(filtered)
    if not verbose:
        os.remove(reads.name)
    if verbose:
        print("Racon output is:" + outputCons)

    return(outputCons)

def assembly(output, cwd, threads, read_mapping, od, fastq, fast5):
    jfc_out = tempfile.NamedTemporaryFile(dir=dirpathTEMP,suffix=".fasta", delete=False)
    print("RUNNING JELLYFISH")
    jfc = JELLYFISH_COUNT % output
    jfd = JELLYFISH_DUMP % jfc_out.name
    print(jfc)
    print(jfd)
    jellyfishc = sb.Popen(jfc, shell=True, cwd=cwd, stdout=sb.PIPE)
    jellyfishd = sb.Popen(jfd, shell=True, cwd=cwd, stdin=jellyfishc.stdout)
    jellyfishd.communicate()
    count = 0
    kmer = dirpathTEMP + uuid.uuid4().hex + ".fasta"
    with open(kmer, "w") as fh:
        for record in SeqIO.parse(jfc_out.name, "fasta"):
            if int(record.id) > 10 and len(record.seq) == 100:
                repetitive = 0
                while 10 >= repetitive:
                    repetitive += 1
                    count += 1
                    record.id = "kmer_" + str(count)
                    SeqIO.write(record, fh, "fasta")
    print(kmer)
    tmp_dir = tempfile.mkdtemp(dir=dirpathTEMP)
    print("RUNNING SPADES")
    spa = SPADES % (kmer, tmp_dir)
    print(spa)
    try:
        spade = sb.Popen(spa, shell=True, cwd=cwd, stderr=sb.PIPE, stdout=sb.PIPE)
        spade.communicate()
        assembled_contigs = os.path.join(tmp_dir, "contigs.fasta")
        bam = tempfile.NamedTemporaryFile(dir=dirpathTEMP,suffix=".bam", delete=False)
        print("RUNNING MINIMAP")
        m = MINIMAP_S % (threads, assembled_contigs, read_mapping)
        ss = SAMSORT % bam.name
        print(m)
        minimap = sb.Popen(m, shell=True, cwd=cwd, stdout=sb.PIPE, stderr=sb.PIPE)
        samtools_sort = sb.Popen(ss, shell=True, cwd=cwd, stdin=minimap.stdout, stdout=sb.PIPE, stderr=sb.PIPE)
        samtools_sort.communicate()
        si = SAMINDEX % bam.name
        samtoos_index = sb.Popen(si, shell=True, cwd=cwd, stdout=sb.PIPE, stderr=sb.PIPE)
        samtoos_index.communicate()
        ni = NANOPOLISHI % (fast5, read_mapping)
        print(ni)
        nanopolish_index = sb.Popen(ni, shell=True, cwd=cwd, stdout=sb.PIPE, stderr=sb.PIPE)
        nanopolish_index.communicate()
        regions = []
        print("RUNNING NANOPOLISH")
        dict_vcf_fasta = []
        for record in SeqIO.parse(assembled_contigs, "fasta"):
            region = record.id + ":0-" + str(len(record.seq))
            vcf = tempfile.NamedTemporaryFile(dir=dirpathTEMP,prefix="polished.", suffix=".vcf", delete=False)
            nv = NANOPOLISHV % (vcf.name, region, read_mapping, bam.name, assembled_contigs, threads)
            print(nv)
            nanopolish_var = sb.Popen(nv, shell=True, cwd=cwd, stdout=sb.PIPE, stderr=sb.PIPE)
            nanopolish_var.communicate()
            if os.path.isfile(vcf.name) and os.path.getsize(vcf.name) > 0:
                regions.append(vcf.name)
            else:
                dict_vcf_fasta.append(record)
        nva = NANOPOLISHVA % (assembled_contigs, " ".join(regions))
        print(nva)
        mlst_done = os.path.join(od, fastq + ".MLST.done.fasta")
        mlst = os.path.join(assembled_contigs + ".MLST.fasta")
        with open(mlst, "w") as fh:
            nanopolish_vcf = sb.Popen(nva, shell=True, cwd=cwd, stdout=fh, stderr=sb.PIPE)
        nanopolish_vcf.communicate()
        bam = tempfile.NamedTemporaryFile(dir=dirpathTEMP,suffix=".bam", delete=False)
        bi = BWAI % mlst
        bwa_index = sb.Popen(bi, shell=True, cwd=dirpathTEMP)
        bwa_index.communicate()
        bm = BWA % (threads, mlst, output)
        #shutil.copyfile(mlst, mlst_done) #os.path.join(cwd,output + ".fasta"))
        ss = SAMSORT % bam.name
        print(bm)
        bwa_mem = sb.Popen(bm, shell=True, cwd=dirpathTEMP, stdout=sb.PIPE, stderr=sb.PIPE)
        samtools_sort = sb.Popen(ss, shell=True, cwd=cwd, stdin=bwa_mem.stdout, stdout=sb.PIPE, stderr=sb.PIPE)
        samtools_sort.communicate()
        si = SAMINDEX % bam.name
        samtoos_index = sb.Popen(si, shell=True, cwd=cwd, stdout=sb.PIPE, stderr=sb.PIPE)
        samtoos_index.communicate()
        vcf = tempfile.NamedTemporaryFile(dir=dirpathTEMP,suffix=".vcf", delete=False)
        bcfm = BCFTOOLS_MP % (mlst, bam.name) #FREEBAYES
        bcfc = BCFTOOLS_CALL % vcf.name
        print(bcfm)
        print(bcfc)
        bcftools_mpile = sb.Popen(bcfm, shell=True, cwd=cwd, stdout=sb.PIPE, stderr=sb.PIPE)
        bcftools_call = sb.Popen(bcfc, shell=True, cwd=cwd, stdin=bcftools_mpile.stdout)
        bcftools_call.communicate()
        bcfi = BCFTOOLS_IN % vcf.name

        #BGZIP = "bgzip %s"
        #TABIX = "tabix -p vcf %s"
        bgzip = BGZIP % (vcf.name)
        vcf_gz = vcf.name + ".gz"
        tabix = TABIX % (vcf_gz)
        bgzip_run = sb.Popen(bgzip, shell=True, cwd=cwd, stdout=sb.PIPE, stderr=sb.PIPE)
        bgzip_run.communicate()
        print(bgzip)
        tabix_run = sb.Popen(tabix, shell=True, cwd=cwd, stdout=sb.PIPE, stderr=sb.PIPE)
        tabix_run.communicate()
        print(tabix)
        print(bcfi)
        bcftools_index = sb.Popen(bcfi, shell=True, cwd=cwd, stdout=sb.PIPE, stderr=sb.PIPE)
        bcftools_index.communicate()
        bcfco = BCFTOOLS_CO % (mlst, mlst_done, vcf_gz)
        print(bcfco)
        #print(mlst_done)
        bcftools_call = sb.Popen(bcfco, shell=True, cwd=cwd, stdout=sb.PIPE, stderr=sb.PIPE)
        bcftools_call.communicate()
        return(mlst_done)

    except:
        mlst_done = ""
        return (mlst_done)

def mlstget(species):

    req = requests.get(mlst_database)
    root = ET.fromstring(req.text)
    mlst_all = []
    for child in root:
        if species in child.text:
            mlst_comb = child.findall("mlst/database/profiles")
            mlst_location = child.findall("mlst/database/loci")
            for item in mlst_comb:
                mlst_pattern_ftp =item.find('url').text
            for item in mlst_location:
                mlst_ftps = item.findall('locus/url')
                for mlst_ftp in mlst_ftps:
                    mlst_seq = requests.get(mlst_ftp.text).text
                    for record in SeqIO.parse(StringIO(mlst_seq), "fasta"):
                        mlst_all.append(record)

    mlst_complete = dirpathTEMP + uuid.uuid4().hex + ".fasta"
    with open (mlst_complete, "w") as fh:
        SeqIO.write(mlst_all, fh, "fasta")
    patternsST=[]
    mlst_pattern = requests.get(mlst_pattern_ftp).text.split("\n")
    for combination in mlst_pattern:
        if combination.startswith("ST"):
            orderST = combination.split("\t")[1:-1]
        else:
            patternsST.append([combination.split("\t")[0],"".join(combination.split("\t")[1:-1])])
    cmd_make = " ".join(["makeblastdb", "-dbtype","nucl","-in" , mlst_complete])
    makedb = sb.Popen(cmd_make, shell=True)
    # print(cline)
    makedb.communicate()
    return(mlst_complete)

def plothystogramSpecies(files_kt_list, blastdb, output, od, dict_species_all, name_table, dict_table_all, figure_file):
    files_kt = " ".join(files_kt_list)
    k = KTIMPORTTAX % (blastdb, output + ".html", files_kt)
    print(k)
    ktimport = sb.Popen(k, shell=True, cwd=od)
    ktimport.communicate()


    species_pd = pd.DataFrame.from_dict(dict_species_all, orient='index')
    species_name_colors = [name for barcode in dict_species_all for name in dict_species_all[barcode]]
    species_name_colors = list(dict.fromkeys(species_name_colors))
    cmap = convert_color(species_name_colors)
    cmap2 = []
    for i in cmap:
        cmap2.append(i[1])
    species_abs = pd.DataFrame.from_dict(dict_table_all, orient='index')
    species_abs.transpose().to_excel(name_table)
    cmap1 = LinearSegmentedColormap.from_list("my_colormap", cmap2)
    species_pd.head()
    species_pd.fillna(0)
    f = plt.figure()
    species_pd.plot(kind="bar", stacked=True, colormap=cmap1, ax=f.gca())
    plt.legend(loc='center left', bbox_to_anchor=(1.0, 0.5))
    plt.title("Abundance species")
    plt.xlabel("barcode")
    plt.ylabel("Relative abundance(%)")
    # plt.show()
    plt.savefig(figure_file)
    return()

def main():
    print("PROGRAM STARTED")
    args = setting()
    global dirpathTEMP
    global dirpathTEMP
    dirpathTEMP = tempfile.mkdtemp(dir=args.temp)
    global verbose
    verbose = args.verbose

    if args.species != "":
        fastamlst = mlstget(args.species)
    cwd = args.folder_fastq
    od = os.getcwd()
    Bio.Entrez.email = args.email
    new_string = re.sub("[^0-9a-zA-Z]+", "_", args.search)
    name_plot = args.output + ".pdf"
    name_table = args.output + ".xlsx"
    figure_file = os.path.join(od, name_plot)
    name_table = os.path.join(od, name_table)
    #os.environ["BLASTDB"] = "/data2/blastdb"
    blastdb = os.environ['BLASTDB']
    outputRacon = ""
    if args.search != "":
        database_name = os.path.join(blastdb, new_string)
        print("Using database %s" % database_name)
        database_name_check = os.path.join(blastdb, new_string + ".nal")
        if not os.path.exists(database_name_check) or args.update:
            handle_n = Bio.Entrez.esearch(db="nucleotide", term=args.search)
            record_n = Bio.Entrez.read(handle_n)
            retmax = args.retmax
            retstart = 0
            repeat = (str(int(record_n['Count'])/retmax)).split(".")[0]
            print("FOUND " + str(record_n['Count']) + " SEQUENCES\n")
            repeat = int(repeat) + 1
            id_download = []
            count = 0
            while (count < repeat):
                count += 1
                id_download.append([args.search,retmax, retstart, args.email])
                retstart = count * retmax
            result_list_tqdm = []
            if repeat > 10:
                threads = 10
            else:
                threads = repeat
            with Pool(processes=threads) as pool:
                for result in tqdm(pool.imap(func=download_fasta, iterable=id_download), total=len(id_download)):
                    result_list_tqdm.append(result)
            records = "\n".join(result_list_tqdm)
            fp = tempfile.NamedTemporaryFile(dir=dirpathTEMP, delete= False)
            fp.write(records.encode())
            cm = "blastdb_aliastool -db nt -dbtype nucl -gilist %s -out %s" % (fp.name, database_name)
            print(cm)
            blastalias = sb.Popen(cm, shell= True)
            blastalias.communicate()
    elif args.database != "":
        fasta_all = []
        if os.path.isabs(args.database):
            if args.database.endswith("fasta"):
                for rec in SeqIO.parse(args.database, "fasta"):
                    desc = rec.description
                    if 200 < len(rec.seq) < 10000:
                        single_elem = desc.split(" ")
                        if not "sp." in single_elem[2]:
                            if len(single_elem) >= 3:
                                species = " ".join([single_elem[1], single_elem[2]])
                                rec.description = species
                                fasta_all.append(rec)
                database_name = os.path.join(blastdb,args.database + "clean.fasta")
                SeqIO.write(fasta_all, database_name , "fasta")
                print("Using database %s" % database_name)

            elif args.database.endswith(".nal") or args.database.endswith(".Nal"):
                database_name = os.path.abspath(args.database)
                interm = database_name.split(".")[:-1]
                database_name = ".".join(interm)
                print("Using database %s" % database_name)
            else:
                print("DATABASE NOT FOUND")
        elif os.path.isfile(os.path.abspath(args.database)):
            database_name_orig = os.path.abspath(args.database)
            if database_name_orig.endswith("fasta"):
                for rec in SeqIO.parse(database_name_orig, "fasta"):
                    desc = rec.description
                    if 200 < len(rec.seq) < 10000:
                        single_elem = desc.split(" ")
                        if not "sp." in single_elem[2]:
                            if len(single_elem) >= 3:
                                species = " ".join([single_elem[1], single_elem[2]])
                                rec.description = species
                                fasta_all.append(rec)
                database_name = os.path.join(blastdb,args.database + "clean.fasta")
                print("Using database %s" % database_name)
                SeqIO.write(fasta_all, database_name , "fasta")
            elif database_name_orig.endswith(".nal") or database_name_orig.endswith(".Nal"):
                interm = database_name_orig.split(".")[:-1]
                database_name = ".".join(interm)
            else:
                print("DATABASE NOT FOUND")
        elif os.path.isfile(os.path.join(blastdb, args.database)):
            database_name_orig = os.path.join(blastdb, args.database)
            if database_name_orig.endswith("fasta"):
                for rec in SeqIO.parse(database_name_orig, "fasta"):
                    desc = rec.description
                    if 200 < len(rec.seq) < 10000:
                        single_elem = desc.split(" ")
                        if not "sp." in single_elem[2]:
                            if len(single_elem) >= 3:
                                species = " ".join([single_elem[1], single_elem[2]])
                                rec.description = species
                                fasta_all.append(rec)
                database_name = os.path.join(blastdb, args.database + "clean.fasta")
                print("Using database %s" % database_name)
                SeqIO.write(fasta_all, database_name, "fasta")
            elif database_name_orig.endswith(".nal") or database_name_orig.endswith(".Nal"):
                interm = database_name_orig.split(".")[:-1]
                database_name = ".".join(interm)
                print("Using database %s" % database_name)
            else:
                print("DATABASE NOT FOUND")
        else:
            print("DATABASE NOT FOUND")
    else:
        database_name = os.path.join(blastdb, "nt")
        print("Using database %s" % database_name)
    if "," in args.barcode:
        barcode_no_space = args.barcode.replace(" ","")
        list_barcode = barcode_no_space.split(",")
    elif "-" in args.barcode and args.barcode.startswith("barcode"):
        list_barcode = []
        barcode_name = args.barcode.replace("barcode","")
        b_min, b_max = barcode_name.split("-")
        list_number = range(int(b_min), int(b_max) + 1)
        for elem in list_number:
            number = (f"{elem:02d}")
            list_barcode.append("barcode" + str(number))

    else:
        list_barcode = [args.barcode]
    dict_species_all = {}
    dict_table_all = {}
    number_species = []
    files_kt_list = []
    mlstCount = {}
    for barcode in list_barcode:
        dir_fastq = os.path.join(cwd, barcode)
        print("EXECUTING BARCODE: %s" % barcode)
        if os.path.exists(dir_fastq):
            if os.path.isfile(dir_fastq):
                read_mapping = dir_fastq
                if verbose:
                    print(read_mapping)
            else:
                for root, dirs, files in os.walk(dir_fastq, topdown=False):
                    reads_fastq = tempfile.NamedTemporaryFile(dir=dirpathTEMP, suffix=".fastq", delete=False)
                    with open(reads_fastq.name, "w") as outfile:
                        for filename in files:
                            if filename.endswith("gz"):
                                filename = os.path.join(cwd, barcode, filename)
                                with gzip.open(filename, 'rb') as infile:
                                    for line in infile:
                                        outfile.write(line.decode())

                read_mapping = reads_fastq.name
            reads_porechop = tempfile.NamedTemporaryFile(dir=dirpathTEMP, suffix=".fastq", delete=False)
            min, max = args.interval.split("-")
            count = 0
            with open(reads_porechop.name, "w") as fh:
                for record in SeqIO.parse(read_mapping, "fastq"):
                    if int(min) < len(record.seq) < int(max):
                        if args.subset > 0:
                            if count < args.subset:
                                count += 1
                                SeqIO.write(record, fh, "fastq")
                        else:
                            SeqIO.write(record, fh, "fastq")
            print("RUNNING PORECHOP")
            reads_porechop_1 = tempfile.NamedTemporaryFile(dir=dirpathTEMP, suffix=".fastq", delete=False)
            porechop_cmd = PORECHOP % (reads_porechop.name, args.threads, reads_porechop_1.name)
            porechop = sb.Popen(porechop_cmd, shell=True, cwd=cwd, stderr=sb.PIPE, stdout=sb.PIPE)
            porechop.communicate()
            result_blast = False
            if args.assemble or args.use_assembled or args.correct:
                if args.folder_fast5 == "" and args.use_assembled:
                    sys.exit("FAST5 NOT PASSED WHILE ASKED FOR ASSEMBLED READS FOR ALIGNEMENT")
                else:
                    if args.correct:
                        fastx_all = []
                        barcode_corr = racon(reads_porechop_1.name, args.threads, cwd, args.num_reads_cor, args.mafft)
                        if verbose:
                            print("USING " + barcode_corr + " FOR BLAST")
                        for record in SeqIO.parse(barcode_corr, "fasta"):
                            if record.seq != "":
                                fp = tempfile.NamedTemporaryFile(dir=dirpathTEMP, suffix=".fasta", mode="w", delete=False)
                                fo = tempfile.NamedTemporaryFile(dir=dirpathTEMP, mode="w", prefix=barcode,
                                                                 suffix=".blastn", delete=False)
                                SeqIO.write(record, fp, "fasta")
                                fastx_all.append([fp.name, database_name, fo.name])
                    else:
                        fastx_all = []
                        if verbose:
                            print("USING " + reads_porechop_1.name + " FOR BLAST")
                        for record in SeqIO.parse(reads_porechop_1.name, "fastq"):
                            if record.seq != "":
                                fp = tempfile.NamedTemporaryFile(dir=dirpathTEMP, suffix=".fasta", mode="w", delete=False)
                                fo = tempfile.NamedTemporaryFile(dir=dirpathTEMP, mode="w", prefix=barcode,
                                                                 suffix=".blastn", delete=False)
                                SeqIO.write(record, fp, "fasta")
                                fastx_all.append([fp.name, database_name, fo.name])
                                if not barcode.endswith("raw"):
                                    if barcode.endswith("/"):
                                        barcode = barcode.split("/")[-2] + "raw"
                                    else:
                                        barcode = barcode.split("/")[-1] + "raw"
                    if args.use_assembled or args.assemble:
                        if not args.correct:
                            barcode_corr = racon(reads_porechop_1.name, args.threads, cwd, args.num_reads_cor, args.mafft)
                        assembled = assembly(barcode_corr, cwd, args.threads, read_mapping, od, barcode,
                                                 args.folder_fast5)
                        if verbose:
                            print("USING " + assembled + " FOR BLAST")
                        if args.use_assembled:
                            fastx_all = []
                            print(assembled)
                            if assembled != "":
                                for record in SeqIO.parse(assembled, "fasta"):
                                    if record.seq != "":
                                        fp = tempfile.NamedTemporaryFile(dir=dirpathTEMP, suffix=".fasta", mode="w", delete=False)
                                        fo = tempfile.NamedTemporaryFile(dir=dirpathTEMP, mode="w", prefix=barcode, suffix=".blastn", delete=False)
                                        SeqIO.write(record, fp, "fasta")
                                        fastx_all.append([fp.name, database_name, fo.name])
                            else:
                                continue
            else:
                if verbose:
                    print("USING " + reads_porechop_1.name + " FOR BLAST")
                for record in SeqIO.parse(reads_porechop_1.name, "fastq"):
                    if record.seq != "":
                        fp = tempfile.NamedTemporaryFile(dir=dirpathTEMP, suffix=".fasta", mode="w", delete=False)
                        fo = tempfile.NamedTemporaryFile(dir=dirpathTEMP, mode="w", prefix=barcode, suffix=".blastn", delete=False)
                        SeqIO.write(record, fp, "fasta")
                        fastx_all.append([fp.name, database_name,fo.name])
                        if not barcode.endswith("raw"):
                            if barcode.endswith("/"):
                                barcode = barcode.split("/")[-2]+ "raw"
                            else:
                                barcode = barcode.split("/")[-1]+ "raw"
            # if arg
            if len(fastx_all) < int(args.min_reads) and not args.use_assembled:
                print("Not analysing barcode " + barcode +". Less than " + str(args.min_reads) + " reads")
                continue
            result_list = []
            result_list_mlst = []
            correct_racon_blast = []
            if args.species and os.path.exists(fastamlst):
                hystogram = []
                for record in SeqIO.parse(outputRacon, "fasta"):
                    if record.seq != "":
                        fp = tempfile.NamedTemporaryFile(dir=dirpathTEMP, suffix=".fasta", mode="w", delete=False)
                        fo = tempfile.NamedTemporaryFile(dir=dirpathTEMP, mode="w", prefix=barcode, suffix=".blastn",
                                                         delete=False)
                        SeqIO.write(record, fp, "fasta")
                        correct_racon_blast.append([fp.name, fastamlst,fo.name])

                with Pool(processes=args.threads) as pool:
                    for result in tqdm(pool.imap(func=blast, iterable=correct_racon_blast), total=len(correct_racon_blast)):
                        result_list_mlst.append(result)
                for result in result_list_mlst:
                    with open(result, "r") as fh:
                        first_line = fh.readline()
                        second_line = fh.readline()
                    if len(first_line.split("\t")) > 3 and len(second_line.split("\t")) > 3 and first_line.split("\t")[2] > second_line.split("\t")[2]:
                        hystogram.append(first_line.split("\t")[1])
                a = dict(Counter(hystogram))
                mlstCount[barcode] = a
            print("RUNNING BLAST USING "  + str(args.threads) + " THREADS")
            with Pool(processes=args.threads) as pool:
                for result in tqdm(pool.imap(func=blast, iterable=fastx_all), total=len(fastx_all)):
                    result_list.append(result)
            result_list_tqdm = []
            name_blast = os.path.join(od, args.output + barcode + ".blast.txt")
            with open(name_blast, "w") as new_file:
                for name in result_list:
                    data_single_file = []
                    with open(name) as f:
                        for line in f:
                            new_file.write(line)
                            data_single_file.append(line)
                        new_file.write("\n")
                    blast_out_signle = "\n".join(data_single_file)
                    result_list_tqdm.append(blast_out_signle)
                    if not args.verbose:
                        os.remove(name)
            for line in result_list_tqdm:
                if line != "" and not result_blast:
                    result_blast = True
            if not result_blast:
                continue
            read_best_hit = {}
            for output in result_list_tqdm:
                if output != "":
                    align_species = {}
                    align = output.split("\n")
                    for single_align in align:
                        align_elm = single_align.split("\t")
                        if len(single_align) > 3:
                            align_species_score = {}
                            read = align_elm[0]
                            score = float(align_elm[2])
                            align_species_score[align_elm[4]] =  [align_elm[3]]
                            if align_elm[3] != "N/A":
                                if score > 200 and float(align_elm[4]) >= args.percentage :
                                    if score in align_species:
                                        for key in  align_species[score]:
                                            if align_elm[4] in align_species_score:
                                                align_species[score][key].append(align_elm[3])
                                            else:
                                                align_species[score][key] = [align_elm[3]]
                                    else:
                                        align_species[score] = align_species_score
                    if align_species:
                        list_align_species = list(align_species.items())
                        list_align_species.sort(reverse=True)
                        align_species_ident = list_align_species[0][1]
                        list_align_species_b = list(align_species_ident.items())
                        list_align_species_b.sort(reverse=True)
                        align_species_ident_b = list_align_species_b[0][1]
                        read_best_hit[read] = align_species_ident_b
            dict_match = []
            for match in read_best_hit:
                dict_match.append([match, read_best_hit[match]])
            result_list_tqdm = []
            with Pool(processes=args.threads) as pool:
                for result in tqdm(pool.imap(func=counpute_count, iterable=dict_match), total=len(dict_match)):
                    result_list_tqdm.append(result)


            kt_barcode = tempfile.NamedTemporaryFile(dir=dirpathTEMP,suffix=".txt", prefix=barcode, delete=False, mode = "w")
            files_kt_list.append(kt_barcode.name)
            for value in result_list_tqdm:
                for key in value[0]:
                    read = value[0][key][0]
                    taxid = ncbi.get_name_translator([key])
                    if bool(taxid):
                        try:
                            species_count = str(taxid[key][0]).split(";")[0]
                            kt_barcode.write(read + "\t" + species_count + "\n")
                        except:
                            print("NO TAXID FOR " + taxid[key][0])
                            continue
            species_dict = {}
            genera_dict = {}
            for result in result_list_tqdm:
                sp = list(result[0].keys())
                gen = list(result[1].keys())
                if len(sp) > 0:
                    if sp[0] in species_dict:
                        species_dict[sp[0]] = species_dict[sp[0]]+ 1
                    else:
                        species_dict[sp[0]] = 1
                if len(gen) > 0:
                    if gen[0] in genera_dict:
                        genera_dict[gen[0]] = [genera_dict[gen[0]][0] + 1]
                    else:
                        genera_dict[gen[0]] = [1]
            total_reads_mapped = sum(species_dict.values())
            min_value = total_reads_mapped * (float(args.min)/100)
            if min_value < 1 :
                if not args.assemble:
                    min_value = 1
            print("minimum number of reads used " + str(min_value))
            species_retain = {}
            specied_final = {}
            for key in species_dict:
                current_value = species_dict[key]
                if current_value >= min_value:#float(args.min):
                    species_retain[key] = species_dict[key]
            for key in species_retain:
                specied_final[key] = species_retain[key] / sum(species_retain.values()) * 100
                number_species.append(key)
            dict_species_all[barcode] = specied_final
            dict_table_all[barcode] = dict(sorted(species_dict.items(), key=lambda item : item[1]))

        else:
            print("PATH FOR BARCODE " + barcode + " DO NOT EXISTS. CHECK BARCODE NAME OR PATH" )
    if len(files_kt_list) > 0:
        plothystogramSpecies(files_kt_list,blastdb,args.output,od,dict_species_all,name_table,dict_table_all,figure_file)
    else:
        print("NOTHING TO PLOT")


if __name__ == '__main__':
    main()
