#! /usr/bin/python3.10
import sys
import uuid
import gzip
import random
import string
from matplotlib import pyplot as plt
import matplotlib
import subprocess as sb
import re
import os
from Bio import Entrez
from Bio import SeqIO
from Bio.Blast.Applications import NcbiblastnCommandline
from multiprocessing import Pool
from difflib import SequenceMatcher
from tqdm import tqdm
import tempfile
import argparse
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap
import csv
from ete3 import NCBITaxa


matplotlib.style.use('ggplot')
KTIMPORTTAX = "ktImportTaxonomy -tax %s -o %s %s"
MINIMAP = "minimap2 -a -x ava-ont -t %s %s %s "
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

def get_desired_ranks(taxid, desired_ranks):
    lineage = ncbi.get_lineage(taxid)
    lineage2ranks = ncbi.get_rank(lineage)
    ranks2lineage = dict((rank, taxid) for (taxid, rank) in lineage2ranks.items())
    return {'{}_id'.format(rank): ranks2lineage.get(rank, '<not present>') for rank in desired_ranks}

def main(taxids, desired_ranks, path):
    with open(path, 'w') as csvfile:
        fieldnames = ['{}_id'.format(rank) for rank in desired_ranks]
        writer = csv.DictWriter(csvfile, delimiter='\t', fieldnames=fieldnames)
        writer.writeheader()
        for taxid in taxids:
            writer.writerow(get_desired_ranks(taxid, desired_ranks))

#@Gooey(optional_cols=2, program_name="Subparser Layout Demo")
def setting():
    parser = argparse.ArgumentParser(prog='quantify and detect pathogens', usage='%(prog)s [options]')
    parser.add_argument("-o","--output", nargs="?", default="output", help="output name")
    parser.add_argument("--barcode", help="barcode to analyse", required=True)
    parser.add_argument("-f5","--folder_fast5", help="folder containing fast5 files",default="/data/sharedata/")
    parser.add_argument("-t", "--threads", nargs="?", type = int, default="10", help="number of thres")
    parser.add_argument("-d", "--database", nargs="?", choices=['Plant_Bacteria_Funghi.nal', 'nt.nal', 'ITS_16S_18S_28S_LSU_SSU.nal',
                                                                'Metazoa_Organism_OR_metazoa_All_Fields_.nal','Xanthomonas_genomes.nal','curatedXylellaDatabase.nal'],
                        default="Plant_Bacteria_Funghi.nal",
                        help="database name; can be a fasta file or a subset of ncbi database; Default is Plant_Bacteria_Funghi")
    parser.add_argument("-w", "--workdir", nargs="?", default="/tmp/")
    parser.add_argument("--folder_fastq", help='folder where fastq file are located', required=True, default="/data/sharedata/")#, widget='DirChooser')
    parser.add_argument("-e", "--email", nargs="?", default="", help="email for blast search", required=True)
    parser.add_argument("-s", "--search", nargs="?", default="", help="ncbi search to retrieve GIs")
    parser.add_argument("-m", "--min",  default="10", help="minimum number of reads to plot a species as fraction of total mappd reads [0-100]")
    parser.add_argument("-p", "--percentage", type = float, default="90", help="minimum identity to consider a valid alignment")
    parser.add_argument("-mr", "--min_reads", type = float, default="100", help="minimum number of reads to sequence to consider a sample in the analysis")
    parser.add_argument("-a", "--assemble", action='store_true', help="assembled-reads", default=False)
    parser.add_argument("-ua", "--use_assembled", action='store_true', help="use assembled reads for discovery", default=False)
    parser.add_argument("-dh", "--database_history", action='store_true', help="", default=False)
    parser.add_argument("-c", "--correct", action='store_true', help="correct or not reads via racon", default=False)
    parser.add_argument("-u", "--update", action='store_true', help="update database", default=False)#nargs="?", default="", required=True)
    parser.add_argument('--verbose', help='be verbose', dest='verbose', action='store_true', default=False)
    args = parser.parse_args()
    return args


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
    blastn_cline = NcbiblastnCommandline(db=elm[1], query=elm[0], evalue=0.001, out=elm[2] ,outfmt = "6 qseqid sseqid bitscore sscinames pident evalue staxids qlen" )
    try:
        blastn_cline()
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
    Entrez.email = email
    handle = Entrez.esearch(db="nucleotide", term=search, retmax=retmax, retstart=retstart)
    record = Entrez.read(handle)
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

def racon(od, fastq, threads, fast5, cwd):
    for root, dirs, files in os.walk(os.path.join(cwd,fastq), topdown=False):
        reads_fastq = tempfile.NamedTemporaryFile(suffix=".fastq", delete=False)
        with open(reads_fastq.name, "w") as outfile:
            for filename in files:
                filename = os.path.join(cwd, fastq, filename)
                with gzip.open(filename, 'rb') as infile:
                    for line in infile:
                        outfile.write(line.decode())

    read_mapping=reads_fastq.name
    print("RUNNING PORECHOP")
    reads_porechop_1 = tempfile.NamedTemporaryFile(suffix=".fastq", delete=False)
    reads_porechop = "/tmp/" + uuid.uuid4().hex + ".fastq"
    porechop_cmd = PORECHOP % (read_mapping, threads, reads_porechop_1.name)
    porechop = sb.Popen(porechop_cmd, shell=True, cwd=cwd, stderr=sb.PIPE, stdout=sb.PIPE)
    porechop.communicate()
    min = 1
    max = 200000
    with open(reads_porechop, "w") as fh:
        for record in SeqIO.parse(reads_porechop_1.name, "fastq"):
            if min < len(record.seq) < max:
                SeqIO.write(record, fh, "fastq")
    sam = tempfile.NamedTemporaryFile(suffix=".sam", delete=False)
    reads = tempfile.NamedTemporaryFile(suffix=".fasta", delete=False)
    print("RUNNING MINIMAP")
    m = MINIMAP % (threads, reads_porechop, reads_porechop)
    print(m)
    minimap = sb.Popen(m, shell=True, cwd=cwd, stdout=sam, stderr=sb.PIPE)
    minimap.communicate()
    print("RUNNING RACON")
    r = RACON % (threads, reads_porechop, sam.name, reads_porechop)
    print(r)
    racon_cmd = sb.Popen(r, shell=True, cwd=cwd, stdout=reads, stderr=sb.PIPE)
    racon_cmd.communicate()
    sam = tempfile.NamedTemporaryFile(suffix=".sam")
    print("RUNNING MINIMAP")
    m = MINIMAP % (threads, reads.name, reads.name)
    print(m)
    minimap = sb.Popen(m, shell=True, cwd=cwd, stdout=sam, stderr=sb.PIPE)
    minimap.communicate()
    print("RUNNING RACON")
    r = RACON % (threads, reads.name, sam.name, reads.name)
    print(r)
    output = tempfile.NamedTemporaryFile(suffix=".fasta", delete=False)
    jfc_out = tempfile.NamedTemporaryFile(suffix=".fasta", delete=False)
    racon_cmd = sb.Popen(r, shell=True, cwd=cwd, stdout=output, stderr=sb.PIPE)
    racon_cmd.communicate()
    print("RUNNING JELLYFISH")
    jfc = JELLYFISH_COUNT % output.name
    jfd = JELLYFISH_DUMP % jfc_out.name
    print(jfc)
    print(jfd)
    jellyfishc = sb.Popen(jfc, shell=True, cwd=cwd, stdout=sb.PIPE)
    jellyfishd = sb.Popen(jfd, shell=True, cwd=cwd, stdin=jellyfishc.stdout)
    jellyfishd.communicate()
    count = 0
    kmer = "/tmp/" + uuid.uuid4().hex + ".fasta"
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
    tmp_dir = tempfile.mkdtemp(dir="/tmp")
    print("RUNNING SPADES")
    spa = SPADES % (kmer, tmp_dir)
    print(spa)
    try:
        spade = sb.Popen(spa, shell=True, cwd=cwd, stderr=sb.PIPE, stdout=sb.PIPE)
        spade.communicate()
        assembled_contigs = os.path.join(tmp_dir, "contigs.fasta")
        bam = tempfile.NamedTemporaryFile(suffix=".bam", delete=False)
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
        ni = NANOPOLISHI % (os.path.join(fast5,fastq), reads_porechop)
        print(ni)
        nanopolish_index = sb.Popen(ni, shell=True, cwd=cwd, stdout=sb.PIPE, stderr=sb.PIPE)
        nanopolish_index.communicate()
        regions = []
        print("RUNNING NANOPOLISH")
        dict_vcf_fasta = []
        for record in SeqIO.parse(assembled_contigs, "fasta"):
            region = record.id + ":0-" + str(len(record.seq))
            vcf = tempfile.NamedTemporaryFile(prefix="polished.", suffix=".vcf", delete=False)
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
        bam = tempfile.NamedTemporaryFile(suffix=".bam", delete=False)
        bi = BWAI % mlst
        bwa_index = sb.Popen(bi, shell=True, cwd="/tmp")
        bwa_index.communicate()
        bm = BWA % (threads, mlst, output.name)
        #shutil.copyfile(mlst, mlst_done) #os.path.join(cwd,output + ".fasta"))
        ss = SAMSORT % bam.name
        print(bm)
        bwa_mem = sb.Popen(bm, shell=True, cwd="/tmp", stdout=sb.PIPE, stderr=sb.PIPE)
        samtools_sort = sb.Popen(ss, shell=True, cwd=cwd, stdin=bwa_mem.stdout, stdout=sb.PIPE, stderr=sb.PIPE)
        samtools_sort.communicate()
        si = SAMINDEX % bam.name
        samtoos_index = sb.Popen(si, shell=True, cwd=cwd, stdout=sb.PIPE, stderr=sb.PIPE)
        samtoos_index.communicate()
        vcf = tempfile.NamedTemporaryFile(suffix=".vcf", delete=False)
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
        return(mlst_done, "A")
    except:
        return(output.name, "B")


def analysis():
    print("PROGRAM STARTED")
    args = setting()
    mem_threads = 0
    cwd = args.folder_fastq
    od = os.getcwd()
    Entrez.email = args.email
    new_string = re.sub("[^0-9a-zA-Z]+", "_", args.search)
    name_plot = args.output + ".pdf"
    name_table = args.output + ".xlsx"
    figure_file = os.path.join(od, name_plot)
    name_table = os.path.join(od, name_table)
    os.environ["BLASTDB"] = "/data2/blastdb"
    blastdb = os.environ['BLASTDB']
    if args.search != "":
        database_name = os.path.join(blastdb, new_string)
        print("Using database %s" % database_name)
        database_name_check = os.path.join(blastdb, new_string + ".nal")
        if not os.path.exists(database_name_check) or args.update:
            handle_n = Entrez.esearch(db="nucleotide", term=args.search)
            record_n = Entrez.read(handle_n)
            retmax = 10000000
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
            fp = tempfile.NamedTemporaryFile(dir="/tmp", delete= False)
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
        list_barcode = args.barcode.split(",")
    else:
        list_barcode = [args.barcode]
    dict_species_all = {}
    dict_table_all = {}
    number_species = []
    files_kt_list = []
    for barcode in list_barcode:
        dir_fastq = os.path.join(cwd, barcode)
        if os.path.exists(dir_fastq):
            result_blast = False
            print("EXECUTING BARCODE: %s" % barcode)
            fastx_all = []
            fastx_all_assembled = []
            if args.assemble or args.use_assembled:
                if args.folder_fast5 == "":
                    print("ASSEMBLY NOT DONE, PLEASE ADD FAST5 PARAMETER - RUNNING ALIGNMENT")
                    if args.use_assembled:
                        sys.exit("FAST5 NOT PASSED WHILE ASKED FOR ASSEMBLED READS FOR ALIGNEMENT")
                else:
                    barcode_corr, spade_t = racon(od, barcode, args.threads, args.folder_fast5, cwd)
                    for record in SeqIO.parse(barcode_corr, "fasta"):
                        fp = tempfile.NamedTemporaryFile(dir="/tmp", suffix=".fasta", mode="w", delete=False)
                        fo = tempfile.NamedTemporaryFile(dir="/tmp", mode="w", prefix=barcode, suffix=".blastn", delete=False)
                        SeqIO.write(record, fp, "fasta")
                        fastx_all_assembled.append([fp.name, database_name, fo.name])
                        # else:
                        #     print(barcode + " BARCODE NOT DONE WITH CONSENSUS SEQUENCES")
                        #     continue
            if args.use_assembled:
                if spade_t == "A":
                    print("USING ASSEMBLED READS FOR ALIGNEMENT TO DATABASE FOR BARCODE " + barcode)
                else:
                    print("USING NON-ASSEMBLED READS FOR ALIGNEMENT TO DATABASE FOR BARCODE " + barcode)
                fastx_all = fastx_all_assembled
            else:
                print("RUNNING PORECHOP")
                reads_porechop = tempfile.NamedTemporaryFile(suffix=".fastq", delete=False)
                porechop_cmd = PORECHOP % (barcode, args.threads, reads_porechop.name)
                porechop = sb.Popen(porechop_cmd, shell=True, cwd=cwd, stderr=sb.PIPE, stdout=sb.PIPE)
                porechop.communicate()
                if args.correct:
                    sam = tempfile.NamedTemporaryFile(suffix=".sam", delete=False)
                    reads = tempfile.NamedTemporaryFile(suffix=".fasta", delete=False)
                    print("RUNNING MINIMAP")
                    m = MINIMAP % (args.threads, reads_porechop.name, reads_porechop.name)
                    print(m)
                    minimap = sb.Popen(m, shell=True, cwd=cwd, stdout=sam, stderr=sb.PIPE)
                    minimap.communicate()
                    print("RUNNING RACON")
                    r = RACON % (args.threads, reads_porechop.name, sam.name, reads_porechop.name)
                    print(r)
                    racon_cmd = sb.Popen(r, shell=True, cwd=cwd, stdout=reads, stderr=sb.PIPE)
                    racon_cmd.communicate()
                    if int(os.path.getsize(reads.name)) > 0:
                        sam = tempfile.NamedTemporaryFile(suffix=".sam")
                        print("RUNNING MINIMAP")
                        m = MINIMAP % (args.threads, reads.name, reads.name)
                        minimap = sb.Popen(m, shell=True, cwd=cwd, stdout=sam, stderr=sb.PIPE)
                        minimap.communicate()
                        print("RUNNING RACON")
                        r = RACON % (args.threads, reads.name, sam.name, reads.name)
                        output = tempfile.NamedTemporaryFile(suffix=".fasta", delete=False)
                        racon_cmd = sb.Popen(r, shell=True, cwd=cwd, stdout=output, stderr=sb.PIPE)
                        racon_cmd.communicate()
                        count = 0
                        for record in SeqIO.parse(output.name, "fasta"):
                            count += 1
                            fp = tempfile.NamedTemporaryFile(dir="/tmp", suffix=".fasta", mode="w", delete=False)
                            fo = tempfile.NamedTemporaryFile(dir="/tmp", mode="w", prefix=barcode, suffix=".blastn", delete=False)
                            SeqIO.write(record, fp, "fasta")
                            fastx_all.append([fp.name, database_name,fo.name])

                    else:
                        count = 0
                        print(reads_porechop.name)
                        for record in SeqIO.parse(reads_porechop.name, "fastq"):
                            count += 1
                            fp = tempfile.NamedTemporaryFile(dir="/tmp", suffix=".fasta", mode="w", delete=False)
                            fo = tempfile.NamedTemporaryFile(dir="/tmp", mode="w", prefix=barcode, suffix=".blastn", delete=False)
                            SeqIO.write(record, fp, "fasta")
                            fastx_all.append([fp.name, database_name, fo.name])
                            if not barcode.endswith("raw"):
                                barcode = barcode + "raw"
                else:
                    print(reads_porechop.name)
                    for record in SeqIO.parse(reads_porechop.name, "fastq"):
                        fp = tempfile.NamedTemporaryFile(dir="/tmp", suffix=".fasta", mode="w", delete=False)
                        fo = tempfile.NamedTemporaryFile(dir="/tmp", mode="w", prefix=barcode, suffix=".blastn", delete=False)
                        SeqIO.write(record, fp, "fasta")
                        fastx_all.append([fp.name, database_name,fo.name])
                        if not barcode.endswith("raw"):
                            if barcode.endswith("/"):
                                barcode = barcode.split("/")[-2]+ "raw"
                            else:
                                barcode = barcode.split("/")[-1]+ "raw"
                if len(fastx_all) < int(args.min_reads):
                    print("Not analysing barcode " + barcode +". Less than 50 reads")
                    continue
            result_list = []
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

            kt_barcode = tempfile.NamedTemporaryFile(suffix=".txt", prefix=barcode, delete=False, mode = "w")
            files_kt_list.append(kt_barcode.name)
            for value in result_list_tqdm:
                for key in value[0]:
                    read = value[0][key][0]
                    taxid = ncbi.get_name_translator([key])
                    #print("before")
                    #print(taxid)
                    #print("after")
                    if bool(taxid):
                        try:
                            #print(str(taxid[key][0]))
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
    # KTIMPORTTAX
        else:
            print("PATH FOR BARCODE " + barcode + " DO NOT EXISTS. CHECK BARCODE NAME OR PATH" )
    if len(files_kt_list) > 0:
        files_kt = " ".join(files_kt_list)
        k = KTIMPORTTAX % (blastdb, args.output + ".html", files_kt)
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
        #plt.show()
        plt.savefig(figure_file)
    else:
        print("NOTHING TO PLOT")




if __name__ == '__main__':
    analysis()
