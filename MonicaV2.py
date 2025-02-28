import argparse
import sys
import time
import uuid
import fileinput
import string
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
import random
from Bio import SeqIO
import seaborn as sns
import matplotlib.pyplot as plt
from multiprocessing import Pool
from io import StringIO
import shutil
from tqdm import tqdm
import tempfile
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap
from simple_slurm import Slurm
from ete4 import NCBITaxa
import requests
import xml.etree.ElementTree as ET
import resource

resource.setrlimit(
    resource.RLIMIT_CORE,
    (resource.RLIM_INFINITY, resource.RLIM_INFINITY))

import warnings
warnings.filterwarnings("ignore")



matplotlib.style.use('ggplot')
KTIMPORTTAX = "ktImportTaxonomy -tax %s -o %s %s"
MINIMAP = "minimap2 -x ava-ont -t %s -N %s %s %s > %s"
MINIMAP_S = "minimap2 -a -x map-ont --secondary=no -t %s %s %s "
MINIMAP_SY = "minimap2 -a -x map-ont -t %s %s %s "
SAMSORT = "samtools sort -o %s -O BAM"
SAMINDEX = "samtools index %s"
RACON = "racon -t %s -u -f %s %s %s"
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
BLASTN = "blastn -db %s -query %s -max_target_seqs 10 -evalue 0.001 -out %s -outfmt \"6 qseqid sseqid bitscore sscinames pident evalue staxids qlen\""

mlst_database = "https://pubmlst.org/data/dbases.xml"

ncbi = NCBITaxa()

parser = argparse.ArgumentParser(prog='quantify and detect pathogens', usage='%(prog)s [options]')
parser.add_argument("-o","--output", nargs="?", default="output", help="output name")
parser.add_argument("-r", "--retmax", nargs="?", type = int, default="100000", help="number of id to retrieve per itaration")
parser.add_argument("--barcode", help="barcode to analyse", required=True)
parser.add_argument("-f5","--folder_fast5", help="folder containing fast5 files",default="/mnt/data3ZEUSPV/SequenceData/")
parser.add_argument("-cwd","--working_directory", help="working directory",default="")
parser.add_argument("-t", "--threads", nargs="?", type = int, default="10", help="number of threads")
parser.add_argument("-d", "--database", nargs="?", choices=['Plant_Bacteria_Funghi.nal', 'nt.nal', 'ITS_16S_18S_28S_LSU_SSU.nal',
                                                            'Metazoa_Organism_OR_metazoa_All_Fields_.nal','Xanthomonas_genomes.nal','curatedXylellaDatabase.nal',
                                                            'customITSreduced.nal','customITSreducedReduced.nal','BacteriaFungiVitis.nal',
                                                            '50_4000RefSeq_biomol_genomic_ITS_16S_18S_28S_LSU_SSU.nal',
                                                            'bacterialGenome2024_01_24.fna.nal','emu','16S_ITS.nal','uniteMonica'],
                    default="",
                    help="database name; can be a fasta file or a subset of ncbi database; Default is Plant_Bacteria_Funghi")
parser.add_argument("-dbe", "--database_external", nargs="?", default="")
parser.add_argument("-tmp", "--temp", nargs="?", default='/mnt/dataQNAP/tmp/')
parser.add_argument("--folder_fastq", help='folder where fastq file are located', required=True,
                    default="/mnt/data3ZEUSPV/SequenceData/")#, widget='DirChooser')
parser.add_argument("-e", "--email", nargs="?", default="", help="email for blast search", required=True)
parser.add_argument("-s", "--search", nargs="?", default="", help="ncbi search to retrieve GIs")
parser.add_argument("-sp", "--species", nargs="?", default="", help="PUBMLST species to use")
parser.add_argument("-m", "--min",  default="1", help="minimum number of reads to plot a species as "
                                                       "fraction of total mappd reads [0-100]", type=float)
parser.add_argument("-mr", "--min_reads", type = float, default="100", help="minimum number of reads to "
                                                                            "sequence to consider a sample in the "
                                                                            "analysis")
parser.add_argument("-int", "--interval", default="100-20000", help="minimum and max read length"
                                                                            " to use in the correction step; default [100-20000]")
parser.add_argument("-su", "--subset", default="0", type = int, help="number of reads to use for each barcode [0]")
parser.add_argument("-mrc", "--num_reads_cor", type=int, default="15",
                    help="minimum number of reads to use in overlapping for reads correction [0]. 0 = no filtering")
parser.add_argument("-a", "--assemble", action='store_true', help="assembled-reads", default=False)
parser.add_argument("-sl", "--slurm", action='store_true', help="activate slurm grid", default=False)
parser.add_argument("-ua", "--use_assembled", action='store_true', help="use assembled reads for "
                                                                        "discovery", default=False)
parser.add_argument("-dh", "--database_history", action='store_true', help="", default=False)
parser.add_argument("-c", "--correct", action='store_true', help="correct or not reads via racon",
                    default=False)
parser.add_argument("-mafft", "--mafft", action='store_true', help="correct or not reads via racon",
                    default=False)
parser.add_argument("-u", "--update", action='store_true', help="update database", default=False)
parser.add_argument("-u_e3", "--update_ete", action='store_true', help="update database ete3", default=False)#nargs="?", default="", required=True)#nargs="?", default="", required=True)
parser.add_argument("-v",'--verbose', help='be verbose', dest='verbose', action='store_true', default=False)
parser.add_argument("-cu", "--cleanup", action='store_true', help="clean up verbose previous run",
                    default=False)
parser.add_argument("-slt", "--slurm_threads", nargs="?", type = int, default="0", help="number of threads for slurm")


def parse_blast(name):
        species = ""
        unambiguous = True
        with open(name) as f:
            for line in f:
                qseqid, sseqid, bitscore, sscinames, pident, evalue, staxids, qlen = line.rstrip().split("\t")
                if staxids != "0" and not ";" in staxids:
                    if species == "" and unambiguous:
                        try:
                            lineage = ncbi.get_lineage(staxids)
                            species = [staxids, bitscore, qseqid, pident, qseqid, sseqid, lineage]
                            staxidsC, bitscoreC, qseqidC, pidentC, qseqidC, sseqidC, lineageC = species
                        except:
                            print("TAXID " + str(staxids) + " NOT FOUND")
                    else:
                        if unambiguous and bitscoreC < bitscore:
                            try:
                                lineage = ncbi.get_lineage(staxids)
                                species = [staxids, bitscore, bitscore, pident, qseqid, sseqid, lineage]
                                staxidsC, bitscoreC, qseqidC, pidentC, qseqidC, sseqidC, lineageC = species
                            except:
                                print("TAXID " + str(staxids) + " NOT FOUND")

                        elif unambiguous and bitscoreC == bitscore and staxidsC != staxids and len(lineage) > 1:
                            if unambiguous and pidentC < pident:
                                try:
                                    lineage = ncbi.get_lineage(staxids)
                                    species = [staxids, bitscore, bitscore, pident, qseqid, sseqid, lineage]
                                    staxidsC, bitscoreC, qseqidC, pidentC, qseqidC, sseqidC, lineageC = species
                                except:
                                    print("TAXID " + str(staxids) + " NOT FOUND")
                            else:
                                lineage_old = "_".join(list(map(str, lineageC)))
                                lineage_new = "_".join(list(map(str, lineage)))
                                both_lineage = {lineage_old, lineage_new}
                                common_taxid = common_start(*both_lineage)
                                if "_" in common_taxid and len(common_taxid.split("_")) >= 2:
                                    common_taxid_uniq = common_taxid.split("_")[-2]
                                    species = [common_taxid_uniq, bitscore, bitscore, pident, qseqid, sseqid,
                                               ncbi.get_lineage(common_taxid_uniq)]
                                else:
                                    species = []
                                    continue
                if staxids == "0" and sscinames == "N/A":
                    if ":" in sseqid:
                        staxids = sseqid.split(":")[0]
                        if species == "" and unambiguous:
                            try:
                                lineage = ncbi.get_lineage(staxids)
                                species = [staxids, bitscore, qseqid, pident, qseqid, sseqid, lineage]
                                staxidsC, bitscoreC, qseqidC, pidentC, qseqidC, sseqidC, lineageC = species
                            except:
                                print("TAXID " + str(staxids) + " NOT FOUND")
                                lineage = []
                        else:
                            if unambiguous and bitscoreC < bitscore:
                                try:
                                    lineage = ncbi.get_lineage(staxids)
                                    species = [staxids, bitscore, qseqid, pident, qseqid, sseqid, lineage]
                                    staxidsC, bitscoreC, qseqidC, pidentC, qseqidC, sseqidC, lineageC = species
                                except:
                                    print("TAXID " + str(staxids) + " NOT FOUND")
                                    lineage = []
                            elif unambiguous and bitscoreC == bitscore and staxidsC != staxids and len(lineage) > 1:
                                if unambiguous and pidentC < pident:
                                    try:
                                        lineage = ncbi.get_lineage(staxids)
                                        species = [staxids, bitscore, qseqid, pident, qseqid, sseqid, lineage]
                                        staxidsC, bitscoreC, qseqidC, pidentC, qseqidC, sseqidC, lineageC = species
                                    except:
                                        print("TAXID " + str(staxids) + " NOT FOUND")
                                        lineage = []
                                else:
                                    lineage_old = "_".join(list(map(str, lineageC)))
                                    lineage_new = "_".join(list(map(str, lineage)))
                                    both_lineage = {lineage_old, lineage_new}
                                    common_taxid = common_start(*both_lineage)
                                    if "_" in common_taxid and len(common_taxid.split("_")) >= 2:
                                        common_taxid_uniq = common_taxid.split("_")[-2]
                                        species = [common_taxid_uniq, bitscore, qseqid, pident, qseqid, sseqid, ncbi.get_lineage(common_taxid_uniq)]
                                    else:
                                        species = []
                                        continue


        if not verbose:
            os.remove(name)
        return (species)

def get_desired_ranks(taxid, desired_ranks):
    value_taxid, score = taxid
    try:
        lineage = ncbi.get_lineage(value_taxid)
        lineage2ranks = ncbi.get_rank(lineage)
        ranks2lineage = dict((rank, value_taxid) for (value_taxid, rank) in lineage2ranks.items())
        ranks2lineage_select = {'{}'.format(rank): ranks2lineage.get(rank, '<not present>') for rank in desired_ranks}
        for key in ranks2lineage_select:
            if ranks2lineage_select[key] != "<not present>":
                ranks2lineage_select[key] = list(ncbi.get_taxid_translator([ranks2lineage_select[key]]).values())[0]
        taxid_dict = {"tax_id": value_taxid,"abundance":score}
        taxid_dict.update(ranks2lineage_select)
        return(taxid_dict)
    except:
        print(str(taxid))
        return ()

def desired_rank_file(taxids, desired_ranks):
    taxid_dict = {}
    for taxid in tqdm(taxids):
        rank = get_desired_ranks(taxid, desired_ranks)
        taxid_dict[taxid[0]] = rank
    return(taxid_dict)

def best_align_match(samfile, iteration):
    if iteration > 0:
        samDict = {}
        count = 0
        filtered = samfile + ".reduced.paf"
        with open(samfile, "r") as sam, open(filtered, "w") as reduced:
            for line in sam:
                #line = line.decode()
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
    with open(file, "r") as sam:
        for line in sam:
            #line = line.decode()
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
        fasta = tempfile.NamedTemporaryFile(dir=dirpathTEMP, suffix=".fasta", delete=False, mode="w")
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
    aln = tempfile.NamedTemporaryFile(dir=dirpathTEMP,  suffix=".fasta", delete=False, mode="w")
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

def blast_slurm(elms):
    
    query, db, out, verbose, cleanup, jobs_id = elms
    if verbose:
        slurm = Slurm(
            ntasks=1,
            cpus_per_task=1,
            export=f'BLASTDB=/mnt/dataQNAP/blastdb/',
        )
    else:
        slurm = Slurm(
            ntasks=1,
            cpus_per_task=1,
            output=f'/dev/null',
            error=f'/dev/null',
            export=f'BLASTDB=/mnt/dataQNAP/blastdb/',
    )

    blastn_cline = BLASTN % (db, query, out)
    #print(blastn_cline)
    slurm.add_arguments(job_name=jobs_id)
    try:
        slurm.sbatch(blastn_cline)
    except:
        print("BLASTN failed " + blastn_cline)
    if verbose:
        print(slurm)

def wait_for_file_to_stabilize(filepath):

    check_interval = 2
    wait_time = 10
    previous_size = -1
    stable_time = 0

    while True:
        # Get the current file size
        current_size = os.path.getsize(filepath)

        # If file size hasn't changed, increment the stable time
        if current_size == previous_size:
            stable_time += check_interval
        else:
            stable_time = 0  # Reset if size changes

        # If the file size is stable for the desired wait time, exit the loop
        if stable_time >= wait_time:
            #print(f"File has stabilized after {stable_time} seconds.")
            break

        # Update the previous size and wait before checking again
        previous_size = current_size
        #print(f"File size: {current_size} bytes, waiting for {check_interval} seconds...")
        time.sleep(check_interval)

def extract_aligned_reads(paf_file, fasta_file, num_reads_cor):
    aligned_read_ids = {}

    # Read the PAF file and extract the read IDs of aligned reads
    with open(paf_file, 'r') as paf:
        for line in paf:
            fields = line.split('\t')
            if fields[0] in aligned_read_ids:
                if len(aligned_read_ids[fields[0]]) < int(num_reads_cor):
                    aligned_read_ids[fields[0]].append(fields[5])
            else:
                aligned_read_ids[fields[0]] = [fields[5]]

    record_dict = SeqIO.to_dict(SeqIO.parse(fasta_file, "fastq"))

    file_names_list = []
    for sequence in aligned_read_ids:
        collection_sequence = []
        for sequence_id in aligned_read_ids[sequence]:
            collection_sequence.append(record_dict[sequence_id])
        file_name = os.path.join(dirpathTEMP, sequence + "out_seq")
        file_names_list.append(file_name)
        with open(file_name, "w") as output_handle:
            SeqIO.write(collection_sequence, output_handle, "fasta")

    return (file_names_list)

def run_cap3(fasta_file):
    keep_fasta = []
    # for fasta in fasta_file:
    cap3 = sb.Popen(['/home/lfaino/Downloads/CAP3/cap3', fasta_file], stdout=sb.PIPE, stderr=sb.PIPE, cwd = dirpathTEMP)
    out, err = cap3.communicate()
    file_contig = fasta_file +".cap.contigs"
    if os.path.exists(file_contig) and os.path.getsize(file_contig) > 10:
        for record in SeqIO.parse(file_contig, 'fasta'):
            fasta = fasta_file.split("/")[-1].split(".")[0].replace("out_seq", "")
            record.id = fasta
            keep_fasta.append(record)

    if len(keep_fasta) == 1:
        return keep_fasta[0]

def slurm_others(command, slurm_threads,verbose,tmp_slurm, email):
    job_id = ''.join(random.choices('0123456789abcdef', k=8))
    software = command.split(" ")[0]
    while True:
        job_master_error = tmp_slurm + f'/{job_id}_{software}.err'
        job_master_out = tmp_slurm + f'/{job_id}_{software}.out'
        slurm = Slurm(
            cpus_per_task=slurm_threads,
            quiet = True
        )
        slurm.add_arguments(job_name=job_id)
        slurm.add_arguments(mail_user=email)
        slurm.add_arguments(output=job_master_out)
        slurm.add_arguments(error=job_master_error)


        slurm.sbatch(command)
        if verbose:
            print(slurm)
        print_out = True
        while True:
            result = sb.run(['squeue'], stdout=sb.PIPE, stderr=sb.PIPE, text=True)
            if print_out:
                print("WAITING FOR SLURM TO FINISH " + software)
                print_out = False
            if not job_id in result.stdout:
                break
            time.sleep(60)
        if "racon" in command:
            output_file = command.split(" ")[-1]
            wait_for_file_to_stabilize(output_file)
        elif "minimap2" in command:
            output_file = command.split(" ")[-1]
            wait_for_file_to_stabilize(output_file)
        else:
            print("DONE" + command)
        if os.path.exists(job_master_error):
            with open(job_master_error, "r") as f:
                lines = f.readlines()
                if lines:
                    last_line = lines[-1].strip()
                    if "Killed" in last_line:
                        slurm_threads -= 10
                        if "racon" in command:
                            command_list = command.split(" ")
                            command_list[2] = str(slurm_threads)
                            command = " ".join(command_list)
                            os.remove(job_master_out)
                            os.remove(job_master_error)
                        elif "minimap2" in command:
                            command_list = command.split(" ")
                            command_list[4] = str(slurm_threads)
                            command = " ".join(command_list)
                            os.remove(job_master_out)
                            os.remove(job_master_error)
                        else:
                            print("I CAN NOT SET NEW THREADS NUMBER")
                    else:
                        break

def blast(elm):

    query, db, out, verbose, cleanup, slurm = elm
    open(out, 'a').close()
    if os.path.isfile(out) and os.path.getsize(out) > 0 and verbose and not cleanup:
        return out
    else:
        #BLASTN = "blastn -db %s query %s -max_target_seqs 10 -evalue 0.001 -out %s -outfmt \"6 qseqid sseqid bitscore sscinames pident evalue staxids qlen\""
        blastn_cline = BLASTN % (db, query, out)
        try:
            blastn = sb.Popen(blastn_cline, shell=True, stderr=sb.PIPE)
            blastn.communicate()
            if not verbose:
                os.remove(query)
        except ValueError:
            print(blastn_cline)
        return(out)

def download_fasta(elm):
    search, retmax, retstart, email = elm
    Bio.Entrez.email = email
    handle = Bio.Entrez.esearch(db="nucleotide", term=search, retmax=retmax, retstart=retstart)
    record = Bio.Entrez.read(handle)
    records = str("\n".join(record["IdList"]))
    return (records)

def cap3(fastq, threads, cwd, num_reads_cor, mafft, verbose, slurm, slurm_threads, tmp_slurm, email):

    paf = tempfile.NamedTemporaryFile(dir=dirpathTEMP,suffix=".paf", delete=False)
    reads_out = tempfile.NamedTemporaryFile(dir=dirpathTEMP,suffix=".fasta", delete=False)
    print("RUNNING MINIMAP")

    if slurm:
        m = MINIMAP % (slurm_threads, num_reads_cor, fastq, fastq, paf.name)
        if num_reads_cor == 0:
            m = m.replace("-N 0", " ")
        slurm_others(m, slurm_threads, verbose, tmp_slurm,email)
    else:
        m = MINIMAP % (threads, num_reads_cor, fastq, fastq, paf.name)
        if num_reads_cor == 0:
            m = m.replace("-N 0", " ")
        minimap = sb.Popen(m, shell=True, cwd=cwd, stderr=sb.PIPE)
        minimap.communicate()
    if verbose:
        print(m)
    fasta_file = extract_aligned_reads(paf.name, fastq, num_reads_cor)
    result_cap3 = []
    print("RUNNING CAP3")
    with Pool(processes=threads) as pool:
        for result in tqdm(pool.imap(func=run_cap3, iterable=fasta_file), total=len(fasta_file)):
            result_cap3.append(result)
    print(reads_out.name)
    cleaned_arr = list(filter(None, result_cap3))
    with open(reads_out.name, 'w') as out_fasta:
        for record in cleaned_arr:
            SeqIO.write(record, out_fasta, 'fasta')
    return (reads_out.name)

def racon_assembly(fastq, threads, cwd, num_reads_cor, verbose):
    sam = tempfile.NamedTemporaryFile(dir=dirpathTEMP,suffix=".paf.gz", delete=False)
    reads = tempfile.NamedTemporaryFile(dir=dirpathTEMP,suffix=".fasta", delete=False)
    print("RUNNING MINIMAP")
    m = MINIMAP % (threads, num_reads_cor, fastq, fastq, sam.name)
    print(m)
    minimap = sb.Popen(m, shell=True, cwd=cwd, stderr=sb.PIPE)
    minimap.communicate()
    filtered = best_align_match(sam.name, num_reads_cor)
    print("RUNNING RACON")
    r = RACON % (threads, fastq, filtered , fastq)
    if verbose:
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
    m = MINIMAP % (threads, num_reads_cor, reads.name, fastq, sam.name)
    if verbose:
        print(m)
    minimap = sb.Popen(m, shell=True, cwd=cwd, stderr=sb.PIPE)
    minimap.communicate()
    filtered = best_align_match(sam.name, num_reads_cor)
    print("RUNNING RACON")
    r = RACON % (threads, fastq, filtered, reads.name)
    if verbose:
        print(r)
    output = tempfile.NamedTemporaryFile(dir=dirpathTEMP,suffix=".fasta", delete=False)
    racon_cmd = sb.Popen(r, shell=True, cwd=cwd, stdout=output, stderr=sb.PIPE)
    racon_cmd.communicate()
    if mafft:
        print("RUNNING MAFFT FOR READ CORRECTION AND ASSEMBLY")
        sam = tempfile.NamedTemporaryFile(dir=dirpathTEMP,suffix=".paf.gz", delete=False)
        print("RUNNING MINIMAP")
        m = MINIMAP % (threads, num_reads_cor, output.name, output.name, sam.name)
        if verbose:
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
    else:
        return (output.name)

def assembly(output, cwd, threads, read_mapping, od, fastq, fast5, verbose):

    print("STARTED TO RUN ASSEMBLY")
    jfc_out = tempfile.NamedTemporaryFile(dir=dirpathTEMP,suffix=".jfc.fasta", delete=False)
    print("RUNNING JELLYFISH")
    jfc = JELLYFISH_COUNT % output
    jfd = JELLYFISH_DUMP % jfc_out.name
    if verbose:
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
                SeqIO.write(record, fh, "fasta")
    if verbose:
        print(kmer)
    tmp_dir = tempfile.mkdtemp(dir=dirpathTEMP)
    print("RUNNING SPADES")
    spa = SPADES % (kmer, tmp_dir)
    if verbose:
        print(spa)
    try:
        spade = sb.Popen(spa, shell=True, cwd=cwd, stderr=sb.PIPE, stdout=sb.PIPE)
        spade.communicate()
        assembled_contigs = os.path.join(tmp_dir, "contigs.fasta")
        bam = tempfile.NamedTemporaryFile(dir=dirpathTEMP,suffix=".bam", delete=False)
        print("RUNNING MINIMAP")
        m = MINIMAP_S % (threads, assembled_contigs, read_mapping)
        ss = SAMSORT % bam.name
        if verbose:
            print(m)
        minimap = sb.Popen(m, shell=True, cwd=cwd, stdout=sb.PIPE, stderr=sb.PIPE)
        samtools_sort = sb.Popen(ss, shell=True, cwd=cwd, stdin=minimap.stdout, stdout=sb.PIPE, stderr=sb.PIPE)
        samtools_sort.communicate()
        si = SAMINDEX % bam.name
        samtoos_index = sb.Popen(si, shell=True, cwd=cwd, stdout=sb.PIPE, stderr=sb.PIPE)
        samtoos_index.communicate()
        ni = NANOPOLISHI % (fast5, read_mapping)
        if verbose:
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
            if verbose:
                print(nv)
            nanopolish_var = sb.Popen(nv, shell=True, cwd=cwd, stdout=sb.PIPE, stderr=sb.PIPE)
            nanopolish_var.communicate()
            if os.path.isfile(vcf.name) and os.path.getsize(vcf.name) > 0:
                regions.append(vcf.name)
            else:
                dict_vcf_fasta.append(record)
        nva = NANOPOLISHVA % (assembled_contigs, " ".join(regions))
        if verbose:
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
        if verbose:
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
        if verbose:
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
        if verbose:
            print(bgzip)
        tabix_run = sb.Popen(tabix, shell=True, cwd=cwd, stdout=sb.PIPE, stderr=sb.PIPE)
        tabix_run.communicate()
        if verbose:
            print(tabix)
            print(bcfi)
        bcftools_index = sb.Popen(bcfi, shell=True, cwd=cwd, stdout=sb.PIPE, stderr=sb.PIPE)
        bcftools_index.communicate()
        bcfco = BCFTOOLS_CO % (mlst, mlst_done, vcf_gz)
        if verbose:
            print(bcfco)
        #print(mlst_done)
        bcftools_call = sb.Popen(bcfco, shell=True, cwd=cwd, stdout=sb.PIPE, stderr=sb.PIPE)
        bcftools_call.communicate()
        print("FINISHED TO RUN ASSEMBLY")

        return(mlst_done)

    except:
        mlst_done = ""
        print("FINISHED TO RUN ASSEMBLY")
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
    makedb.communicate()
    return(mlst_complete)

def plothystogramSpecies(files_kt_list, blastdb, output, od, verbose):
    files_kt = " ".join(files_kt_list)
    k = KTIMPORTTAX % (blastdb, output + ".html", files_kt)
    print(k)
    ktimport = sb.Popen(k, shell=True, cwd=od, stderr = sb.PIPE, stdout=sb.PIPE)
    out, err = ktimport.communicate()
    if verbose:
        print(err)
        print(out)

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

def histograplot(dict_species_all, figure_file):
    species_pd = pd.DataFrame.from_dict(dict_species_all, orient='index')
    species_name_colors = [name for barcode in dict_species_all for name in dict_species_all[barcode]]
    species_name_colors = list(dict.fromkeys(species_name_colors))
    cmap = convert_color(species_name_colors)
    cmap2 = []
    for i in cmap:
        cmap2.append(i[1])
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

def heatmap_clustering(dict_species_all, figure_file):


    df = pd.DataFrame.from_dict(dict_species_all)
    # df1 = df.transpose()
    df.fillna(0, inplace=True)
    df_percentage = df.div(df.sum(axis=0), axis=1) * 100
    df_filtered = df_percentage[(df_percentage <= 1).any(axis=1)]

    plt.figure(figsize=(10, 8))
    sns.heatmap(df_filtered, cmap='viridis', annot=True, fmt=".1f")

    plt.title('Heatmap delle prime 30 specie più abbondanti nei campioni')
    plt.xlabel('Campioni')
    plt.ylabel('Microorganismi')

    cluster_map = sns.clustermap(df_filtered, cmap='viridis', annot=True, fmt=".1f", metric="euclidean", method="average",
                                 figsize=(12, 10))
    plt.title('Clustermap delle prime 30 specie più abbondanti nei campioni')
    plt.savefig(figure_file)

def main():
    print("PROGRAM STARTED")
    # args = setting()
    args = parser.parse_args()
    job_id = ''.join(random.choices('0123456789abcdef', k=8))
    global dirpathTEMP
    dirpathTEMP = tempfile.mkdtemp(dir=args.temp)
    global temp
    temp = args.temp
    if args.slurm:
        temp = "/mnt/dataQNAP/tmp"
    global verbose
    verbose = args.verbose
    if args.update_ete:
        ncbi.update_taxonomy_database()
    if args.species != "":
        fastamlst = mlstget(args.species)

    global cwd
    cwd = args.folder_fastq
    global od
    od = os.getcwd()

    Bio.Entrez.email = args.email
    new_string = re.sub("[^0-9a-zA-Z]+", "_", args.search)
    name_plot_hist = args.output + ".histogram.pdf"
    name_plot_heat= args.output + ".heatmap.pdf"
    name_table = args.output + ".xlsx"
    figure_file_hist = os.path.join(od, name_plot_hist)
    figure_file_heat = os.path.join(od, name_plot_heat)
    name_table = os.path.join(od, name_table)
    os.environ["BLASTDB"] = "/mnt/dataQNAP/blastdb/"
    blastdb = os.environ['BLASTDB']
    outputRacon = ""
    blastdbcmd_cmd = "blastdbcmd -db " + args.database_external + " -info"
    blastdbcmd = sb.Popen(blastdbcmd_cmd, shell=True, stderr=sb.PIPE, stdout=sb.PIPE)
    out, err = blastdbcmd.communicate()
    if err.decode() == "":
        database_name = args.database_external
    elif args.search != "":
        database_name = os.path.join(blastdb, new_string)
        print("USING DATABASE NAME %s" % database_name)
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
            #if repeat > 10:
            #    threads = 10
            #else:
            #    threads = repeat
            for subset in tqdm(id_download):
                result = download_fasta(subset)
                result_list_tqdm.append(result)
                time.sleep(10)
            records = "\n".join(result_list_tqdm)
            fp = tempfile.NamedTemporaryFile(dir=dirpathTEMP, delete= False)
            fp.write(records.encode())
            cm = "blastdb_aliastool -db nt -dbtype nucl -gilist %s -out %s" % (fp.name, database_name)
            print(cm)
            blastalias = sb.Popen(cm, shell= True)
            blastalias.communicate()
    elif args.database != "":
        #fasta_all = []
        if os.path.isfile(os.path.join(blastdb, args.database)):
            if args.database.endswith(".nal") or args.database.endswith(".Nal"):
                database_name = os.path.join(blastdb, args.database)
                interm = database_name.split(".")[:-1]
                database_name = ".".join(interm)
    elif args.database == "" and args.database_external != "":
        fasta_all = []
        if args.database_external.endswith("fasta"):
            for rec in SeqIO.parse(args.database, "fasta"):
                desc = rec.description
                if 200 < len(rec.seq) < 10000:
                    single_elem = desc.split(" ")
                    if not "sp." in single_elem[2]:
                        if len(single_elem) >= 3:
                            species = " ".join([single_elem[1], single_elem[2]])
                            rec.description = species
                            fasta_all.append(rec)
            database_name = os.path.join(blastdb, args.database + "clean.fasta")
            SeqIO.write(fasta_all, database_name, "fasta")
            print("Using database %s" % database_name)
        elif args.database_external.endswith(".nal") or args.database_external.endswith(".Nal"):
            database_name_abs = os.path.abspath(args.database_external)
            database_name_rel = os.path.join(blastdb, args.database_external)
            if os.path.isfile(database_name_abs):
                database_name = database_name_abs
                interm = database_name.split(".")[:-1]
                database_name = ".".join(interm)
                print("Using database %s" % database_name)
            elif os.path.isfile(database_name_rel):
                database_name = database_name_rel
                interm = database_name.split(".")[:-1]
                database_name = ".".join(interm)
                print("Using database %s" % database_name)
            else:
                interm = args.database_external.split(".")[:-1][0].split("/")[-1]
                blastdbcmd_cmd = "blastdbcmd -db " + interm + " -info"
                blastdbcmd = sb.Popen(blastdbcmd_cmd, shell=True, stderr=sb.PIPE, stdout=sb.PIPE)
                out, err = blastdbcmd.communicate()
                if err.decode() == "":
                    database_name = os.path.join(blastdb, interm)
                else:
                    sys.exit("DATABASE NOT FOUND")

        else:
            blastdbcmd_cmd = "blastdbcmd -db " + args.database_external + " -info"
            blastdbcmd = sb.Popen(blastdbcmd_cmd, shell=True, stderr=sb.PIPE, stdout=sb.PIPE)
            out, err = blastdbcmd.communicate()
            if err.decode() == "":
                database_name = os.path.join(blastdb, args.database_external)
            else:
                sys.exit("DATABASE NOT FOUND")
    else:
        database_name = os.path.join(blastdb, "nt")
        print("USING DATABASE NAME %s" % database_name)
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
    dict_table_all = {}
    files_kt_list = []
    mlstCount = {}
    min, max = args.interval.split("-")
    for barcode in list_barcode:
        fastx_all = []
        list_species_verbose = []
        rank_path = os.path.join(od, barcode + "_rank.txt")
        count = 0
        dir_fastq = os.path.join(cwd, barcode)
        print("EXECUTING BARCODE: %s" % barcode)
        if os.path.exists(dir_fastq):
            if os.path.isfile(dir_fastq):
                reads_fastq = tempfile.NamedTemporaryFile(dir=dirpathTEMP, suffix=".fastq", delete=False, mode="w")
                if dir_fastq.endswith("gz"):
                    with gzip.open(dir_fastq, 'rt') as f:
                        if args.subset > 0:
                            for record in SeqIO.parse(f, "fastq"):
                                if int(min) < len(record.seq) < int(max):
                                    if count < args.subset:
                                        count += 1
                                        SeqIO.write(record, reads_fastq, "fastq")
                        else:
                            for record in SeqIO.parse(f, "fastq"):
                                if int(min) < len(record.seq) < int(max):
                                    SeqIO.write(record, reads_fastq, "fastq")
                elif dir_fastq.endswith("fastq"):
                    if args.subset > 0:
                        for record in SeqIO.parse(dir_fastq, "fastq"):
                            if int(min) < len(record.seq) < int(max):
                                if count < args.subset:
                                    count += 1
                                    SeqIO.write(record, reads_fastq, "fastq")
                    else:
                        for record in SeqIO.parse(dir_fastq, "fastq"):
                            if int(min) < len(record.seq) < int(max):
                                SeqIO.write(record, reads_fastq, "fastq")
                else:
                    print("UNKNOWN BARCODE TYPE. ACCEPTED FILES OR FOLDERS")
                read_mapping = reads_fastq.name
                if verbose:
                    print(read_mapping)
            elif os.path.isdir(dir_fastq):
                for root, dirs, files in os.walk(dir_fastq, topdown=False):
                    reads_fastq = tempfile.NamedTemporaryFile(dir=dirpathTEMP, suffix=".fastq", delete=False)
                    with open(reads_fastq.name, "w") as outfile:
                        for filename in files:
                            if filename.endswith("gz"):
                                filename = os.path.join(cwd, barcode, filename)
                                with gzip.open(filename, 'rb') as infile:
                                    for line in infile:
                                        outfile.write(line.decode())
                            elif filename.endswith("fastq"):
                                filename = os.path.join(cwd, barcode, filename)
                                with open(filename, 'rb') as infile:
                                    for line in infile:
                                        outfile.write(line.decode())
                reads_fastq_select = tempfile.NamedTemporaryFile(dir=dirpathTEMP, suffix=".fastq", delete=False, mode="w")
                if args.subset > 0:
                    for record in SeqIO.parse(reads_fastq.name, "fastq"):
                        if int(min) < len(record.seq) < int(max):
                            if count < args.subset:
                                count += 1
                                SeqIO.write(record, reads_fastq_select, "fastq")
                else:
                    for record in SeqIO.parse(reads_fastq.name, "fastq"):
                        if int(min) < len(record.seq) < int(max):
                            SeqIO.write(record, reads_fastq_select, "fastq")

                read_mapping = reads_fastq_select.name
            else:
                print("UNKNOWN BARCODE TYPE. ACCEPTED FILES OR FOLDERS")
            count = 0
            print("RUNNING PORECHOP")
            reads_porechop_1 = tempfile.NamedTemporaryFile(dir=dirpathTEMP, suffix=".fastq", delete=False)
            porechop_cmd = PORECHOP % (read_mapping, args.threads, reads_porechop_1.name)
            if verbose:
                print(porechop_cmd)
            porechop = sb.Popen(porechop_cmd, shell=True, cwd=cwd, stderr=sb.PIPE, stdout=sb.PIPE)
            porechop.communicate()
            result_blast = False
            if args.assemble or args.use_assembled or args.correct:
                if args.folder_fast5 == "" and args.use_assembled:
                    sys.exit("FAST5 NOT PASSED WHILE ASKED FOR ASSEMBLED READS FOR ALIGNEMENT")
                else:
                    if args.correct:
                        barcode_corr = cap3(reads_porechop_1.name, args.threads, cwd, args.num_reads_cor, args.mafft, args.verbose,
                                             args.slurm, args.slurm_threads, temp, args.email)
                        if verbose:
                            print("USING " + barcode_corr + " FOR BLAST")
                        fasta_uncorrect = SeqIO.to_dict(SeqIO.parse(reads_porechop_1.name, "fastq"))
                        fasta_correct = SeqIO.to_dict(SeqIO.parse(barcode_corr, "fasta"))
                        print("WE CORRECTED " + str(len(fasta_correct)) + " OUT OF " + str(len(fasta_uncorrect)))
                        uncorrected_records = [key for key in fasta_uncorrect if not key in fasta_correct]
                        for id in uncorrected_records:
                            fasta_correct[id] = fasta_uncorrect[id]
                        for key in fasta_correct:
                            record = fasta_correct[key]
                            if record.seq != "":
                                fp = os.path.join(dirpathTEMP,
                                                  record.id + ".fasta")  # tempfile.NamedTemporaryFile(dir=dirpathTEMP, suffix=".fasta", mode="w", delete=False)
                                fo = os.path.join(dirpathTEMP,
                                                  record.id + ".blastn")
                                SeqIO.write(record, fp, "fasta")
                                fastx_all.append([fp, database_name, fo, args.verbose, args.cleanup,job_id])
                    if args.use_assembled or args.assemble:
                        # barcode_corr_ass = racon_assembly(reads_porechop_1.name, args.threads, cwd, args.num_reads_cor, verbose)
                        assembled = assembly(reads_porechop_1.name, cwd, args.threads, read_mapping, od, barcode,
                                                 args.folder_fast5, args.verbose )

                        if args.use_assembled:
                            if verbose:
                                print("USING " + assembled + " FOR BLAST")
                            if args.verbose:
                                print(assembled)
                            if assembled != "":
                                for record in SeqIO.parse(assembled, "fasta"):
                                    if record.seq != "":
                                        fp = os.path.join(dirpathTEMP,
                                                          record.id + ".fasta")  # tempfile.NamedTemporaryFile(dir=dirpathTEMP, suffix=".fasta", mode="w", delete=False)
                                        fo = os.path.join(dirpathTEMP,
                                                          record.id + ".blastn")  # tempfile.NamedTemporaryFile(dir=dirpathTEMP, mode="w", prefix=barcode, suffix=".blastn",
                                        SeqIO.write(record, fp, "fasta")
                                        fastx_all.append([fp, database_name, fo, args.verbose, args.cleanup,job_id])
                            else:
                                continue
                    if len(fastx_all) == 0:
                        if verbose:
                            print("USING UNCORRECTED READS FOR BLAST")
                        for record in SeqIO.parse(reads_porechop_1.name, "fastq"):
                            if record.seq != "":
                                fp = os.path.join(dirpathTEMP,
                                                  record.id + ".fasta")  # tempfile.NamedTemporaryFile(dir=dirpathTEMP, suffix=".fasta", mode="w", delete=False)
                                fo = os.path.join(dirpathTEMP,
                                                  record.id + ".blastn")  # tempfile.NamedTemporaryFile(dir=dirpathTEMP, mode="w", prefix=barcode, suffix=".blastn",
                                SeqIO.write(record, fp, "fasta")
                                fastx_all.append([fp, database_name, fo, args.verbose, args.cleanup, job_id])
                                if not barcode.endswith("raw"):
                                    if barcode.endswith("/"):
                                        barcode = barcode.split("/")[-2] + "raw"
                                    else:
                                        barcode = barcode.split("/")[-1] + "raw"
            else:
                if verbose:
                    print("USING " + reads_porechop_1.name + " FOR BLAST")
                if 'fastx_all' in locals():
                    for record in SeqIO.parse(reads_porechop_1.name, "fastq"):
                        if record.seq != "":
                            fp = os.path.join(dirpathTEMP, record.id + ".fasta")# tempfile.NamedTemporaryFile(dir=dirpathTEMP, suffix=".fasta", mode="w", delete=False)
                            fo = os.path.join(dirpathTEMP, record.id + ".blastn") #tempfile.NamedTemporaryFile(dir=dirpathTEMP, mode="w", prefix=barcode, suffix=".blastn",
                            SeqIO.write(record, fp, "fasta")
                            fastx_all.append([fp, database_name,fo, args.verbose, args.cleanup, job_id])
                            if not barcode.endswith("raw"):
                                if barcode.endswith("/"):
                                    barcode = barcode.split("/")[-2]+ "raw"
                                else:
                                    barcode = barcode.split("/")[-1]+ "raw"

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
                        fp = os.path.join(dirpathTEMP, record.id + ".fasta")# tempfile.NamedTemporaryFile(dir=dirpathTEMP, suffix=".fasta", mode="w", delete=False)
                        fo = os.path.join(dirpathTEMP, record.id + ".blastn") #tempfile.NamedTemporaryFile(dir=dirpathTEMP, mode="w", prefix=barcode, suffix=".blastn",
                                                        # delete=False)
                        SeqIO.write(record, fp, "fasta")
                        correct_racon_blast.append([fp, fastamlst,fo, args.verbose, args.cleanup, job_id])

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
            if args.slurm:
                result_list = []
                for index, line in enumerate(fastx_all):
                    blast_slurm(line)
                    result_list.append(line[2])
                    if args.slurm_threads > 0:
                        if (index + 1) % args.slurm_threads == 0:
                            time.sleep(90)  # Pause the loop for 5 seconds

                while True:
                    result = sb.run(['squeue'], stdout=sb.PIPE, stderr=sb.PIPE, text=True)
                    print("WAITING FOR SLURM TO FINISH")
                    if not job_id in result.stdout:
                        break
                    time.sleep(15)
                # for p in Path(od).glob("slurm-*.out"):
                #     p.unlink()
            else:
                cores = args.threads
                print("RUNNING BLAST USING "  + str(cores) + " THREADS")
                with Pool(processes=cores) as pool:
                    for result in tqdm(pool.imap(func=blast, iterable=fastx_all), total=len(fastx_all)):
                        result_list.append(result)
            existing_files = [f for f in result_list if os.path.exists(f)]
            dict_species_new ={}
            barcode_all_excel =  os.path.join(od, barcode + "_blastn_results.txt")
            with open(barcode_all_excel, 'w') as file:
                input_lines = fileinput.input(existing_files)
                file.writelines(input_lines)


            blastn_results = []
            with Pool(processes=args.threads) as pool:
                for result in tqdm(pool.imap(func=parse_blast, iterable=existing_files), total=len(result_list)):
                    blastn_results.append(result)
            # for fasta_single in fastx_all:
            #     os.remove(fasta_single[0])

            count_blast_hit = 0
            for species in blastn_results:
                if len(species) == 7:
                    count_blast_hit += 1
                    if verbose:
                        list_species_verbose.append(species)
                    if species[0] in dict_species_new:
                        dict_species_new[species[0]] = dict_species_new[species[0]] + 1
                    else:
                        dict_species_new[species[0]] = 1
            if count_blast_hit == 0:
                continue
            print("NUMBER OF READS CLASSIFIED " + str(count_blast_hit) + " RAPPRESENTING " + str((count_blast_hit/len(result_list)) *100) + " OF THE TOTAL NUMBER OF READS" )

            if args.verbose:
                name_ranks = os.path.join(od, args.output + barcode + ".full_ranks.txt")
                with open(name_ranks, "w") as fh:
                    for line in list_species_verbose:
                        try:
                            taxname = ncbi.get_taxid_translator([str(line[0])])
                            seqid = line[4]
                            fh.write("\t".join([(list(taxname.values()))[0],seqid + "\n"]))
                        except:
                            print("RANK NOT FOUND")
            dict_name = {}
            desired_ranks = ["species","genus","family","order","class","phylum","clade","superkingdom","subspecies","species subgroup","species group"]
            for key in dict_species_new:
                try:
                    name_species = ncbi.get_taxid_translator([key])
                    speciesname = (list(name_species.values()))[0]
                except:
                    speciesname = key
                dict_name[speciesname] = dict_species_new[key]
            dict_table_all[barcode] = dict_name
            kt_barcode = tempfile.NamedTemporaryFile(dir=dirpathTEMP, suffix=".txt", prefix=barcode, delete=False)
            files_kt_list.append(kt_barcode.name)
            seq_count = 0
            total = 0
            for key in dict_species_new:
                total +=  dict_species_new[key]
            min_value_to_use = (int(total) * int(args.min) / 100)
            with open(kt_barcode.name, "w") as fh:
                if args.verbose:
                    print(kt_barcode.name)
                for fastq_line in blastn_results:
                    if len(fastq_line) > 3:
                        output_kt = fastq_line[2] + "\t" + fastq_line[0] + "\n"
                        fh.write(output_kt)
            taxid_list = [(key, dict_species_new[key]) for key in dict_species_new]
            taxid_dict_rank = desired_rank_file(taxid_list, desired_ranks)#, rank_path)
            taxid_pd_rank = pd.DataFrame.from_dict(taxid_dict_rank, orient='index')
            taxid_pd_rank['abundance']  = (taxid_pd_rank['abundance'] / taxid_pd_rank['abundance'].sum())
            taxid_pd_rank[taxid_pd_rank['abundance']  < args.min / 100] = None
            taxid_pd_rank.dropna(subset=['abundance'], inplace=True)
            taxid_pd_rank.to_csv(rank_path,sep="\t",index=False)
        else:
            print("PATH FOR BARCODE " + barcode + " DO NOT EXISTS. CHECK BARCODE NAME OR PATH" )



    species_abs = pd.DataFrame.from_dict(dict_table_all, orient='index')
    species_abs.fillna(0)
    species_abs.transpose().to_excel(name_table)
    for key_barcode in dict_table_all:
        dict_single_barcode = dict_table_all[key_barcode]
        total_sum = sum(dict_single_barcode.values())
        keys_to_remove = [key for key, value in dict_single_barcode.items() if value < (total_sum * (args.min/100))]
        for key in keys_to_remove:
            del dict_single_barcode[key]
        dict_table_all[key_barcode] = dict_single_barcode

    if len(files_kt_list) > 0:
        plothystogramSpecies(files_kt_list,blastdb,args.output,od, args.verbose)
        histograplot(dict_table_all, figure_file_hist)
        heatmap_clustering(dict_table_all, figure_file_heat)
    else:
        print("NOTHING TO PLOT")

    # if not args.verbose:
    #     shutil.rmtree(dirpathTEMP)

if __name__ == '__main__':
    main()
