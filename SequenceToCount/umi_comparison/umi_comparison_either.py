
import argparse
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet.IUPAC import unambiguous_dna, ambiguous_dna
from collections import namedtuple
import csv
#from __future__ import division
import fileinput
import glob
import gzip
import itertools as IT
from multiprocessing import Pool, Lock
import os
import random
import re
import sys
import subprocess

###########################################################################
## Struct Definitions
###########################################################################
  
sampleStruct = namedtuple("sampleStruct","Sample FilePrefix intMultiBCArray")
pairedFastqFileStruct = namedtuple("pairedFastqFileStruct","FwdFastq RevFastq MatchPrefix")

###########################################################################
## Global Variables
###########################################################################

sampleArray = []
indexToSampleMap = {}
templateSeqLengthsDict = {}
templateSeqArray = [] #one to two inner arrays, describing the template corresponding to each read
validFeatureTypes = {"X" : "BC", "D" : "multiplexBC", "U" : "UMI"}
validFeatureCounts = {"X" : 0, "D" : 0, "U" : 0}
usedFastqFiles = []
constantRegionFastaFilenames = [] # same dimensions as templateSeqArray
constantRegionNames = [] # same dimensions as templateSeqArray
allConstRegionsFileName = "allCRs.fasta"
maxNInReads = 3
fileBufferSize = 10000 #make sure we aren't constantly accessing the disk but small enough to reasonably keep in memory
readsPerSampleForErrors = 10000
expectedBCLength = 0

#Blast parameters
constantRegionsBlastParams = ["-word_size", "6","-outfmt","6","-evalue","1E0"]


cmdLineArgParser = argparse.ArgumentParser(description="Process Illumina PCR Amplicon Fastq files to count barcode frequencies. Must have blastn and dada2 (R Bioconductor  package) installed to run. Make sure there are no spaces in any file names or directory paths, the program will not work otherwise. Make sure that files are already split by their illumina indices (N700 / S500 index). This code will not function if there are multiple barcodes or internal multiplexing barcodes within a single read. Please use the readLength flag to solve this issue. \n\nGenerated Files:\n\nnumReadsFoundPerSample: One file per input fastq file, it details the assignment of each read to a combination of inline indices if they exist, and the subsequent filtering of reads. The file has 2 columns, the first being a unique identifier for the fastq file and inline index combination identified, and the second being an array of 6 numbers binning the reads into the following categories: \n\t1. correct reads used for subsequent mapping. \n\t2. There are too many Ns in the read.\n\t3. UMIs have been requested and not identified in the read/\n\t4. There is no barcode found in the read.\n\t5. The combination of indices do not match a sample specified in the sampleFile.\n\t6. A second check for finding a valid barcode in the read.\n\n For each sample in SampleFile with at least one barcode found a number of files are generated.\n R1(R2).fastq - raw (possibly truncated depending on the readLength parameter) reads associated with this sample\n barcode.fastq - the portion of the reads associated with all barcode sequences concatenated together, in the same order as the R1/R2.fastq files\n UMISeqs.tab - tab delmited UMI sequences if they exist and are being used, in the same order as the R1/R2.fastq files.\n readBarcodeID.txt - the ID number of the barcode in the clusteredBCs.fasta file that each read was mapped to, in the same order as the R1/R2.fastq files.\nbarcodeCalls.tab - Count of each barcode in the sample, removing UMI duplicates if asked for. Line 1 contains the counts for barcode 1 (defined in clusteredBCs.fasta), line 2 for barcode 2 etc.\n\nallBarcodeCalls.tab - tab delimited concatenation of all of the barcodeCalls.tab files, with a header row identifying which column comes from which sample.")
cmdLineArgParser.add_argument("-fastqDir", dest="fastqDir", help="directory location of fastq files",required=True)
cmdLineArgParser.add_argument("-outputDir", dest="outputDir", help="location of output directory",required=True)
cmdLineArgParser.add_argument("-templateSeq", dest="templateSeqFile", help="Template sequence of amplicon locus. This file contains a single line with standard DNA bases. UMI (unique molecular identifier) sequences are coded as U, multiplexing indices are coded as D and barcode loci coded as X. If these features have different lengths between samples, define the template using the longest possible length of each feature. Every feature annotated must be covered by the sequencing data, and no feature can span the exact middle of the template sequence when using paired end data.",required=True)
cmdLineArgParser.add_argument("-sample", dest="sampleFile", help="File defining the samples present in the sequencing data. This file is tab delimited, with no header line. The column values are: Sample Name\t File Prefix\t internal multiplexing barcode 1\t internal multiplexing barcode 2... The internal multiplexing barcode columns must correspond to the names of the sequences in the multiBCFasta file. Do not use spaces in any of the columns for file / directory names, as this tends to behave poorly.",required=True)
cmdLineArgParser.add_argument("-multiBCFasta", dest="multiBCFastaFile", help="A multi-line fasta file defining multiplexing tag sequences. Required if there are multiplexing tags within the amplicon sequence as defined by the templateSeq file.")
cmdLineArgParser.add_argument("-barcodeList", dest="barcodeListFile", help="Optional fasta file specifying the barcodes present in the sample. If file is not supplied, unique barcodes will be identified de novo. The name for each sequence must be unique. If the file is being generated manually, the sequence must be a simple concatenation of all barcode regions as defined in the template sequence in 5'-3' order.")

cmdLineArgParser.add_argument("-barcode5PrimeTrimLength", dest="barcode5PrimeTrimLength", default=0,  help="Number of bp to trim from the 5' end of each barcode. ")
cmdLineArgParser.add_argument("-barcode3PrimeTrimLength", dest="barcode3PrimeTrimLength", default=0,  help="Number of bp to trim from the 3' end of each barcode. ")
cmdLineArgParser.add_argument("-bcNGapLength", dest="bcNGapLength", default=0,  help="Number of bp of Ns to put between barcodes from forward and reverse reads. Use if the sequence is not covering the entire barcode and you you have provided the list of valid barcodes using the barcodeList argument.")
cmdLineArgParser.add_argument("-blastPath", dest="blastPath", help="BLAST installation directory if it is not in the Path already", default="")
cmdLineArgParser.add_argument("-useBowtie2", dest="useBowtie2", help="Flag to use Bowtie2 instead of the default BWA mem for barcode mapping. ", action="store_true")
cmdLineArgParser.add_argument("-bowtie2Path", dest="bowtie2Path", help="Bowtie2 installation directory if it is not in the Path already", default="")
cmdLineArgParser.add_argument("-bwaPath", dest="bwaPath", help="BWA installation directory if it is not in the Path already", default="")
cmdLineArgParser.add_argument("-demultiplexOnly", dest="demultiplexOnly", action="store_true",  help="Use flag if you want to only split the raw fastq files and not continue with the rest of the barcode counting. This is useful when distributing demultiplexing across several machines, i.e. in a cluster.")
cmdLineArgParser.add_argument("-DNAclustPath", dest="DNAclustPath", help="DNAclust installation directory if it is not in the Path already", default="")
cmdLineArgParser.add_argument("-numThreads", dest="numThreads", default=1,  help="Number of threads to be used for computation.")
cmdLineArgParser.add_argument("-pairedEnd", dest="pairedEnd", action="store_true",  help="Use if sequencing data is paired end")
cmdLineArgParser.add_argument("-readLength", dest="readLength", default=100,type=int,  help="Expected length of each read from sequencing machine. Default = 100. Reduce this number from the true read length if necessary such that non-constant regions of the barcode locus are not shared between reads. This does not modify the input fastq files, but effectively truncates reads before processing")
cmdLineArgParser.add_argument("-remapBarcodes", dest="remapBarcodes", action="store_true",  help="Set to True if you want to remap barcodes even if the files already exist")
cmdLineArgParser.add_argument("-skipSplitFastq", dest="skipSplitFastq", action="store_true",  help="Use flag if you want to skip the splitting of the raw fastq files (i.e. if you have already done this and do not want to redo it).")
cmdLineArgParser.add_argument("-useUMI", dest="UMI", action="store_true",  help="Use flag if you want to remove PCR duplicate reads using UMI data")
cmdLineArgParser.add_argument("-skipBowtieBuild", dest="skipBowtieBuild", action="store_true",  help="Use flag to skip building the bowtie database (helpful if you are running parallel processes that use this database)")

umi_methods = ['either']


def get_connected_group(node,graph, already_seen):
        result = []
        nodes = set([node])
        while nodes:
            node = nodes.pop()
            already_seen.add(node)
            nodes = nodes or graph[node] - already_seen
            result.append(node)
        return result, already_seen

def get_all_connected_groups(graph):
    already_seen = set()
    result = []
    for node in graph:
        if node not in already_seen:
            connected_group, already_seen = get_connected_group(node, graph, already_seen)
            result.append(connected_group)
    return result
            
def connected_components(neighbors):
    seen = set()
    def component(node):
        nodes = set([node])
        while nodes:
            node = nodes.pop()
            seen.add(node)
            nodes |= neighbors[node] - seen
            yield node
    for node in neighbors:
        if node not in seen:
            yield component(node)

def n_components(graph):
    new_graph = {node: set(edge for edge in edges)
             for node, edges in graph.items()}
    components = []
    for component in connected_components(new_graph):
        c = set(component)
        components.append([edge for edges in graph.values()
                                for edge in edges
                                if c.intersection(edge)])

    return len(components)

def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

def parseTemplateSeq():
    sequence = ""
    with open(args.templateSeqFile) as infile:
        sequence = infile.readline().strip()
    if(len(re.sub("[ACGTacgtXUDNn]","",sequence))>0): # if the sequence has disallowed characters, quit
        eprint("Template sequence has illegal characters!\n")
        sys.exit(0)
    if(args.multiBCFastaFile=="" and sequence.count("D")>0): #if there is no multiplexing file but we find multiplexing loci in the template, quit
        eprint("Multiplexing indices exist but no multiBCFastaFile provided!\n")
        sys.exit(0)
    #if it is paired end, split template using readSize, 
    templateArray = []
    if(args.pairedEnd):
        templateArray.append(sequence[0:int(args.readLength)])
        templateArray.append(sequence[int(len(sequence)-args.readLength):(len(sequence)+1)])
    else:
        templateArray.append(sequence)
        
    expectedBarcodeLength = sequence.count("X")
    #collapse non-constant features into single character and add to global variable
    global templateSeqArray
    for seq in templateArray:
        seq = re.sub("D+","\tD\t",seq)
        seq = re.sub("U+","\tU\t",seq)
        seq = re.sub("X+","\tX\t",seq)
        seq = re.sub("\t\t","\t",seq)
        seq = re.sub("\t\t","\t",seq)
        seq = re.sub("\t\t","\t",seq)
        seq = re.sub("^\t","",seq)
        seq = re.sub("\t$","",seq)
        templateSeqArray.append(seq.split("\t"))

def parseSampleFile():
    global sampleArray
    global indexToSampleMap
    numInlineIndices = 0
    for seqAr in templateSeqArray:
        for seq in seqAr:
            if(seq == "D"):
                numInlineIndices+=1
    with open(args.sampleFile) as infile: # go through every line in the sample file
        for line in infile:
            if(line.strip() == ""):
                continue
            lineSplit = line.strip().split("\t") #separate values in each line of the sample file
            if(len(lineSplit) != 2 + numInlineIndices): # if the line has not enough columns, error and quit
                eprint("Sample file line:\n"+line+"has incorrect number of columns!\n")
                eprint(str(len(lineSplit)))
                eprint(str(numInlineIndices))
                sys.exit(0)
            sampleIndexArray = "_".join(lineSplit[1:(len(lineSplit)+1)])
            sampleIdentityArray = lineSplit[0]
            
            if sampleIndexArray in indexToSampleMap: #if this particular file / index combination has been seen already, quit
                eprint("Sample file line:\n"+line+"has been assigned to a sample already!\n")
                sys.exit(0)
            
            #else catalog the sample
            
            sampleArray.append(sampleStruct(lineSplit[0],lineSplit[1],lineSplit[2:(len(lineSplit)+1)]))
            indexToSampleMap[sampleIndexArray] = sampleIdentityArray

def identifyUsedFastqFiles():
    usedFastqFileDict = {}
    for sample in sampleArray:
        patternString = re.compile(".*"+sample.FilePrefix+".*.fastq(.gz)?")
        readFiles = glob.glob(args.fastqDir+"*")
        readFiles2 = list(filter(patternString.match,readFiles))
        readFiles2.sort()
        myfwd = None
        myrev = None
        if(len(readFiles2)==0 or (len(readFiles2)==1 and args.pairedEnd) or (len(readFiles2)==2 and not args.pairedEnd) or len(readFiles2)>2): #check that we found the expected number of fastq files (based on single or paired end data expected)
            eprint(readFiles2)
            eprint("Incorrect number of matching fastq files found for "+sample.FilePrefix)
            sys.exit(1)
        myfwd = readFiles2[0]
        if(args.pairedEnd):
            myrev = readFiles2[1]
        usedFastqFileDict[pairedFastqFileStruct(FwdFastq=myfwd,RevFastq=myrev,MatchPrefix=sample.FilePrefix)]=1
    global usedFastqFiles
    usedFastqFiles = usedFastqFileDict.keys()

def mapBarcodes(mySamp):
    #only run on this sample if the output file doesn't exist or flag has been set
    indexString = mySamp
    if (os.path.isfile(args.outputDir+indexString+"_barcode.fastq") and (not os.path.isfile(args.outputDir+indexString+"_readBarcodeID.txt") or args.remapBarcodes)):
    
        #make a dictionary to map barcode names in the barcode list fasta file to consecutive numbers for indexing in a vector.
        BCNameToIdxDict = {}
        totalNumBCs = 1
        with open(args.barcodeListFile,"r") as infile:
           
            for record in SeqIO.parse(infile,"fasta"):
                if(record.id in BCNameToIdxDict): #we have found a duplicate entry in the barcode list. quit!
                    eprint("Duplicate entry "+record.id+" found in input barcode list!")
                    sys.exit(1)
                BCNameToIdxDict[record.id]=totalNumBCs
                totalNumBCs +=1
        
        bcFastqFile = args.outputDir+indexString+"_barcode.fastq"
        bcSamFile = args.outputDir+indexString+"_barcode.sam"
        bcIDFile = args.outputDir+indexString+"_readBarcodeID.txt"
        mapQualFile = args.outputDir+indexString+"_readMappingQuality.txt"
        
        # if(args.useBowtie2): #bowtie2 call if flagged
        #     subprocess.call([args.bowtie2Path+"bowtie2","-L 10","-q","--very-sensitive-local","-x "+args.barcodeListFile,"-U"+bcFastqFile,"-S"+bcSamFile])
        # else: #bwa mem call
        #     with open(bcSamFile, "w") as outfile:
        #         subprocess.call([args.bwaPath+"bwa","mem","-k 10","-y 12",args.barcodeListFile,bcFastqFile], stdout = outfile)
        
        # #get the barcode match for each read and put into a single column output file. do the same for mapping quality
        # with open(bcIDFile,"w") as outfile:
        #     subprocess.call("grep -v '^@' "+bcSamFile+" | cut -f 3",stdout=outfile, shell=True)
        # with open(mapQualFile,"w") as outfile:
        #     subprocess.call("grep -v '^@' "+bcSamFile+" | cut -f 5",stdout=outfile, shell=True)
        

        BCUMIMap = {}
        BCCountList = {}
        BCUMIDupCountList = {}
        totalUnmappedReads = 0

        bcIDFileHandle = open(bcIDFile,"r")
        mapQFileHandle = open(mapQualFile,"r")
        UMIFileHandle = open(args.outputDir+indexString+"_UMISeqs.tab","r")

        for umi_method in umi_methods:
            BCUMIMap[umi_method] = {}
            BCCountList[umi_method] = [0]*int(totalNumBCs-1)
            BCUMIDupCountList[umi_method] = [0]*int(totalNumBCs-1)
            # totalUnmappedReads[umi_method] = 0
        
        unique_tracker_either = {} # initialize unique tracker

        #for each read bc / umi pair
        for bcID, UMIstring, mapQ in zip(bcIDFileHandle, UMIFileHandle, mapQFileHandle):
            bcID = bcID.strip()
            mapQ = mapQ.strip()

            # f_umi = UMIstring.split('\t')[0]
            # r_umi = UMIstring.split('\t')[1]

            if '\t' in UMIstring:
                f_umi = UMIstring.split('\t')[0]
                r_umi = UMIstring.split('\t')[1]
            else:
                f_umi = ''
                r_umi = ''
                # just skip ones with only one or the other
           


            if(bcID == "*" or bcID == "" or mapQ == "*" or int(mapQ) <= 20):
                totalUnmappedReads = totalUnmappedReads + 1
                continue
            else:
                bcID = BCNameToIdxDict[bcID] #get the internal index corresponding to the matched barcode

                for umi_method in umi_methods:

                    if umi_method == 'either':

                        bc = str(bcID)

                        f_umi = 'F'+f_umi
                        r_umi = 'R'+r_umi

                        if bc not in unique_tracker_either.keys():
                            unique_tracker_either[bc] = {}
                            unique_tracker_either[bc][f_umi] = set([r_umi])
                            unique_tracker_either[bc][r_umi] = set([f_umi])
                        else:
                            if f_umi not in unique_tracker_either[bc].keys():
                                unique_tracker_either[bc][f_umi] = set([r_umi])
                            else:
                                unique_tracker_either[bc][f_umi].add(r_umi)
                            if r_umi not in unique_tracker_either[bc].keys():
                                unique_tracker_either[bc][r_umi] = set([f_umi])
                            else:
                                unique_tracker_either[bc][r_umi].add(f_umi)

        for umi_method in umi_methods:
            for bc in unique_tracker_either.keys():
                this_bc_components = n_components(unique_tracker_either[bc])
                BCCountList[umi_method][int(bc)-1] = this_bc_components
                BCUMIDupCountList[umi_method][int(bc)-1] = len(unique_tracker_either[bc])-this_bc_components
        
        # the file is parsed, close filehandles
        bcIDFileHandle.close()
        UMIFileHandle.close()
        mapQFileHandle.close()

        
        for umi_method in umi_methods:

            #write total count data to file
            with open(args.outputDir+'umi_comparison/'+f'{umi_method}/'+indexString+"_barcodeCounts.tab","w") as outFileHandle:
                for countVal in BCCountList[umi_method]:
                    outFileHandle.write(str(countVal)+"\n")
            with open(args.outputDir+'umi_comparison/'+f'{umi_method}/'+indexString+"_numUnmappedReads.txt","w") as outFileHandle:
                outFileHandle.write(str(totalUnmappedReads))
            if(args.UMI):
                with open(args.outputDir+'umi_comparison/'+f'{umi_method}/'+indexString+"_barcodeUMIDupCounts.tab","w") as outFileHandle:
                    for countVal in BCUMIDupCountList:
                        outFileHandle.write(str(countVal)+"\n")

args = cmdLineArgParser.parse_args()
args.bcNGapLength = int(args.bcNGapLength)
parseTemplateSeq()
parseSampleFile()
identifyUsedFastqFiles()
# createConstRegionFasta()

uniqueSamples = {}
for sample in sampleArray:
    uniqueSamples[sample.Sample] = 1

with Pool(processes = int(args.numThreads)) as pool:
    pool.map(mapBarcodes, uniqueSamples.keys())

