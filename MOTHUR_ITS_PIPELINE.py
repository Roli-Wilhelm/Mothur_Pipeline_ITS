#!/usr/bin/python
import sys, os, re, getopt, glob, numpy as np
import timeit
import itertools

start = timeit.default_timer()

Usage = """
Usage:  ./MOTHUR_ITS_PIPELINE.py -o FUNGI_0.07 -n FUNGI_0.07 -i IIKFCBR02.fasta -f IIKFCBR02.oligos -d 11 -l 200 -p 8

REQUIRED ARGUMENTS:
		-o	output directory
	
                -n      the name given to the project

                -i      the input file (accepted: .fasta, .sff, .shhh.fasta. MUST SPECIFY if the latter two)
			(if providing .fasta, one must provide a .qual file)

		-f	the .oligos file

		-d	max number of base-pair differences used in clustering
			(For example, if your sequences will be chopped to 200bp, giving a value of 10 will yield a dissimilarity of 0.05)

		-l	maximum length of sequence. ALL sequences will be chopped to this uniform length.

		-p	the number of processors available to use

DEPENDENCY:
		- mothur v.1.32.1 (written using, probably valide for earlier and later versions)
		- UNITE db v6 
		- crunchclust
		- deunique_mothurlist.pl

OPTIONAL:
		-S <Y>	Input is .sff file
			(Must provide "Y")

		-R <Y>	SELECT IF you've already run .shhh flows and do not want to re-run (obviously b/c it takes forever), choose this option.
			(Required: *.shhh.names, *.shhh.groups, *.shhh.fasta)

		-C <Y>	SELECT IF you've already run this script on the same data and you wish to just re-do the OTU clustering step.
			(Required: *.trim.unique.pick.chop.precluster.unique.fasta, *.trim.unique.pick.chop.precluster.names, *.pick.chop.groups)

		-O 	If you get an error regarding your "flow.files", there is a chance you must specify a different option for flow file to use (ex. B). Do so here.
	
                -b <V>  Enter DEBUG MODE - saves debugging information to "error.log" 
                        (Must provide "V")

UTILITY:
This script will performs the general steps in the mothur pipeline for fungal ITS sequences. This script will require tailoring for different datasets, but provides a simple archetype as well as a quick implementation of whatever you decide. The script does do "summary.seqs" at various points along the way, but the user will have to manually examine the outputs to see how the cleaning went. There is one point where manual intervention is required and the script will wait until one does so."


NOTES: 
1) Minflows is set to 180 (i.e. the minimum length of sequence will be 180bp)
2) Maxflows is set to 450 (i.e. the maximum length of sequence will be 450bp)

"""

if len(sys.argv)<3:
        print Usage
        sys.exit()


# Store input and output file names
NAME=''
INPUT=''
OLIGOS=''
DIFF=''
RERUN=''
LENGTH=''
SFF=''
CRUNCH=''
ORDER=''
PROCESSORS=''
OUTPUT=''
Debug=''

## Read command line args
myopts, args = getopt.getopt(sys.argv[1:],"n:i:o:S:b:d:R:l:f:C:O:p:")

###############################
# o == option
# a == argument passed to the o
###############################
for o, a in myopts:
    if o == '-o':
        OUTPUT= a
    if o == '-n':
        NAME= a
    if o == '-i':
        INPUT= a
    if o == '-d':
        DIFF= a
    if o == '-l':
        LENGTH= a
    if o == '-f':
        OLIGOS= a
    if o == '-S':
        SFF= a
    if o == '-R':
        RERUN= a
    if o == '-C':
        CRUNCH= a
    if o == '-O':
        ORDER= a
    if o == '-p':
        PROCESSORS= a
    if o == '-b':
        Debug= a

## Print Debug Info (not really used much in this script)
if Debug:
        print "You are in Debugging Mode"
        error = open("error.log", "w")

## Create Output Directory
if len(OUTPUT)>0:
        if os.path.exists('./' + OUTPUT):
                print "\nOutput Folder Exists - Caution: Files May Be Re-Written"
        else:
                os.mkdir(OUTPUT)

## Get Basename of Input File to Use in Naming Throughout the Script
if SFF:
	INPUT = re.sub(".sff", "", INPUT)
	print INPUT

	if Debug:
        	print "The basename of your input file is:\n"
	        error.write(INPUT+"\n")

elif RERUN:
	INPUT = re.sub(".shhh.fasta", "", INPUT)
	print INPUT

	if Debug:
        	print "The basename of your input file is:\n"
	        error.write(INPUT+"\n")
else:
	INPUT = re.sub(".fasta", "", INPUT)
	print INPUT

	if Debug:
        	print "The basename of your input file is:\n"
	        error.write(INPUT+"\n")

## Same for oligos file
OLIGOS = re.sub(".oligos", "", OLIGOS)
print OLIGOS

## Unpack SFF and Perform Quality Filtering with SHHH Flows
if SFF:
	## The chemistry CAN be different for various sequencing runs, MOTHUR provides for this based on the "ORDER" of flows are read. This is sometimes necessary to specify.
	if ORDER:
		os.system(' '.join([
			"for n in",
			INPUT+".sff; do mothur \"# sffinfo(sff=$n, flow=T)\"; done"
		]))

		os.system(' '.join([
			"for n in",
			INPUT+".flow; do mothur \"# trim.flows(flow=$n,",
			"oligos="+OLIGOS+".oligos,"
			"pdiffs=2, bdiffs=1, minflows=180, maxflows=450,",
			"order="+ORDER+",",
			"processors="+PROCESSORS+")\"; done"
		]))

		os.system(' '.join([
			"for n in",
			INPUT+".flow.files; do mothur \"# shhh.flows(file=$n,",
			"order="+ORDER+",",
			"processors="+PROCESSORS+")\"; done"
		]))

	else:
		os.system(' '.join([
			"for n in",
			INPUT+".sff; do mothur \"# sffinfo(sff=$n, flow=T)\"; done"
		]))

		os.system(' '.join([
			"for n in",
			INPUT+".flow; do mothur \"# trim.flows(flow=$n,",
			"oligos="+OLIGOS+".oligos,"
			"pdiffs=2, bdiffs=1, minflows=180, maxflows=450,",
			"processors="+PROCESSORS+")\"; done"
		]))

		os.system(' '.join([
			"for n in",
			INPUT+".flow.files; do mothur \"# shhh.flows(file=$n,",
			"processors="+PROCESSORS+")\"; done"
		]))


## Perform regular pipeline on data from .sff or from data previously cleaned using shhh.flows
if SFF or RERUN:
	## If user specified to skip to Crunch Clust step, do it
	if CRUNCH:
		os.system(' '.join([
			"crunchclust",
			"--diff "+DIFF+" --in",
			INPUT+".shhh.trim.unique.pick.chop.precluster.unique.fasta",
			"--out ./crunch.clstr --d_all --endgaps"
		]))
	
		## Get user to do some manual manipulatios. I tried coding these, but the first had too many special charaters that were causing headaches and since we're already stopping it makes it easier to have the user do the second.
		print "Open an alternate window (keep it open for the next two steps).\nAt this point you will have to manually enter the command located in the file entitled \"post_crunchclust.add\". Follow the instructions given in that file.\n"
		print "When you are done, hit enter"
		go = raw_input()

		print "At this point you will have to manually make changes to the \"Fungi_temp.list\" file. Please do the following:\n"
		print "Using a text editor like vim or nano in a separate window, add the dissimilarity cutoff followed by the number of sequences in your file separated by a tab.\n For example, \"0.05	10340\""
		print "When you are done, hit enter"
		go = raw_input()

		## Use Deunique_mothurlist.pl to prepare the names file based on crunchclust output
		os.system(' '.join([
			"perl deunique_mothurlist.pl",
			"-n "+INPUT+".shhh.trim.unqiue.pick.chop.precluster.unique.names",
			"-m Fungi_temp.list",
			"-o "+NAME+"_final.list"
		]))

		## Simplify All Names to Final Versions
		os.system(' '.join([
			"cp",
			INPUT+".shhh.trim.unique.pick.chop.precluster.unique.names",
			NAME+"_final.names"
		]))
	
		os.system(' '.join([
			"cp",
			INPUT+".shhh.trim.unique.pick.chop.precluster.unique.fasta",
			NAME+"_final.fasta"
		]))

		os.system(' '.join([
			"cp",
			INPUT+".shhh.pick.chop.groups",
			NAME+"_final.groups"
		]))

		os.system(' '.join([
			"mv",
			"Fungi_temp.list",
			NAME+"_final.an.list"
		]))

		## Do Classification Using UNITE.db
		os.system(' '.join([
			"for n in",
			NAME+"_final.fasta; do mothur \"# classify.seqs(fasta=$n,",
			"name="+NAME+"_final.names,",
			"group="+NAME+"_final.groups,",
			"template=~/Phylogenetic_Gene_Databases/UNITE_ITS/UNITEv6_sh_dynamic.fasta,",
			"taxonomy=~/Phylogenetic_Gene_Databases/UNITE_ITS/UNITEv6_sh_dynamic.tax, cutoff=50,",
			"processors="+PROCESSORS+")\"; done"
		]))

		os.system(' '.join([
			"mv",
			NAME+"_final.UNITEv6_sh_dynamic.wang.taxonomy",
			"./"+OUTPUT+"/"+NAME+"_final.taxonomy"
		]))

		## Move all Final Files to OUTPUT directory
		os.system(' '.join([
			"mv",
			NAME+"_final.*",
			"./"+OUTPUT+"/"
		]))

	else:
		## Start from the earliest point in the MOTHUR pipeline

		## Make Record of Starting Stastics of the Library
		os.system(' '.join([
			"for n in",
			INPUT+".shhh.fasta; do mothur \"# summary.seqs(fasta=$n,",
			"name="+INPUT+".shhh.names,",
			"processors="+PROCESSORS+")\"; done"
		]))

		## More QC
		os.system(' '.join([
			"for n in",
			INPUT+".shhh.fasta; do mothur \"# trim.seqs(fasta=$n,",
			"oligos="+OLIGOS+".oligos, name="+INPUT+".shhh.names,",
			"maxambig=0, maxhomop=8, bdiffs=1, pdiffs=2, minlength=180,",
			"processors="+PROCESSORS+")\"; done"
		]))

		os.system(' '.join([
			"for n in",
			INPUT+".shhh.trim.fasta; do mothur \"# unique.seqs(fasta=$n,",
			"name="+INPUT+".shhh.names)\"; done"
		]))

		os.system(' '.join([
			"for n in",
			INPUT+".shhh.trim.unique.fasta; do mothur \"# summary.seqs(fasta=$n,",
			"name="+INPUT+".shhh.trim.names,",
			"processors="+PROCESSORS+")\"; done"
		]))

		os.system(' '.join([
			"for n in",
			INPUT+".shhh.trim.unique.fasta; do mothur \"# chimera.uchime(fasta=$n,",
			"name="+INPUT+".shhh.trim.names,",
			"group="+INPUT+".shhh.groups,",
			"processors="+PROCESSORS+")\"; done"
		]))

		os.system(' '.join([
			"for n in",
			INPUT+".shhh.trim.unique.uchime.accnos; do mothur \"# remove.seqs(accnos=$n,",
			"fasta="+INPUT+".shhh.trim.unique.fasta,",
			"name="+INPUT+".shhh.trim.names,",
			"group="+INPUT+".shhh.groups)\"; done"
		]))

		## Chop Sequences to Specified Length
		os.system(' '.join([
			"for n in",
			INPUT+".shhh.trim.unique.pick.fasta; do mothur \"# chop.seqs(fasta=$n,",
			"name="+INPUT+".shhh.trim.pick.names,",
			"group="+INPUT+".shhh.pick.groups,",
			"numbases="+LENGTH+", keep=front,",
			"processors="+PROCESSORS+")\"; done"
		]))
	
		os.system(' '.join([
			"for n in",
			INPUT+".shhh.trim.unique.pick.chop.fasta; do mothur \"# summary.seqs(fasta=$n,",
			"name="+INPUT+".shhh.trim.pick.chop.names,",
			"processors="+PROCESSORS+")\"; done"
		]))

		## Reduce the number of sequences prior to Crunch Clust
		os.system(' '.join([
			"for n in",
			INPUT+".shhh.trim.unique.pick.chop.fasta; do mothur \"# pre.cluster(fasta=$n,",
			"name="+INPUT+".shhh.trim.pick.chop.names,",
			"group="+INPUT+".shhh.pick.chop.groups, diffs=2,",
			"processors="+PROCESSORS+")\"; done"
		]))

		os.system(' '.join([
			"for n in",
			INPUT+".shhh.trim.unique.pick.chop.precluster.fasta; do mothur \"# unique.seqs(fasta=$n,",
			"name="+INPUT+".shhh.trim.unique.pick.chop.precluster.names)\"; done"
		]))

		## Perform Crunch Clust as above
		os.system(' '.join([
			"crunchclust",
			"--diff "+DIFF+" --in",
			INPUT+".shhh.trim.unique.pick.chop.precluster.unique.fasta",
			"--out ./crunch.clstr --d_all --endgaps"
		]))

		print "Open an alternate window (keep it open for the next two steps).\nAt this point you will have to manually enter the command located in the file entitled \"post_crunchclust.add\". Follow the instructions given in that file.\n"
		print "When you are done, hit enter"
		go = raw_input()

		print "At this point you will have to manually make changes to the \"Fungi_temp.list\" file. Please do the following:\n"
		print "Using a text editor like vim or nano in a separate window, add the dissimilarity cutoff followed by the number of sequences in your file separated by a tab.\n For example, \"0.05	10340\""
		print "When you are done, hit enter"
		go = raw_input()

		os.system(' '.join([
			"perl deunique_mothurlist.pl",
			"-n "+INPUT+".shhh.trim.unqiue.pick.chop.precluster.unique.names",
			"-m Fungi_temp.list",
			"-o "+NAME+"_final.list"
		]))

		os.system(' '.join([
			"cp",
			INPUT+".shhh.trim.unique.pick.chop.precluster.unique.names",
			NAME+"_final.names"
		]))
	
		os.system(' '.join([
			"cp",
			INPUT+".shhh.trim.unique.pick.chop.precluster.unique.fasta",
			NAME+"_final.fasta"
		]))

		os.system(' '.join([
			"cp",
			INPUT+".shhh.pick.chop.groups",
			NAME+"_final.groups"
		]))

		os.system(' '.join([
			"mv",
			"Fungi_temp.list",
			NAME+"_final.an.list"
		]))

		os.system(' '.join([
			"for n in",
			NAME+"_final.fasta; do mothur \"# classify.seqs(fasta=$n,",
			"name="+NAME+"_final.names,",
			"group="+NAME+"_final.groups,",
			"template=~/Phylogenetic_Gene_Databases/UNITE_ITS/UNITEv6_sh_dynamic.fasta,",
			"taxonomy=~/Phylogenetic_Gene_Databases/UNITE_ITS/UNITEv6_sh_dynamic.tax, cutoff=50,",
			"processors="+PROCESSORS+")\"; done"
		]))

		os.system(' '.join([
			"mv",
			NAME+"_final.UNITEv6_sh_dynamic.wang.taxonomy",
			"./"+OUTPUT+"/"+NAME+"_final.taxonomy"
		]))

		os.system(' '.join([
			"mv",
			NAME+"_final.*",
			"./"+OUTPUT+"/"
		]))

## If user wants to provide only a .fasta and .qual file for analysis.
else:		
	if CRUNCH:
		os.system(' '.join([
			"crunchclust",
			"--diff "+DIFF+" --in",
			INPUT+".trim.unique.pick.chop.precluster.unique.fasta",
			"--out ./crunch.clstr --d_all --endgaps"
		]))

		print "Open an alternate window (keep it open for the next two steps).\nAt this point you will have to manually enter the command located in the file entitled \"post_crunchclust.add\". Follow the instructions given in that file.\n"
		print "When you are done, hit enter"
		go = raw_input()

		print "At this point you will have to manually make changes to the \"Fungi_temp.list\" file. Please do the following:\n"
		print "Using a text editor like vim or nano in a separate window, add the dissimilarity cutoff followed by the number of sequences in your file separated by a tab.\n For example, \"0.05	10340\""
		print "When you are done, hit enter"
		go = raw_input()

		os.system(' '.join([
			"perl deunique_mothurlist.pl",
			"-n "+INPUT+".trim.unique.pick.chop.precluster.unique.names",
			"-m Fungi_temp.list",
			"-o "+NAME+"_final.list"
		]))

		os.system(' '.join([
			"cp",
			INPUT+".trim.unique.pick.chop.precluster.unique.names",
			NAME+"_final.names"
		]))

		os.system(' '.join([
			"cp",
			INPUT+".trim.unique.pick.chop.precluster.unique.fasta",
			NAME+"_final.fasta"
		]))

		os.system(' '.join([
			"cp",
			INPUT+".pick.chop.groups",
			NAME+"_final.groups"
		]))

		os.system(' '.join([
			"mv",
			"Fungi_temp.list",
			NAME+"_final.an.list"
		]))

		os.system(' '.join([
			"for n in",
			NAME+"_final.fasta; do mothur \"# classify.seqs(fasta=$n,",
			"name="+NAME+"_final.names,",
			"group="+NAME+"_final.groups,",
			"template=~/Phylogenetic_Gene_Databases/UNITE_ITS/UNITEv6_sh_dynamic.fasta,",
			"taxonomy=~/Phylogenetic_Gene_Databases/UNITE_ITS/UNITEv6_sh_dynamic.tax, cutoff=50,",
			"processors="+PROCESSORS+")\"; done"
		]))

		os.system(' '.join([
			"mv",
			NAME+"_final.UNITEv6_sh_dynamic.wang.taxonomy",
			"./"+OUTPUT+"/"+NAME+"_final.taxonomy"
		]))

		os.system(' '.join([
			"mv",
			NAME+"_final.*",
			"./"+OUTPUT+"/"
		]))

	else:
		os.system(' '.join([
			"for n in",
			INPUT+".fasta; do mothur \"# summary.seqs(fasta=$n,",
			"processors="+PROCESSORS+")\"; done"
		]))
	
		os.system(' '.join([
			"for n in",
			INPUT+".fasta; do mothur \"# trim.seqs(fasta=$n, oligos="+OLIGOS+".oligos, qfile="+INPUT+".qual, maxambig=0, maxhomop=8, bdiffs=1, pdiffs=2, qwindowaverage=35, qwindowsize=50, minlength=180,",
			"processors="+PROCESSORS+")\"; done"
		]))

		os.system(' '.join([
			"for n in",
			INPUT+".trim.fasta; do mothur \"# summary.seqs(fasta=$n,",
			"processors="+PROCESSORS+")\"; done"
		]))

		os.system(' '.join([
			"for n in",
			INPUT+".trim.fasta; do mothur \"# unique.seqs(fasta=$n,",
			"name="+INPUT+".names)\"; done"
		]))

		os.system(' '.join([
			"for n in",
			INPUT+".trim.unique.fasta; do mothur \"# summary.seqs(fasta=$n,",
			"name="+INPUT+".trim.names,",
			"processors="+PROCESSORS+")\"; done"
		]))

		os.system(' '.join([
			"for n in",
			INPUT+".trim.unique.fasta; do mothur \"# chimera.uchime(fasta=$n,",
			"name="+INPUT+".trim.names,",
			"group="+INPUT+".groups,",
			"processors="+PROCESSORS+")\"; done"
		]))

		os.system(' '.join([
			"for n in",
			INPUT+".trim.unique.uchime.accnos; do mothur \"# remove.seqs(accnos=$n,",
			"fasta="+INPUT+".trim.unique.fasta,",
			"name="+INPUT+".trim.names,",
			"group="+INPUT+".groups)\"; done"
		]))

		os.system(' '.join([
			"for n in",
			INPUT+".trim.unique.pick.fasta; do mothur \"# chop.seqs(fasta=$n,",
			"name="+INPUT+".trim.pick.names,",
			"group="+INPUT+".pick.groups,",
			"processors="+PROCESSORS+")\"; done"
		]))

		os.system(' '.join([
			"for n in",
			INPUT+".trim.unique.pick.chop.fasta; do mothur \"# summary.seqs(fasta=$n,",
			"name="+INPUT+".trim.pick.chop.names,",
			"numbases="+LENGTH+", keep=front,",
			"processors="+PROCESSORS+")\"; done"
		]))

		os.system(' '.join([
			"for n in",
			INPUT+".trim.unique.pick.chop.fasta; do mothur \"# pre.cluster(fasta=$n,",
			"name="+INPUT+".trim.pick.chop.names,",
			"group="+INPUT+".pick.chop.groups, diffs=2,",
			"processors="+PROCESSORS+")\"; done"
		]))

		os.system(' '.join([
			"for n in",
			INPUT+".trim.unique.pick.chop.precluster.fasta; do mothur \"# unique.seqs(fasta=$n,",
			"name="+INPUT+".trim.unique.pick.chop.precluster.names)\"; done"
		]))

		os.system(' '.join([
			"crunchclust",
			"--diff "+DIFF+" --in",
			INPUT+".trim.unique.pick.chop.precluster.unique.fasta",
			"--out ./crunch.clstr --d_all --endgaps"
		]))

		print "Open an alternate window (keep it open for the next two steps).\nAt this point you will have to manually enter the command located in the file entitled \"post_crunchclust.add\". Follow the instructions given in that file.\n"
		print "When you are done, hit enter"
		go = raw_input()

		print "At this point you will have to manually make changes to the \"Fungi_temp.list\" file. Please do the following:\n"
		print "Using a text editor like vim or nano in a separate window, add the dissimilarity cutoff followed by the number of sequences in your file separated by a tab.\n For example, \"0.05	10340\""
		print "When you are done, hit enter"
		go = raw_input()

		os.system(' '.join([
			"perl deunique_mothurlist.pl",
			"-n "+INPUT+".trim.unique.pick.chop.precluster.unique.names",
			"-m Fungi_temp.list",
			"-o "+NAME+"_final.list"
		]))
	
		os.system(' '.join([
			"cp",
			INPUT+".trim.unique.pick.chop.precluster.unique.names",
			NAME+"_final.names"
		]))

		os.system(' '.join([
			"cp",
			INPUT+".trim.unique.pick.chop.precluster.unique.fasta",
			NAME+"_final.fasta"
		]))

		os.system(' '.join([
			"cp",
			INPUT+".pick.chop.groups",
			NAME+"_final.groups"
		]))

		os.system(' '.join([
			"mv",
			"Fungi_temp.list",
			NAME+"_final.an.list"
		]))

		os.system(' '.join([
			"for n in",
			NAME+"_final.fasta; do mothur \"# classify.seqs(fasta=$n,",
			"name="+NAME+"_final.names,",
			"group="+NAME+"_final.groups,",
			"template=~/Phylogenetic_Gene_Databases/UNITE_ITS/UNITEv6_sh_dynamic.fasta,",
			"taxonomy=~/Phylogenetic_Gene_Databases/UNITE_ITS/UNITEv6_sh_dynamic.tax, cutoff=50,",
			"processors="+PROCESSORS+")\"; done"
		]))

		os.system(' '.join([
			"mv",
			NAME+"_final.UNITEv6_sh_dynamic.wang.taxonomy",
			"./"+OUTPUT+"/"+NAME+"_final.taxonomy"
		]))

		os.system(' '.join([
			"mv",
			NAME+"_final.*",
			"./"+OUTPUT+"/"
		]))

stop = timeit.default_timer()

print stop - start

