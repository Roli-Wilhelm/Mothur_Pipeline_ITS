#!/usr/bin/python
import sys, os, re, getopt, glob, numpy as np, random
import timeit
import itertools

start = timeit.default_timer()

Usage = """
Usage:  ./MOTHUR_ITS_PIPELINE.py -o FUNGI_0.055 -n FUNGI_0.055 -i IIKFCBR02.fasta -f IIKFCBR02.oligos -d 11 -l 200 -p 8

REQUIRED ARGUMENTS:
		-o	output directory
	
                -n      the name given to the project

		EITHER
                -i      the input file (accepted: .fasta, .sff, .shhh.fasta. MUST SPECIFY if the latter two)
			(if providing .fasta, one must provide a .qual file)

		&&

		-f	the .oligos file

                OR JUST:
                -m      specify directory if you would like to combine multiple .sff files from a specific directory.
                        (the directory must also contain all the associated .oligos files for each run)

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
1) PROCESSORS ARE SET TO 6. MANUALLY CHANGE IF DESIRED.
2) Minflows is set to 180 (i.e. the minimum length of sequence will be 180bp)
3) Maxflows is set to 450 (i.e. the maximum length of sequence will be 450bp)
4) If you are merging multiple runs AND if there are overlapping barcodes, the script CAN handle this, HOWEVER
it makes the assumption that the first 9 nucleotides of your primer are not degenerate (i.e. the same for all sequences).
This is b/c the barcodes will be made unique for all of your reads and to ensure short barcode length sequences
within your sequences are not erroneously changed, the script searches for barcode + 9 nucleotides of primer.

Usage:  ./MOTHUR_ITS_PIPELINE.py -o FUNGI_0.055 -n FUNGI_0.055 -i IIKFCBR02.fasta -f IIKFCBR02.oligos -d 11 -l 200 -p 8

or

Usage:  ./MOTHUR_ITS_PIPELINE.py -o FUNGI_0.055 -n FUNGI_0.055 -m ~/FUNGI_SFF/ -d 11 -l 200 -p 8

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
MIX=''
Debug=''
NO_SHHH=''

## Read command line args
myopts, args = getopt.getopt(sys.argv[1:],"n:i:o:S:b:d:R:l:f:C:O:p:m:")

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
    if o == '-m':
        MIX= a
    if o == '-b':
        Debug= a


####
# Function to Find Duplicates in List
####
def list_duplicates(seq):
        seen = set()
        # adds all elements it doesn't know yet to seen and all other to seen_twice
        seen_add = seen.add

        # turn the set into a list (as requested)
        seen_twice = set( x for x in seq if x in seen or seen_add(x) )

        return list(seen_twice)

def id_generator(size, chars):
        return ''.join(random.choice(chars) for _ in range(size))


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
if MIX:
        INPUT = []
        if Debug:
                print "You chose to merge multiple files for input.\n"

        if SFF:
                for FILE in glob.glob("./"+MIX+"/*.sff"):
                        FILE = re.sub(MIX, "", FILE)
                        FILE = re.sub("./", "", FILE)
                        INPUT.append(re.sub(".sff", "", FILE))

                        if Debug:
                                error.write("The base filename of your input file is:"+re.sub(".sff", "", FILE)+"\n")
        elif RERUN:
                for FILE in glob.glob("./"+MIX+"/*.shhh.fasta"):
                        FILE = re.sub(MIX, "", FILE)
                        FILE = re.sub("./", "", FILE)
                        INPUT.append(re.sub(".shhh.fasta", "", FILE))

                        if Debug:
                                error.write("The base filename of your input file is:"+re.sub(".shhh.fasta", "", FILE)+"\n")

	else:
                for FILE in glob.glob("./"+MIX+"/*.fasta"):
                        FILE = re.sub(MIX, "", FILE)
                        FILE = re.sub("./", "", FILE)
                        INPUT.append(re.sub(".fasta", "", FILE))

                        if Debug:
                                error.write("The base filename of your input file is:"+re.sub(".fasta", "", FILE)+"\n")
        if Debug:
                error.write(str(INPUT)+"\n")

else:
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

        OLIGOS = re.sub(".oligos", "", OLIGOS)
        print OLIGOS

#Process and merge output from .sff files
if MIX and SFF:
        for FILE in INPUT:
                if Debug:
                        error.write("You are now processing: "+FILE+"\n\n")
                ## The chemistry CAN be different for various sequencing runs, MOTHUR provides for this based on the "ORDER" of flows are read. This is sometimes necessary to specify.
                if ORDER:
                        os.system(' '.join([
                                "for n in",
                                FILE+".sff; do mothur \"# sffinfo(sff=$n, flow=T)\"; done"
                        ]))

                        if Debug:
                                error.write(' '.join([
                                        "for n in",
                                FILE+".sff; do mothur \"# sffinfo(sff=$n, flow=T)\"; done"
                        ]))

                        if Debug:
                                error.write(' '.join([
                                        "for n in",
                                        FILE+".sff; do mothur \"# sffinfo(sff=$n, flow=T)\"; done"
                                ])+"\n")


                        os.system(' '.join([
                                "for n in",
                                FILE+".flow; do mothur \"# trim.flows(flow=$n,",
                                "oligos="+FILE+".oligos,",
                                "pdiffs=2, bdiffs=1, minflows=180, maxflows=450,",
                                "order="+ORDER+",",
                                "processors="+PROCESSORS+")\"; done"
                        ]))

                        if Debug:
                                error.write(' '.join([
                                        "for n in",
                                        FILE+".flow; do mothur \"# trim.flows(flow=$n,",
                                        "oligos="+FILE+".oligos,"
                                        "pdiffs=2, bdiffs=1, minflows=180, maxflows=450,",
                                        "order="+ORDER+",",
                        ]))

                        if Debug:
                                error.write(' '.join([
                                        "for n in",
                                        FILE+".flow; do mothur \"# trim.flows(flow=$n,",
                                        "oligos="+FILE+".oligos,"
                                        "pdiffs=2, bdiffs=1, minflows=180, maxflows=450,",
                                        "order="+ORDER+",",
                                        "processors="+PROCESSORS+")\"; done"
                                ])+"\n")

                        os.system(' '.join([
                                "for n in",
                                FILE+".flow.files; do mothur \"# shhh.flows(file=$n,",
                                "order="+ORDER+",",
                                "processors="+PROCESSORS+")\"; done"
                        ]))

                        if Debug:
                                error.write(' '.join([
                                        "for n in",
                                        FILE+".flow.files; do mothur \"# shhh.flows(file=$n,",
                                        "order="+ORDER+",",
                                        "processors="+PROCESSORS+")\"; done"
                                ])+"\n")

                else:
                        os.system(' '.join([
                                "for n in",
                                FILE+".sff; do mothur \"# sffinfo(sff=$n, flow=T)\"; done"
                        ]))

                        if Debug:
                                error.write(' '.join([
                                        "for n in",
                                        FILE+".sff; do mothur \"# sffinfo(sff=$n, flow=T)\"; done"
                                ])+"\n")

                        os.system(' '.join([
                                "for n in",
                                FILE+".flow; do mothur \"# trim.flows(flow=$n,",
                                "oligos="+FILE+".oligos,",
                                "pdiffs=2, bdiffs=1, minflows=180, maxflows=450,",
                                "processors="+PROCESSORS+")\"; done"
                        ]))

                        if Debug:
                                error.write(' '.join([
                                        "for n in",
                                        FILE+".flow; do mothur \"# trim.flows(flow=$n,",
                                        "oligos="+FILE+".oligos,"
                                        "pdiffs=2, bdiffs=1, minflows=180, maxflows=450,",
                                        "processors="+PROCESSORS+")\"; done"
                                ])+"\n")

                        os.system(' '.join([
                                "for n in",
                                FILE+".flow.files; do mothur \"# shhh.flows(file=$n,",
                                "processors="+PROCESSORS+")\"; done"
                        ]))

                        if Debug:
                                error.write(' '.join([
                                        "for n in",
                                        FILE+".flow.files; do mothur \"# shhh.flows(file=$n,",
                                        "processors="+PROCESSORS+")\"; done"
                                ])+"\n")

if MIX and not SFF and not RERUN:
        NO_SHHH = "TRUE"
        print "NO_SHHH set to TRUE"

################################################################################
# Deal with the fact that the same barcode sequences may be used in multiple runs
################################################################################

#Strategy is to get a list of duplicates, and run through them by changing the oligo's file and then the .fasta sequence with randomly generated sequences
if MIX:

        #Concatenate OLIGOS Files
        SAMPLE_LOCATION_DICT = {}

        #Do concatenation
        COMBO_NAME = NAME+"_Combined_Libraries"

        if Debug:
                error.write("Your files hae been concatenated with the base name: "+COMBO_NAME+"\n")

        count = 0

	## Write Combined Oligos File and Make Dictionary of Which 
        with open(MIX+'/'+COMBO_NAME+'.oligos', 'w') as outfile:
                for fname in INPUT:
                        if not re.search("Combined_Libraries", fname):
                                with open(MIX+"/"+fname+".oligos") as infile:
                                        if count != 0:
                                                next(infile)

                                        for line in infile:
                                                line = line.strip("\r\n")
                                                #Write to new concatenated oligos file
                                                outfile.write(line+"\n")
                                                line = line.split()

                                                #Store which pyrotag sample came from which oligos file
                                                if np.size(line)>2:
                                                        fname = re.sub("\\./","",fname)
                                                        SAMPLE_LOCATION_DICT[line[2]] = fname
                                        count = count + 1

        ## Check for Duplicates
        DUPLICATE_DICT = {}
        DUPLICATE_LIST = []
        OLIGOS = MIX+'/'+COMBO_NAME

        with open(OLIGOS+".oligos") as infile:
                #Grab the first 9 nucleotides of the primer for future use
                PRIMER = infile.readline()
                PRIMER = PRIMER.strip()
                SPLIT = PRIMER.split()
                PRIMER_8 = SPLIT[1][0:8]

                #Read in all barcodes & sample IDs
                for line in infile:
                        line = line.strip()
                        line = line.split()

                        #Make list of all barcodes
                        DUPLICATE_LIST.append(line[1])

                        #Make dictionary of sample name and barcode
                        if not DUPLICATE_DICT.has_key(line[1]):
                                DUPLICATE_DICT[line[1]] = [line[2]]
                        else:
                                DUPLICATE_DICT[line[1]].append(line[2])

        #Get list of duplicates
        DUPLICATE_BARCODES = list_duplicates(DUPLICATE_LIST)

        #Change sequence barcodes
        mock_oligos = open("TEMP.oligos", "w")
        mock_oligos.write(PRIMER+"\n")

        with open(OLIGOS+".oligos") as infile:
                next(infile)

                for line in infile:
                        DUPLICATE_SEQ = line.strip("\r\n")
                        DUPLICATE_SEQ = DUPLICATE_SEQ.split()
                        DUPLICATE_SEQ = DUPLICATE_SEQ[1]

                        #See if sequence is duplicated
                        if DUPLICATE_SEQ in DUPLICATE_BARCODES:

                                if np.size(DUPLICATE_DICT[DUPLICATE_SEQ]) > 1:
                                        DUPLICATE_ID = DUPLICATE_DICT[DUPLICATE_SEQ][0]

                                        print "Now Substituting Barcodes Found in Sample: "+DUPLICATE_ID+" due to overlap with another sample.\n"

                                        #Only substitute barcode if there is > 1 instance
                                        if re.search(DUPLICATE_ID, line):
                                                BARCODE_LENGTH = len(DUPLICATE_SEQ)
                                                NEW_BARCODE = id_generator(BARCODE_LENGTH, "TCGA")

                                                line = re.sub(DUPLICATE_SEQ, NEW_BARCODE, line)
                                                mock_oligos.write(line)

                                                #Remove element just in case there are more than two duplications (this would mean the duplicate list
                                                #would contain multiples of hte same DUPLICATE_SEQ and you'll cycle through until that list is exhausted
                                                del DUPLICATE_DICT[DUPLICATE_SEQ][0]

                                                if Debug:
                                                        error.write("The Length of NO_SHHH is: "+str(len(NO_SHHH))+"\n")

                                                ## Find Correct NAME file
                                                CORRECT_NAME = SAMPLE_LOCATION_DICT[DUPLICATE_ID]

                                                if len(NO_SHHH) > 1:
                                                        #Replace all instances in the fasta file with sed
                                                        os.system(' '.join([
                                                                "sed",
                                                                "-i",
                                                                "\'s/"+DUPLICATE_SEQ+PRIMER_8+"/"+NEW_BARCODE+PRIMER_8+"/g\'",
                                                                MIX+'/'+CORRECT_NAME+".fasta"
                                                        ]))

                                                        if Debug:
                                                                error.write(' '.join([
                                                                        "sed",
                                                                        "-i",
                                                                        "\'s/"+DUPLICATE_SEQ+PRIMER_8+"/"+NEW_BARCODE+PRIMER_8+"/g\'",
                                                                        MIX+'/'+CORRECT_NAME+".fasta"+"\n"
                                                                ]))

                                                else:
                                                        os.system(' '.join([
                                                                "sed",
                                                                "-i",
                                                                "\'s/"+DUPLICATE_SEQ+PRIMER_8+"/"+NEW_BARCODE+PRIMER_8+"/g\'",
                                                                MIX+'/'+CORRECT_NAME+".shhh.fasta"
                                                        ]))

                                                        if Debug:
                                                                error.write(' '.join([
                                                                        "sed",
                                                                        "-i",
                                                                        "\'s/"+DUPLICATE_SEQ+PRIMER_8+"/"+NEW_BARCODE+PRIMER_8+"/g\'",

                                                                        MIX+'/'+CORRECT_NAME+".shhh.fasta"+"\n"
                                                                ]))

			                        if Debug:
							error.write("Finished Processing: "+DUPLICATE_ID+"\n")
                                        else:
                                                mock_oligos.write(line)
                                else:
                                        mock_oligos.write(line)

                        #Write oligos that have nothing to do with duplicates
                        else:
                                mock_oligos.write(line)

                mock_oligos.close()

                #Re-name old oligos file and new oligos file
                os.system(' '.join([
                        "mv",
                        MIX+'/'+COMBO_NAME+".oligos",
                        MIX+'/'+COMBO_NAME+".original.oligos"
                ]))

                os.system(' '.join([
                        "mv",
                        "TEMP.oligos",
                        MIX+'/'+COMBO_NAME+".oligos"
                ]))

## Concatenate multiple files into one
if MIX and SFF or MIX and RERUN:

        #FASTA
        with open(MIX+'/'+COMBO_NAME+'.shhh.fasta', 'w') as outfile:
            for fname in INPUT:
                if not re.search("Combined_Libraries", fname):
                        with open(fname+".shhh.fasta") as infile:
                            for line in infile:
                                line = line.strip("\r\n")
                                outfile.write(line+"\n")

        #NAMES
        with open(MIX+'/'+COMBO_NAME+'.shhh.names', 'w') as outfile:
            for fname in INPUT:
                if not re.search("Combined_Libraries", fname):
                        with open(fname+".shhh.names") as infile:
                            for line in infile:
                                line = line.strip("\r\n")
                                outfile.write(line+"\n")

        OLIGOS = MIX+'/'+COMBO_NAME
        INPUT = MIX+'/'+COMBO_NAME

        if Debug:
                error.write(str(OLIGOS))

# Do concatenation of other files
elif MIX and not SFF and not RERUN:

        #Do concatenation
        COMBO_NAME = NAME+"_Combined_Libraries"

        #FASTA
        with open(MIX+'/'+COMBO_NAME+'.fasta', 'w') as outfile:
            for fname in INPUT:
                with open(fname+".fasta") as infile:
                    for line in infile:
                        line = line.strip("\r\n")
                        outfile.write(line+"\n")

        #NAMES
        with open(MIX+'/'+COMBO_NAME+'.names', 'w') as outfile:
            for fname in INPUT:
                with open(fname+".names") as infile:
                    for line in infile:
                        line = line.strip("\r\n")
                        outfile.write(line+"\n")

        OLIGOS = MIX+'/'+COMBO_NAME
        INPUT = MIX+'/'+COMBO_NAME

        if Debug:
                error.write(str(OLIGOS))

if not MIX and SFF:
        ## The chemistry CAN be different for various sequencing runs, MOTHUR provides for this based on the "ORDER" of flows are read. This is sometimes necessary to specify.
        if ORDER:
                os.system(' '.join([
                        "for n in",
                        INPUT+".sff; do mothur \"# sffinfo(sff=$n, flow=T)\"; done"
                ]))

                if Debug:
                        error.write(' '.join([
                                "for n in",
                                INPUT+".sff; do mothur \"# sffinfo(sff=$n, flow=T)\"; done\n"
                        ]))

                os.system(' '.join([
                        "for n in",
                        INPUT+".flow; do mothur \"# trim.flows(flow=$n,",
                        "oligos="+OLIGOS+".oligos,"
                        "pdiffs=2, bdiffs=1, minflows=180, maxflows=450,",
                        "order="+ORDER+",",
                        "processors="+PROCESSORS+")\"; done"
                ]))

                if Debug:
                        error.write(' '.join([
                                "for n in",
                                INPUT+".flow; do mothur \"# trim.flows(flow=$n,",
                                "oligos="+OLIGOS+".oligos,"
                                "pdiffs=2, bdiffs=1, minflows=180, maxflows=450,",
                                "order="+ORDER+",",
                                "processors="+PROCESSORS+")\"; done\n"
                        ]))

                os.system(' '.join([
                        "for n in",
                        INPUT+".flow.files; do mothur \"# shhh.flows(file=$n,",
                        "order="+ORDER+",",
                        "processors="+PROCESSORS+")\"; done"
                ]))

                if Debug:
                        error.write(' '.join([
                                "for n in",
                                INPUT+".flow.files; do mothur \"# shhh.flows(file=$n,",
                                "order="+ORDER+",",
                                "processors="+PROCESSORS+")\"; done\n"
                        ]))

        else:
                os.system(' '.join([
                        "for n in",
                        INPUT+".sff; do mothur \"# sffinfo(sff=$n, flow=T)\"; done"
                ]))

                if Debug:
                        error.write(' '.join([
                             "for n in",
                                INPUT+".sff; do mothur \"# sffinfo(sff=$n, flow=T)\"; done\n"
                        ]))

                os.system(' '.join([
                        "for n in",
                        INPUT+".flow; do mothur \"# trim.flows(flow=$n,",
                        "oligos="+OLIGOS+".oligos,"
                        "pdiffs=2, bdiffs=1, minflows=180, maxflows=450,",
                        "processors="+PROCESSORS+")\"; done"
                ]))

                if Debug:
                        error.write(' '.join([
                                "for n in",
                                INPUT+".flow; do mothur \"# trim.flows(flow=$n,",
                                "oligos="+OLIGOS+".oligos,"
                                "pdiffs=2, bdiffs=1, minflows=180, maxflows=450,",
                                "processors="+PROCESSORS+")\"; done\n"
                        ]))

                os.system(' '.join([
                       "for n in",
                        INPUT+".flow.files; do mothur \"# shhh.flows(file=$n,",
                        "processors="+PROCESSORS+")\"; done"
                ]))

                if Debug:
                        error.write(' '.join([
                                "for n in",
                                INPUT+".flow.files; do mothur \"# shhh.flows(file=$n,",
                                "processors="+PROCESSORS+")\"; done\n"

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

	        ## Provide a version of the taxonomy file acceptable for importing into R
	        os.system(' '.join([
        	        "cp",
                	"./"+OUTPUT+"/"+NAME+"_final.taxonomy",
	                "./"+OUTPUT+"/"+NAME+"_final.R.taxonomy"
        	]))

	        os.system(' '.join([
        	        "sed -i 's/\t/;/g'",
                	"./"+OUTPUT+"/"+NAME+"_final.R.taxonomy"
	        ]))

        	os.system(' '.join([
                	"sed -i 's/;$//g'",
	                "./"+OUTPUT+"/"+NAME+"_final.R.taxonomy"
        	]))

		## Move all Final Files to OUTPUT directory
		os.system(' '.join([
			"mv",
			NAME+"_final.*",
			"./"+OUTPUT+"/"
		]))

		## Cat all logfiles in order of creation and move
                os.system(' '.join([
                        "cat",
                        "$(ls -t mothur.*)",
                        ">",
			"./"+OUTPUT+"/"+NAME+".mothur.logfiles"
                ]))

	else:
		## Start from the earliest point in the MOTHUR pipeline post-shhh.flows()
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
			"maxhomop=8, bdiffs=1, pdiffs=2, minlength=180,",
			"processors="+PROCESSORS+")\"; done"
		]))

		os.system(' '.join([
			"for n in",
			INPUT+".shhh.trim.fasta; do mothur \"# unique.seqs(fasta=$n,",
			"name="+INPUT+".shhh.trim.names)\"; done"
		]))

		os.system(' '.join([
			"for n in",
			INPUT+".shhh.trim.unique.fasta; do mothur \"# summary.seqs(fasta=$n,",
			"name="+INPUT+".shhh.trim.unique.names,",
			"processors="+PROCESSORS+")\"; done"
		]))

		os.system(' '.join([
			"for n in",
			INPUT+".shhh.trim.unique.fasta; do mothur \"# chimera.uchime(fasta=$n,",
			"name="+INPUT+".shhh.trim.unique.names,",
			"group="+INPUT+".shhh.groups,",
			"processors="+PROCESSORS+")\"; done"
		]))

		os.system(' '.join([
			"for n in",
			INPUT+".shhh.trim.unique.uchime.accnos; do mothur \"# remove.seqs(accnos=$n,",
			"fasta="+INPUT+".shhh.trim.unique.fasta,",
			"name="+INPUT+".shhh.trim.unique.names,",
			"group="+INPUT+".shhh.groups)\"; done"
		]))

		## Chop Sequences to Specified Length
		os.system(' '.join([
			"for n in",
			INPUT+".shhh.trim.unique.pick.fasta; do mothur \"# chop.seqs(fasta=$n,",
			"name="+INPUT+".shhh.trim.unique.pick.names,",
			"group="+INPUT+".shhh.pick.groups,",
			"numbases="+LENGTH+", keep=front,",
			"processors="+PROCESSORS+")\"; done"
		]))
	
		os.system(' '.join([
			"for n in",
			INPUT+".shhh.trim.unique.pick.chop.fasta; do mothur \"# summary.seqs(fasta=$n,",
			"name="+INPUT+".shhh.trim.unique.pick.chop.names,",
			"processors="+PROCESSORS+")\"; done"
		]))

		## Reduce the number of sequences prior to Crunch Clust
		os.system(' '.join([
			"for n in",
			INPUT+".shhh.trim.unique.pick.chop.fasta; do mothur \"# pre.cluster(fasta=$n,",
			"name="+INPUT+".shhh.trim.unique.pick.chop.names,",
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

	        ## Provide a version of the taxonomy file acceptable for importing into R
	        os.system(' '.join([
        	        "cp",
                	"./"+OUTPUT+"/"+NAME+"_final.taxonomy",
	                "./"+OUTPUT+"/"+NAME+"_final.R.taxonomy"
        	]))

	        os.system(' '.join([
        	        "sed -i 's/\t/;/g'",
                	"./"+OUTPUT+"/"+NAME+"_final.R.taxonomy"
	        ]))

        	os.system(' '.join([
                	"sed -i 's/;$//g'",
	                "./"+OUTPUT+"/"+NAME+"_final.R.taxonomy"
        	]))

		os.system(' '.join([
			"mv",
			NAME+"_final.*",
			"./"+OUTPUT+"/"
		]))

		## Cat all logfiles in order of creation and move
                os.system(' '.join([
                        "cat",
                        "$(ls -t mothur.*)",
                        ">",
			"./"+OUTPUT+"/"+NAME+".mothur.logfiles"
                ]))

## If user wants to provide only a .fasta and .qual file for analysis.
else:		
	## Same as above, but with the .shhh removed from all commands
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

	        ## Provide a version of the taxonomy file acceptable for importing into R
	        os.system(' '.join([
        	        "cp",
                	"./"+OUTPUT+"/"+NAME+"_final.taxonomy",
	                "./"+OUTPUT+"/"+NAME+"_final.R.taxonomy"
        	]))

	        os.system(' '.join([
        	        "sed -i 's/\t/;/g'",
                	"./"+OUTPUT+"/"+NAME+"_final.R.taxonomy"
	        ]))

        	os.system(' '.join([
                	"sed -i 's/;$//g'",
	                "./"+OUTPUT+"/"+NAME+"_final.R.taxonomy"
        	]))

		os.system(' '.join([
			"mv",
			NAME+"_final.*",
			"./"+OUTPUT+"/"
		]))

		## Cat all logfiles in order of creation and move
                os.system(' '.join([
                        "cat",
                        "$(ls -t mothur.*)",
                        ">",
			"./"+OUTPUT+"/"+NAME+".mothur.logfiles"
                ]))

	else:
		## Same as .shhh pipeline (except .shhh removed from commmand) EXCEPT trim.seqs(is done according to a moving window)		
		## Make Record of Starting Stastics of the Library
		os.system(' '.join([
			"for n in",
			INPUT+".fasta; do mothur \"# summary.seqs(fasta=$n,",
			"name="+INPUT+".names,",
			"processors="+PROCESSORS+")\"; done"
		]))

		## More QC
		os.system(' '.join([
			"for n in",
			INPUT+".fasta; do mothur \"# trim.seqs(fasta=$n,",
			"oligos="+OLIGOS+".oligos,",
			"qfile="+INPUT+".qual,",
			"maxhomop=8, bdiffs=1, pdiffs=2, qwindowaverage=35, qwindowsize=50, minlength=180,",
			"processors="+PROCESSORS+")\"; done"
		]))

		os.system(' '.join([
			"for n in",
			INPUT+".trim.fasta; do mothur \"# unique.seqs(fasta=$n,",
			"name="+INPUT+".trim.names)\"; done"
		]))

		os.system(' '.join([
			"for n in",
			INPUT+".trim.unique.fasta; do mothur \"# summary.seqs(fasta=$n,",
			"name="+INPUT+".trim.unique.names,",
			"processors="+PROCESSORS+")\"; done"
		]))

		os.system(' '.join([
			"for n in",
			INPUT+".trim.unique.fasta; do mothur \"# chimera.uchime(fasta=$n,",
			"name="+INPUT+".trim.unique.names,",
			"group="+INPUT+".groups,",
			"processors="+PROCESSORS+")\"; done"
		]))

		os.system(' '.join([
			"for n in",
			INPUT+".trim.unique.uchime.accnos; do mothur \"# remove.seqs(accnos=$n,",
			"fasta="+INPUT+".trim.unique.fasta,",
			"name="+INPUT+".trim.unique.names,",
			"group="+INPUT+".groups)\"; done"
		]))

		## Chop Sequences to Specified Length
		os.system(' '.join([
			"for n in",
			INPUT+".trim.unique.pick.fasta; do mothur \"# chop.seqs(fasta=$n,",
			"name="+INPUT+".trim.unique.pick.names,",
			"group="+INPUT+".pick.groups,",
			"numbases="+LENGTH+", keep=front,",
			"processors="+PROCESSORS+")\"; done"
		]))
	
		os.system(' '.join([
			"for n in",
			INPUT+".trim.unique.pick.chop.fasta; do mothur \"# summary.seqs(fasta=$n,",
			"name="+INPUT+".trim.unique.pick.chop.names,",
			"processors="+PROCESSORS+")\"; done"
		]))

		## Reduce the number of sequences prior to Crunch Clust
		os.system(' '.join([
			"for n in",
			INPUT+".trim.unique.pick.chop.fasta; do mothur \"# pre.cluster(fasta=$n,",
			"name="+INPUT+".trim.unique.pick.chop.names,",
			"group="+INPUT+".pick.chop.groups, diffs=2,",
			"processors="+PROCESSORS+")\"; done"
		]))

		os.system(' '.join([
			"for n in",
			INPUT+".trim.unique.pick.chop.precluster.fasta; do mothur \"# unique.seqs(fasta=$n,",
			"name="+INPUT+".trim.unique.pick.chop.precluster.names)\"; done"
		]))

		## Perform Crunch Clust as above
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
			"-n "+INPUT+".trim.unqiue.pick.chop.precluster.unique.names",
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

	        ## Provide a version of the taxonomy file acceptable for importing into R
	        os.system(' '.join([
        	        "cp",
                	"./"+OUTPUT+"/"+NAME+"_final.taxonomy",
	                "./"+OUTPUT+"/"+NAME+"_final.R.taxonomy"
        	]))

	        os.system(' '.join([
        	        "sed -i 's/\t/;/g'",
                	"./"+OUTPUT+"/"+NAME+"_final.R.taxonomy"
	        ]))

        	os.system(' '.join([
                	"sed -i 's/;$//g'",
	                "./"+OUTPUT+"/"+NAME+"_final.R.taxonomy"
        	]))

		os.system(' '.join([
			"mv",
			NAME+"_final.*",
			"./"+OUTPUT+"/"
		]))

		## Cat all logfiles in order of creation and move
                os.system(' '.join([
                        "cat",
                        "$(ls -t mothur.*)",
                        ">",
			"./"+OUTPUT+"/"+NAME+".mothur.logfiles"
                ]))

stop = timeit.default_timer()

print stop - start

