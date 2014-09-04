Mothur_Pipeline_ITS
===================

Standard Pipeline for Cleaning and Annotating Fungal ITS region phylogenetic gene libraries

**Provides ALL necessary files to use the R script**: "Mothur_to_PhyloSeq_to_Differential_Abundance_Tools"

(in one case, a specially formatted taxonomy file)

This script is not meant to be a completely hands-off method for processing ITS region data. It IS fully automated, but only in so far as it executes all the pre-set instructions. If your library requires different processing parameters, you will have to manually alter the script. It will do a decent job on your data as is, based on the fact multiple colleagues in my lab have used an nearly identical set of commands on large datasets (500+ libraries).

Ths script is meant to offer a simple archetype. There is a different pipeline offered through QIIME (<http://nbviewer.ipython.org/github/qiime/qiime/blob/1.8.0/examples/ipynb/Fungal-ITS-analysis.ipynb>), which I have previously used successfully. However, I have not compared the two. The mothur GUI would be a cumbersome way of performing this script b/c it is more customized in its use of crunchclust for OTU generation. 

Dependencies:

    - mothur (built with v.1.32.1)

    - UNITE db v6 (formatted for mothur) (script expects it to be located in ~/Phylogenetic_Databases/UNITE_ITS/*)
      (available here: <http://unite.ut.ee/repository.php>)

    - crunchclust (available here: <https://code.google.com/p/crunchclust/downloads/list>)
    
    - Perl script "deunique_mothurlist.pl" (included in this repository)
    
    - Info file post_crunchclust.add (included in this repository)


This script will require mothur AND crunchclust to be in your PATH. 
