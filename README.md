# usher_analyse_tool
Tools for analysing usher results. 

``des.json`` is lineage designations and is 
downloaded from https://nextstrain.org/charon/getDataset?prefix=staging/nextclade/sars-cov-2/ .
It is used to tune some of the ill-defined branches on usher tree(usually caused by deletions) to their correct form.

``reference_seq.txt`` is SARS-Cov-2 reference seq.

To analyse, run 

``python analyse.py --usher usher.json --important-threshold 2 ``

``usher.json`` is the json file downloaded from usher (upload fasta to https://genome.ucsc.edu/cgi-bin/hgPhyloPlace to get it)

``important-threshold`` is the threshold of uploaded seqs for important branches to display.

The program will automatically filter out all seqs of the 2 forms:
1：singlet seqs with >5 reversions compared with a designated lineage
2：undesignated lineages with (spike/Orf9b/new stop codon/start codon removed).


Packages required:
``json``
``argparse``
``copy``

Add highlighted lineage features to highlight samples with too many reversions, branches with additional undesignated Spike or Orf9b mutation, or branches with early stop signal/start codon destroyed. 

Link to view important seqs on the most recent uploads (date replaced to the most recent date and num refers to number of seqs (1:0-1000, 2:1000-2000 etc)): 
Currently updating every 3 days. 

https://nextstrain.org/fetch/raw.githubusercontent.com/aviczhl2/usher_analyse_tool/main/date-num.json?branchLabel=Spike%20mutations&f_userOrOld=highlighted%20sample

Example for seqs 0-1000 at 7-10(referring to 7-6~7-9 seqs on GISAID)：
https://nextstrain.org/fetch/raw.githubusercontent.com/aviczhl2/usher_analyse_tool/main/7-10-1.json?branchLabel=Spike%20mutations&f_userOrOld=highlighted%20sample






