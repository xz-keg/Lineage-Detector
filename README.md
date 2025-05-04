# Lineage Detector
Tools for analysing usher results. 

``des.json`` is lineage designations and is 
downloaded from https://nextstrain.org/charon/getDataset?prefix=staging/nextclade/sars-cov-2/ .
``curl -o des.json "https://nextstrain.org/charon/getDataset?prefix=staging/nextclade/sars-cov-2/"``
It is used to tune some of the ill-defined branches on usher tree(usually caused by deletions) to their correct form.

``reference_seq.txt`` is SARS-Cov-2 reference seq.

To analyse, run 

``python analyse.py --usher usher.json --important-threshold 2 ``

``usher.json`` is the json file downloaded from usher (upload fasta to https://genome.ucsc.edu/cgi-bin/hgPhyloPlace to get it)
Update on 2024-11-19: Usher now uses the .gz form for those files, please add a .gz to their filename and unzip them before processing. 

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
Currently updating every 3 days. (2024-7-1)

Update on 2024-11-19: Change analysing period to 4 days due to reduced global sequencing.

Update on 2025-1-20: Add Orf3a as one of the default highlight proteins. 

Update on 2025-3-18: Change analysing period to 5 days due to further reduced global sequencing.

Update on 2025-4-20: Change analysing period to 7 days due to further reduced global sequencing.

https://nextstrain.org/fetch/raw.githubusercontent.com/xz-keg/Lineage-Detector/main/date-num.json?branchLabel=Spike%20mutations&f_userOrOld=highlighted%20sample

Example for seqs 0-1000 at 7-10(referring to 7-6~7-9 seqs on GISAID)：
https://nextstrain.org/fetch/raw.githubusercontent.com/xz-keg/Lineage-Detector/main/2024-7-10-1.json?branchLabel=Spike%20mutations&f_userOrOld=highlighted%20sample

Tech Report:
Lineage Detector: Efficient Tool for Detecting New SARS-Cov-2 Lineages
https://www.biorxiv.org/content/10.1101/2024.11.01.621557v1

@article{zou2024lineage,
  title={Lineage Detector: Efficient Tool for Detecting New SARS-Cov-2 Lineages},
  author={Zou, Xu},
  journal={bioRxiv},
  pages={2024--11},
  year={2024},
  publisher={Cold Spring Harbor Laboratory}
}


