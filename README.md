# CodonOptimizer 
Repository that contains simple python script for DNA analysis, specifically:
* DNA translation
* DNA codon optimizing
* Detailed table with more INFO about the process of optimizing codon 
## Functions:
* DNA_translation - returns Protein sequence
* download_codons_table- returns codons table in a CSV format representing amino_acid,codon, relative_frequency  
* csv_data_to_codons_dict- returns a dictionary to represent codons table  
* UpdatedDic- returns a dictionary containing only highest codon frequencies for each Amino Acid 
* translateOptimizedSeq - returns the new or the optimized sequence based on high frequency feature
* changesTable - returns a table figure to show more details about optimizing process
## Development
* Functions are being adapted to read FASTA format.
* We usually use codons tables from http://www.kazusa.or.jp/
## Usage
* Just type Fasta file name and codon table URL 
* attached an example for fasta files called sequence.fasta 


