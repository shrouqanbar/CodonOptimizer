import re
import urllib.request
from decimal import Decimal, getcontext
import plotly.graph_objects as go
#from prettytable import PrettyTable
import argparse
import sys


#Parse Fasta File:
def parse_fasta(file):
    # The first line in a FASTA record is the title line.
    first_line = file.readline()
    # Double-check that it's a FASTA file by looking for the '>'.
    if not first_line.startswith(">"):
        raise TypeError("Not a FASTA file: %r" % first_line)
    # Read the sequence lines
    sequence_lines = []
    while 1:
        line = file.readline().rstrip()
        # Reached the end of the record or end of the file
        if line == "":
            break
        sequence_lines.append(line)
    sequence ="".join(sequence_lines)
    return sequence

#Translate original sequence to PROTEIN
def DNA_translation(sequence):
    table = {
        'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
        'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
        'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
        'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
        'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
        'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
        'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
        'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
        'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
        'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
        'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
        'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
        'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
        'TAC': 'Y', 'TAT': 'Y', 'TAA': '_', 'TAG': '_',
        'TGC': 'C', 'TGT': 'C', 'TGA': '_', 'TGG': 'W',
    }
    protein = ''
    if len(sequence) % 3 == 0:
        for i in range(0, len(sequence), 3):
            codon = sequence[i:i + 3]
            protein += table[codon]
    return protein

#Transform a CSV string of a codon table to a dict.
def csv_data_to_codons_dict(csv_string):
    result = {}
    for line in csv_string.split("\n")[1:]:
        aa, codon, freq = line.split(',')
        if aa not in result:  # new amino acid if not existed
            result[aa] = {}
        result[aa][codon] = float(freq)
    return result

#download table & convert it to a CSV file
urlopen = urllib.request.urlopen
def download_codons_table(fileURL):
    codonRegexpr = r"([ATGCU]{3}) ([A-Z]|\*) (\d.\d+)"
    html_content = urlopen(fileURL).read().decode().replace("\n", " ")  # all the page in HTML
    # print(html_content)
    csv_data = "\n".join(["amino_acid,codon,frequency"] + sorted([
        "%s,%s,%s" % (aa, codon, freq)
        for codon, aa, freq in re.findall(codonRegexpr, html_content)

    ]))

    return csv_data_to_codons_dict(csv_data)

#To replace U by T
def U_replaced_by_T(table):
    return {
        aa: {
            codon.replace('U', 'T'): freq
            for codon, freq in aa_data.items()
        }
        for aa, aa_data in table.items()
    }
#new Dictionary wirh heighest frequencies
def UpdatedDic(table):
    dic = {}
    for d_aa, d_codon in table.items():
        s = max(d_codon.values())
        for (key, value) in d_codon.items():
            if value == s:
                dic.update({d_aa: key})
    return dic

def translateOptimizedSeq(protein,table):
    OPSeq = ""
    for i in range(0, len(protein)):
        amino = protein[i]
        OPSeq += UpdatedDic(table)[amino]
    return OPSeq



#Print changes Table in console
'''def changesTable(sequence,OPSeq,table):
    listOfOldCodons = []
    listOfOPcodons = []
    listOfChangedAA = []  # changed amino acids
    fractionsDic = {}  # fractions
    changedSeqDic = {}  # changed seq codons
    changedOpDic = {}  # changed optimized codons
    listOfIndexes = []  # index of optimized codons
    counter = 0
    for i, j in zip(range(0, len(sequence), 3), range(0, len(OPSeq), 3)):
        originalCodon = sequence[i:i + 3]
        Optimizedcodon = OPSeq[j:j + 3]
        if originalCodon != Optimizedcodon:
            listOfIndexes.append(i + 1)
            counter += 1
            listOfOldCodons.append(originalCodon)
            listOfOPcodons.append(Optimizedcodon)
            for aa, aa_data in table.items():
                for k, v in aa_data.items():
                    if k == originalCodon:
                        listOfChangedAA.append(aa)
                        fractionsDic.update({aa: aa_data})
                        changedSeqDic.update({originalCodon: v})
                    elif k == Optimizedcodon:
                        changedOpDic.update({Optimizedcodon: v})
    freqDiff = []
    getcontext().prec = 2
    for i, j in zip(range(0, len(listOfOldCodons)), range(0, len(listOfOPcodons))):
        freqDiff.append(Decimal(changedOpDic[listOfOPcodons[j]]) - Decimal(changedSeqDic[listOfOldCodons[i]]))
    x = PrettyTable()
    column_names = ["#", "POS", "Old", "New","diff","amino","Fractions"]
    x.add_column(column_names[0], [i for i in range(1, counter+1)])
    x.add_column(column_names[1], [listOfIndexes [i] for i in range(0, len(listOfIndexes ))])
    x.add_column(column_names[2], [listOfOldCodons[i] for i in range(0, len(listOfOldCodons))])
    x.add_column(column_names[3], [listOfOPcodons[i] for i in range(0, len(listOfOPcodons))])
    x.add_column(column_names[4], [freqDiff[i] for i in range(0, len(freqDiff))])
    x.add_column(column_names[5], [listOfChangedAA[i] for i in range(0, len(listOfChangedAA))])
    x.add_column(column_names[6], [str(v) for i in listOfChangedAA for k, v in updatedTable.items() if i == k])
    print(x)
#html output
#print(x.get_html_string())'''

#Another way to draw Changes Table in details
def changesTable(sequence,OPSeq,table):
    headerColor = 'white'
    listOfOldCodons = []
    listOfOPcodons = []
    listOfChangedAA = []  # changed amino acids
    fractionsDic = {}  # fractions
    changedSeqDic = {}  # changed seq codons
    changedOpDic = {}  # changed optimized codons
    listOfIndexes = []  # index of optimized codons
    counter = 0
    for i, j in zip(range(0, len(sequence), 3), range(0, len(OPSeq), 3)):
        originalCodon = sequence[i:i + 3]
        Optimizedcodon = OPSeq[j:j + 3]
        if originalCodon != Optimizedcodon:
            listOfIndexes.append(i + 1)
            counter += 1
            listOfOldCodons.append(originalCodon)
            listOfOPcodons.append(Optimizedcodon)
            for aa, aa_data in table.items():
                for k, v in aa_data.items():
                    if k == originalCodon:
                        listOfChangedAA.append(aa)
                        fractionsDic.update({aa: aa_data})
                        changedSeqDic.update({originalCodon: v})
                    elif k == Optimizedcodon:
                        changedOpDic.update({Optimizedcodon: v})
    freqDiff = []
    getcontext().prec = 2
    for i, j in zip(range(0, len(listOfOldCodons)), range(0, len(listOfOPcodons))):
        freqDiff.append(Decimal(changedOpDic[listOfOPcodons[j]]) - Decimal(changedSeqDic[listOfOldCodons[i]]))
    fig = go.Figure(data=[go.Table(
        header=dict(
            values=['<b>#</b>', '<b>POS</b>', '<b>Old</b>', '<b>New</b>', '<b>diff</b>', '<b>amino</b>',
                    '<b>Fractions</b>'],
            line_color='lightgrey',
            fill_color=headerColor,
            align=['center', 'center'],
            font=dict(color='black', size=11)
        ),
        cells=dict(
            values=[
                [i for i in range(1, counter + 1)],
                [listOfIndexes [i] for i in range(0, len(listOfIndexes ))],
                [listOfOldCodons[i] for i in range(0, len(listOfOldCodons))],
                [listOfOPcodons[i] for i in range(0, len(listOfOPcodons))],
                [freqDiff[i] for i in range(0, len(freqDiff))],
                [listOfChangedAA[i] for i in range(0, len(listOfChangedAA))],
                [str(v) for i in listOfChangedAA for k, v in table.items() if i == k]],
            line_color='lightgrey',
            # 2-D list of colors for alternating rows
            # fill_color = [[rowOddColor,rowEvenColor,rowOddColor, rowEvenColor,rowOddColor]*5],
            align=['center', 'center'],
            font=dict(color='darkslategray', size=11)
        ))
    ])
    return fig.show()

def Main():
    # to enter a list "nargs='+'"
    parser = argparse.ArgumentParser(description='Codon Optimizing Process')
    parser.add_argument("fasta_file_path", nargs='?', help="Enter Fasta file path... ", type=argparse.FileType('r'),
                        default=sys.stdin)
    parser.add_argument("usage_bias_url", help="Paste Codon Usage Bias Table (Standard Format)", type=str)
    parser.add_argument("-o", "--output", help="Output optimized sequence to a file. ", action="store_true")
    args = parser.parse_args()
    sequence = parse_fasta(args.fasta_file_path)
    proteinSeq = DNA_translation(sequence)
    print("Sequences before and after optimization: \n", "Protein Sequence: \n", proteinSeq)
    print("Original Sequence [FASTA]: \n", sequence)
    codonsTable = download_codons_table(args.usage_bias_url)
    updatedTable = U_replaced_by_T(codonsTable)
    OptimizedCodon = translateOptimizedSeq(proteinSeq,updatedTable)
    print("Optimized Sequence [FASTA]: \n", OptimizedCodon)
    if args.output:
        f = open("optimizedSequence.txt", "a")
        f.write(OptimizedCodon)
    #print("Changes Table in details: \n")
    changesTable(sequence, OptimizedCodon,updatedTable)
if __name__ == '__main__':
    Main()
