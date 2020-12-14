

# At this step you import modules. Modules in Python are like "Apps" that enable you to perform some complex tasks without having to write thousands of line of code
from Bio import SeqIO
from Bio import SeqUtils
from Bio.SeqUtils.ProtParam import ProteinAnalysis


#
file_1 ='protein.fasta'

#
file_2 ='protein_fragments.fasta'



##########################################################
# Create a function to compute molecular weight
##########################################################


def function_molecular_weight(sequence_text):
    # Here is a Python dictionary in which I stored the molecular weights listed in the above URL
    weights = {
        'A': 71.0788,
        'R': 156.1875,
        'N': 114.1038,
        'D': 115.0886,
        'C': 103.1388,
        'E': 129.1155,
        'Q': 128.1307,
        'G': 57.0519,
        'H': 137.1411,
        'I': 113.1594,
        'L': 113.1594,
        'K': 128.1741,
        'M': 131.1926,
        'F': 147.1766,
        'P': 97.1167,
        'S': 87.0782,
        'T': 101.1051,
        'W': 186.2132,
        'Y': 163.1760,
        'V': 99.1326
    }

    # Here we compute the mass, using the weights of the dictionary
    weight = sum(weights[p] for p in sequence_text)

    return weight

##########################################################
# Treat File 1
##########################################################

#
for seq_record in SeqIO.parse(open(file_1, mode='r'), 'fasta'):


    # Compute Molecular weight for the protein sequence
    protein_molecular_weight = function_molecular_weight(str(seq_record.seq))


# Here we print the result -> The result seems far from the mass on the teachers' notes
print("The protein is : " + str(seq_record.seq))
print("The molecular weight of that protein is " + str(protein_molecular_weight))


##########################################################
# Treat File 2 : Insert sequence parts in a List
##########################################################


# Define an empty list in which to insert Sequence parts
list_sequence_parts = []

# Loop through all sequences in the 2nd input file to count the number of sequences in it
for seq_record in SeqIO.parse(open(file_2, mode='r'), 'fasta'):

    # For each sequence found, insert it in the lsit
    list_sequence_parts.insert(0, str(seq_record.seq))

print('Total of sequence parts that could be used to re-build it  : ' + str(len(list_sequence_parts)))

print('List of sequences :')
print(list_sequence_parts)



##########################################################
# Build a list containing all possible combinations of sequences
##########################################################


import itertools

# Create an empty list in which all combinations will be stored
tested_sequences_combinations_FULL = []

# This will fill the list with all possible combinations of the Proteins
for i in range(1, len(list_sequence_parts)+1):

    els = [list(x) for x in itertools.combinations(list_sequence_parts, i)]
    tested_sequences_combinations_FULL.extend(els)



print('All Combinations list : ')
print(tested_sequences_combinations_FULL)



##########################################################
# Compute the molecular weight of each combination
##########################################################

# For each combiantion
for combination in tested_sequences_combinations_FULL :

    # Give 0 as default value for total weight
    Total_Weight = 0

    # For each sequence inside the combination
    for sequence in combination :

        # Add the weight of the sequence to the total weight of the combination
        Total_Weight += function_molecular_weight(str(sequence))


    print('Molecular weight of combination : ' + str(combination) + ' - is = ' + str(Total_Weight))



