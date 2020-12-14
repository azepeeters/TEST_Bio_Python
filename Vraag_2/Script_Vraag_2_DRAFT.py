

# At this step you import modules. Modules in Python are like "Apps" that enable you to perform some complex tasks without having to write thousands of line of code
from Bio import SeqIO
from Bio import SeqUtils
from Bio.SeqUtils.ProtParam import ProteinAnalysis


#
file_1 ='protein.fasta'

#
file_2 ='protein_fragments.fasta'


##########################################################
# Treat File 1
##########################################################

#
for seq_record in SeqIO.parse(open(file_1, mode='r'), 'fasta'):


    # Put the sequence in a "Protein Analysis" function that will enable us to analyse that protein sequence afterwards
    X = ProteinAnalysis(str(seq_record.seq))

    # Compute Molecular weight for the protein sequence
    seq_record.molecular_weight = X.molecular_weight()



# Here we print the result -> The result seems far from the mass on the teachers' notes
print("The protein is : " + seq_record.seq)
print("The molecular weight of that protein is " + str(seq_record.molecular_weight))


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

        # Put the sequence in a "Protein Analysis" function that will enable us to analyse that protein sequence afterwards
        X = ProteinAnalysis(sequence)

        # Add the weight of the sequence to the total weight of the combination
        Total_Weight += X.molecular_weight()


    print('Molecular weight of combination : ' + str(combination) + ' - is = ' + str(Total_Weight))



