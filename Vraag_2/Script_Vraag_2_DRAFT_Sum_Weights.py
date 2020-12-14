

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

    return round(weight, 4)

##########################################################
# Treat File 1
##########################################################

#
for seq_record in SeqIO.parse(open(file_1, mode='r'), 'fasta'):


    protein_sequence = str(seq_record.seq)

    # Compute Molecular weight for the protein sequence
    protein_molecular_weight = function_molecular_weight(protein_sequence)


# Here we print the result -> The result seems far from the mass on the teachers' notes
print("The protein is : " + protein_sequence)
print("The molecular weight of that protein is " + str(protein_molecular_weight))


##########################################################
# Treat File 2 : Insert all fragments in a List (only if they are relevant for the protein)
##########################################################


# Define an empty list in which to insert relevant fragments
list_relevant_fragments = []

# Loop through all fragments found in the 2nd input file
for seq_record in SeqIO.parse(open(file_2, mode='r'), 'fasta'):

    # Check if fragment is relevant by checking if it is found somewhere in the protein
    if protein_sequence.find(str(seq_record.seq)) >= 0:

        # For each relevant fragments found, insert it in the list
        list_relevant_fragments.insert(0, str(seq_record.seq))


print('Total of unique fragments that could be used to re-build the protein  : ' + str(len(list_relevant_fragments)))

print('The list of these fragments is :')
print(list_relevant_fragments)


##########################################################
# Convert list iof fragments to list of weights
##########################################################


# Define an empty list in which to insert relevant fragments
list_relevant_fragments_weights = []

for fragment in list_relevant_fragments :

    # For each fragments, compute the weight and insert it in the list
    list_relevant_fragments_weights.insert(0, function_molecular_weight(str(fragment)))


print('The weight of these fragments is :')
print(list_relevant_fragments_weights)



##########################################################
# Generate all possible combinations of the fragments, and compute the weight of each combinations
#   For each combination, assess if its weight matches the total proteins weight
##########################################################


import itertools


# This will make a loop from 1 to X, with X being the total number of different fragments in the 2nd file
for i in range(1, len(list_relevant_fragments_weights)+1):

    # printing which size of combinations is being tested here
    print('Testing combinations of ' + str(i) + ' fragments')

    # This will create a list with all possible combinations of fragments, with a max of i fragments (i is defined here above)
    fragment_combination_weight = [list(x) for x in itertools.combinations(list_relevant_fragments_weights, i)]

    # Compute the weight of each combination
    for combination in fragment_combination_weight:

        #print('Combination = ')
        #print(combination)

        # Give 0 as default value for total weight
        Total_Weight = 0

        # Add the weight of each fragment in the combination to the Total_Weight
        for fragment in combination:

            #print('Fragment = ')
            #print(fragment)

            # Add the weight of the sequence to the total weight of the combination
            Total_Weight += fragment

        # Test if the Total_Weight of that combination equals the weight of the complete protein (rounded at 4 decimal digits)
        if round(Total_Weight, 4) == round(protein_molecular_weight, 4) :

            print('This combination''s weight matches the protein weight : ')
            print(combination)
            print('With as weight : ' + str(round(protein_molecular_weight, 4)))

