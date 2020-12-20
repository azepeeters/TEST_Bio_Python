

# At this step you import modules. Modules in Python are like "Apps" that enable you to perform some complex tasks without having to write thousands of line of code
from Bio import SeqIO
from Bio import SeqUtils
from Bio.SeqUtils.ProtParam import ProteinAnalysis


#
file_1 ='short_protein.fasta'

#
file_2 ='short_protein_fragments.fasta'



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


print('VRAAG 1 :')

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


print('VRAAG 2 :')


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
# Generate all possible combinations of the fragments, and compute the weight of each combinations
#   For each combination, assess if its weight matches the total proteins weight
##########################################################


import itertools


# Make an empty list in which to store matching combinations
combinations_with_same_weight = []

# This will make a loop from 1 to X, with X being the total number of different fragments in the 2nd file
for i in range(1, len(list_relevant_fragments)+1):

    # printing which size of combinations is being tested here
    print('Testing combinations of ' + str(i) + ' fragments')

    # This will create a list with all possible combinations of fragments, with a max of i fragments (i is defined here above)
    fragment_combination = [list(x) for x in itertools.combinations(list_relevant_fragments, i)]

    # Compute the weight of each combination
    for combination in fragment_combination:

        # Give 0 as default value for total weight
        Total_Weight = 0

        # Add the weight of each fragment in the combination to the Total_Weight
        for fragment in combination:

            # Add the weight of the sequence to the total weight of the combination
            Total_Weight += function_molecular_weight(str(fragment))

        # Test if the Total_Weight of that combination equals the weight of the complete protein (rounded at 4 decimal digits)
        if round(Total_Weight, 4) == round(protein_molecular_weight, 4) :

            print('This combination''s weight matches the protein weight : ')
            print(combination)
            print('With as weight : ' + str(round(protein_molecular_weight, 4)))

            # For each relevant fragments found, insert it in the list
            combinations_with_same_weight.insert(0, combination)



print('We found ' + str(len(combinations_with_same_weight)) + ' combinations with the same weight than the protein')


##########################################################
# Try to rebuild the protein using the combinations that match based on weight
##########################################################

print('VRAAG 3 :')


#--------------------------------
#   Create function to try to rebuild protein based on fragments
#--------------------------------


def check_protein_match(protein, fragment_list, matching_fragment_list):

    # Remove fragments already matched
    fragment_list = [i for i in fragment_list if i not in matching_fragment_list]

    # Stop if the protein was completely rebuilt
    if len(protein) == 0 :

        print('This combination can be used to rebuild the protein : ' + str(matching_fragment_list))

    # Move on if still any part of the protein to find
    else:

        # Test for all the fragments in the list
        for fragment in fragment_list:

            # Does the protein start with that fragment?
            if protein.startswith(str(fragment)) == True:

                # Add the matching fragment to the list
                matching_fragment_list.insert(len(matching_fragment_list), str(fragment))

                # Remove the fragment found from the protein-characters
                protein = protein[len(str(fragment)):]

                # Apply the same function to the remaining parts
                check_protein_match(protein, fragment_list, matching_fragment_list)


#--------------------------------
#   Test all weight-matching combinations with the above created function
#--------------------------------


for combination in combinations_with_same_weight:

    # Make an empty list to start with
    matching_fragment_list = []

    # Test the combination
    check_protein_match(protein_sequence, combination, matching_fragment_list)
