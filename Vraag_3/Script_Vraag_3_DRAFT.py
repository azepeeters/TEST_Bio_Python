

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



def check_protein_match(protein, fragment_list, matching_fragment_list):

    # Remove fragments already matched
    fragment_list = [i for i in fragment_list if i not in matching_fragment_list]

    # Continue if still any parts of the protein to find
    if len(protein) == 0 :

        print('Done looping, a possible combination is : ')
        print(matching_fragment_list)

    # Move on if still any part of the protein to find
    else:

        # Test for all the fragments in the list
        for fragment in fragment_list:

            # Does the protein start with that fragment?
            if protein.startswith(str(fragment)) == True:

                #print('Matching fragment : ' + str(fragment))

                # Add the matching fragment to the list
                matching_fragment_list.insert(len(matching_fragment_list), str(fragment))

                #print('Its length is : ' + str(len(str(fragment))))

                # Remove the fragment found from the protein list
                protein = protein[len(str(fragment)):]

                #print('Remaining protein = ' + str(protein))

                # Remove the fragment from the list of fragments to be tested
                #fragment_list.remove(fragment)

                check_protein_match(protein, fragment_list, matching_fragment_list)

                #if len(protein) > 0 and len(fragment_list) > 0 :
                    #check_protein_match(protein, fragment_list, matching_fragment_list)

                #else :
                    #print('Done looping, a possible combination is : ')
                    #print(matching_fragment_list)


matching_fragment_list = []
test_protein_sequence = protein_sequence
fragment_list = list_relevant_fragments

check_protein_match(test_protein_sequence, fragment_list, matching_fragment_list)

