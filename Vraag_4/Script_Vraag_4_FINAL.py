

from Bio import SeqIO




print('VRAAG 4 :')


#
file_1 ='short_protein.fasta'


##########################################################
# Create a table listing all possible tRNA sequences per protein
##########################################################


table={
        "F":["UUU", "UUC"],
        "L":["UUA", "UUG", "CUU", "CUC", "CUA", "CUG"],
        "S":["UCU", "UCC", "UCA", "UCG", "AGU", "AGC"],
        "Y":["UAU", "UAC"],
        "C":["UGU", "UGC"],
        "W":["UGG"],
        "P":["CCU", "CCC", "CCA", "CCG"],
        "H":["CAU", "CAC"],
        "Q":["CAA", "CAG"],
        "R":["CGU", "CGC", "CGA", "CGG", "AGA", "AGG"],
        "I":["AUU", "AUC", "AUA"],
        "M":["AUG"],
        "T":["ACU", "ACC", "ACA", "ACG"],
        "N":["AAU", "AAC"],
        "K":["AAA", "AAG"],
        "V":["GUU", "GUC", "GUA", "GUG"],
        "A":["GCU", "GCC", "GCA", "GCG"],
        "D":["GAU", "GAC"],
        "E":["GAA", "GAG"],
        "G":["GGU", "GGC", "GGA", "GGG"]
    }


##########################################################
# For each PROTEIN in the input file : list all tRNA combinations (only if less then 5000 combinations)
##########################################################


#
for seq_record in SeqIO.parse(open(file_1, mode='r'), 'fasta'):


    protein_sequence = str(seq_record.seq)

    # Give start value of 1
    possible_combinations = 1

    # Compute the total numer of possible combinations of tRNA for the protein sequence
    for p in protein_sequence :

        # Multiply the value of each loop with the output value of the previous loop
        possible_combinations *= len(table[p])


    print("The protein sequence is : " + protein_sequence)
    print("The number of possible tRNA combinations is : " + str(possible_combinations))


    ##########################################################
    # Write error message if more then 5000 combinations
    ##########################################################


    #
    if possible_combinations > 5000 :

        print('More then 5000 tRNA combinations : we will not go further')

    else :

        print('Protein sequence with at most 5000 tRNA combinations : we WILL go further to list these tRNA combinations')

        ##########################################################
        # List all possible tRNA combinations
        ##########################################################

        # Dit module gebruiken om alle mogelijke combinaties te genereren
        from itertools import product

        # Combute the tRNA-list version of the protein sequence
        res = [table[p] for p in protein_sequence]

        print('The "tRNA-list version" of the protein sequence is : ')
        print(res)

        print('The list of tRNA combinations for the protein sequence is : ')
        print(list(product(*res)))



