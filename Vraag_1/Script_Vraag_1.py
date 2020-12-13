

# At this step you import modules. Modules in Python are like "Apps" that enable you to perform some complex tasks without having to write thousands of line of code
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.SeqUtils.ProtParam import ProteinAnalysis

# We here define the name of the input file (= file that we will treat)
file_in ='Fasta_test.fasta'

# Here we define of the output file (= file that we will here create, containing the final results)
file_out='Gene_seq_out_20201208.fasta'


# We here open the output file (to be able to write data inot it later in this script)
with open(file_out, 'w') as f_out:

    # We here open the input file and for each protein sequence in the input file, we will perform the tasks in the belows "indentation block"
    for seq_record in SeqIO.parse(open(file_in, mode='r'), 'fasta'):

        # Put the sequence in a "Protein Analysis" function that will enable us to analyse that protein sequence afterwards
        X = ProteinAnalysis(str(seq_record.seq))

        # Compute Molecular weight for the protein sequence
        seq_record.molecular_weight = X.molecular_weight()

        # Add the molecular weight to the description
        seq_record.description= str(seq_record.description) + ' - Molec weight = ' + str(seq_record.molecular_weight)

        # Print (= show on screen) some information found about that
        print('SequenceID = '  + seq_record.id)
        print('Molecular Weight = ' + str(seq_record.molecular_weight))
        print('New Description = ' + seq_record.description)
        print('Sequence = ' + seq_record.seq)

        # Print a line with a blank value to print a space between the protein sequences
        print(' ')

        # Save the above treated data in a record
        record_to_write = SeqRecord(seq_record.seq, id= seq_record.id, description=seq_record.description)

        # Write the record to output FASTA file
        r=SeqIO.write(record_to_write, f_out, 'fasta')

