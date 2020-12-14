

# At this step you import modules. Modules in Python are like "Apps" that enable you to perform some complex tasks without having to write thousands of line of code
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.SeqUtils.ProtParam import ProteinAnalysis

# We here define the name of the input file (= file that we will treat)
file_in ='protein.fasta'


# We here open the input file and for each protein sequence in the input file, we will perform the tasks in the belows "indentation block"
for seq_record in SeqIO.parse(open(file_in, mode='r'), 'fasta'):

    # Put the sequence in a "Protein Analysis" function that will enable us to analyse that protein sequence afterwards
    X = ProteinAnalysis(str(seq_record.seq))

    # Compute Molecular weight for the protein sequence
    seq_record.molecular_weight = X.molecular_weight()


    # Here we print the result -> The result seems far from the mass on the teachers' notes
    print("The protein is : " + seq_record.seq )
    print("The molecular weight of that protein is " + str(seq_record.molecular_weight))

