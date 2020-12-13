
# Here is the module to use for computing the molecular mass
from Bio.SeqUtils.ProtParam import ProteinAnalysis

# Here is the protein sequence to weight
Protein = 'DHPFWKQTACKHV'

# Here we compute the mass, using the ProteinAnalysis module without any parameter
weight = ProteinAnalysis(Protein).molecular_weight()

# Here we print the result -> The result seems far from the mass on the teachers' notes
print ("The molecular weight of protein " + Protein + " is " + str(weight) )

# Here we compute the mass, using the ProteinAnalysis module with the monoisotopic=True parameter
weight = ProteinAnalysis(Protein, monoisotopic=True).molecular_weight()

# Here we print the result -> The result seems far from the mass on the teachers' notes
print ("The molecular weight of protein " + Protein + " is " + str(weight) )

