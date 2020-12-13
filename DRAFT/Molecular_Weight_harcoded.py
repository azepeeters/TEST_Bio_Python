

# Source of amino acid masses (average))
# http://education.expasy.org/student_projects/isotopident/htdocs/aa-list.html

# Here is a Python dictionary in which I stored the molecular weights listed in the above URL
weights = {
'A':71.0788,
'R':156.1875,
'N':114.1038,
'D':115.0886,
'C':103.1388,
'E':129.1155,
'Q':128.1307,
'G':57.0519,
'H':137.1411,
'I':113.1594,
'L':113.1594,
'K':128.1741,
'M':131.1926,
'F':147.1766,
'P':97.1167,
'S':87.0782,
'T':101.1051,
'W':186.2132,
'Y':163.1760,
'V':99.1326
}


# Here is the protein sequence to weight
Protein = 'DHPFWKQTACKHV'

# Here we compute the mass, using the weights of the dictionary
weight = sum(weights[p] for p in Protein)

# Here we print the result -> The result seems close to the mass on the teachers' notes
print ("The molecular weight of protein " + Protein + " is " + str(weight) )


