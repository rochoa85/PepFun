#!/usr/bin/python

"""
Pepfun: bioinformatics and cheminformatics protocols for peptide-related computational analysis

From publication "Pepfun: bioinformatics and cheminformatics protocols for peptide-related computational analysis"
Journal of Cheminformatics 
Authors: Rodrigo Ochoa, Roman Laskowski, Pilar Cossio
Year: 2020

Script to test the PepFun functionalities
"""

########################################################################################
# Authorship
########################################################################################

__author__ = "Rodrigo Ochoa"
__credits__ = ["Rodrigo Ochoa", "Roman Laskowski", "Pilar Cossio"]
__license__ = "MIT"
__version__ = "1.0"
__email__ = "rodrigo.ochoa@udea.edu.co"

########################################################################################
# Import PepFun functions
########################################################################################

from peptide_functions import *

########################################################################################
# Main File
########################################################################################

if __name__ == '__main__':
    
    ############################################################
    # Example input data
    ############################################################
    
    sequence="KMGRLFR"
    pdb_file="auxiliar/example_structure.pdb"
    chain="C"
    pep_conformation="linear"
    contact_threshold=4.0
    
    ############################################################
    # Test peptide sequence functions
    ############################################################
    
    print("### TEST PEPTIDE SEQUENCE FUNCTIONS ###")
    pep=peptide_sequence(sequence)
    print("The main peptide sequence is: {}".format(pep.sequence))
    
    pep.compute_peptide_charges()
    print("Net charge at pH 7: {}".format(pep.netCharge))
    
    pep.calculate_properties_from_mol()
    print("Molecular weight: {}".format(pep.mol_weight))
    
    pep.calculate_properties_from_sequence()
    print("Average Hydrophobicity: {}".format(pep.avg_hydro))
    print("Isoelectric Point: {}".format(pep.isoelectric_point))
    
    matrix=matlist.structure
    pep.align_position_matrix("KAGRSFR",matrix)
    print("Alignment score by position with peptide KAGRSFR is: {}".format(pep.score_matrix))
    
    pep.align_position_local("KAGRSFR")
    print("The number of dismatches with peptide KAGRSFR are: {}".format(pep.dismatch))
    
    similarity=pep.similarity_pair("KMGRLFR","KAGRSFR",matrix)
    print("The similarity between peptides KMGRLFR and KAGRSFR is: {}".format(similarity))
    
    pep.similar_smiles("KAGRSFR")
    print("The SMILES similarity is: {}".format(pep.smiles_similarity))
    
    pep.solubility_rules()
    print("{} solubility rules failed from 5".format(pep.solubility_rules_failed))
    
    pep.synthesis_rules()
    print("{} synthesis rules failed from 5".format(pep.synthesis_rules_failed))
    
    ############################################################
    # Test peptide structure functions
    ############################################################
    
    print("### TEST PEPTIDE STRUCTURE FUNCTIONS ###")
    
    pepStr=peptide_structure(pdb_file,chain)
    print("Peptide sequence based on the PDB file is: {}".format(pepStr.sequence))
    
    pepStr.get_secondary_structure()
    print("The predicted secondary structure is: {}".format(pepStr.total_dssp))
    
    pepStr.get_hydrogen_bonds()
    pepStr.plot_hydrogen_bonds(pep_conformation)
    pepStr.get_heavy_atom_contacts(contact_threshold)
    
    print("The total number of contacts are: {}".format(pepStr.total_contacts))
    print("The total number of hydrogen bonds are: {}".format(pepStr.number_hydrogen_bonds))
    print("The following are the details per amino acid in the peptide:")
    print(pepStr.positions)
    
    ############################################################
    # Test additional functions
    ############################################################
    
    print("### TEST ADDITIONAL FUNCTIONS ###")
    
    # Select one of the following methods:
    
    #list_peptides=generate_sequences(2,"natural")
    #list_peptides=combinatorial_library(10,"natural")
    list_peptides=generate_peptide_pattern("XERTX")
    
    # Print the peptides and store them in an external file
    print("The number of peptides generated in the library are: {}".format(len(list_peptides)))
    library_report=open("auxiliar/library_generated.txt","w")
    for pep in list_peptides: library_report.write("{}\n".format(pep))
    library_report.close()