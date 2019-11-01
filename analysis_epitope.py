#!/usr/bin/python

# For peptide_match


# Own module
import peptide_functions as pepfun

epitope="DYDVVYLKPLAGMYK"

pep=pepfun.peptide_sequence(epitope)
pep.compute_peptide_charges()
print "The net charge of %s is:" %epitope,pep.netCharge
pep.helical_wheel()
pep.generate_structure_modpep()
print pep.smiles
pep.generate_conformer()
# pep.calculate_properties_from_mol()
# print pep.mol_logp
pep.similar_smiles("DYDVVYLKPLAGKYK")
print pep.smiles_similarity
#pep.blast_online()
# Check if the peptide is in a local version of UniProt
#pepfun.peptide_match(epitope)

# peptide_lib=pepfun.combinatorial_library(8)
# print len(peptide_lib)
# count_dict=pepfun.frequencies_library(peptide_lib)
# print peptide_lib
# 
# library_file=open("list_combinatorial_library.txt","wb")
# for seq in peptide_lib:
#     library_file.write("%s\n" %seq)
# library_file.close()

# # Example structure
pdb_file="example_structure.pdb"
chain="C"
pepStr=pepfun.peptide_structure(pdb_file,chain)
print pepStr.sequence
pepStr.get_secondary_structure()
print pepStr.total_dssp
pepStr.get_hydrogen_bonds()

# pepStr.get_secondary_structure_PROSS()
# print pepStr.total_pross
# contact_threshold=4.0
# pepStr.get_heavy_atom_contacts(contact_threshold)
# print pepStr.total_contacts
# pepStr.calculate_dihedrals()
# pepStr.peptide_minimization()
# print pepStr.positions
