#!/usr/bin/python

"""
Pepfun: bioinformatics and cheminformatics protocols for peptide-related computational analysis

From publication "Pepfun: bioinformatics and cheminformatics protocols for peptide-related computational analysis"
Journal of Cheminformatics 
Authors: Rodrigo Ochoa, Lucy Jimenez, Roman Laskowski, ..., Pilar Cossio
Year: 2019

Third-party tools required:

BioPython: https://biopython.org/wiki/Download - Ubuntu package: python-rdkit
RDKit: https://github.com/rdkit/rdkit/releases - Ubuntu package: python-biopython
"""

########################################################################################
# Authorship
########################################################################################

__author__ = "Rodrigo Ochoa"
__credits__ = ["Rodrigo Ochoa", "Lucy Jimenez","Roman Laskowski", "...", "Pilar Cossio"]
__license__ = "MIT"
__version__ = "1.0"
__email__ = "rodrigo.ochoa@udea.edu.co"

########################################################################################
# Modules to import
########################################################################################

# External modules
import math
import itertools
import subprocess
import argparse
import re
import os
import numpy as np
from random import randint
from igraph import *

# BioPython
from Bio import pairwise2
from Bio.SubsMat import MatrixInfo as matlist
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio.PDB import *
from Bio import SeqIO

# RDKit
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import DataStructs
from rdkit.Chem import Crippen
from rdkit.Chem import Lipinski
from rdkit.Chem import Descriptors

# Local scripts
# import conformers


########################################################################################
# General Functions
########################################################################################

def aminoacidSMILES(amino):
    """
    Obtain the SMILES representation of a particular amino acid
    
    Arguments:
    amino -- One-letter code of the amino acid
    
    Return:
    smiles -- 1D Chemical representation of the amino acid 
    """
    
    # Dictionary with the SMILES per amino acid
    aminoacids = {'G':{'SMILES': 'NCC(=O)O'},
                  'A':{'SMILES': 'N[C@@]([H])(C)C(=O)O'},
                  'R':{'SMILES': 'N[C@@]([H])(CCCNC(=N)N)C(=O)O'},
                  'N': {'SMILES': 'N[C@@]([H])(CC(=O)N)C(=O)O'},
                  'D': {'SMILES': 'N[C@@]([H])(CC(=O)O)C(=O)O'},
                  'C': {'SMILES': 'N[C@@]([H])(CS)C(=O)O'},
                  'E': {'SMILES': 'N[C@@]([H])(CCC(=O)O)C(=O)O'},
                  'Q': {'SMILES': 'N[C@@]([H])(CCC(=O)N)C(=O)O'},
                  'H': {'SMILES': 'N[C@@]([H])(CC1=CN=C-N1)C(=O)O'},
                  'I': {'SMILES': 'N[C@@]([H])(C(CC)C)C(=O)O'},
                  'L': {'SMILES': 'N[C@@]([H])(CC(C)C)C(=O)O'},
                  'K': {'SMILES': 'N[C@@]([H])(CCCCN)C(=O)O'},
                  'M': {'SMILES': 'N[C@@]([H])(CCSC)C(=O)O'},
                  'F': {'SMILES': 'N[C@@]([H])(Cc1ccccc1)C(=O)O'},
                  'P': {'SMILES': 'N1[C@@]([H])(CCC1)C(=O)O'},
                  'S': {'SMILES': 'N[C@@]([H])(CO)C(=O)O'},
                  'T': {'SMILES': 'N[C@@]([H])(C(O)C)C(=O)O'},
                  'W': {'SMILES': 'N[C@@]([H])(CC(=CN2)C1=C2C=CC=C1)C(=O)O'},
                  'Y': {'SMILES': 'N[C@@]([H])(Cc1ccc(O)cc1)C(=O)O'},
                  'V': {'SMILES': 'N[C@@]([H])(C(C)C)C(=O)O'}}
    
    # Store the SMILES in a variable
    smiles=aminoacids[amino]['SMILES']
    return smiles

########################################################################################

def generate_sequences(long_peptide,aa_type="natural"):
    """
    Method to generate a combinatorial library of peptide sequences using natural and d-amino acids
    
    Arguments:
    long_peptide -- Lenght of the peptides that are going to be generated
    aa_type -- Nature of the amino acids that can be included. Options: natural (default), d_aminoacids, combined
    
    Return:
    frags -- List with the sequences generated
    """
    
    # List of the amino acids that can be included
    natural=["A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y"]
    d_aminoacids=["[dA]","[dC]","[dD]","[dE]","[dF]","[dH]","[dI]","[dK]","[dL]","[dM]","[dN]","[dP]","[dQ]","[dR]","[dT]","[dV]","[dW]"]
    combined=["A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y",
              "[dA]","[dC]","[dD]","[dE]","[dF]","[dH]","[dI]","[dK]","[dL]","[dM]","[dN]","[dP]","[dQ]","[dR]","[dT]","[dV]","[dW]"]
    
    # List where the sequnces will be stored
    frags=[]
    if aa_type=="natural":
        # Generate the combinatorial library
        for i in itertools.product(natural, repeat=long_peptide):
            frags.append(''.join(map(str, i)))
    elif aa_type=="d_aminoacids":
        # Generate the combinatorial library
        for i in itertools.product(d_aminoacids, repeat=long_peptide):
            frags.append(''.join(map(str, i)))
    elif aa_type=="combined":
        # Generate the combinatorial library
        for i in itertools.product(combined, repeat=long_peptide):
            if "[" in f: frags.append(''.join(map(str, i)))
    else:
        print "The type of amino acid is invalid"
        
    # Return the list of peptide sequences
    return frags

########################################################################################

def generate_peptide_pattern(pattern):
    """
    Method to generate a peptide sequences following a pattern using only natural amino acids
    
    Arguments:
    pattern -- Pattern based on XX to generate the list of sequences
    
    Return:
    list_peptides -- List with the sequences generated
    """
    
    # List of natural amino acids that can be included
    natural=["A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y"]
    
    list_peptides=[]
    
    # Combination of possible sequences following the given pattern
    for pep in itertools.product(natural, repeat=pattern.count('X')):
        pep = list(pep)
        mut_pep = []
        for res in pattern:
            # check if any amino acid can be included in the sequence
            if  res != 'X': mut_pep.append(res)
            else: mut_pep.append(pep.pop(0))
        list_peptides.append("".join(mut_pep))
    
    # Return the list of peptide sequences
    return list_peptides

########################################################################################

def combinatorial_library(long_peptide,aa_type="natural"):
    """
    Method to generate a combinatorial library of peptide sequences but limiting the frequency of amino acids in certain positions
    
    Arguments:
    long_peptide -- Lenght of the peptides that are going to be generated
    aa_type -- Nature of the amino acids that can be included. Options: natural (default), d_aminoacids, combined
    
    Return:
    peptides -- List with the sequences generated
    """
    
    # List of the amino acids that can be included
    natural=["A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y"]
    d_aminoacids=["[dA]","[dC]","[dD]","[dE]","[dF]","[dH]","[dI]","[dK]","[dL]","[dM]","[dN]","[dP]","[dQ]","[dR]","[dT]","[dV]","[dW]"]
    combined=["A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y",
              "[dA]","[dC]","[dD]","[dE]","[dF]","[dH]","[dI]","[dK]","[dL]","[dM]","[dN]","[dP]","[dQ]","[dR]","[dT]","[dV]","[dW]"]
    
    # List where the sequences will be stored
    peptides=[]
    if aa_type=="natural":
        # Counter of the peptide length
        for i in range(0,long_peptide):
            for aa in natural:
                # Run the process 3 times per amino acids
                for j in range(0,3):
                    sequence=[None]*long_peptide
                    for k in range(0,long_peptide):
                        if k==i:
                            sequence[k]=aa
                        else:
                            # Loop to check that the sequence does not have more than 3 equal amino acids
                            accept_change=0
                            while accept_change==0:
                                posRand=randint(0,len(natural)-1)
                                new_aa=natural[posRand]
                                if sequence.count(new_aa)<=3:
                                    accept_change=1
                            # Add the amino acid after checking the restraint
                            sequence[k]=new_aa
                    
                    # Join the amino acids and append the new sequence
                    seq_string=''.join(sequence)
                    peptides.append(seq_string)
        
    # Return the list
    return peptides

########################################################################################

def frequencies_library(peptides):
    """
    Calculate the frequency of each amino acid in a particular peptide library
    
    Arguments:
    peptides -- List of peptides in the library
    
    Return:
    count_dict -- Dictionary with the numbers per each natural amino acid
    """
    
    # Create counter dictionary
    count_dict={}
    for i in range(1,len(peptides[0])+1):
        count_dict[i]={"A":0,"R":0,"N":0,"D":0,"C":0,"Q":0,"E":0,"G":0,"H":0,"I":0,"L":0,"K":0,"M":0,"F":0,"P":0,"S":0,"T":0,"W":0,"Y":0,"V":0}
    
    # Read the sequences and count each amino acid    
    for sequence in peptides:
        for pos,aa in enumerate(sequence):
            count_dict[pos+1][aa]+=1
    
    # Print the dictionary
    for key in count_dict:
        print key,count_dict[key]
    
    # Return the dictionary
    return count_dict

########################################################################################
# Classes and Functions
########################################################################################

class peptide_sequence:
    """
    Class with functions to perform different type of analysis using a peptide sequence as an object
    """
    
    def __init__(self,sequence):
        """
        Inititalize the class calculating some basic properties
        
        Arguments:
        sequence -- Peptide sequence
        
        Return
        Based on the sequence, the class start with some counters, a neutral pH, the peptide length and the SMILES representation
        """
        self.sequence=sequence
        self.pH=7
        self.solubility_rules_failed=0
        self.length_peptide=len(self.sequence)
        self.synthesis_rules_failed=0
        
        # Loop to create the SMILES
        connect_smiles='O'
        for res in sequence:
             connect_smiles=connect_smiles[:-1]
             smiles=aminoacidSMILES(res)          
             connect_smiles=connect_smiles+smiles
        self.smiles=connect_smiles
    
    ############################################################################    
    def align_position_matrix(self,peptide_to_match,matrix):
        """
        Align position by position a peptide of reference with another one
        
        Arguments:
        peptide_to_match -- String with the sequence of a peptide that will be compared
        matrix -- A dictionary of biopython with substitution scores.
        
        Return:
        score -- A numerical value to rank the alignment
        """
        self.score_matrix=0
        
        for i,p in enumerate(self.sequence):
            # Generate the tuples with the pair of amino acids to compare
            pair1=(p,peptide_to_match[i])
            pair2=(peptide_to_match[i],p)
            
            # Obtain the score from the matrix and sum up the values
            if pair1 in matrix: value=matrix[pair1]
            else: value=matrix[pair2]
            self.score_matrix+=value
    
    ############################################################################
    def align_position_local(self,peptide_to_match):
        """
        Align position by position a peptide of reference with another one
        
        Arguments:
        peptide_to_match -- String with the sequence of a peptide that will be compared
        
        Return:
        dismatch - 
        """
        self.dismatch=0
        
        for i,p in enumerate(self.sequence):
            # Generate the tuples with the pair of amino acids to compare
            if p!=peptide_to_match[i]:
                self.dismatch+=1
        
        self.match=len(self.sequence)-self.dismatch
        
    ############################################################################
    def similarity_pair(self,peptide1,peptide2,matrix):
        """
        Function to calculate similarity between two peptide sequences
        
        Arguments:
        peptide1 -- Sequence of one of the input peptides
        peptide2 -- Sequence of the second input peptide
        matrix -- BioPython matrix chosen for running the analysis
        
        Return:
        sim_val -- similarity based on the aligments between the peptides and themselves - Max value is 1
        """
        # Alignment between peptide1 and peptide2
        pep=peptide_sequence(peptide1)
        pep.align_position_matrix(peptide2,matrix)
        score1_2=pep.score_matrix
        
        # Alignment between peptide1 with itself
        pep.align_position_matrix(peptide1,matrix)
        score1_1=pep.score_matrix
        
        # Alignment between peptide2 with itself
        pep2=peptide_sequence(peptide2)
        pep2.align_position_matrix(peptide2,matrix)
        score2_2=pep2.score_matrix
        
        # Calculate similarity value
        sim_val=float(score1_2)/math.sqrt(float(score1_1*score2_2))
        
        # Return similarity
        return sim_val
    
    ############################################################################    
    def similar_smiles(self,peptide_to_match):
        """
        Calculate similarity but using SMILES representations of the peptides
        
        Arguments:
        peptide_to_match -- peptide sequence that will be compared
        
        Return:
        SMILES similarity based on Morgan Fingerprints and Tanimoto coefficient
        """
        
        # Generate molecule from sequence
        mol1 = Chem.MolFromSmiles(self.smiles)
        mol1.SetProp("_Name",self.sequence)
        
        connect_smiles='O'
        for res in peptide_to_match:
             connect_smiles=connect_smiles[:-1]
             smiles=aminoacidSMILES(res)          
             connect_smiles=connect_smiles+smiles

        mol2 = Chem.MolFromSmiles(connect_smiles)
        mol2.SetProp("_Name",peptide_to_match)
        
        # Calculate the fingerprints and the similarity
        fp1=AllChem.GetMorganFingerprintAsBitVect(mol1,2,2048)
        fp2=AllChem.GetMorganFingerprintAsBitVect(mol2,2,2048)
        
        self.smiles_similarity=DataStructs.TanimotoSimilarity(fp1,fp2)

    ############################################################################
    def compute_peptide_charges(self,pH_internal=7):
        """
        Function to calculate the average net charge based on pka values
        
        Arguments:
        pH_internal -- By default is 7
        
        Return:
        The net charge based on reported pka values
        """
        # Set the general variables and pka terms
        self.pH=pH_internal
        self.netCharge=0.0
        pka_alpha_amino={'G':9.60,'A':9.69,'V':9.62,'L':9.60,'I':9.68,'M':9.21,'F':9.13,'W':9.39,'P':10.60,'S':9.15,
                         'T':9.10,'C':10.78,'Y':9.11,'N':8.84,'Q':9.13,'D':9.82,'E':9.67,'K':8.95,'R':9.04,'H':9.17}
        pka_alpha_carboxy={'G':2.34,'A':2.34,'V':2.32,'L':2.36,'I':2.36,'M':2.28,'F':1.83,'W':2.38,'P':1.99,'S':2.21,
                           'T':2.63,'C':1.71,'Y':2.2,'N':2.02,'Q':2.17,'D':2.09,'E':2.19,'K':2.18,'R':2.17,'H':1.82}
        pka_sidechain_positive={'K':10.79,'R':12.48,'H':6.04}
        pka_sidechain_negative={'D':3.86,'E':4.25,'C':8.33,'Y':10.07}
        
        # Calculate the net charge for the extreme groups (without modifications)
        amino=self.sequence[0]
        carboxy=self.sequence[-1]
        self.netCharge=self.netCharge+(math.pow(10,pka_alpha_amino[amino])/(math.pow(10,pka_alpha_amino[amino])+math.pow(10,self.pH)))
        self.netCharge=self.netCharge-(math.pow(10,self.pH)/(math.pow(10,pka_alpha_carboxy[carboxy])+math.pow(10,self.pH)))
        
        # Calculate the net charge for the charged amino acid side chains
        for aa in self.sequence:
            if aa in pka_sidechain_positive:
                self.netCharge=self.netCharge+(math.pow(10,pka_sidechain_positive[aa])/(math.pow(10,pka_sidechain_positive[aa])+math.pow(10,self.pH)))
            if aa in pka_sidechain_negative:
                self.netCharge=self.netCharge-(math.pow(10,self.pH)/(math.pow(10,pka_sidechain_negative[aa])+math.pow(10,self.pH)))
    
    ############################################################################           
    def calculate_properties_from_mol(self):
        """
        Function to calculate some molecular properties based on RDKit functionalities
        
        Return:
        Static physico-chemical properties: molecular weight, crippen logP, number of hydrogen bond acceptors and donors
        """
        
        # Generate molecule from sequence
        mol = Chem.MolFromSmiles(self.smiles)
        mol.SetProp("_Name",self.sequence)
        
        # Calculate the descriptors
        self.num_hdonors = Lipinski.NumHDonors(mol)
        self.num_hacceptors = Lipinski.NumHAcceptors(mol)
        self.mol_weight = Descriptors.MolWt(mol)
        self.mol_logp = Crippen.MolLogP(mol)
    

    
    ############################################################################           
    def calculate_properties_from_sequence(self):
        """
        Function to calculate some molecular properties based on RDKit functionalities
        
        Arguments:
        Sequence - amino acid sequence of the peptide
        
        Return:
        Average Eisenberg hydrophobicity
        ProtParam parameters: Isolectric point, aromaticity, instability index, amino acid percentage
        """
        
        # Hydrophobicity -> Eisenberg scale
        hydrophobicity = {'A':  0.620,'R': -2.530,'N': -0.780,'D': -0.900,'C':  0.290,'Q': -0.850,'E': -0.740,'G':  0.480,'H': -0.400,'Y':  0.260,
                          'I':  1.380,'L':  1.060,'K': -1.500,'M':  0.640,'F':  1.190,'P':  0.120,'S': -0.180,'T': -0.050,'W':  0.810,'V':  1.080}
        self.avg_hydro=sum([hydrophobicity[resi] for resi in self.sequence])
        
        # ProParam properties
        prot_parameters=ProteinAnalysis(self.sequence)
        self.aromaticity=prot_parameters.aromaticity()
        self.aa_percent=prot_parameters.get_amino_acids_percent()
        self.instability_index=prot_parameters.instability_index()
        self.isoelectric_point=prot_parameters.isoelectric_point()
        
    ############################################################################  
    def solubility_rules(self):
        """
        Function to calculate some solubility rules based on recommendations of http://bioserv.rpbs.univ-paris-diderot.fr/services/SolyPep/
        
        Output:
        solubility_rules_failed - return the number of rules faild based on the criteria
        """
        # Rule N1. Number of hydrophobic or charged residues
        hydro_residues=['V','I','L','M','F','W','C']
        charged_residues=['H','R','K','D','E']
                
        
        count_hydro_charged=0
        for aa in self.sequence:
            if aa in hydro_residues or aa in charged_residues: count_hydro_charged+=1
        
        # This condition should change depending on the sequence length
        hydro_char_threshold=float(self.length_peptide)*0.45
        if count_hydro_charged > hydro_char_threshold: self.solubility_rules_failed+=1
        
        # Rule N2. Computed peptide charge
        charge_threshold=1
        self.compute_peptide_charges()
        if self.netCharge > 1:
            self.solubility_rules_failed+=1
            
        # Rule N3. Glycine or Proline content in the sequence
        count_gly_pro=0
        for aa in self.sequence:
            if aa == "G" or aa=="P": count_gly_pro+=1
        # Check threshold
        if count_gly_pro > 1: self.solubility_rules_failed+=1
        
        # Rule N4. First or last amino acid charged
        count_charge=0
        if self.sequence[0] in charged_residues:
            count_charge+=1
        if self.sequence[-1] in charged_residues:
            count_charge+=1
        # Check threshold
        if count_charge > 0: self.solubility_rules_failed+=1
        
        # Rule N5. Any amino acid represent more than 25% of the total sequence
        prot_parameters=ProteinAnalysis(self.sequence)
        aa_content=prot_parameters.get_amino_acids_percent()
        for aa in aa_content:
            if aa_content[aa]>=0.3:
                self.solubility_rules_failed+=1
                break
    
    ############################################################################
    def synthesis_rules(self):
        """
        Function to check some synthesis rules based on empirical recommendations
        
        Return:
        synthesis_rules_failed - return the number of rules faild based on the criteria
        """
        # Presence of forbiden motifs
        forbidden_motifs = {'2-prolines':r'[P]{3,}','DG-DP':r'D[GP]','N-Q-Nterminal':r'^[NQ]',}
        for motif in forbidden_motifs:
            if re.search(forbidden_motifs[motif],self.sequence):
                self.synthesis_rules_failed+=1
                
        # test if there are charged residues every 5 amino acids
        charged_residues=['H','R','K','D','E']
        counter_charged = 0
        for residue in self.sequence:
            counter_charged += 1
            if residue in charged_residues:
                counter_charged = 0
            if counter_charged >= 5:
                self.synthesis_rules_failed+=1
                
        # Check if there are oxidation-sensitive amino acids
        aa_oxidation=['M','C','W']
        for aa in self.sequence:
            if aa in aa_oxidation:
                self.synthesis_rules_failed+=1
                break
        
    ############################################################################
    def generate_conformer(self):
        """
        Function to generate basic conformer using RDKit
        
        Return:
        Structure predicted in sdf format
        NOTE: using OpenBabel the structure can be converted to PDB format using:
        babel -isdf structure_{peptide}.sdf -opdb structure_{peptide}.pdb
        """
        # Generate molecule from smiles
        mol = Chem.MolFromSmiles(self.smiles)
        mol.SetProp("_Name",self.sequence)
        
        print "Generating the basic conformer for peptide {}".format(self.sequence)
        
        # Use UFF force field
        AllChem.EmbedMolecule(mol)
        AllChem.UFFOptimizeMolecule(mol)
        writer = AllChem.SDWriter("structure_%s.sdf" %self.sequence)
        writer.write(mol)
                
########################################################################################
    
class peptide_structure:
    
    """
    Class with functions to perform different type of analysis using a peptide structure alone or in complex with a protein
    """

    def __init__(self,pdb_file,chain):
        """
        Inititalize the class calculating some basic properties
        
        Arguments:
        pdb_file -- PDB file with the required information
        chain -- chain containing the peptide in the file
        
        Return:
        Based on the structure, the class start with some counters, the sequence derived from the structure, its lenght and a dictionary with the aa positions
        """
        
        self.pdb_file=pdb_file
        self.chain=chain
        self.aminoacids_back={"ALA":"A","ASP":"D","GLU":"E","PHE":"F","HIS":"H","ILE":"I","LYS":"K","LEU":"L","MET":"M","GLY":"G",
                              "ASN":"N","PRO":"P","GLN":"Q","ARG":"R","SER":"S","THR":"T","VAL":"V","TRP":"W","TYR":"Y","CYS":"C"}
        self.aminoacids={"A":"ALA","D":"ASP","E":"GLU","F":"PHE","H":"HIS","I":"ILE","K":"LYS","L":"LEU","M":"MET","G":"GLY",
                         "N":"ASN","P":"PRO","Q":"GLN","R":"ARG","S":"SER","T":"THR","V":"VAL","W":"TRP","Y":"TYR","C":"CYS"}
        self.sequence=""
        self.initial_id_residue=0
        self.positions={}
        
        # Read the structure in BioPython format
        parser = PDBParser()
        self.reference = parser.get_structure('REF',pdb_file)
        for ch in self.reference[0]:
            if ch.get_id()==chain:
                for i,residue in enumerate(ch):
                    seq=self.aminoacids_back[residue.get_resname()]
                    # Save the sequence
                    self.sequence=self.sequence+str(seq)
                    if i==0: self.initial_id_residue=residue.get_full_id()[3][1]
                    # Save the positions
                    self.positions[i+1]={"aa":str(seq),"full_aa":residue.get_resname()}
        
        # Store sequence lenght
        self.len_sequence=len(self.sequence)

        
    ############################################################################
    
    def get_secondary_structure(self):
        """
        Function to calculate the secondary structure and the accessible surface area using the auxiliary program mkdssp
        NOTE: mkdssp can be downloaded from the website:
        
        Return:
        The dictionary positions will store the calculated data of dssp and the asa values
        total_dssp will contain the complete predicted secondary structure with the following conventions:
        B -
        H - alpha helix
        E - beta sheets
        """
        
        model=self.reference[0]
        dssp = DSSP(model,self.pdb_file,dssp='auxiliar/mkdssp')
        self.total_dssp=""
        
        # Loop over the keys from the dssp response to store ss and asa values
        for keys in list(dssp.keys()):
            if keys[0]==self.chain:
                position=keys[1][1]
                self.positions[position]["dssp"]=dssp[keys][2]
                self.positions[position]["asa"]=dssp[keys][3]
                self.total_dssp=self.total_dssp+dssp[keys][2]
    
    ############################################################################
    
    def get_hydrogen_bonds(self):
        """
        Function to calculate hydrogen bonds with the other chains
        
        Return:
        Dictionary containing the peptide amino acids and the protein amino acids forming hydrogen bonds
        File with the predicted hydrogen bonds
        """
        
        # Read the protein structure
        model=self.reference[0]
        dssp = DSSP(model,self.pdb_file,dssp='auxiliar/mkdssp')
        self.total_dssp=""
        
        # Loop over the keys from the dssp response to store ss and asa values
        total_index={}
        list_chains=[]
        reference_position=0
        list_chains.append(list(dssp.keys())[0][0])
        self.hbonds_peptide={}
        
        # iterate over the dssp keys to store the corresponding hydrogen bonds based on the atoms positions
        for keys in list(dssp.keys()):
            if keys[0]==list_chains[-1]:
                total_index[reference_position+keys[1][1]]=(keys[0],dssp[keys][1],keys[1][1])
                last_position=reference_position+keys[1][1]
            else:
                list_chains.append(keys[0])
                reference_position=last_position
                total_index[reference_position+keys[1][1]]=(keys[0],dssp[keys][1],keys[1][1])
                last_position=reference_position+keys[1][1]
                
            if keys[0]==self.chain:
                position=keys[1][1]
                amino_name=dssp[keys][1]+str(keys[1][1])
                if amino_name not in self.hbonds_peptide: self.hbonds_peptide[amino_name]=[]
                if dssp[keys][6]<-5:
                    interactions_pos=last_position+dssp[keys][6]
                    self.hbonds_peptide[amino_name].append(total_index[interactions_pos])
                if dssp[keys][8]<-5:
                    interactions_pos=last_position+dssp[keys][8]
                    self.hbonds_peptide[amino_name].append(total_index[interactions_pos])
                if dssp[keys][10]<-5:
                    interactions_pos=last_position+dssp[keys][10]
                    self.hbonds_peptide[amino_name].append(total_index[interactions_pos])
                if dssp[keys][12]<-5:
                    interactions_pos=last_position+dssp[keys][12]
                    self.hbonds_peptide[amino_name].append(total_index[interactions_pos])
        
        # Iterate over the peptide residues to show the hydrogen bonds
        self.number_hydrogen_bonds=0
        output_hydrogen_bonds=open("auxiliar/predicted_hbs_{}.txt".format(self.sequence),"w")
        print "These are the hydrogen bonds detected:"
        for residue in self.hbonds_peptide:
            for partners in self.hbonds_peptide[residue]:
                print "{} interacts with residue {}{} from chain {}".format(residue,partners[1],partners[2],partners[0])
                output_hydrogen_bonds.write("{} interacts with residue {}{} from chain {}\n".format(residue,partners[1],partners[2],partners[0]))
                self.number_hydrogen_bonds+=1
        output_hydrogen_bonds.close()
    
    ############################################################################
    
    def plot_hydrogen_bonds(self):       
        """
        Function to plot the hydrogen bonds using the igraph module of Python.
        NOTE: Installation instructions: https://igraph.org/python/
        
        Return:
        PNG file with the graph of the hydrogen bonds
        """
        
        # Generate fragments
        fragment_chains={}
        for amino in self.hbonds_peptide:
            for receptor in self.hbonds_peptide[amino]:
                chain=receptor[0]
                if chain not in fragment_chains: fragment_chains[chain]=[]
                if receptor not in fragment_chains[chain]:
                    fragment_chains[chain].append(receptor)
        
        for ch in fragment_chains:
            fragment_chains[ch].sort(key=lambda x: x[2])
        
        
        # Get graph order and properties of the nodes and edges
        names=[0]*len(self.hbonds_peptide)
        pep_interactions=[]
        chain_nodes=[]
        chain_colors={"PEP":"red"}
        edge_colors={"5":"black","4":"orange","8":"orange"}
        list_colors=["cyan","green","magenta","blue"]
        type_interactions=[]
        
        # Iterate over the predicted hydrogen bonds
        for i,amino in enumerate(self.hbonds_peptide):
            pos=int(amino[1:])
            names[pos-1]=amino
            chain_nodes.append("PEP")
            
            if i!=len(self.hbonds_peptide)-1:
                pep_interactions.append((i,i+1))
                type_interactions.append(5)
        
        # Iterate over the fragments to assign the properties and types of interactions
        reference=len(self.hbonds_peptide)
        for num_chain,chains in enumerate(fragment_chains):
            chain_colors[chains]=list_colors[num_chain]
            for count,residues in enumerate(fragment_chains[chains]):
                chain_nodes.append(chains)
                names.append(residues[0]+"-"+residues[1]+str(residues[2]))
                number_res=reference+count
                
                for pep_amino in self.hbonds_peptide:
                    if residues in self.hbonds_peptide[pep_amino]:
                        pep_interactions.append((names.index(pep_amino),number_res))
                        type_interactions.append(self.hbonds_peptide[pep_amino].count(residues)*4)
                        
            reference=len(names)
            
        # Create graph
        g = Graph(pep_interactions)
        g.vs["name"] = names
        g.vs["chains"] = chain_nodes
        g.es["count_int"] = type_interactions
        color_per_chain=[chain_colors[chColor] for chColor in g.vs["chains"]]
        color_per_edge=[edge_colors[str(edColor)] for edColor in type_interactions]
    
        # Use a default layout 
        layout =g.layout_fruchterman_reingold()
        #layout = [(0,0),(0,1),(0,2),(0,3),(0,4),(0,5),(0,6),(0,7),(0,8),(0,9),(0,10)]
        visual_style = {}
        visual_style["vertex_size"] = 50
        visual_style["vertex_color"] = color_per_chain
        visual_style["edge_color"] = color_per_edge
        visual_style["vertex_label"] = g.vs["name"]
        visual_style["edge_width"] = type_interactions
        visual_style["layout"] = layout
        visual_style["bbox"] = (800, 800)
        visual_style["margin"] = 80
        plot(g,"auxiliar/plot_hbs_{}.png".format(self.sequence), **visual_style)
        
    ############################################################################
    
    def get_heavy_atom_contacts(self,contact_threshold):
        """
        Function to count heavy atom contacts per residue
        
        Arguments:
        contact_threshold -- theshold to define when a contact is created
        
        Return:
        Add to the amino acids dictionary the number of contacts per amino acid
        total_contacts -- total number of contacts calculated for the full peptide
        """
        
        # Get the other chains present in the PDB file
        chains=[]
        for ch in self.reference[0]:
            chLetter=ch.get_id()
            chains.append(chLetter)
        
        # Store the count per amino acid in the peptide
        countDict={}
        
        # Loop over the chains
        for chValues in chains:
            if chValues!=self.chain:
                for residue in self.reference[0][self.chain]:
                    # Store the distances and a counter
                    counter=0
                    distances=[]
                    for residueA in self.reference[0][chValues]:
                        for atomA in residueA:
                            idAtomA=atomA.get_id()
                            if idAtomA[0] != "H" and idAtomA[0].isdigit() == False:
                                for atom in residue:
                                    idAtom = atom.get_id()
                                    if idAtom[0] != "H" and idAtom[0].isdigit() == False:
                                        # Get the distance differences
                                        diff = atom.coord - atomA.coord
                                        diffValue=np.sqrt(np.sum(diff * diff))
                                        distances.append(diffValue)
                    # Count the number of contacts based on the defined threshold
                    for d in distances:
                        if d < contact_threshold: counter+=1
                    
                    # Get the information per residue
                    res=residue.get_resname()
                    resPos=str(residue.get_full_id()[3][1])
                    if res+resPos not in countDict:
                        countDict[res+resPos]=counter
                    else:
                        countDict[res+resPos]+=counter
        
        # Calculate the final contacts per aa and in total
        self.total_contacts=0
        for i,aa in enumerate(self.sequence):
            residue=self.aminoacids[aa]
            position=i+1
            if residue+str(position) in countDict:
                self.positions[position]["contacts"]=countDict[residue+str(position)]
                self.total_contacts+=countDict[residue+str(position)]
    

########################################################################################
########################################################################################
########################################################################################
# Main execution
########################################################################################
########################################################################################
########################################################################################
if __name__ == '__main__':
    
    # Script arguments
    parser = argparse.ArgumentParser(description='Pepfun: bioinformatics and cheminformatics protocols for peptide-related computational analysis')
    parser.add_argument('-s', dest='pep_seq', action='store',required=True,
                        help='Sequence of peptide to be analyzed')
    parser.add_argument('-p', dest='pep_str', action='store', default="auxiliar/example_structure.pdb",
                        help='Structure that will be used for the analysis')
    # Pending add more
    args = parser.parse_args()
    
    ####################################################################################
    # Assignment of parameters
    sequence=args.pep_seq
    # sequence="KMGRLFR"
    
    ############################################################
    # Test peptide sequence functions
    ############################################################
    print
    print "### TEST PEPTIDE SEQUENCE FUNCTIONS ###"
    pep=peptide_sequence(sequence)
    print "The main peptide sequence is: ",pep.sequence
    
    pep.compute_peptide_charges()
    print "Net charge at pH 7: ",pep.netCharge
    
    pep.calculate_properties_from_mol()
    print "Molecular weight: ",pep.mol_weight
    
    pep.calculate_properties_from_sequence()
    print "Average Hydrophobicity: ",pep.avg_hydro
    print "Isoelectric Point: ",pep.isoelectric_point
    
    matrix=matlist.structure
    pep.align_position_matrix("KAGRSFR",matrix)
    print "Alignment score by position with peptide KAGRSFR is: ",pep.score_matrix
    
    pep.align_position_local("KAGRSFR")
    print "The number of dismatches with peptide KAGRSFR are: ",pep.dismatch
    
    similarity=pep.similarity_pair("KMGRLFR","KAGRSFR",matrix)
    print "The similarity between peptides KMGRLFR and KAGRSFR is: ",similarity
    
    pep.similar_smiles("KAGRSFR")
    print "The SMILES similarity is: ",pep.smiles_similarity
    
    pep.solubility_rules()
    print "{} solubility rules failed from 5".format(pep.solubility_rules_failed)
    
    pep.synthesis_rules()
    print "{} synthesis rules failed from 5".format(pep.synthesis_rules_failed)
    
    ############################################################
    # Test peptide structure functions
    ############################################################
    print
    print "### TEST PEPTIDE STRUCTURE FUNCTIONS ###"
    pdb_file="auxiliar/example_structure.pdb"
    chain="C"
    pepStr=peptide_structure(pdb_file,chain)
    print "Peptide sequence bsaed on the PDB file is: ",pepStr.sequence
    
    pepStr.get_secondary_structure()
    print "The predicted secondary structure is: ",pepStr.total_dssp
    
    pepStr.get_hydrogen_bonds()
    pepStr.plot_hydrogen_bonds()
    
    contact_threshold=4.0
    pepStr.get_heavy_atom_contacts(contact_threshold)
    
    print "The total number of contacts are: ",pepStr.total_contacts
    print "The following are the details per amino acid in the peptide:"
    print pepStr.positions
    
    ############################################################
    # Test additional functions
    ############################################################
    print
    print "### TEST ADDITIONAL FUNCTIONS ###"
    list_peptides=generate_sequences(2,"natural")
    print "The number of peptides generated in the library are: ",len(list_peptides)
    library_report=open("auxiliar/library_generated.txt","w")
    for pep in list_peptides: library_report.write("{}\n".format(pep))
    library_report.close()
     
