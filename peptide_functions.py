#!/usr/bin/python

"""
Pepfun: bioinformatics and cheminformatics protocols for peptide-related computational analysis

From publication "Pepfun: bioinformatics and cheminformatics protocols for peptide-related computational analysis"
Journal of Cheminformatics 
Authors: Rodrigo Ochoa, Lucy Jimenez, Roman Laskowski, ..., Pilar Cossio
Year: 2019

Third-party tools required:

BioPython: https://biopython.org/wiki/Download
RDKit: https://github.com/rdkit/rdkit/releases
"""

########################################################################################
# Authorship
########################################################################################

__author__ = "Rodrigo Ochoa"
__credits__ = ["Rodrigo Ochoa", "Lucy Jimenez","Roman Laskowski", "...", "Pilar Cossio"]
__license__ = "GPL"
__version__ = "1.0"
__email__ = "rodrigo.ochoa@udea.edu.co"

########################################################################################
# Modules to import
########################################################################################

# External modules
import math
import itertools
import subprocess
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
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio import SeqIO

# RDKit
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import DataStructs
from rdkit.Chem import Crippen
from rdkit.Chem import Lipinski
from rdkit.Chem import Descriptors

# # PyRosetta
# from pyrosetta import *
# from pyrosetta.teaching import *
# init()

# Local scripts
# import conformers

# PeptideMatch
#import swagger_client 
#from swagger_client.rest import ApiException
#from pprint import pprint

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
            f='.'.join(map(str, i))
            frags.append('.'.join(map(str, i)))
    elif aa_type=="d_aminoacids":
        # Generate the combinatorial library
        for i in itertools.product(d_aminoacids, repeat=long_peptide):
            f='.'.join(map(str, i))
            frags.append('.'.join(map(str, i)))
    elif aa_type=="combined":
        # Generate the combinatorial library
        for i in itertools.product(combined, repeat=long_peptide):
            f='.'.join(map(str, i))
            if "[" in f: frags.append('.'.join(map(str, i)))
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
    
    # PENDIENTE
    # Calcular protein scales y dar la opcion de compararlas segun la ventana
    # Adaptar funciones de PepLibGen - Synth Rules
    # Adaptar la generacion de secuencias aleatorias
    # Agregar la helical wheel
    
    def __init__(self,sequence):
        self.sequence=sequence
        self.pH=7
        self.solubility_rules_failed=0
        self.length_peptide=len(self.sequence)
        self.synthesis_rules_failed=0
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
    def pairwise_alignment(self,peptide_to_match):
        """
        """
        alignment = pairwise2.align.globalxx(self.sequence,peptide_to_match)
        print alignment
        # TODO: Fix the response
    
    ############################################################################
    def blast_online(self):
        """
        """
        # Create a temporal fasta file with the sequence
        fasta_file=open("%s.fasta" %self.sequence,"wb")
        fasta_file.write(">sequence\n%s" %self.sequence)
        fasta_file.close()
        
        record = SeqIO.read("%s.fasta" %self.sequence, format="fasta")
        result_handle = NCBIWWW.qblast("blastp", "nr", record.format("fasta"), word_size=2, expect=20000.0, matrix_name="PAM30", gapcosts="9 1")
        os.system("rm %s.fasta" %self.sequence)
        b_record = NCBIXML.read(result_handle)
        
        # Parse the results
        for alignment in b_record.alignments:
            for hsp in alignment.hsps:
                #if hsp.expect < E_VALUE_THRESH:
                print hsp.identities
                print hsp.positives
                print hsp.gaps
                print hsp.align_length
                print hsp.strand
                print hsp.query_start
                print '****Alignment****' 
                print 'sequence: %s' % alignment.title
                print 'length: %i' % alignment.length
                print 'e value: %f' % hsp.expect
                print hsp.query[0:75] + '...'
                print hsp.match[0:75] + '...'
                print hsp.sbjct[0:75] + '...'
                break
        #TODO: Fix the response with percentages
    ############################################################################
    def compute_peptide_charges(self,pH_internal=7):
        """
        Function to calculate the average net charge based on pka values
        
        Arguments:
        Sequence - amino acid sequence of the peptide
        pH - By default is 7
        
        Return:
        The net charge based on pka values recorded from ...
        PENDING - Get a plot of charge vs time
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
        
        Arguments:
        Sequence - amino acid sequence of the peptide
        
        Return:
        Static physico-chemical properties: molecular weigth, crippen logP, number of hydrogen bond acceptors and donors
        SMILES - 2D representation of the molecule
        """
        
        # Generate molecule from sequence
        p = subprocess.Popen(['java', '-jar','auxiliar/seq_smiles_peptides.jar','-i',self.sequence], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        self.smiles, err = p.communicate()
        mol = Chem.MolFromSmiles(self.smiles)
        mol.SetProp("_Name",self.sequence)
        
        # # Create an rdkit object with the peptide sequence using HELM functionality
        # pep=".".join(list(self.sequence))
        # mol = Chem.MolFromHELM("PEPTIDE1{%s}$$$$" %pep)
        
        # Generate the smiles
        #self.smiles=Chem.MolToSmiles(mol)
        self.num_hdonors = Lipinski.NumHDonors(mol)
        self.num_hacceptors = Lipinski.NumHAcceptors(mol)
        self.mol_weight = Descriptors.MolWt(mol)
        self.mol_logp = Crippen.MolLogP(mol)
    
    ############################################################################    
    def similar_smiles(self,peptide_to_match):
        """
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
        
        # # Create an rdkit object with the peptide sequence using HELM functionality
        # pep1=".".join(list(self.sequence))
        # pep2=".".join(list(peptide_to_match))
        # mol1 = Chem.MolFromHELM("PEPTIDE1{%s}$$$$" %pep1)
        # mol2 = Chem.MolFromHELM("PEPTIDE1{%s}$$$$" %pep2)
        # # Add the sequence as a name in the mol object
        # sequence1=Chem.MolToSequence(mol1)
        # mol1.SetProp("_Name",sequence1)
        # sequence2=Chem.MolToSequence(mol2)
        # mol2.SetProp("_Name",sequence2)
        
        # Calculate the fingerprints and the similarity
        fp1=AllChem.GetMorganFingerprintAsBitVect(mol1,2,2048)
        fp2=AllChem.GetMorganFingerprintAsBitVect(mol2,2,2048)
        
        self.smiles_similarity=DataStructs.TanimotoSimilarity(fp1,fp2)
    
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
        
        Arguments:
        sequence - amino acid sequence of the peptide
        
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
        Function to check some synthesis rules based on previous recommendations
        
        Arguments:
        Sequence - amino acid sequence of the peptide
        
        Return:
        
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
    def search_matches(self):
        """
        Function to check some synthesis rules based on previous recommendations
        
        Arguments:
        Sequence - amino acid sequence of the peptide
        
        Return:
        
        """
        # Presence of forbiden motifs
        forbidden_motifs = {'2-prolines':r'[P]{3,}','DG-DP':r'D[GP]','N-Q-Nterminal':r'^[NQ]',}
        for motif in forbidden_motifs:
            if re.search(forbidden_motifs[motif],self.sequence):
                self.synthesis_rules_failed+=1
    
    ############################################################################
    def helical_wheel(self):
        """
        """
        
        # Create a temporal fasta file with the sequence
        fasta_file=open("%s.fasta" %self.sequence,"wb")
        fasta_file.write(">sequence\n%s" %self.sequence)
        fasta_file.close()
        
        # Call the program
        os.system("auxiliar/pepwheel -sequence %s.fasta -graph png" %self.sequence)
        os.system("mv pepwheel.1.png {peptide}_wheel.png; rm {peptide}.fasta".format(peptide=self.sequence))
        
    
    ############################################################################
    def generate_conformer(self):
        """
        """
        # Generate molecule from smiles
        mol = Chem.MolFromSmiles(self.smiles)
        mol.SetProp("_Name",self.sequence)
        
        print "Generating the basic conformer for peptide %s" %self.sequence
        
        AllChem.EmbedMolecule(mol)
        AllChem.UFFOptimizeMolecule(mol)
        writer = AllChem.SDWriter("structure_%s.sdf" %self.sequence)
        writer.write(mol)
        
        # Convert to PDB file
        os.system("babel -isdf structure_{peptide}.sdf -opdb structure_{peptide}.pdb".format(peptide=self.sequence))
        
        # Generation mmff - for more than 7 rotatable bonds  
        #engine = conformers.ConformerGenerator(pool_multiplier=200, rmsd_threshold=0.3, force_field='mmff')
        # engine = conformers.ConformerGenerator(rmsd_threshold=0.3, force_field='mmff')
        # mol2 = engine.generate_conformers(mol) 
        # pdb_structure=open("structure_%s.mol" %self.sequence,"wb")
        # pdb_structure.write(Chem.MolToMolBlock(mol2))
        # pdb_structure.close()
        
    ################################################################
    
    def generate_structure_modpep(self):
        """
        """
        # Create a temporal fasta file with the sequence
        fasta_file=open("%s.fasta" %self.sequence,"wb")
        fasta_file.write(">sequence\n%s" %self.sequence)
        fasta_file.close()
        
        os.system("auxiliar/modpep {peptide}.fasta structure_{peptide}.pdb -n 1 -L auxiliar/".format(peptide=self.sequence))
        os.system("sed -i 's# A  # C  #g' structure_{peptide}.pdb; rm {peptide}.fasta".format(peptide=self.sequence))
                
    ############################################################################
    

class peptide_structure:
    
    # Initializer
    def __init__(self,pdb_file,chain):
        self.pdb_file=pdb_file
        self.chain=chain
        self.aminoacids_back={"ALA":"A","ASP":"D","GLU":"E","PHE":"F","HIS":"H","ILE":"I","LYS":"K","LEU":"L","MET":"M","GLY":"G",
                              "ASN":"N","PRO":"P","GLN":"Q","ARG":"R","SER":"S","THR":"T","VAL":"V","TRP":"W","TYR":"Y","CYS":"C"}
        self.aminoacids={"A":"ALA","D":"ASP","E":"GLU","F":"PHE","H":"HIS","I":"ILE","K":"LYS","L":"LEU","M":"MET","G":"GLY",
                         "N":"ASN","P":"PRO","Q":"GLN","R":"ARG","S":"SER","T":"THR","V":"VAL","W":"TRP","Y":"TYR","C":"CYS"}
        self.sequence=""
        self.initial_id_residue=0
        self.positions={}
        
        # Read the structure
        parser = PDBParser()
        self.reference = parser.get_structure('REF',pdb_file)
        for ch in self.reference[0]:
            if ch.get_id()==chain:
                for i,residue in enumerate(ch):
                    seq=self.aminoacids_back[residue.get_resname()]
                    self.sequence=self.sequence+str(seq)
                    if i==0: self.initial_id_residue=residue.get_full_id()[3][1]
                    self.positions[i+1]={"aa":str(seq),"full_aa":residue.get_resname()}
        
        self.len_sequence=len(self.sequence)

        
    ############################################################################
    
    # Function to calculate the secondary structure
    def get_secondary_structure(self):
        """
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
    
    # Function to calculate hydrogen bonds with the other chains
    def get_hydrogen_bonds(self):
        """
        """
        model=self.reference[0]
        dssp = DSSP(model,self.pdb_file,dssp='auxiliar/mkdssp')
        self.total_dssp=""
        
        # Loop over the keys from the dssp response to store ss and asa values
        total_index={}
        list_chains=[]
        reference_position=0
        list_chains.append(list(dssp.keys())[0][0])
        hbonds_peptide={}
        for keys in list(dssp.keys()):
            #print keys
            # if keys[0] not in list_chains:
            #     list_chains.append(keys[0])
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
                #print position
                #print dssp[keys]
                amino_name=dssp[keys][1]+str(keys[1][1])
                if amino_name not in hbonds_peptide: hbonds_peptide[amino_name]=[]
                if dssp[keys][6]<-5:
                    interactions_pos=last_position+dssp[keys][6]
                    #print interactions_pos,total_index[interactions_pos]
                    hbonds_peptide[amino_name].append(total_index[interactions_pos])
                if dssp[keys][8]<-5:
                    interactions_pos=last_position+dssp[keys][8]
                    #print interactions_pos,total_index[interactions_pos]
                    hbonds_peptide[amino_name].append(total_index[interactions_pos])
                if dssp[keys][10]<-5:
                    interactions_pos=last_position+dssp[keys][10]
                    #print interactions_pos,total_index[interactions_pos]
                    hbonds_peptide[amino_name].append(total_index[interactions_pos])
                if dssp[keys][12]<-5:
                    interactions_pos=last_position+dssp[keys][12]
                    #print interactions_pos,total_index[interactions_pos]
                    hbonds_peptide[amino_name].append(total_index[interactions_pos])
        print hbonds_peptide
        
        # Generate fragments
        fragment_chains={}
        for amino in hbonds_peptide:
            for receptor in hbonds_peptide[amino]:
                chain=receptor[0]
                if chain not in fragment_chains: fragment_chains[chain]=[]
                if receptor not in fragment_chains[chain]:
                    fragment_chains[chain].append(receptor)
        
        for ch in fragment_chains:
            fragment_chains[ch].sort(key=lambda x: x[2])
        
        # Generate string
        fragments_final=[]
        for ch in fragment_chains:
            for i,element in enumerate(fragment_chains[ch]):
                if i==0:
                    frag=element[1]
                    pos=element[2]
                else:
                    if pos+1==element[2]:
                        frag+=element[1]
                        pos=element[2]
                    else:
                        fragments_final.append(frag)
                        frag=element[1]
                        pos=element[2]
            fragments_final.append(frag)
        
        print fragment_chains
        print fragments_final
        helm_string=""
        for number,fr in enumerate(fragments_final):
            helm_string+="PEPTIDE%d{%s}|" %(number+1,".".join(fr))
        print helm_string[:-1]
        
        ###################################################
        # Get graph order
        names=[0]*len(hbonds_peptide)
        pep_interactions=[]
        chain_nodes=[]
        chain_colors={"PEP":"red"}
        edge_colors={"5":"black","4":"orange","8":"orange"}
        
        list_colors=["cyan","green","magenta","blue"]
        type_interactions=[]
        
        for i,amino in enumerate(hbonds_peptide):
            pos=int(amino[1:])
            names[pos-1]=amino
            chain_nodes.append("PEP")
            
            if i!=len(hbonds_peptide)-1:
                pep_interactions.append((i,i+1))
                type_interactions.append(5)
        
        reference=len(hbonds_peptide)
        for num_chain,chains in enumerate(fragment_chains):
            chain_colors[chains]=list_colors[num_chain]
            for count,residues in enumerate(fragment_chains[chains]):
                chain_nodes.append(chains)
                names.append(residues[0]+"-"+residues[1]+str(residues[2]))
                number_res=reference+count
                
                for pep_amino in hbonds_peptide:
                    if residues in hbonds_peptide[pep_amino]:
                        pep_interactions.append((names.index(pep_amino),number_res))
                        type_interactions.append(hbonds_peptide[pep_amino].count(residues)*4)
                        
            reference=len(names)
            
        print len(names)
        print len(pep_interactions)
        print type_interactions
        print len(chain_nodes)
        print chain_colors
        
        
        
        g = Graph(pep_interactions)
        g.vs["name"] = names
        g.vs["chains"] = chain_nodes
        g.es["count_int"] = type_interactions
        color_per_chain=[chain_colors[chColor] for chColor in g.vs["chains"]]
        color_per_edge=[edge_colors[str(edColor)] for edColor in type_interactions]
        print color_per_chain
        
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
        
        plot(g, **visual_style)
        ###################################################
            
        
    ############################################################################
    
    # Function to calculate the secondary structure with PROSS
    def get_secondary_structure_PROSS(self):
        """
        """
        bash="python auxiliar/PROSS.py %s > output" %self.pdb_file
        os.system("%s" %bash)
        
        os.system("csplit -s output /Chain:/ {*}; tail -n +3 xx03 > pp2.txt; rm xx*")
        structure=[x.strip() for x in open("pp2.txt")]
        
        # Store a string with the full prediction of PROSS
        self.total_pross=""
        
        for s in structure:
            info=s.split()
            res=info[1]
            resPos=info[0]
            if info[2]=="P":
                self.positions[int(resPos)]["pross"]="P"
                self.total_pross=self.total_pross+"P"
            else:
                if info[3] in ("Dk","Dl","Ek","El"):
                    self.positions[int(resPos)]["pross"]="P"
                    self.total_pross=self.total_pross+"P"
                else:
                    self.positions[int(resPos)]["pross"]=info[2]
                    self.total_pross=self.total_pross+info[2]
        
        os.system("rm output pp2.txt")
    
    ############################################################################
    
    # Function to count heavy atom contacts per residue
    def get_heavy_atom_contacts(self,contact_threshold):
        """
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
    
    ############################################################################
    
    def calculate_dihedrals(self):
        """
        """
        aaChi2=["R","N","D","Q","H","I","L","M","F","P","Y","E","K"]
        
        pose = pose_from_pdb(self.pdb_file)
        
        # Loop over the residues of the chain
        position=1
        for i in range(self.initial_id_residue,self.len_sequence+1):
            residue=pose.pdb_info().pdb2pose(self.chain, i)
            # Backbone dihedrals
            self.positions[position]["phi"]=pose.phi(residue)
            self.positions[position]["psi"]=pose.psi(residue)
            # Print Chi1
            self.positions[position]["chi1"]=pose.chi(1,residue)
            # Print Chi2
            if self.positions[position]["aa"] in aaChi2:
                self.positions[position]["chi2"]=pose.chi(2,residue)
            else:
                self.positions[position]["chi2"]="-"
            position+=1
    


    ################################################################
    
    def peptide_minimization(self):
        """
        """
        pose = pose_from_pdb(self.pdb_file)
         
        scorefxn = get_fa_scorefxn()
        mm = MoveMap()
        mm.set_bb(True)
        mm.set_chi(True) ## For side chain minimization
        #mm.set_bb_true_range(1,4) ## Here range means 1 to 4 or the minimzation is applied only to residues 1 and4 ??
        
        minmover = MinMover(mm, scorefxn, 'dfpmin', 0.00001, True) ## I don't know the meaning of dfpmin,10,True. I saw it somewhere and used it
        minmover.apply(pose)
        
        pose.dump_pdb("%s_minimized.pdb" %self.sequence)
        os.system("csplit {peptide}_minimized.pdb /All/; mv xx00 {peptide}_minimized.pdb; rm xx01".format(peptide=self.sequence))


############################################################################
############################################################################
############################################################################
    
if __name__ == '__main__':
    
    pep=peptide_sequence("KMGRLFR")
    print "Sequence: ",pep.sequence
    pep.compute_peptide_charges()
    print "Net charge at pH 7: ",pep.netCharge
    pep.calculate_properties_from_mol()
    print "Molecular weight: ",pep.mol_weight
    pep.calculate_properties_from_sequence()
    print "Average Hydrophobicity: ",pep.avg_hydro
    print "Isoelectric Point: ",pep.isoelectric_point
    
    matrix=matlist.structure
    pep.align_position_matrix("KAGRSFR",matrix)
    print "Alignment pos_by_pos score: ",pep.score_pos_pos
    pep.solubility_rules()
    print "%d solubility rules failed from 5" %pep.solubility_rules_failed
    pep.synthesis_rules()
    print "%d synthesis rules failed from 5" %pep.synthesis_rules_failed
    
    pep.align_position_local("KAGRSFR")
    print "The number of dismatches are %d" %pep.dismatch
    
    # list_peptides=generate_sequences(2,"natural")
    # print len(list_peptides)
