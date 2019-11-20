#!/usr/bin/python

"""
Pepfun: bioinformatics and cheminformatics protocols for peptide-related computational analysis
Fragment docking protocol to predict peptide bound conformations
NOTE: The protocol requires of auxiliary programas and Unix system commands - Tested on Ubuntu 16.04

From publication "Pepfun: bioinformatics and cheminformatics protocols for peptide-related computational analysis"
Journal of Cheminformatics 
Authors: Rodrigo Ochoa, Lucy Jimenez, Roman Laskowski, ..., Pilar Cossio
Year: 2019

Third-party tools required:

BioPython: https://biopython.org/wiki/Download - Ubuntu package: python-rdkit
RDKit: https://github.com/rdkit/rdkit/releases - Ubuntu package: python-biopython
Modeller: https://salilab.org/modeller/download_installation.html
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

import numpy as np
import unittest
import subprocess
import itertools
import os
import multiprocessing

# RDKit
from rdkit import Chem
from rdkit.Chem import AllChem

# BioPython
from Bio.PDB import *

# Pyrosetta
# from pyrosetta import *
# from pyrosetta.teaching import *
# from pyrosetta.toolbox import cleanATOM

# Modeller
import modeller
from modeller.automodel import *

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

################################################################

def modelling_fragments(pepT,pepM,step,receptor,number_chains):
    """
    Function to increasily model fragment of the desired peptide
    
    Arguments:
    pepT -- Template peptide sequence
    pepM -- Peptide sequence to model
    step -- Step of the modelling process during the growing of the sequence
    receptor -- Structure used to dock the peptide
    number_chains -- Number of chains in the protein structure (max 3)
    
    Return
    PDB structure of the peptide after growing flanking amino acids
    """
    
    # Copy the PDB files that are going to be used in the modelling
    os.system('cp output/step%d/model1_step%d.pdb .' %(step,step))
    if number_chains==1: os.system("sed -i 's/ A  / B  /g' model1_step%d.pdb" %step)
    if number_chains==2: os.system("sed -i 's/ A  / C  /g' model1_step%d.pdb" %step)
    os.system("cat %s.pdb model1_step%d.pdb > complex1_step%d.pdb" %(receptor,step,step))
    
    # Start the Modeller environment
    code = 'complex1_step%d' %step
    e = modeller.environ()
    m = modeller.model(e, file=code)
    aln = modeller.alignment(e)
    aln.append_model(m, align_codes=code)
    aln.write(file=code+'.seq')
   
    # Edit the information of the sequence to store the sequences in the Modeller format
    infoSeq=[x.strip() for x in open(code+'.seq')]
    header=[]
    sequenceLine=''
    for info in infoSeq:    
        if ">" not in info and ":" not in info:
            if info:
                sequenceLine+=info
        else:
            if info: header.append(info)
   
    # Store the sequences in variables according to Modeller format
    last_line=sequenceLine.split("/")
    if len(last_line)==2:
        sequenceTemp=last_line[0]+"/"+pepT+"*"
        sequenceMod=last_line[0]+"/"+pepM+"*"
    if len(last_line)==3:
        sequenceTemp=last_line[0]+"/"+last_line[1]+"/"+pepT+"*"
        sequenceMod=last_line[0]+"/"+last_line[1]+"/"+pepM+"*"

    seqTempList = [sequenceTemp[i:i+75] for i in range(0, len(sequenceTemp), 75)]
    seqModList = [sequenceMod[i:i+75] for i in range(0, len(sequenceMod), 75)]
   
    # Create the alignment file
    alignmentFile=open("alignment.ali","w")
    for h in header: alignmentFile.write(h+"\n")
    for s in seqTempList: alignmentFile.write(s+"\n")
    alignmentFile.write("\n>P1;complex1_step{}_fill\nsequence:::::::::\n".format(step))
    for s in seqModList: alignmentFile.write(s+"\n")
    alignmentFile.close()
    
    # Directories for input atom files
    e.io.atom_files_directory = ['.', '../atom_files']
    a = automodel(e, alnfile='alignment.ali', knowns='complex1_step{}'.format(step), sequence='complex1_step{}_fill'.format(step))
    a.starting_model= 1
    a.ending_model  = 1
    a.make()
   
    # Read the structure obtained from Modeller
    parser = PDBParser()
    structure = parser.get_structure('REF',"complex1_step{}_fill.B99990001.pdb".format(step))
    model = structure[0]
    if number_chains==1: ch_selected="B"
    if number_chains==2: ch_selected="C"
    out_structure=model[ch_selected]
    numbers_to_change=[]
    
    # Store the numbers of the residues that should be changed
    for num_aa,residue in enumerate(out_structure):
        numbers_to_change.append(int(residue.get_full_id()[3][1]))

    # Save a temportal PDB file
    io = PDBIO()
    io.set_structure(out_structure)
    io.save("post-modelled.pdb")
    
    # Replace the atoms and residues
    for newN,oldN in enumerate(numbers_to_change):
        diffN=len(str(oldN))-len(str(newN+1))
        if diffN==1: os.system("sed -i 's/ {} {} / A {}  /g' post-modelled.pdb".format(ch_selected,oldN,newN+1))
        if diffN==2: os.system("sed -i 's/ {} {} / A {}   /g' post-modelled.pdb".format(ch_selected,oldN,newN+1))
        if diffN==3: os.system("sed -i 's/ {} {} / A {}    /g' post-modelled.pdb".format(ch_selected,oldN,newN+1))
    
    # Delete temporal files and save result
    os.system("rm complex1_step{s}* model1_step{s}* alignment.ali".format(s=step))
    os.system("mv post-modelled.pdb pepM_step{}.pdb".format(step+1))

##########################################################################################

def generateConformer(sequence):
    """
    Function to generate a basic conformer based on the sequence
    
    Arguments:
    sequence -- Peptide fragment sequence
    
    Return:
    Peptide structure in PDB format
    """
    
    # Generate SMILES
    connect_smiles='O'
    for res in sequence:
         connect_smiles=connect_smiles[:-1]
         smiles=aminoacidSMILES(res)          
         connect_smiles=connect_smiles+smiles
    final_smiles=connect_smiles
    
    # Generate molecule from smiles
    mol = Chem.MolFromSmiles(final_smiles)
    mol.SetProp("_Name",sequence)
    print "Generating the basic conformer for peptide {}".format(sequence)
    
    # Generate the conformer with the UFF force field
    AllChem.EmbedMolecule(mol)
    AllChem.UFFOptimizeMolecule(mol)
    writer = AllChem.SDWriter("pepM_{}.sdf".format(sequence))
    writer.write(mol)
    
    # Convert to PDB file
    os.system("babel -isdf pepM_{peptide}.sdf -opdb pepM_{peptide}.pdb".format(peptide=sequence))
    os.system("rm pepM_{}.sdf".format(sequence))
    
################################################################    

def protonation(molecule,pH):
    """
    Function to protonate the molecules based on a particular pH
    
    Arguments:
    molecule -- Molecule that is going to be parameterized
    pH -- Definition of the pH to protonate
    
    Return:
    File after doing the protonation (.pqr format)
    """
    
    os.system("./scripts/pdb2pqr-linux/pdb2pqr --with-ph={} --ph-calc-method=propka --drop-water --apbs-input --ff=amber --verbose --chain --summary {mol}.pdb {mol}.pqr".format(pH,mol=molecule))
    os.system("rm {mol}.in {mol}-input.p {mol}.summary {mol}.propka".format(mol=molecule))

################################################################

def generate_box_initial(sequence,center_x,center_y,center_z,pdb):
    """
    Function to generate the initial box for the docking
    
    Arguments:
    sequence -- Sequence of the peptide fragment that will be docking
    center_x -- center in x of the defined box
    center_y -- center in y of the defined box
    center_z -- center in z of the defined box
    pdb -- Structure file
    
    Return:
    Configuration file with the coordinates
    size_x -- Initial size of the cubic box
    """
    
    # Peptide lenght
    number_amino=len(sequence)
    
    # Read the structure
    parser = PDBParser()
    structure = parser.get_structure('PEP', pdb)
    model = structure[0]
    
    # Get the difference in coordinates
    diff=model["A"][1]["CA"].coord-model["A"][number_amino]["CA"].coord
    diffValue=np.sqrt(np.sum(diff * diff))
    size_x=abs(diffValue*2.5)
    size_y=abs(diffValue*2.5)
    size_z=abs(diffValue*2.5)
    
    # Exhaustiveness. This will depend on the number of cores
    exhaustiveness=int(multiprocessing.cpu_count())
    
    # Open the config file to write the information
    config=open("config.txt","wb")
    config.write("center_x={}\n".format(center_x))
    config.write("center_y={}\n".format(center_y))
    config.write("center_z={}\n".format(center_z))
    config.write("size_x={}\n".format(size_x))
    config.write("size_y={}\n".format(size_y))
    config.write("size_z={}\n".format(size_z))
    config.write("exhaustiveness={}".format(exhaustiveness))
    config.close()
    
    # Return the cubic size
    return size_x

################################################################

def generate_box(sequence,center_x,center_y,center_z,pdb,initial_size):
    """
    Function to generate the subsequent boxes for the docking
    
    Arguments:
    sequence -- Sequence of the peptide fragment that will be docking
    center_x -- center in x of the defined box
    center_y -- center in y of the defined box
    center_z -- center in z of the defined box
    pdb -- Structure file
    initial_size -- size of the previous box
    
    Return:
    Configuration file with the coordinates
    """
    
    # Peptide lenght
    number_amino=len(sequence)
    
    # Read the pdb structure
    parser = PDBParser()
    structure = parser.get_structure('PEP', pdb)
    model = structure[0]
    
    # Get the difference in coordinates    
    diff_x=model["A"][1]["CA"].coord[0]-model["A"][number_amino]["CA"].coord[0]
    diff_y=model["A"][1]["CA"].coord[1]-model["A"][number_amino]["CA"].coord[1]
    diff_z=model["A"][1]["CA"].coord[2]-model["A"][number_amino]["CA"].coord[2]
    
    # Assign the sizes based on the growing direction of the peptide
    if abs(diff_x*2.5)>initial_size:
        size_x=abs(diff_x*2.5)
    else:
        size_x=initial_size
        
    if abs(diff_y*2.5)>initial_size:
        size_y=abs(diff_y*2.5)
    else:
        size_y=initial_size
    
    if abs(diff_z*2.5)>initial_size:
        size_z=abs(diff_z*2.5)
    else:
        size_z=initial_size
    
    # Exhaustiveness. This will depend on the number of cores
    exhaustiveness=int(multiprocessing.cpu_count())
    
    # Create the configuration file
    config=open("config.txt","wb")
    config.write("center_x={}\n".format(center_x))
    config.write("center_y={}\n".format(center_y))
    config.write("center_z={}\n".format(center_z))
    config.write("size_x={}\n".format(size_x))
    config.write("size_y={}\n".format(size_y))
    config.write("size_z={}\n".format(size_z))
    config.write("exhaustiveness={}".format(exhaustiveness))
    config.close()

################################################################

def docking(receptor,peptide,sequence,stepNumber):
    """
    Function to execute the docking
    
    Arguments:
    receptor -- Protein structure
    peptide -- Fragment that will be docked
    sequence -- Sequence of the peptide
    stepNumber -- Number of the step during the fragment growing
    """
    
    # List of amino acids
    aminoacids={"A":"ALA","D":"ASP","E":"GLU","F":"PHE","H":"HIS","I":"ILE","K":"LYS","L":"LEU","M":"MET","G":"GLY",
                "N":"ASN","P":"PRO","Q":"GLN","R":"ARG","S":"SER","T":"THR","V":"VAL","W":"TRP","Y":"TYR","C":"CYS"}
    
    # Check if this is the frist step
    if stepNumber==0:
        # Convert the molecules in PDBQT format
        os.system("./scripts/pythonsh scripts/prepare_receptor4.py -r {rec}.pqr -o {rec}.pdbqt -C -U waters".format(rec=receptor))
        os.system("./scripts/pythonsh scripts/prepare_ligand4.py -l {pep}.pqr  -C -U '' -B -o {pep}.pdbqt".format(pep=peptide))
        
        # Run Autodock Vina- NOTE: Install using: sudo apt-get install autodock-vina
        print "Docking step number {} ...".format(stepNumber)
        os.system("vina --receptor {}.pdbqt --ligand {}.pdbqt --log score.log --out out.pdbqt --config config.txt".format(receptor,peptide))
        
        # Split and get the first model
        os.system("csplit out.pdbqt /MODEL/ {{*}}; mv xx01 model1_step{}.pdbqt".format(stepNumber))
        os.system("rm xx* score.log out.pdbqt")
        
        # Count the number of models
        bash="ls | grep '_step{}.pdbqt' | grep 'model' | wc -l".format(stepNumber)
        numberModels = subprocess.check_output(['bash','-c', bash])
        return int(numberModels)
        
    else:
        # Select the amino acids that will be assigned as rigid
        rigidAA=[]
        for i in range(1,len(sequence)+1):
            if i!=1 and i!=2 and i!=len(sequence)-1 and i!=len(sequence):
                rigidAA.append(i)
        
        # Assign the rigid atoms
        rigidAtoms=[]
        for rigid in rigidAA:
            if rigid <=9:
                os.system("grep {} {}.pqr | grep 'A   {}' | awk '{{print $2}}' > atom.temp".format(aminoacids[sequence[rigid-1]],peptide,rigid))
            else:
                os.system("grep {} {}.pqr | grep 'A  {}' | awk '{{print $2}}' > atom.temp".format(aminoacids[sequence[rigid-1]],peptide,rigid))
            num=[x.strip() for x in open("atom.temp")]
            rigidAtoms=rigidAtoms+num
            os.system("rm atom.temp")
        
        # Get the atoms that will be inactivated
        inactive=[]
        os.system("./scripts/pythonsh scripts/prepare_ligand4.py -s -l {pep}.pqr -C -U '' -B -o {pep}.pdbqt".format(pep=peptide))
        os.system("grep REMARK {}.pdbqt | grep '  A  ' | awk '{{print $6\"\t\"$8}}' | sed 's/_/ /g' | awk '{{print $2\"_\"$4}}' > active.temp".format(peptide))
        bonds=[x.strip() for x in open("active.temp")]
        for b in bonds:
            data=b.split("_")
            if data[0] in rigidAtoms or data[1] in rigidAtoms:
                inactive.append("{}_{}".format(str(int(data[0])-1),str(int(data[1])-1)))
        os.system("rm active.temp")
        
        commandInactive="_".join(inactive)
        
        # Prepare the ligand based on the inactivated atoms
        os.system("./scripts/pythonsh scripts/prepare_ligand4.py -s -l {}.pqr -C -U '' -B -o {}.pdbqt -I {}".format(peptide,peptide,commandInactive))
        
        # Run Autodock Vina
        print "Docking step number {} ...".format(stepNumber)
        os.system("vina --receptor {}.pdbqt --ligand {}.pdbqt --log score.log --out out.pdbqt --config config.txt".format(receptor,peptide))
        
        # Get the first model
        os.system("csplit out.pdbqt /MODEL/ {{*}}; mv xx01 model1_step{}.pdbqt".format(stepNumber))
        os.system("rm xx* score.log out.pdbqt")
        
        # Count number models
        bash="ls | grep '_step{}.pdbqt' | grep 'model' | wc -l".format(stepNumber)
        numberModels = subprocess.check_output(['bash','-c', bash])
        return int(numberModels)

################################################################

def new_coordinates(base,model_pdbqt):
    """
    Function to annotate the coordinates from the pqr file
    
    Arguments:
    base -- Structure protonated
    model_pdbqt -- Model requiring the information
    
    Return:
    PDB file with the new coordinates
    """
    
    # Read the structure
    print "Processing protein %s ..." %model_pdbqt
    parser = PDBParser()
    structure = parser.get_structure('PEP', "%s.pqr" %base)
    model = structure[0]
    counter=1
    batoms=["C","N","O","H","S"]
    
    # Iterate over the residues
    for residue in model["A"]:
        resName=residue.get_resname()
        resNumber=residue.get_full_id()[3][1]
        
        for atom in residue:
            idAtom = atom.get_id()
            # Check the spaces of the PDB file based on the residue number
            if resNumber<=9:
                if idAtom in batoms:
                    bash="grep {} {}.pdbqt | grep 'A   {}' | grep '  {}  ' | awk '{{print $7}}'".format(resName,model_pdbqt,str(resNumber),idAtom)
                    new_x = subprocess.check_output(['bash','-c', bash])
                    bash="grep {} {}.pdbqt | grep 'A   {}' | grep '  {}  ' | awk '{{print $8}}'".format(resName,model_pdbqt,str(resNumber),idAtom)
                    new_y = subprocess.check_output(['bash','-c', bash])
                    bash="grep {} {}.pdbqt | grep 'A   {}' | grep '  {}  ' | awk '{{print $9}}'".format(resName,model_pdbqt,str(resNumber),idAtom)
                    new_z = subprocess.check_output(['bash','-c', bash])
                    atom.set_coord([float(new_x),float(new_y),float(new_z)])
                else:
                    if len(idAtom)==4: idAtom=idAtom[-1]+idAtom[:3]
                        
                    bash="grep {} {}.pdbqt | grep 'A   {}' | grep ' {} ' | awk '{{print $7}}'".format(resName,model_pdbqt,str(resNumber),idAtom)
                    new_x = subprocess.check_output(['bash','-c', bash])
                    bash="grep {} {}.pdbqt | grep 'A   {}' | grep ' {} ' | awk '{{print $8}}'".format(resName,model_pdbqt,str(resNumber),idAtom)
                    new_y = subprocess.check_output(['bash','-c', bash])
                    bash="grep {} {}.pdbqt | grep 'A   {}' | grep ' {} ' | awk '{{print $9}}'".format(resName,model_pdbqt,str(resNumber),idAtom)
                    new_z = subprocess.check_output(['bash','-c', bash])
                    atom.set_coord([float(new_x),float(new_y),float(new_z)])
            else:
                if idAtom in batoms:
                    bash="grep {} {}.pdbqt | grep 'A  {}' | grep '  {}  ' | awk '{{print $7}}'".format(resName,model_pdbqt,str(resNumber),idAtom)
                    new_x = subprocess.check_output(['bash','-c', bash])
                    bash="grep {} {}.pdbqt | grep 'A  {}' | grep '  {}  ' | awk '{{print $8}}'".format(resName,model_pdbqt,str(resNumber),idAtom)
                    new_y = subprocess.check_output(['bash','-c', bash])
                    bash="grep {} {}.pdbqt | grep 'A  {}' | grep '  {}  ' | awk '{{print $9}}'".format(resName,model_pdbqt,str(resNumber),idAtom)
                    new_z = subprocess.check_output(['bash','-c', bash])
                    atom.set_coord([float(new_x),float(new_y),float(new_z)])
                else:
                    if len(idAtom)==4: idAtom=idAtom[-1]+idAtom[:3]
                        
                    bash="grep {} {}.pdbqt | grep 'A  {}' | grep ' {} ' | awk '{{print $7}}'".format(resName,model_pdbqt,str(resNumber),idAtom)
                    new_x = subprocess.check_output(['bash','-c', bash])
                    bash="grep {} {}.pdbqt | grep 'A  {}' | grep ' {} ' | awk '{{print $8}}'".format(resName,model_pdbqt,str(resNumber),idAtom)
                    new_y = subprocess.check_output(['bash','-c', bash])
                    bash="grep {} {}.pdbqt | grep 'A  {}' | grep ' {} ' | awk '{{print $9}}'".format(resName,model_pdbqt,str(resNumber),idAtom)
                    new_z = subprocess.check_output(['bash','-c', bash])
                    atom.set_coord([float(new_x),float(new_y),float(new_z)])
            
            counter+=1
    
    # Save the final PDB structure
    io = PDBIO()
    io.set_structure(structure)
    io.save('%s.pdb' %model_pdbqt)
    
########################################################################################
########################################################################################
########################################################################################
# Main execution
########################################################################################
########################################################################################
########################################################################################
if __name__ == '__main__':
    
    # # Script arguments
    # parser = argparse.ArgumentParser(description='Fragment docking protocol to predict peptide bound conformations')
    # parser.add_argument('-s', dest='pep_seq', action='store',required=True,
    #                     help='Sequence of peptide to be analyzed')
    # parser.add_argument('-p', dest='pep_ph', action='store', default=5.3,
    #                     help='pH of the system')
    # # Pending add more
    # args = parser.parse_args()
    
    ####################################################################################
    # Assignment of parameters
    # sequence=args.pep_seq
    
    ####################################################################################
    # Input data
    ####################################################################################
    ref_peptide="ENPVVHFFKNIVTPR"
    peptide="FFK"
    receptor="1BX2-ABP"
    number_chains=2
    pH=5.3
    # Center of reference selected by the user
    center_x=-5.256
    center_y=26.343
    center_z=-2.416
    
    ####################################################################################
    # Start functions
    ####################################################################################
    
    # Predict peptide conformer
    generateConformer(peptide)
     
    # Protonation
    protonation(receptor,pH)
    protonation("pepM_{}".format(peptide),pH)
    
    # Docking
    initial_size=generate_box_initial(peptide,center_x,center_y,center_z,"pepM_{}.pdb".format(peptide))
    num_models=docking(receptor,"pepM_{}".format(peptide),peptide,0)
    new_coordinates("pepM_{}".format(peptide),"model1_step0")
    
    # Move to a folder with step0
    os.system("mkdir output/step0; mv pepM_{}.* *_step0.pdbqt *_step0.pdb config.txt output/step0".format(peptide))
     
    # list of future fragments for growing the peptide
    list_fragments=[]
    ref_limit=ref_peptide.index(peptide)
    count=0
    limit=0
    while limit==0:
        count+=1
        low_limit=ref_limit-count
        if low_limit<0: low_limit=0
        up_limit=ref_limit+len(peptide)+count
        if up_limit>len(ref_peptide): up_limit=len(ref_peptide)
        list_fragments.append(ref_peptide[low_limit:up_limit])
        if low_limit==0 and up_limit==len(ref_peptide):limit=1
    
    # Print the list of fragments
    print list_fragments
    
    # Iterate over the fragments
    final_step=0
    pep_reference=peptide
    for i,frag in enumerate(list_fragments):
        # Get the index of the fragment in the sequence
        index_sub=frag.index(pep_reference)
        pepTemplate=""
        control_size=0
        for pos,aa in enumerate(frag):
            if pos<index_sub:
                pepTemplate=pepTemplate+"-"
            else:
                control_size+=1
                if control_size > len(pep_reference):
                    pepTemplate=pepTemplate+"-"
                else:
                    pepTemplate=pepTemplate+aa
        pep_reference=frag
        
        # Model the fragment
        modelling_fragments(pepTemplate,frag,i,receptor,number_chains)
        protonation("pepM_step{}".format(i+1),pH)
        
        # Generate the box
        generate_box(frag,center_x,center_y,center_z,"pepM_step{}.pdb".format(i+1),initial_size)
        num_models=docking(receptor,"pepM_step{}".format(i+1),frag,i+1)
        
        # Check if no models were produced to exit
        if num_models==0:
           print "Docking error ... exiting the protocol"
           break
        
        # Update the coordinates and save the results
        new_coordinates("pepM_step{}".format(i+1),"model1_step{}".format(i+1))
        os.system("mkdir output/step{num_step}; mv pepM_step{num_step}.* *_step{num_step}.pdb *_step{num_step}.pdbqt *_step{num_step}.pqr config.txt output/step{num_step}".format(num_step=str(i+1)))
        final_step=i+1
        
    # Save the final complex
    os.system('cp output/step{final}/model1_step{final}.pdb .'.format(final=final_step))
    if number_chains==1: os.system("sed -i 's/ A  / B  /g' model1_step{}.pdb".format(final_step))
    if number_chains==2: os.system("sed -i 's/ A  / C  /g' model1_step{}.pdb".format(final_step))
    os.system("cat {}.pdb model1_step{}.pdb > output/final_complex_{}.pdb".format(receptor,final_step,ref_peptide))
    os.system("rm model1_step{}.pdb".format(final_step))
