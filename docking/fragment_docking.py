#!/usr/bin/python2.7

# Import modules
import numpy as np
import unittest
import subprocess
import itertools
import os
import multiprocessing

# Chemical modules
from rdkit import Chem
from rdkit.Chem import AllChem
from Bio.PDB import *

# Pyrosetta
from pyrosetta import *
from pyrosetta.teaching import *
from pyrosetta.toolbox import cleanATOM

# Modeller
import modeller
from modeller.automodel import *

def aminoacidSMILES(amino):
    """
    """
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
    return aminoacids[amino]['SMILES']

################################################################

def modelling_fragments(pepT,pepM,step,receptor,number_chains):
    """
    """
    
    # Get the sequence of the 1qg8 PDB file, and write to an alignment file
    os.system('cp step%d/model1_step%d.pdb .' %(step,step))
    if number_chains==1: os.system("sed -i 's/ A  / B  /g' model1_step%d.pdb" %step)
    if number_chains==2: os.system("sed -i 's/ A  / C  /g' model1_step%d.pdb" %step)
    os.system("cat %s.pdb model1_step%d.pdb > complex1_step%d.pdb" %(receptor,step,step))
   
    code = 'complex1_step%d' %step
    e = modeller.environ()
    m = modeller.model(e, file=code)
    aln = modeller.alignment(e)
    aln.append_model(m, align_codes=code)
    aln.write(file=code+'.seq')
   
    infoSeq=[x.strip() for x in open(code+'.seq')]
    header=[]
    sequenceLine=''
    for info in infoSeq:    
        if ">" not in info and ":" not in info:
            if info:
                sequenceLine+=info
        else:
            if info: header.append(info)
   
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
    alignmentFile=open("alignment.ali","wb")
    for h in header: alignmentFile.write(h+"\n")
    for s in seqTempList: alignmentFile.write(s+"\n")
    alignmentFile.write("\n>P1;complex1_step%d_fill\nsequence:::::::::\n" %step)
    for s in seqModList: alignmentFile.write(s+"\n")
    alignmentFile.close()
    # # directories for input atom files
    e.io.atom_files_directory = ['.', '../atom_files']
    a = automodel(e, alnfile='alignment.ali', knowns='complex1_step%d' %step, sequence='complex1_step%d_fill' %step)
    a.starting_model= 1
    a.ending_model  = 1
    a.make()
   
    
    parser = PDBParser()
    structure = parser.get_structure('REF',"complex1_step%d_fill.B99990001.pdb" %step)
    model = structure[0]
    if number_chains==1: ch_selected="B"
    if number_chains==2: ch_selected="C"
    
    out_structure=model[ch_selected]
    numbers_to_change=[]
    
    for num_aa,residue in enumerate(out_structure):
        numbers_to_change.append(int(residue.get_full_id()[3][1]))

    io = PDBIO()
    io.set_structure(out_structure)
    io.save("post-modelled.pdb")
    

    for newN,oldN in enumerate(numbers_to_change):
        diffN=len(str(oldN))-len(str(newN+1))
        if diffN==1: os.system("sed -i 's/ %s %d / A %d  /g' post-modelled.pdb" %(ch_selected,oldN,newN+1))
        if diffN==2: os.system("sed -i 's/ %s %d / A %d   /g' post-modelled.pdb" %(ch_selected,oldN,newN+1))
        if diffN==3: os.system("sed -i 's/ %s %d / A %d    /g' post-modelled.pdb" %(ch_selected,oldN,newN+1))
    
    os.system("rm complex1_step%d* model1_step%d* alignment.ali" %(step,step))
    os.system("mv post-modelled.pdb pepM_step%d.pdb" %(step+1))

##########################################################################################

def generateConformer(sequence):
    """
    """
    
    connect_smiles='O'
    for res in sequence:
         connect_smiles=connect_smiles[:-1]
         smiles=aminoacidSMILES(res)          
         connect_smiles=connect_smiles+smiles
    final_smiles=connect_smiles
    # Generate molecule from smiles
    mol = Chem.MolFromSmiles(final_smiles)
    mol.SetProp("_Name",sequence)
    
    print "Generating the basic conformer for peptide %s" %sequence
    
    AllChem.EmbedMolecule(mol)
    AllChem.UFFOptimizeMolecule(mol)
    writer = AllChem.SDWriter("pepM_%s.sdf" %sequence)
    writer.write(mol)
    
    # Convert to PDB file
    os.system("babel -isdf pepM_{peptide}.sdf -opdb pepM_{peptide}.pdb".format(peptide=sequence))
    os.system("rm pepM_{peptide}.sdf".format(peptide=sequence))
    
################################################################    

def protonation(molecule):
    os.system("./pdb2pqr-linux/pdb2pqr --with-ph=5.3 --ph-calc-method=propka --drop-water --apbs-input --ff=amber --verbose --chain --summary %s.pdb %s.pqr" %(molecule,molecule))
    os.system("rm {mol}.in {mol}-input.p {mol}.summary {mol}.propka".format(mol=molecule))

################################################################

def generate_box_initial(sequence,center_x,center_y,center_z,pdb):
    
    number_amino=len(sequence)
    
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
    
    config=open("config.txt","wb")
    config.write("center_x=%0.3f\n" %center_x)
    config.write("center_y=%0.3f\n" %center_y)
    config.write("center_z=%0.3f\n" %center_z)
    config.write("size_x=%0.3f\n" %size_x)
    config.write("size_y=%0.3f\n" %size_y)
    config.write("size_z=%0.3f\n" %size_z)
    config.write("exhaustiveness=%d" %exhaustiveness)
    config.close()
    
    return size_x

################################################################

def generate_box(sequence,center_x,center_y,center_z,pdb,initial_size):
    
    number_amino=len(sequence)
    
    parser = PDBParser()
    structure = parser.get_structure('PEP', pdb)
    model = structure[0]
    
    # Get the difference in coordinates    
    diff_x=model["A"][1]["CA"].coord[0]-model["A"][number_amino]["CA"].coord[0]
    diff_y=model["A"][1]["CA"].coord[1]-model["A"][number_amino]["CA"].coord[1]
    diff_z=model["A"][1]["CA"].coord[2]-model["A"][number_amino]["CA"].coord[2]
    
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
    
    config=open("config.txt","wb")
    config.write("center_x=%0.3f\n" %center_x)
    config.write("center_y=%0.3f\n" %center_y)
    config.write("center_z=%0.3f\n" %center_z)
    config.write("size_x=%0.3f\n" %size_x)
    config.write("size_y=%0.3f\n" %size_y)
    config.write("size_z=%0.3f\n" %size_z)
    config.write("exhaustiveness=%d" %exhaustiveness)
    config.close()

################################################################

def docking(receptor,peptide,sequence,stepNumber):
    
    aminoacids={"A":"ALA","D":"ASP","E":"GLU","F":"PHE","H":"HIS","I":"ILE","K":"LYS","L":"LEU","M":"MET","G":"GLY",
                "N":"ASN","P":"PRO","Q":"GLN","R":"ARG","S":"SER","T":"THR","V":"VAL","W":"TRP","Y":"TYR","C":"CYS"}
    
    if stepNumber==0:
        os.system("./scripts/pythonsh scripts/prepare_receptor4.py -r %s.pqr -o %s.pdbqt -C -U waters" %(receptor,receptor))
        os.system("./scripts/pythonsh scripts/prepare_ligand4.py -l %s.pqr  -C -U '' -B -o %s.pdbqt" %(peptide,peptide))
        
        # KEY
        print "Docking step number %d ..." %stepNumber
        os.system("vina --receptor %s.pdbqt --ligand %s.pdbqt --log score.log --out out.pdbqt --config config.txt" %(receptor,peptide))
        #
        
        os.system("csplit out.pdbqt /MODEL/ {*}; mv xx01 model1_step%d.pdbqt" %stepNumber)
        os.system("rm xx* score.log out.pdbqt")
        # Count the number of models
        
        bash="ls | grep '_step%d.pdbqt' | grep 'model' | wc -l" %stepNumber
        numberModels = subprocess.check_output(['bash','-c', bash])
        return int(numberModels)
        
    else:
        rigidAA=[]
        for i in range(1,len(sequence)+1):
            if i!=1 and i!=2 and i!=len(sequence)-1 and i!=len(sequence):
                rigidAA.append(i)
        
        rigidAtoms=[]
        for rigid in rigidAA:
            if rigid <=9:
                os.system("grep %s %s.pqr | grep 'A   %d' | awk '{print $2}' > atom.temp" %(aminoacids[sequence[rigid-1]],peptide,rigid))
            else:
                os.system("grep %s %s.pqr | grep 'A  %d' | awk '{print $2}' > atom.temp" %(aminoacids[sequence[rigid-1]],peptide,rigid))
            num=[x.strip() for x in open("atom.temp")]
            rigidAtoms=rigidAtoms+num
            os.system("rm atom.temp")
        
        # Get the atoms that will be inactivated
        inactive=[]
        os.system("./scripts/pythonsh scripts/prepare_ligand4.py -s -l %s.pqr -C -U '' -B -o %s.pdbqt" %(peptide,peptide))
        os.system("grep REMARK %s.pdbqt | grep '  A  ' | awk '{print $6\"\t\"$8}' | sed 's/_/ /g' | awk '{print $2\"_\"$4}' > active.temp" %peptide)
        bonds=[x.strip() for x in open("active.temp")]
        for b in bonds:
            data=b.split("_")
            if data[0] in rigidAtoms or data[1] in rigidAtoms:
                inactive.append("%s_%s" %(str(int(data[0])-1),str(int(data[1])-1)))
        os.system("rm active.temp")
        
        commandInactive="_".join(inactive)
        
        
        os.system("./scripts/pythonsh scripts/prepare_ligand4.py -s -l %s.pqr -C -U '' -B -o %s.pdbqt -I %s" %(peptide,peptide,commandInactive))
        # KEY
        print "Docking step number %d ..." %stepNumber
        os.system("vina --receptor %s.pdbqt --ligand %s.pdbqt --log score.log --out out.pdbqt --config config.txt" %(receptor,peptide))
        #
        
        os.system("csplit out.pdbqt /MODEL/ {*}; mv xx01 model1_step%d.pdbqt" %stepNumber)
        os.system("rm xx* score.log out.pdbqt")
        
        # Count number models
        bash="ls | grep '_step%d.pdbqt' | grep 'model' | wc -l" %stepNumber
        numberModels = subprocess.check_output(['bash','-c', bash])
        return int(numberModels)

################################################################

def new_coordinates(base,model_pdbqt):
    print "Processing protein %s ..." %model_pdbqt
    parser = PDBParser()
    structure = parser.get_structure('PEP', "%s.pqr" %base)
    model = structure[0]
    counter=1
    batoms=["C","N","O","H","S"]
    
    for residue in model["A"]:
        resName=residue.get_resname()
        resNumber=residue.get_full_id()[3][1]
        
        for atom in residue:
            idAtom = atom.get_id()

            if resNumber<=9:
                if idAtom in batoms:
                    bash="grep %s %s.pdbqt | grep 'A   %s' | grep '  %s  ' | awk '{print $7}'" %(resName,model_pdbqt,str(resNumber),idAtom)
                    new_x = subprocess.check_output(['bash','-c', bash])
                    bash="grep %s %s.pdbqt | grep 'A   %s' | grep '  %s  ' | awk '{print $8}'" %(resName,model_pdbqt,str(resNumber),idAtom)
                    new_y = subprocess.check_output(['bash','-c', bash])
                    bash="grep %s %s.pdbqt | grep 'A   %s' | grep '  %s  ' | awk '{print $9}'" %(resName,model_pdbqt,str(resNumber),idAtom)
                    new_z = subprocess.check_output(['bash','-c', bash])
                    atom.set_coord([float(new_x),float(new_y),float(new_z)])
                else:
                    if len(idAtom)==4: idAtom=idAtom[-1]+idAtom[:3]
                        
                    bash="grep %s %s.pdbqt | grep 'A   %s' | grep ' %s ' | awk '{print $7}'" %(resName,model_pdbqt,str(resNumber),idAtom)
                    new_x = subprocess.check_output(['bash','-c', bash])
                    bash="grep %s %s.pdbqt | grep 'A   %s' | grep ' %s ' | awk '{print $8}'" %(resName,model_pdbqt,str(resNumber),idAtom)
                    new_y = subprocess.check_output(['bash','-c', bash])
                    bash="grep %s %s.pdbqt | grep 'A   %s' | grep ' %s ' | awk '{print $9}'" %(resName,model_pdbqt,str(resNumber),idAtom)
                    new_z = subprocess.check_output(['bash','-c', bash])
                    atom.set_coord([float(new_x),float(new_y),float(new_z)])
            else:
                if idAtom in batoms:
                    bash="grep %s %s.pdbqt | grep 'A  %s' | grep '  %s  ' | awk '{print $7}'" %(resName,model_pdbqt,str(resNumber),idAtom)
                    new_x = subprocess.check_output(['bash','-c', bash])
                    bash="grep %s %s.pdbqt | grep 'A  %s' | grep '  %s  ' | awk '{print $8}'" %(resName,model_pdbqt,str(resNumber),idAtom)
                    new_y = subprocess.check_output(['bash','-c', bash])
                    bash="grep %s %s.pdbqt | grep 'A  %s' | grep '  %s  ' | awk '{print $9}'" %(resName,model_pdbqt,str(resNumber),idAtom)
                    new_z = subprocess.check_output(['bash','-c', bash])
                    atom.set_coord([float(new_x),float(new_y),float(new_z)])
                else:
                    if len(idAtom)==4: idAtom=idAtom[-1]+idAtom[:3]
                        
                    bash="grep %s %s.pdbqt | grep 'A  %s' | grep ' %s ' | awk '{print $7}'" %(resName,model_pdbqt,str(resNumber),idAtom)
                    new_x = subprocess.check_output(['bash','-c', bash])
                    bash="grep %s %s.pdbqt | grep 'A  %s' | grep ' %s ' | awk '{print $8}'" %(resName,model_pdbqt,str(resNumber),idAtom)
                    new_y = subprocess.check_output(['bash','-c', bash])
                    bash="grep %s %s.pdbqt | grep 'A  %s' | grep ' %s ' | awk '{print $9}'" %(resName,model_pdbqt,str(resNumber),idAtom)
                    new_z = subprocess.check_output(['bash','-c', bash])
                    atom.set_coord([float(new_x),float(new_y),float(new_z)])
            
            counter+=1
    
    io = PDBIO()
    io.set_structure(structure)
    io.save('%s.pdb' %model_pdbqt)
    

################################################################        
# Main execution
if __name__ == '__main__':
    # Input
    ref_peptide="ENPVVHFFKNIVTPR"
    peptide="FFK"
    receptor="1BX2-ABP"
    number_chains=2
    
    # Predict peptide conformer
    generateConformer(peptide)
     
    # Protonation
    protonation(receptor)
    protonation("pepM_%s" %peptide)
    
    # Docking
    # Center of reference selected by the user
    center_x=-5.256
    center_y=26.343
    center_z=-2.416
    initial_size=generate_box_initial(peptide,center_x,center_y,center_z,"pepM_%s.pdb" %peptide)
    
    num_models=docking(receptor,"pepM_%s" %peptide,peptide,0)
    new_coordinates("pepM_%s" %peptide,"model1_step0")
    
    # Move to a folder with step0
    os.system("mkdir step0; mv pepM_%s.* *_step0.pdbqt *_step0.pdb config.txt step0" %peptide)
    # 
    # list of future fragments
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
    
    print list_fragments
    
    
    #flagEnd=0
    final_step=0
    pep_reference=peptide
    for i,frag in enumerate(list_fragments):
        
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
        
        modelling_fragments(pepTemplate,frag,i,receptor,number_chains)
        protonation("pepM_step%d" %(i+1))
        
        generate_box(frag,center_x,center_y,center_z,"pepM_step%s.pdb" %(i+1),initial_size)
        num_models=docking(receptor,"pepM_step%d" %(i+1),frag,i+1)
        if num_models==0:
           print "Docking error ... exiting the protocol"
           break
        
        new_coordinates("pepM_step%d" %(i+1),"model1_step%d" %(i+1))
        os.system("mkdir step{num_step}; mv pepM_step{num_step}.* *_step{num_step}.pdb *_step{num_step}.pdbqt *_step{num_step}.pqr config.txt step{num_step}".format(num_step=str(i+1)))
        final_step=i+1
        
    # Save the final complex
    os.system('cp step%d/model1_step%d.pdb .' %(final_step,final_step))
    if number_chains==1: os.system("sed -i 's/ A  / B  /g' model1_step%d.pdb" %final_step)
    if number_chains==2: os.system("sed -i 's/ A  / C  /g' model1_step%d.pdb" %final_step)
    os.system("cat %s.pdb model1_step%d.pdb > final_complex_%s.pdb" %(receptor,final_step,ref_peptide))
    os.system("rm model1_step%d.pdb" %final_step)
