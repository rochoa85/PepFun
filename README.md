# PepFun project

## Bioinformatics and cheminformatics protocols for peptide analysis

* From publication "PepFun: open source protocols for peptide-related computational analysis". Link: https://www.mdpi.com/1420-3049/26/6/1664
* Molecules, 2021
* Authors: Rodrigo Ochoa, Pilar Cossio

## Purpose

PepFun is a compilation of bioinformatics and chemoinformatics functionalities that are easy to implement and personalize for studying peptides at different levels: sequence, structure and large datasets. The package has been created under the python scripting language based on built-in functions and methods available in the open source projects BioPython and RDKit. Some of the prediction and characterization tools were tested with two datasets of peptide binders of known protein systems, the MHC class II and the Granzyme B protease.

## Third-party tools required:

- BioPython: https://biopython.org/wiki/Download (recommended version 1.76)
- RDKit: https://github.com/rdkit/rdkit/releases

To execute PepFun with python3 and the required dependencies, the best option is to generate a conda virtual environment. A guide to install Conda is available here: https://docs.conda.io/projects/conda/en/latest/user-guide/install/. After the installation, the virtual environment required for PepFun can be created using the following command:

`conda create -c rdkit -n pepfun-env rdkit biopython matplotlib scipy pip pycairo nb_conda_kernels`

To activate the environment you can use the command:

`source activate pepfun-env`

After entering the virtual environment, you can install the igraph module for python3.6 using pip:

`python3.6 -m pip install python-igraph`

**EXTRA**: An addition was included to predict peptide conformers using MODELLER and PSIPRED. To install modeller with conda use the following command:

`conda install -c salilab modeller`

After that, add the academic license key in the file suggested. The license can be acquired in the Modeller website: https://salilab.org/modeller/registration.html

**NOTES:**
1. If the script is called as a module from a different folder, the path where the pepfun.py script is can be added using the following commands:
```
import sys
sys.path.append('<PATH-TO-PEPFUN>')
```
2. The package was build under a Unix environment, but it can be used under any other OS based on the provided paths.

3. To run the Jupyter tutorial, please verify that Jupyter Notebook is installed in your environment.

## How to run the script

PepFun can be called to run some basic sequence- and structure-based analysis with peptides using the following syntax:

`pepfun.py [-h] -m MODE [-s PEP_SEQ] [-p PEP_STR]
                            [-c PEP_CHAIN] [-b PEP_CONFORMATION]
                            [-d DSSP_ROUTE] [-t CONTACT_THRESHOLD]`
                                       
where the arguments are:

```
arguments:
  -h, --help            show this help message and exit
  -m MODE               Choose a mode to run the script from two options: 1)
                        sequence, 2) structure.
  -s PEP_SEQ            Sequence of peptide to be analyzed
  -p PEP_STR            Structure that will be used for the analysis
  -c PEP_CHAIN          Chain of the peptide in the structure
  -r CONFORMER_MODE     Mode to generate the basic peptide conformer. 
                        Options: rdkit (default), modeller (with SS restraints).
  -b PEP_CONFORMATION   Conformation of the peptide. Options: linear
                        (default), cyclic
  -d DSSP_ROUTE         Route where the mkdssp program is located. By default
                        it is located in the auxiliar folder
  -t CONTACT_THRESHOLD  Threshold (in angstroms) to count contacts between the
                        peptide and the protein
 ```
 
The main argument required is the mode, which is *sequence* or *structure* depending on the analysis. If the *sequence* mode is selected, an amino acid string should be provided. If the *structure* mode is selected, various arguments should be provided. These include the path to the PDB file having the peptide alone or in complex with a protein, the chain ID of the peptide in the structure, the conformation of the peptide (by default it is *linear*, or can be *cyclic*), and a contact threshold to count the interactions of the peptide with the protein chains (the default is *4.0* (angstroms)).

PepFun can be also ran as a module within another script to use its functionalities and combine the output with other tools. See the Jupyter Tutorial for examples (tutorial_PepFun.ipynb). 

To test that PepFun is working, the following command can be run to calculate some sequence and structural data from an example peptide:

`python test.py`

The output files can be checked to confirm that PepFun is ready to go.

## User examples

### Calculate sequence properties

As an example to run the script, we can calculate a set of properties for a peptide sequence. For that purpose, the script can be called as:

`python pepfun.py -m sequence -s GYTRTEGSDF`

**NOTE: Remember to activate first the conda virtual environment as explained previously.**

The output is a file with the name `sequence_analysis_[sequence].txt`, where multiple properties for the sequence are reported. An example of this and the explanation of the parameters is:

```
###########################################
The main peptide sequence is: GYTRTEGSDF
###########################################

Calculated properties:
Net charge at pH 7: -1.0008535231234645
Molecular weight: 1132.152
Average Hydrophobicity: -2.04
Isoelectric Point: 4.37030029296875
Instability index: 3.740000000000001
Number of hydrogen bond acceptors: 18
Number of hydrogen bond donors: 20
Crippen LogP: -7.4280300000000254
1 solubility rules failed from 5. The rules violated are the number(s): 3. *For details of the rules see the last part of the report.
0 synthesis rules failed from 5. The rules violated are the number(s): None. **For details of the rules see the last part of the report.
###########################################

Explanation of the calculated parameters:
- Net charge: average net charge based on pka values of each amino acid. By default a pH=7 is used.
- Molecular weight: calculated in g/mol using the SMILES representation of the peptide.
- Average Hydrophobicity: calculated by averaging the values of each amino acid hydrophobicity value from the Eisenberg scale.
- Isoelectric point: obtained from the ProtParam package of the expasy server.
- Instability Index: from ProtParam. It provides an estimate of the stability of the peptide in a test tube. Values smaller than 40 is predicted as stable, a value above 40 predicts as unstable.
- Number of hydrogen bond acceptors: calculated using the SMILES representation of the peptide.
- Number of hydrogen bond donors: calculated using the SMILES representation of the peptide.
- Crippen LogP: estimation of the octanol/water partition coefficient using the Ghose/Crippen approach available in the RDKit project.
###########################################

The last two results are the number of solubility and synthesis rules violated. The higher the number of rules violated, the lower the probability to be solubilized or synthesized experimentally (https://bioserv.rpbs.univ-paris-diderot.fr/services/SolyPep/).

*List of solubility rules violations:
1. Discard if the number of charged and/or of hydrophobic amino acids exceeds 45%
2. Discard if the absolute total peptide charge at pH 7 is more than +1
3. Discard if the number of glycine or proline is more than one in the sequence
4. Discard if the first or the last amino acid is charged
5. Discard if any amino acid represents more than 25% of the total sequence

**List of synthesis rules violations:
1. Discard if 2 prolines are consecutive
2. Discard if the motifs DG and DP are present in the sequence
3. Discard if the sequences ends with N or Q residues
4. Discard if there are charged residues every 5 amino acids
5. Discard if there are oxidation-sensitive amino acids (M, C or W)
```

In addition, a file named `structure_[sequence].pdb` will contain the peptide structure predicted by the conformer option in RDKit (see main publication), which will have the correct numeration and order of the amino acids based on the input sequence. With the `-r modeller` option it is possible to predict the conformer using Modeller together with secondary structure restraints to be more accurate about the prediction. 


### Calculate structural information of a protein-peptide complex

We can use PepFun to obtain some information about a structure containing a peptide alone or in complex with a protein. This is an example for using a structure provided in the auxiliar folder of the code:

`python pepfun.py -m structure -p auxiliar/example_structure.pdb -c C -b linear -t 4.0 -d auxiliar/mkdssp`

Here, we are analyzing a peptide bound to the MHC class II protein, which chain ID is C and has a linear bound conformation. To count the number of contacts, we selected as threshold a value of 4.0 (angstroms). After running the script, we obtain the file `structure_analysis_[sequence].txt` with the report of some interaction observables:

```
Peptide sequence based on the PDB file is: NPVVHFFKNIVTPRTPPPSQ
The total number of contacts are: 181
The total number of hydrogen bonds are: 25
```

In addition, we report the details of the hydrogen bonds detected between the peptide and the protein, and a plot of the interactions based on the selected conformation: linear or cyclic. The file name is `plot_hbs_[sequence].png`.

```
These are the hydrogen bonds detected:
V4 interacts with residue R50 from chain A
H5 interacts with residue F51 from chain A
H5 interacts with residue G49 from chain A
F6 interacts with residue R50 from chain A
F6 interacts with residue A52 from chain A
F7 interacts with residue F51 from chain A
F7 interacts with residue S53 from chain A
K8 interacts with residue T77 from chain B
K8 interacts with residue H81 from chain B
N9 interacts with residue S53 from chain A
N9 interacts with residue S53 from chain A
I10 interacts with residue K12 from chain B
I10 interacts with residue K12 from chain B
V11 interacts with residue Q57 from chain A
V11 interacts with residue A59 from chain A
T12 interacts with residue Q10 from chain B
T12 interacts with residue I63 from chain A
P13 interacts with residue K67 from chain A
R14 interacts with residue A64 from chain A
R14 interacts with residue I63 from chain A
T15 interacts with residue I63 from chain A
T15 interacts with residue L70 from chain A
P16 interacts with residue E59 from chain B
P17 interacts with residue L70 from chain A
P18 interacts with residue L70 from chain A

```
### Generation of Libraries 

To generate different types of peptide libraries, it is required to call a function from pepfun.py. This means that calling pepfun.py as a stand-alone tool is not suitable for this purpose. An example of how to generate a library is available in the Jupyter Notebook tutorial (tutorial_PepFun.ipynb).

### Specialized tutorial for developers

If the user wants to explore in detail the functions and applications using massive datasets, a Jupyter Notebook is provided to run a set of operations with PepFun modules. *To run the example please verify that the conda environment was created successfully.*

## Support

For support, please contact us to the email: rodrigo.ochoa@udea.edu.co. Please cite: "Ochoa & Cossio. PepFun: open source protocols for peptide-related computational analysis. Molecules, 2021." (to be published)

