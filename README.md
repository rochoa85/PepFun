# PepFun project

## Package with bioinformatics and cheminformatics protocols for peptide analysis

* From publication "Pepfun: bioinformatics and cheminformatics protocols for peptide-related computational analysis"
* Journal of Cheminformatics, 2020
* Authors: Rodrigo Ochoa, Roman Laskowski, Pilar Cossio

## Purpose

PepFun is a compilation of bioinformatics and chemoinformatics functionalities that are easy to implement and personalize for studying peptides at different levels: sequence, structure and their interactions. The package has been created under the python scripting language based on built-in functions and methods available in the open source projects BioPython and RDKit. Some of the prediction and characterization tools were tested with two datasets of peptide binders of known protein systems, the MHC class II and the Granzyme B protease. In addition, a fragment-growing peptide docking protocol is provided to predict bound conformations of flexible peptides, and their interactions with relevant protein targets

## Third-party tools required:

- BioPython: https://biopython.org/wiki/Download
- RDKit: https://github.com/rdkit/rdkit/releases

To allow the execution of PepFun with python3 and the required dependencies, the best option is to generate a conda virtual environment. For that purpose, the following command can be used:

`conda create -c rdkit -n pepfun-env rdkit biopython matplotlib scipy pip pycairo nb_conda_kernels`

To activate the environment you can use the command:

`conda activate pepfun-env`

After entering the virtual environment, we can install the igraph module for python3 using pip:

`pip install python-igraph`

## How to run the script

PepFun can be called as a module into another script to use their functionalities and combine the output with other tools. Examples of this mode are provided in the *test.py* script. However, the script by itself can be called to run some basic sequence- and structure-based analysis with peptides using the following syntax:

`peptide_functions.py [-h] -m MODE [-s PEP_SEQ] [-p PEP_STR]
                            [-c PEP_CHAIN] [-b PEP_CONFORMATION]
                            [-t CONTACT_THRESHOLD]`
                                       
where the arguments are:

```
arguments:
  -h, --help            show this help message and exit
  -m MODE               Choose a mode to run the script from two options: 
                        1) sequence, 2) structure.
  -s PEP_SEQ            Sequence of peptide to be analyzed
  -p PEP_STR            Structure that will be used for the analysis
  -c PEP_CHAIN          Chain of the peptide in the structure
  -b PEP_CONFORMATION   Conformation of the peptide in the structure that will
                        be used to plot hydrogen bonds
  -t CONTACT_THRESHOLD  Threshold to count contacts between the peptide and
                        the protein
 ```
 
The main required argument is the mode, which is *sequence* or *structure* depending on the required analysis. If *sequence* mode is selected, an amino acid string should be provided. If *structure* mode is selected, various arguments should be provided. These include the path to the PDB file having the peptide alone or in complex with a protein, the chain ID of the peptide in the structure, the conformation of the peptide (by default is *linear*, or can be *cyclic*), and a contact threshold defined to count the interactions of the peptide with the protein chains (by default is *4.0*).

## User examples

### Calculate sequence properties

As an example to run the script as an user, we can calculate a set of properties. For that purpose, the script can be called as:

`python peptide_functions.py -m sequence -s GYTRTEGSDF`

The output is a file with the name `sequence\_analysis\_[sequence].txt`, where multiple properties for the required sequence are reported. An example of this is:

```
The main peptide sequence is: GYTRTEGSDF
Net charge at pH 7: -1.0008535231234645
Molecular weight: 1132.152
Average Hydrophobicity: -2.04
Isoelectric Point: 4.37030029296875
Instability index: 3.740000000000001
Number of hydrogen bond acceptors: 18
Number of hydrogen bond donors: 20
Crippen LogP: -7.4280300000000254
1 solubility rules failed from 5
0 synthesis rules failed from 5
```

In addition, a file named `structure_[sequence].sdf` will contain the peptide structure predicted by the conformer option in RDKit. These can be converted to a PDB file using OpenBabel to allow the correct numeration and order of the amino acids in the structure file.


### Calculate structural information of a protein-peptide complex

For the second case, we can call the script to obtain some information of a structure containing a peptide alone or in complex with a protein. This is an example using a structure provided in the auxiliar folder of the code:

`python peptide_functions.py -m structure -p auxiliar/example_structure.pdb -c C -b linear -t 4.0`

Here we are analyzing a peptide bound to the MHC class II protein, which chain ID is C and has a linear bound conformation. To count the number of contacts, we selected as threshold a value of 4.0. After running the script, we obtain the file `structure\_analysis\_[sequence].txt` with the report of some interaction observables:

```
Peptide sequence based on the PDB file is: NPVVHFFKNIVTPRTPPPSQ
The total number of contacts are: 181
The total number of hydrogen bonds are: 25
```

In addition, we report the details of the hydrogen bonds detected between the peptide and the protein, and a plot of the interactions based on the selected conformation: linear or cyclic.

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

### Specialized tutorial for developers

In case the user want to explore in detail the functions and applications using massive datasets, a Jupyter notebook is provided to run a set of operations with PepFun modules. *To run successfully the example please verify that the conda environment was created successfully.*

### How to access the fragment-docking script

Because the fragment-docking protocol depends on multiple additional resources and a Unix environment, we provide the code in the form of a docker with all the required dependencies. You can find the docker container here: XXX

An example of the script syntax is as follows::

`fragment_docking.py [-h] -s PEP_SEQ -f PEP_FRAG -t TARGET -n NUM_CHAINS
                           -x CENTER_X -y CENTER_Y -z CENTER_Z [-p PEP_PH]`
                                       
where the arguments are:

```
arguments:
  -h, --help     show this help message and exit
  -s PEP_SEQ     Sequence of peptide to be analyzed
  -f PEP_FRAG    Fragment that will be used as initial anchor
  -t TARGET      Name of the PDB structure used as target
  -n NUM_CHAINS  Number of chains present in the target structure
  -x CENTER_X    Coordinate in X of the box center
  -y CENTER_Y    Coordinate in Y of the box center
  -z CENTER_Z    Coordinate in Z of the box center
  -p PEP_PH      pH of the system
 ```

In addition of the script, the folder require of the target PDB structure file, a folder with a set of necessary scripts, and an output folder where the docking results step by step will be stored.

```
[target].pdb output scripts
```
An example to run the protocol script using the structure provided in the docker folder `/home/docking` is here:

`
python3 fragment_docking.py -s ENPVVHFFKNIVTPR -f FFK -t 1BX2-ABP -n 2 -x "-7.1" -y "24.5" -z "-2.41"
`

The method will start the modelling of the initial script, and the docking of the fragment after each growing step, until the peptide obtain the final desired size. The code can be modified to modifiy box sizes, as well as verify if the conformation of the growing ligand is according to previous findings of the biological system. An example of the output docked result, and the configuration each docking step is available in the output folder as follows:

`final_complex_ENPVVHFFKNIVTPR.pdb  step0  step1  step2	step3  step4  step5  step6`


## Support

In case the protocol is useful for other research projects and require some advice, please contact us to the email: rodrigo.ochoa@udea.edu.co

