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

To allow the execution of PepFun with python3 and the required dependencies, the best option is generate a conda virtual environment. For that purpose, the following command can be used:

`conda create -c rdkit -n pepfun-env rdkit biopython matplotlib scipy pip pycairo nb_conda_kernels`

To activate the environment you can use the command:

`conda activate pepfun-env`

After entering the virtual environment, this notebook can be called with jupyter and selecting the kernel provided by this conda instance. Then, we can install the igraph module for python3 using pip:

`pip install python-igraph`

## How to run the script

The basic command line to run the script is:

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
 
The required arguments are information to access the mySQL MEROP database. For the other parameters the script has default values, including a PDB of reference available in the annotated proteases, or just run the models for all the families and structures present in the dataset. **However, please be aware of changing the Rosetta version through the flag or directly in the script by default.**

## User examples

### Calculate sequence properties

The modelling of substrates bound to proteases has two main modes. One involves the selection of proteases bound to peptides composed of natural amino acids (called *ready* in the script). The second mode is when the substrate has one non-natural amino acids (NNAA) to modify. The following is a command example for the first case:

`python3.5 model_substrate_protease.py -f serine -n 1 -m ready -u user -p password -d merops`

In this case, we are allowing the modelling of only one substrate per structure available in the dataset bound to peptides with natural amino acids. For this script two families are available: serine and cysteine proteases. Here we selected serine as the reference family. Finally the MEROPS credentials are provided to look for reported substrates

After that, the models will be stored in the **models/model_ready** folder with the name `[structure]\_[substrate sequence].pdb`, where the structure is the PDB id and the substrate sequence is the peptide that was fully modelled using the protocol. A report of the model is provided in the following form:
```
pdb,pepTemplate,pepModel,chain,merops,old_aa,new_aa,pos_aa,uniprot
1smf,-CAKSI--,SCAKSIIG,I,S01.151,THR,ALA,3,O60256
...
```
Here the fragment CAKSI was modelled into the 8-mer peptide SCAKSIIG, which is part of the substrate protein with UniProtKB id O60256. The model was performed using the protease structure with PDB id 1smf that belongs to the MEROPS family S01.151 (trypsin). In addition, an amino acid from the original structure had to be changed by another one reported in the substrate, in this case a threonine by an alanine in the position 3 of the peptide.

### Calculate structural information of a protein-peptide complex

For the second case, we require the modification of a NNAA for a natural amino acid based on a similarity threshold defined by the user (called *complete* in the script). It means that we can model a substrate depending on how similar we want the template and new amino acids be. The following is an example based on a particular PDB structure available in the dataset:

`python3.5 model_substrate_protease.py -s 1tps -f serine -n 1 -m complete -t 0.4 -u user -p password -d merops`

Here we are modelling a substrate based on the template present in the structure with PDB id 1tps. The natural amino acid that will replace the present NNAA require to be at least 40\% similar based on the Tanimoto comparison of the amino acid side chains. We selected serine proteases as the reference family, and the MEROPS credentials are also provided to look for reported substrates.

After that, the model is stored in the **models/model_complete** folder with the name `[structure]\_[substrate sequence]\_modelled.pdb`, where the structure is the PDB id and the substrate sequence is the peptide that was fully modelled using the protocol. A report of the model is provided in the following form:

```
pdb,pepTemplate,pepModel,chain,merops,old_aa,new_aa,pos_aa,uniprot,sim
1tps,-LTREL--,FLTRELAE,B,S01.151,DLE,LEU,247,P23396,1.0
```

Here the fragment LTREL was modelled into the 8-mer peptide FLTRELAE, which is part of the substrate protein with UniProtKB id P23396. The model was performed using the protease structure with PDB id 1tps that belongs to the MEROPS family S01.151 (trypsin). In addition, a NNAA from the original structure had to be changed by another one reported in the substrate, in this case the residue DLE by an L-Leucine, with a similarity of 100\%.

## Developer modifications

Tutorial

### How to run the sampling script

After having modelled 8-mer peptides in protease reference structures, it is possible to call the second script for modelling any peptide of interest, run a simulation of the system using the backrub method from Rosetta and calculate structural descriptors from the trajectory. The basic command line to run the script is:

`run_dynamic_proteases.py [-h] -p PATH -s SEQUENCE -c CHAIN [-r ROSETTA]`
                                       
where the arguments are:

```
arguments:
  -h, --help   show this help message and exit
  -p PATH      Path of the structure that will be used to run the analysis
  -s SEQUENCE  Sequence of the peptide that will be modelled and sampled
  -c CHAIN     Chain identifier of the peptide in the structure
  -r ROSETTA   Version of Rosetta that will be implemented
 ```

The required arguments are the path of the model that we want to use as template, the sequence of the peptide that will be modelled, and the chain of the peptide in the structure of reference. **Please be aware of changing the Rosetta version through the flag or directly in the script by default.**

### Use a docker container with all the required dependencies

ref_peptide="vim fr"
    # peptide="FFK"
    # receptor="1BX2-ABP"
    # number_chains=2
    # pH=5.3
    # # Center of reference selected by the user
    # center_x=-7.1
    # center_y=24.5
    # center_z=-2.41

To run the dynamic analysis of a peptide of reference, we can call the script as:

`python3.5 run_dynamic_proteases.py -p models/model_ready/1smf_SCAKSIIG_modelled.pdb -s TGYHKLPR -c B`



## Support

In case the protocol is useful for other research projects and require some advice, please contact us to the email: rodrigo.ochoa@udea.edu.co

