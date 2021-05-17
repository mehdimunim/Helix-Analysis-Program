# Helix Analysis Program


## Presentation

Python code for analyzing protein structure through dynamic studies to obtain dynamic information on proteins by using the model of deformable molecular cylinders for helices, defining geometric descriptors that best describe their deformations and their positioning to one another.


Work done in a team of 4 for second year biological engineering project.

## Workflow

The program has the following workflow:
- DSSP
- Inertia Axis
- Turn Per Residue
- Helices Length
- Center of Mass

## Folders

Our final code is in the folder src.

Our trials are dispatched between other folders.

## General Information

The program follows the following steps:

1- Determination of helices in a protein: an algorithm based on the DSSP principle [https://swift.cmbi.umcn.nl/gv/dssp/]  which determines the secondary structures of a protein by identifying the hydrogen bonds between the carbonyl and amide groups of the main chain.

2- Calculation of the principal axis of each helix from the alpha carbon atoms of the output of our DSSP algorithm.

3- Compute 4 descriptors that describe the helix trajectory of a protein:
- Tilting that define the orientation of each helix from the angles between the axes of each helical structure
- Compression and Extension that the helix undergoes along its trajectory from the length of the helix.
- Displacement of the helix from the position vector of the center of gravity of the atoms defining the helix.
- The rotation of the residues :  the turn angle per residu (TPR) 

## Technologies used

Project is created with:
* Python 2.7 or 3.6
* 
...

## Setup

To run this project :

```
$ cd ../l
...
```


## Flowchart

A flowchart of the whole program can be found :
https://github.com/mehdimunim/Interdisciplinary-Project/blob/main/flowcharts/flowchart%20program.png

## Sources

This program is inspired by Mihaly Meze analysis tool SIMULAID a simulation facilitator (2010) that performs a large number of simulation-related tasks: interconversion and modification of structure and trajectory files, optimization of orientation, and a large variety of analysis functions [https://mezeim01.u.hpc.mssm.edu/simulaid/simulaid.html].


## Contact

Created by :
- Mehdi Munim @mehdimunim  Role: Bioinformatics Support
- Hassna Boudalil @hassnabdl Role: Bioinformatics Support
- Romain Fourel  @Romain-Fourel  Role:  IT Support
- Naima Oulhourchemt @noulhourchemt   Role:  IT Support for Physics

Feel free to contact us !


