# Helix Analysis Program

## Presentation


Proteins are essential components for life. There 
Python code for analyzing molecules, using the model of deformable molecular cylinders.

Loosely based on the exhaustive and polyvalent program Simulaid by Mihaly Mezei (*Simulaid: a simulation facilitator and analysis program, 2010*).

Work done in a team of four for second year biological engineering project.

## Requirements

python 3

numpy

matplotlib

module cogent3

MDAnalysis

biopython=1.68 


## How it works

You have to launch main.py in the folder **src**.

Command line

python.exe main.py filename

filename must be either a static pdb file or a trajectory file (concatenation of pdb files). The program detects if it is a trajectory or static file.

Please be patient with trajectory file.

The output are several graphs along with command line output.

The graphs are saved in same folder

## Examples

python.exe main.py glut1.pdb 

python.exe main.py TPSO_traj.pdb

## Functions

The following functions are available:

- DSSP
- Inertia Axis
- Turn Per Residue
- Helices Length
- Center of Mass
