# C GDY-SAC project

## Usage

1. Construct adsorption models for C1 pathway intermediates. _Not First goal_
1. Construct adsorption models for C2 pathway intermediates.

## Goals/Structures

1. Place the adsorbate molecule into the base model.

2. Write castep input files.
   1. Create folder named by the structure.
   1. .cell file
      - Format:
        1. Lattice vectors
        2. Fractional coords of atoms, sorted in order of element atomic number.
           - Requires sorting existing elements?
        3. Miscellaneous informations.
   1. .param file
   1. dos .cell file
   1. dos .param file
   1. .kptaux
   1. .trjaux
   1. copy SMCastep_Extensions
   1. copy .xsd and .msi files into folder.
3. Write pbs execution scripts.

## Design of codes

To achieve the generation of an adsorbate-attached SAC_GDY_Metal_Ads_Pos.msi
file, we need to do the following:

1. Design proper data structure to handle the lattice and adsorbate
   1. Basic unit: a data structure for atom
      - Element
      - ACL Label (**read from file for the time being** or generate by myself?)
      - XYZ coordinate
      - Atom Id
      - Id in msi tree
   2. Second level unit: a data structure for adsorbate and a structure for
      lattice
      - Adsorbate
        - name
        - array of `struct Atom`
        - coord atom number
        - array of coord atom id (check a hardcoded table for info)
        - array of stem atom id (to determine the direction of the molecule)
        - array of plane atom id (to determine the desired "plane" of the
          molecule to rotate it)
        - **Specify if it is symmetric**
      - Lattice
        - name
        - array of `struct Atom`
        - array of carbon sites ids (known)
1. Parse the current base-lattice file and adsorbate file into memory
1. Prepare the adsorbate molecule to be ready for placement
1. Assemble the adsorbate and lattice together
