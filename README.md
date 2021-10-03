# C GDY-SAC project

## Usage

1. Construct adsorption models for C1 pathway intermediates.
1. Construct adsorption models for C2 pathway intermediates.

## Goals/Structures

1. Place the adsorbate molecule into the base model.
   1. Parse .msi file
      1. molecule
         1. struct ATOM
         ```C
         typedef struct {
            char *element;
            int itemId;
            double x;
            double y;
            double z;
            char *text;
         }
         ```
         2. The coordinates should be reset to set the first element as the origin.
      1. base model
         1. parse all atoms with struct ATOM.
         2. Get the last atom itemId.
   2. Compute coordinates
      1. align the coordinate site of adsorbate to the adsorption site.
      2. Rotate the molecule to align with the selected
         two adsorption sites (for C2).
   3. Attach adsorbate molecule atoms into base model.
      1. Calculate the next itemId for the inserted adsorbate molecule atoms.
      2. Insert the text into .msi files.
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
