#pragma once
#include "atom.h"

struct _Molecule;
typedef struct _Molecule Molecule;

// Memory Management
Molecule *createMolecule(char *name, int atomNum, Atom **atom_arr,
                         int coordAtomNum, int *coordAtomIds, int *stemAtomIds,
                         int *planeAtomIds);
void destroyMolecule(Molecule *molPtr);

// Methods
Atom *Molecule_get_Atom_by_Id(Molecule *, int);
Matrix **Molecule_get_coords(Molecule *);
void Molecule_update_Atom_coords(Molecule *, Matrix **);
Matrix *Molecule_get_vector_ab(Molecule *, int a, int b);
Matrix *Molecule_get_stem_vector(Molecule *);
Matrix *Molecule_get_plane_normal(Molecule *);
Matrix *Molecule_get_centroid(Molecule *);
