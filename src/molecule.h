#pragma once
#include "atom.h"

struct _Molecule
{
    char *name;
    int atomNum;
    Atom **atom_arr;
    int coordAtomNum;
    int *coordAtomIds;
    int stemAtomIds[2];
    int planeAtomIds[3];
    struct Molecule_vtable *vtable;
};
typedef struct _Molecule Molecule;
struct Molecule_vtable
{
    Atom *(*get_atom_by_Id)(Molecule *self, int atomId);
    Matrix *(*get_mol_coords)(Molecule *self);
    void (*update_atom_coords)(Molecule *self, Matrix *MolCoords);
    Matrix *(*get_vector_ab)(Molecule *self, int aId, int bId);
    Matrix *(*get_stem_vector)(Molecule *self);
    Matrix *(*get_plane_normal)(Molecule *self);
    void (*make_upright)(Molecule *self);
    void (*destroy)(Molecule *self);
};

// Memory Management
Molecule *createMolecule(char *name, int atomNum, Atom **atom_arr,
                         int coordAtomNum, int *coordAtomIds, int *stemAtomIds,
                         int *planeAtomIds);
void destroyMolecule(Molecule *molPtr);

// Methods
Atom *Molecule_get_Atom_by_Id(Molecule *, int);
Matrix *Molecule_get_coords(Molecule *);
void Molecule_update_Atom_coords(Molecule *, Matrix *);
Matrix *Molecule_get_vector_ab(Molecule *, int a, int b);
Matrix *Molecule_get_stem_vector(Molecule *);
Matrix *Molecule_get_plane_normal(Molecule *);
void Molecule_make_upright(Molecule *mol);
