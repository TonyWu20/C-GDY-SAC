#pragma once
#include "atom.h"

struct _Molecule
{
    char *name;
    int atomNum;
    Atom **atom_arr;
    struct Molecule_vtable *vtable;
};
typedef struct _Molecule Molecule;

typedef struct
{
    Molecule *_mol;
    int coordAtomNum;
    int *coordAtomIds;
    int stemAtomIds[2];
    int planeAtomIds[3];
    struct Adsorbate_vtable *ads_vtable;
} Adsorbate;

struct Molecule_vtable
{
    Atom *(*get_atom_by_Id)(Molecule *self, int atomId);
    Matrix *(*get_mol_coords)(Molecule *self);
    void (*update_atom_coords)(Molecule *self, Matrix *MolCoords);
    Matrix *(*get_vector_ab)(Molecule *self, int aId, int bId);
};

struct Adsorbate_vtable
{
    Matrix *(*get_stem_vector)(Adsorbate *self);
    Matrix *(*get_plane_normal)(Adsorbate *self);
    void (*make_upright)(Adsorbate *self);
    void (*destroy)(Adsorbate *self);
};

// Memory Management
Molecule *createMolecule(char *name, int atomNum, Atom **atom_arr);
Adsorbate *createAdsorbate(Molecule *newMol, int coordAtomNum,
                           int *coordAtomIds, int *stemAtomIds,
                           int *planeAtomIds);
void destroyAdsorbate(Adsorbate *self);

// Methods
Atom *Molecule_get_Atom_by_Id(Molecule *, int);
Matrix *Molecule_get_coords(Molecule *);
void Molecule_update_Atom_coords(Molecule *, Matrix *);
Matrix *Molecule_get_vector_ab(Molecule *, int a, int b);
Matrix *Adsorbate_get_stem_vector(Adsorbate *);
Matrix *Adsorbate_get_plane_normal(Adsorbate *);
void Adsorbate_make_upright(Adsorbate *mol);
