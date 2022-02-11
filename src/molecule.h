#pragma once
#include "atom.h"

struct _Molecule
{
    char *name;
    int atomNum;
    Atom **atom_arr;
    double *coordMatrix;
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
    int bSym;
    int upperAtomId;
    struct Adsorbate_vtable *ads_vtable;
    struct taskTable *taskLists;
} Adsorbate;

struct taskTable
{
    int taskNum;
    int **tasks;
};
struct taskTable *createTasks(Adsorbate *self);

struct Molecule_vtable
{
    Atom *(*get_atom_by_Id)(Molecule *self, int atomId);
    double *(*get_mol_coords)(Molecule *self);
    void (*update_atom_coords)(Molecule *self, double *MolCoords);
    double *(*get_vector_ab)(Molecule *self, int aId, int bId);
    double *(*get_centroid_ab)(Molecule *self, int aId, int bId);
    void (*apply_transformation)(Molecule *self, double *trans_m,
                                 void (*trans_func)(double *trans_m,
                                                    double *coords,
                                                    double **result));
    char **(*export_text)(Molecule *);
    Molecule *(*duplicate)(Molecule *);
    void (*destroy)(Molecule *self);
};

struct Adsorbate_vtable
{
    double *(*get_stem_vector)(Adsorbate *self);
    double *(*get_plane_normal)(Adsorbate *self);
    void (*make_upright)(Adsorbate *self);
    void (*export_msi)(Adsorbate *self, char *dest);
    Adsorbate *(*duplicate)(Adsorbate *self);
    void (*destroy)(Adsorbate *self);
};

// Memory Management
Molecule *createMolecule(char *name, int atomNum, Atom **atom_arr);
/* Duplicate the molecule */
Molecule *Molecule_duplicate(Molecule *self);
// Memory Management
Adsorbate *createAdsorbate(Molecule *newMol, int coordAtomNum,
                           int *coordAtomIds, int *stemAtomIds,
                           int *planeAtomIds, int bSym, int upperAtomId);
/* Duplicate the adsorbate */
Adsorbate *Adsorbate_duplicate(Adsorbate *self);

// Memory Management
void destroyMolecule(Molecule *self);
// Memory Management
void destroyAdsorbate(Adsorbate *self);

// Methods

// Returns an atom by AtomId
Atom *Molecule_get_Atom_by_Id(Molecule *, int);
// Returns the molecule coords in 4xN matrix
double *Molecule_get_coords(Molecule *);
// Updates the coords to atoms from the input matrix
void Molecule_update_Atom_coords(Molecule *, double *);
// Returns a vector from a to b by Id
double *Molecule_get_vector_ab(Molecule *, int a, int b);
double *Molecule_get_centroid_ab(Molecule *, int a, int b);
double *Adsorbate_get_stem_vector(Adsorbate *);
double *Adsorbate_get_plane_normal(Adsorbate *);
void Adsorbate_make_upright(Adsorbate *mol);
void Molecule_apply_transformation(Molecule *mPtr, double *trans_mat,
                                   void (*transform_func)(double *coords,
                                                          double *trans_mat,
                                                          double **result));
char **Molecule_textblock(Molecule *);

/* Format adsorbate data and write to .msi.
 * param:
 *		dest (char *): can be NULL for default export location
 */
void Adsorbate_export_MSI(Adsorbate *, char *dest);
