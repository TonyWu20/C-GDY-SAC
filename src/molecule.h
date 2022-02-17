#pragma once
#include "ads_database.h"
#include "atom.h"
#include <stdbool.h>

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
    bool bSym;
    bool bVer;
    int upperAtomId;
    char *pathName;
    struct Adsorbate_vtable *ads_vtable;
} Adsorbate;

struct Molecule_vtable
{
    Atom *(*get_atom_by_Id)(Molecule *self, int atomId);
    Matrix *(*get_vector_ab)(Molecule *self, int aId, int bId);
    double *(*get_centroid_ab)(Molecule *self, int aId, int bId);
    void (*apply_transformation)(Molecule *self, Matrix *trans_m,
                                 void (*trans_func)(Matrix *trans_m,
                                                    Matrix *coords,
                                                    Matrix **result));
    char **(*export_text)(Molecule *);
    Molecule *(*duplicate)(Molecule *);
    void (*destroy)(Molecule *self);
};

struct Adsorbate_vtable
{
    Matrix *(*get_stem_vector)(Adsorbate *self);
    Matrix *(*get_plane_normal)(Adsorbate *self);
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
                           int *planeAtomIds, bool bSym, bool bVer,
                           int upperAtomId, char *pathName);
/* Duplicate the adsorbate */
Adsorbate *Adsorbate_duplicate(Adsorbate *self);

// Memory Management
void destroyMolecule(Molecule *self);
// Memory Management
void destroyAdsorbate(Adsorbate *self);

// Methods

// Returns an atom by AtomId
Atom *Molecule_get_Atom_by_Id(Molecule *, int);
// Returns a vector from a to b by Id
Matrix *Molecule_get_vector_ab(Molecule *, int a, int b);
double *Molecule_get_centroid_ab(Molecule *, int a, int b);
Matrix *Adsorbate_get_stem_vector(Adsorbate *);
Matrix *Adsorbate_get_plane_normal(Adsorbate *);
void Adsorbate_make_upright(Adsorbate *mol);
void Molecule_apply_transformation(Molecule *mPtr, Matrix *trans_mat,
                                   void (*transform_func)(Matrix *coords,
                                                          Matrix *trans_mat,
                                                          Matrix **result));
char **Molecule_textblock(Molecule *);

/* Format adsorbate data and write to .msi.
 * param:
 *		dest (char *): can be NULL for default export location
 */
void Adsorbate_export_MSI(Adsorbate *, char *dest);
