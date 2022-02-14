#pragma once
#include "atom.h"
#include "my_maths.h"
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
    Molecule *mol;
    int coordAtomNum;
    int *coordAtomIds;
    int stemAtomIds[2];
    int planeAtomIds[3];
    bool bVertical;
    bool bSym;
    int upperAtomId;
    struct Adsorbate_vtable *vtable;
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
    simd_double3 (*get_vector_ab)(Molecule *self, int aId, int bId);
    simd_double3 (*get_centroid_ab)(Molecule *self, int aId, int bId);
    void (*rotateMol)(Molecule *self, simd_quatd rotation);
    char **(*export_text)(Molecule *);
    Molecule *(*duplicate)(Molecule *);
    void (*destroy)(Molecule *self);
};

struct Adsorbate_vtable
{
    simd_double3 (*get_stem_vector)(Adsorbate *self);
    simd_double3 (*get_plane_normal)(Adsorbate *self);
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
                           int *planeAtomIds, bool bVer, bool bSym,
                           int upperAtomId);
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
simd_double3 Molecule_get_vector_ab(Molecule *, int a, int b);
simd_double3 Molecule_get_centroid_ab(Molecule *, int a, int b);
void Molecule_apply_rotation(Molecule *, simd_quatd rotation);
simd_double3 Adsorbate_get_stem_vector(Adsorbate *);
simd_double3 Adsorbate_get_plane_normal(Adsorbate *);
void Adsorbate_make_upright(Adsorbate *mol);
char **Molecule_textblock(Molecule *);

/* Format adsorbate data and write to .msi.
 * param:
 *		dest (char *): can be NULL for default export location
 */
void Adsorbate_export_MSI(Adsorbate *, char *dest);
