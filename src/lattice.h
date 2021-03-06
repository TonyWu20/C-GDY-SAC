#pragma once
#include "atom.h"
#include "molecule.h"

struct carbon_site
{
    char *name;
    int id;
};
struct _Lattice
{
    Molecule *mol;
    matrix_double3x3 lattice_vectors;
    int metal_site_id;
    char *attached_adsName;
    char *pathName;
    struct Lattice_vtable *vtable;
};

typedef struct _Lattice Lattice;

struct Lattice_vtable
{
    vec_double3 (*get_carbon_chain_vector)(Lattice *self);
    vec_double3 (*get_carbon_metal_vector)(Lattice *self, int);
    Lattice *(*attach_molecule)(Lattice *self, Adsorbate *ads, char *newName);
    /* Rotate lattice to standard orientation "C along Z, A in XZ plane" */
    void (*rotate_to_standard_orientation)(Lattice *self);
    void (*modify_metal)(Lattice *self, const char *metalSymbol, int elementId);
    char *(*exportDir)(Lattice *self);
    void (*export_msi)(Lattice *self);
    void (*destroy)(Lattice *self);
};

/* Create a pointer to Lattice struct and malloc memory for it. Inherit from
 * Molecule object, with adding lattice_vectors (3x3 column-major matrix)
 */
Lattice *createLattice(Molecule *mol, matrix_double3x3 lattice_vectors);

/* Change the element of the metal atom */
void lattice_modify_metal_element(Lattice *self, const char *metal_symbol,
                                  int elementId);

/* Fill metal info after initialization */
void lattice_metal_info(Lattice *self);

/* Free memory occupied by the Lattice struct pointer
 */
void destroyLattice(Lattice *self);

/* Returns a pointer to Atom struct with given index from the atom_arr
 * param:
 *		Lattice *self
 *		int Id: index in atom_arr
 * returns: pointer to Atom
 */
Atom *lattice_get_atom_by_Id(Lattice *self, int);

/* Returns a vector of the carbon chain */
static inline vec_double3 lattice_get_carbon_chain_vector(Lattice *self)
{
    return self->mol->vtable->get_vector_ab(self->mol, 41, 42);
}
/* Returns a vector point from carbon to metal atom */
static inline vec_double3 lattice_get_carbon_metal_vector(Lattice *self,
                                                          int cId)
{
    return self->mol->vtable->get_vector_ab(self->mol, cId,
                                            self->metal_site_id);
}
/* Add a Molecule into the Lattice. The atomId and treeId
 * of the atoms in mol will be updated to follow the order in current Lattice.
 * Returns a new Lattice struct pointer for future exports
 */
Lattice *lattice_attach_molecule(Lattice *self, Adsorbate *ads, char *newName);

/* Rotate lattice to standard orientation "C along Z, A in XZ plane"
 */
void lattice_rotate_to_standard_orientation(Lattice *self);

char *get_carbon_site_name(int siteId);

/* Format and export lattice contents to .msi
 * if dest == NULL, use default
 */
void lattice_export_MSI(Lattice *self);

/* Generate export destination */
char *lattice_export_dest(Lattice *self);
