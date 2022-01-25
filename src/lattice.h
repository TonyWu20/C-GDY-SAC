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
    Molecule *_mol;
    Matrix *lattice_vectors;
    struct carbon_site carbon_sites[7];
    int metal_site_id;
    int metal_order;
    char *metal_family;
    char *metal_symbol;
    struct Lattice_vtable *vtable;
};

typedef struct _Lattice Lattice;

struct Lattice_vtable
{
    Matrix *(*get_carbon_chain_vector)(Lattice *self);
    Matrix *(*get_carbon_metal_vector)(Lattice *self, int);
    Lattice *(*attach_molecule)(Lattice *self, Adsorbate *ads, char *newName);
    void (*export_msi)(Lattice *self, char *pathName);
    void (*destroy)(Lattice *self);
};

/* Create a pointer to Lattice struct and malloc memory for it. Inherit from
 * Molecule object, with adding lattice_vectors (3x3 column-major matrix)
 */
Lattice *createLattice(Molecule *mol, Matrix *lattice_vectors);

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

/* Returns a 4xN matrix containing all atoms' coordinates
 */
Matrix *lattice_get_coords(Lattice *self);

/* Update all atoms' coordinates from the given 4xN matrix
 */
void lattice_update_Atom_coords(Lattice *self, Matrix *);
/* Returns a vector from two given atoms' coordinates */
Matrix *lattice_get_vector_ab(Lattice *self, int a, int b);
/* Returns a vector of the carbon chain */
Matrix *lattice_get_carbon_chain_vector(Lattice *self);
/* Returns a vector point from carbon to metal atom */
Matrix *lattice_get_carbon_metal_vector(Lattice *self, int);
/* Add a Molecule into the Lattice. The atomId and treeId
 * of the atoms in mol will be updated to follow the order in current Lattice.
 * Returns a new Lattice struct pointer for future exports
 */
Lattice *lattice_attach_molecule(Lattice *self, Adsorbate *ads, char *newName);

char *get_carbon_site_name(int siteId);

/* Format and export lattice contents to .msi
 * if dest == NULL, use default
 */
void lattice_export_MSI(Lattice *self, char *pathName);
char *lattice_export_dest(Lattice *self, char *pathName);
