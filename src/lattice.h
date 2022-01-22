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
    char *name;
    int atomNum;
    Atom **atom_arr;
    Matrix *lattice_vectors;
    struct carbon_site carbon_sites[6];
    int metal_site_id;
    struct Lattice_vtable *vtable;
};

typedef struct _Lattice Lattice;

struct Lattice_vtable
{
    Atom *(*get_atom_by_Id)(Lattice *self, int);
    Matrix *(*get_lattice_coords)(Lattice *self);
    void (*update_atom_coords)(Lattice *self, Matrix *);
    Matrix *(*get_vector_ab)(Lattice *self, int, int);
    Matrix *(*get_carbon_chain_vector)(Lattice *self, int, int);
    Lattice *(*attach_molecule)(Lattice *self, Molecule *mol);
    void (*destroy)(Lattice *self);
};

/* Create a pointer to Lattice struct and malloc memory for it.
 * Initiate it with name (char *), atomNum (int), array of atoms (Atom **),
 * lattice_vectors (3x3 column-major matrix)
 */
Lattice *createLattice(char *name, int atomNum, Atom **atom_arr,
                       Matrix *lattice_vectors);

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
/* Returns a vector from two given atoms' coordinates */
Matrix *lattice_get_carbon_chain_vector(Lattice *self, int a, int b);
/* Add a Molecule into the Lattice. The atomId and treeId
 * of the atoms in mol will be updated to follow the order in current Lattice.
 * Returns a new Lattice struct pointer for future exports
 */
Lattice *lattice_attach_molecule(Lattice *self, Molecule *mol);
