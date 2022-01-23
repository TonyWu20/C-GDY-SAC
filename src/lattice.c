#include "lattice.h"
#include "molecule.h"
#include <matrix.h>
#include <stdlib.h>
#include <string.h>

struct carbon_site siteDict[] = {{"c1", 41},       {"c2", 42},
                                 {"c3", 54},       {"c4", 43},
                                 {"far_ring", 52}, {"near_ring", 40}};
struct Lattice_vtable lat_vtable = {lattice_get_carbon_chain_vector,
                                    lattice_get_carbon_metal_vector,
                                    lattice_attach_molecule, destroyLattice};
Lattice *createLattice(Molecule *mol, Matrix *lattice_vectors)
{
    Lattice *new = malloc(sizeof(Lattice));
    new->_mol = mol;
    new->lattice_vectors = lattice_vectors;
    memcpy(new->carbon_sites, siteDict, sizeof(struct carbon_site) * 6);
    new->metal_site_id = 73;
    new->vtable = &lat_vtable;
    return new;
}

void destroyLattice(Lattice *self)
{
    Molecule *mPtr = self->_mol;
    mPtr->vtable->destroy(mPtr);
    destroy_matrix(self->lattice_vectors);
    free(self->lattice_vectors);
    free(self);
}

Matrix *lattice_get_carbon_chain_vector(Lattice *self)
{
    return self->_mol->vtable->get_vector_ab(self->_mol, 41, 42);
}

Matrix *lattice_get_carbon_metal_vector(Lattice *self, int cId)
{
    return self->_mol->vtable->get_vector_ab(self->_mol, cId, 73);
}

Lattice *lattice_attach_molecule(Lattice *self, Adsorbate *ads)
{
    Molecule *lat_mol = self->_mol, *ads_mol = ads->_mol;
    Atom **cur_arr = lat_mol->atom_arr;
    int new_atomNum = lat_mol->atomNum + ads_mol->atomNum;
    cur_arr = realloc(cur_arr, sizeof(Atom *) * new_atomNum);
    Atom *lastAtom = lat_mol->vtable->get_atom_by_Id(lat_mol, lat_mol->atomNum);
    int lastId = lastAtom->vtable->get_atomId(lastAtom);
    for (int i = 0; i < ads_mol->atomNum; ++i)
    {
        Atom *cur = ads_mol->atom_arr[i];
        cur_arr[i + lat_mol->atomNum] = malloc(sizeof(Atom));
        memcpy(cur_arr[i + lat_mol->atomNum], cur, sizeof(*cur));
        cur_arr[i + lat_mol->atomNum]->atomId += lastId;
        cur_arr[i + lat_mol->atomNum]->treeId += lastId;
    }
    int molNameLen = strlen(ads_mol->name);
    int newNameLen = strlen(lat_mol->name) + molNameLen + 1;
    lat_mol->name = realloc(lat_mol->name, newNameLen);
    lat_mol->name = strcat(lat_mol->name, ads_mol->name);
    lat_mol->atom_arr = cur_arr;
    Lattice *new = createLattice(lat_mol, self->lattice_vectors);
    return new;
}
