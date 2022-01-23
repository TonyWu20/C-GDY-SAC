#include "lattice.h"
#include "molecule.h"
#include <matrix.h>
#include <stdio.h>
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

Lattice *lattice_attach_molecule(Lattice *self, Adsorbate *ads, char *newName)
{
    Molecule *lat_mol = self->_mol, *ads_mol = ads->_mol;
    Atom **cur_arr = lat_mol->atom_arr;
    int new_atomNum = lat_mol->atomNum + ads_mol->atomNum;
    Atom *lastAtom = lat_mol->vtable->get_atom_by_Id(lat_mol, lat_mol->atomNum);
    int lastId = lastAtom->vtable->get_atomId(lastAtom);
    int i;
    Atom **new_arr = malloc(sizeof(Atom *) * new_atomNum);
    for (i = 0; i < lat_mol->atomNum; ++i)
    {
        Atom *cur = cur_arr[i];
        new_arr[i] = cur->vtable->dupAtom(cur);
    }
    for (i = 0; i < ads_mol->atomNum; ++i)
    {
        Atom *cur = ads_mol->atom_arr[i];
        new_arr[i + lat_mol->atomNum] = cur->vtable->dupAtom(cur);
        new_arr[i + lat_mol->atomNum]->atomId += lastId;
        new_arr[i + lat_mol->atomNum]->treeId += lastId;
    }
    Molecule *resMol = createMolecule(newName, new_atomNum, new_arr);
    Matrix *lattice_vectors = create_matrix(3, 3);
    copy_matrix(self->lattice_vectors, &lattice_vectors);
    Lattice *new = createLattice(resMol, lattice_vectors);
    return new;
}

char *get_carbon_site_name(int siteId)
{
    switch (siteId)
    {
    case 41:
        return siteDict[0].name;
    case 42:
        return siteDict[1].name;
    case 54:
        return siteDict[2].name;
    case 43:
        return siteDict[3].name;
    case 52:
        return siteDict[4].name;
    case 40:
        return siteDict[5].name;
    default:
        printf("Incorrect carbon site Id %d\n", siteId);
        return NULL;
    }
}
