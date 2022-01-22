#include "lattice.h"
#include <matrix.h>
#include <stdlib.h>
#include <string.h>

struct carbon_site siteDict[] = {{"c1", 41},       {"c2", 42},
                                 {"c3", 54},       {"c4", 43},
                                 {"far_ring", 52}, {"near_ring", 40}};
struct Lattice_vtable lat_vtable = {lattice_get_atom_by_Id,
                                    lattice_get_coords,
                                    lattice_update_Atom_coords,
                                    lattice_get_vector_ab,
                                    lattice_get_carbon_chain_vector,
                                    lattice_attach_molecule,
                                    destroyLattice};
Lattice *createLattice(char *name, int atomNum, Atom **atom_arr,
                       Matrix *lattice_vectors)
{
    Lattice *new = malloc(sizeof(Lattice));
    new->name = strdup(name);
    new->atomNum = atomNum;
    new->atom_arr = atom_arr;
    new->lattice_vectors = lattice_vectors;
    memcpy(new->carbon_sites, siteDict, sizeof(struct carbon_site) * 6);
    new->metal_site_id = 73;
    new->vtable = &lat_vtable;
    return new;
}

void destroyLattice(Lattice *self)
{
    free(self->name);
    for (int i = 0; i < self->atomNum; ++i)
    {
        Atom *cur = self->atom_arr[i];
        cur->vtable->destroy(cur);
    }
    destroy_matrix(self->lattice_vectors);
    free(self->lattice_vectors);
    free(self);
}

Atom *lattice_get_atom_by_Id(Lattice *self, int atomId)
{
    return self->atom_arr[atomId - 1];
}

Matrix *lattice_get_coords(Lattice *self)
{
    Matrix *latCoords = create_matrix(4, self->atomNum);
    for (int i = 0; i < self->atomNum; ++i)
    {
        Atom *cur = self->atom_arr[i];
        for (int j = 0; j < 4; ++j)
        {
            latCoords->value[j][i] = cur->vtable->get_coord(cur)->value[j][0];
        }
    }
    return latCoords;
}

void lattice_update_Atom_coords(Lattice *self, Matrix *new)
{
    for (int i = 0; i < self->atomNum; ++i)
    {
        Atom *cur = self->atom_arr[i];
        cur->vtable->update_coord(cur, new->value[0][i], new->value[1][i],
                                  new->value[2][i]);
    }
    return;
}
Matrix *lattice_get_vector_ab(Lattice *self, int aId, int bId)
{
    Atom *a = self->vtable->get_atom_by_Id(self, aId);
    Atom *b = self->vtable->get_atom_by_Id(self, bId);
    Matrix *a_coord = a->vtable->get_coord(a);
    Matrix *b_coord = b->vtable->get_coord(b);
    Matrix *minus_b = create_matrix(b_coord->lines, b_coord->columns);
    copy_matrix(b_coord, &minus_b);
    multiply_matrix_with_scalar(minus_b, -1.0);
    Matrix *res;
    add_matrices(a_coord, minus_b, &res);
    destroy_matrix(minus_b);
    free(minus_b);
    res->value[3][0] = 1;
    return res;
}

Matrix *lattice_get_carbon_chain_vector(Lattice *self, int a, int b)
{
    return self->vtable->get_vector_ab(self, a, b);
}

Lattice *lattice_attach_molecule(Lattice *self, Molecule *mol)
{
    Atom **cur_arr = self->atom_arr;
    int new_atomNum = self->atomNum + mol->atomNum;
    cur_arr = realloc(cur_arr, sizeof(Atom *) * new_atomNum);
    Atom *lastAtom = self->vtable->get_atom_by_Id(self, self->atomNum);
    int lastId = lastAtom->vtable->get_atomId(lastAtom);
    for (int i = 0; i < mol->atomNum; ++i)
    {
        Atom *cur = mol->atom_arr[i];
        cur_arr[i + self->atomNum] = malloc(sizeof(Atom));
        memcpy(cur_arr[i + self->atomNum], cur, sizeof(*cur));
        cur_arr[i + self->atomNum]->atomId += lastId;
        cur_arr[i + self->atomNum]->treeId += lastId;
    }
    Lattice *new = createLattice(self->name, self->atomNum, cur_arr,
                                 self->lattice_vectors);
    return new;
}
