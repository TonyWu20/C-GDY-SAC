#include "molecule.h"
#include "atom.h"
#include "matrix.h"
#include "my_maths.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// Implementation of Molecule struct

struct Molecule_vtable vtable = {
    Molecule_get_Atom_by_Id,     Molecule_get_coords,
    Molecule_update_Atom_coords, Molecule_get_vector_ab,
    Molecule_get_stem_vector,    Molecule_get_plane_normal,
    Molecule_make_upright,       destroyMolecule};

Molecule *createMolecule(char *name, int atomNum, Atom **atom_arr,
                         int coordAtomNum, int *coordAtomIds, int *stemAtomIds,
                         int *planeAtomIds)
{
    Molecule *newMol = malloc(sizeof(Molecule));
    newMol->name = strdup(name);
    newMol->atomNum = atomNum;
    newMol->atom_arr = atom_arr;
    newMol->coordAtomNum = coordAtomNum;
    newMol->coordAtomIds = malloc(coordAtomNum * sizeof(int));
    memcpy(newMol->coordAtomIds, coordAtomIds, sizeof(int) * coordAtomNum);
    memcpy(newMol->stemAtomIds, stemAtomIds, 2 * sizeof(int));
    memcpy(newMol->planeAtomIds, planeAtomIds, 3 * sizeof(int));
    newMol->vtable = &vtable;
    return newMol;
}

void destroyMolecule(Molecule *molPtr)
{
    free(molPtr->name);
    for (int i = 0; i < molPtr->atomNum; ++i)
    {
        Atom *cur = molPtr->atom_arr[i];
        cur->vtable->destroy(cur);
    }
    free(molPtr->atom_arr);
    free(molPtr->coordAtomIds);
    free(molPtr);
}

// Methods
Atom *Molecule_get_Atom_by_Id(Molecule *mPtr, int atomId)
{
    printf("%d\n", atomId);
    return mPtr->atom_arr[atomId - 1];
}

Matrix *Molecule_get_coords(Molecule *mPtr)
{
    Matrix *MolCoords = create_matrix(4, mPtr->atomNum);
    for (int i = 0; i < mPtr->atomNum; ++i)
    {
        Atom *cur = mPtr->atom_arr[i];
        for (int j = 0; j < 4; ++j)
        {
            MolCoords->value[j][i] = cur->vtable->get_coord(cur)->value[j][0];
        }
    }
    return MolCoords;
}

void Molecule_update_Atom_coords(Molecule *mPtr, Matrix *MolCoords)
{
    for (int i = 0; i < mPtr->atomNum; ++i)
    {
        Atom *cur = mPtr->atom_arr[i];
        cur->vtable->update_coord(cur, MolCoords->value[0][i],
                                  MolCoords->value[1][i],
                                  MolCoords->value[2][i]);
    }
}

Matrix *Molecule_get_vector_ab(Molecule *mPtr, int aId, int bId)
{
    Atom *a = mPtr->vtable->get_atom_by_Id(mPtr, aId);
    Atom *b = mPtr->vtable->get_atom_by_Id(mPtr, bId);
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

Matrix *Molecule_get_stem_vector(Molecule *mPtr)
{
    return mPtr->vtable->get_vector_ab(mPtr, mPtr->stemAtomIds[0],
                                       mPtr->stemAtomIds[1]);
}

Matrix *Molecule_get_plane_normal(Molecule *mPtr)
{
    Matrix *ba = mPtr->vtable->get_vector_ab(mPtr, mPtr->planeAtomIds[0],
                                             mPtr->planeAtomIds[1]);
    Matrix *ca = mPtr->vtable->get_vector_ab(mPtr, mPtr->planeAtomIds[0],
                                             mPtr->planeAtomIds[2]);
    double y_axis[] = {0, 1, 0, 1};
    Matrix *y_base = col_vector_view_array((double *)y_axis, 4);
    Matrix *normal = cross_product(ba, ca);
    destroy_matrix(ba);
    destroy_matrix(ca);
    destroy_matrix(y_base);
    free(ba);
    free(ca);
    free(y_base);
    return normal;
}

void Molecule_make_upright(Molecule *mol)
{
    Matrix *plane_normal = mol->vtable->get_plane_normal(mol); // malloced
    double y_axis[] = {0, 1, 0, 1};
    Matrix *y_base = col_vector_view_array(y_axis, 4);
    double rot_angle = vector_angle(plane_normal, y_base);
    Matrix *mol_coords = mol->vtable->get_mol_coords(mol);
    Matrix *new_coords;
    rotate_around_origin(mol_coords, rot_angle, 'X', &new_coords);
    mol->vtable->update_atom_coords(mol, new_coords);
    destroy_matrix(mol_coords);
    destroy_matrix(y_base);
    destroy_matrix(plane_normal);
    destroy_matrix(new_coords);
    free(mol_coords);
    free(y_base);
    free(plane_normal);
    free(new_coords);
}
