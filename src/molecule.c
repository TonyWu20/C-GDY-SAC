#include "molecule.h"
#include "atom.h"
#include "matrix.h"
#include "my_maths.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// Implementation of Molecule struct

struct Molecule_vtable vtable = {Molecule_get_Atom_by_Id,
                                 Molecule_get_coords,
                                 Molecule_update_Atom_coords,
                                 Molecule_get_vector_ab,
                                 Molecule_get_centroid_ab,
                                 Molecule_apply_transformation,
                                 destroyMolecule};
struct Adsorbate_vtable ads_vtable = {Adsorbate_get_stem_vector,
                                      Adsorbate_get_plane_normal,
                                      Adsorbate_make_upright, destroyAdsorbate};

Molecule *createMolecule(char *name, int atomNum, Atom **atom_arr)
{
    Molecule *newMol = malloc(sizeof(Molecule));
    newMol->name = strdup(name);
    newMol->atomNum = atomNum;
    newMol->atom_arr = atom_arr;
    newMol->vtable = &vtable;
    return newMol;
}

Adsorbate *createAdsorbate(Molecule *newMol, int coordAtomNum,
                           int *coordAtomIds, int *stemAtomIds,
                           int *planeAtomIds)
{
    Adsorbate *ads = malloc(sizeof(Adsorbate));
    ads->_mol = newMol;
    ads->coordAtomNum = coordAtomNum;
    ads->coordAtomIds = malloc(coordAtomNum * sizeof(int));
    memcpy(ads->coordAtomIds, coordAtomIds, sizeof(int) * coordAtomNum);
    memcpy(ads->stemAtomIds, stemAtomIds, 2 * sizeof(int));
    memcpy(ads->planeAtomIds, planeAtomIds, 3 * sizeof(int));
    ads->ads_vtable = &ads_vtable;
    return ads;
}

void destroyMolecule(Molecule *self)
{
    free(self->name);
    for (int i = 0; i < self->atomNum; ++i)
    {
        Atom *cur = self->atom_arr[i];
        cur->vtable->destroy(cur);
    }
    free(self->atom_arr);
    free(self);
}
void destroyAdsorbate(Adsorbate *ads)
{
    ads->_mol->vtable->destroy(ads->_mol);
    free(ads->coordAtomIds);
    free(ads);
}

// Methods
Atom *Molecule_get_Atom_by_Id(Molecule *mPtr, int atomId)
{
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

double *Molecule_get_centroid_ab(Molecule *mPtr, int aId, int bId)
{
    Atom *a = mPtr->vtable->get_atom_by_Id(mPtr, aId);
    Atom *b = mPtr->vtable->get_atom_by_Id(mPtr, bId);
    Matrix *ab = create_matrix(3, 2);
    for (int i = 0; i < 3; ++i)
    {
        ab->value[i][0] = a->coord->value[i][0];
        ab->value[i][1] = b->coord->value[i][0];
    }
    double *c_ab = centroid_of_points(ab);
    destroy_matrix(ab);
    free(ab);
    return c_ab;
}

void Molecule_apply_transformation(Molecule *mPtr, Matrix *trans_mat,
                                   void (*transform_func)(Matrix *trans_mat,
                                                          Matrix *coords,
                                                          Matrix **result))
{
    Matrix *mol_coords = mPtr->vtable->get_mol_coords(mPtr);
    Matrix *result;
    transform_func(trans_mat, mol_coords, &result);
    mPtr->vtable->update_atom_coords(mPtr, result);
    destroy_matrix(mol_coords);
    destroy_matrix(result);
    free(result);
    free(mol_coords);
}

Matrix *Adsorbate_get_stem_vector(Adsorbate *adsPtr)
{
    Molecule *mPtr = adsPtr->_mol;
    return mPtr->vtable->get_vector_ab(mPtr, adsPtr->stemAtomIds[0],
                                       adsPtr->stemAtomIds[1]);
}

Matrix *Adsorbate_get_plane_normal(Adsorbate *adsPtr)
{
    Molecule *mPtr = adsPtr->_mol;
    Matrix *ba = mPtr->vtable->get_vector_ab(mPtr, adsPtr->planeAtomIds[0],
                                             adsPtr->planeAtomIds[1]);
    Matrix *ca = mPtr->vtable->get_vector_ab(mPtr, adsPtr->planeAtomIds[0],
                                             adsPtr->planeAtomIds[2]);
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

void Adsorbate_make_upright(Adsorbate *adsPtr)
{
    Molecule *mPtr = adsPtr->_mol;
    Matrix *plane_normal =
        adsPtr->ads_vtable->get_plane_normal(adsPtr); // malloced
    double y_axis[] = {0, 1, 0, 1};
    Matrix *y_base = col_vector_view_array(y_axis, 4);
    double rot_angle = vector_angle(plane_normal, y_base);
    Matrix *rot_mat = rotationMatrix(rot_angle, 'X');
    mPtr->vtable->apply_transformation(mPtr, rot_mat, rotate_around_origin);
    destroy_matrix(y_base);
    destroy_matrix(plane_normal);
    destroy_matrix(rot_mat);
    free(y_base);
    free(plane_normal);
    free(rot_mat);
}
