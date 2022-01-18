#include "molecule.h"
#include "matrix.h"
#include <stdlib.h>
#include <string.h>

// Implementation of Molecule struct
struct _Molecule
{
    char *name;
    int atomNum;
    Atom **atom_arr;
    int coordAtomNum;
    int *coordAtomIds;
    int *stemAtomIds;
    int *planeAtomIds;
};

Molecule *createMolecule(char *name, int atomNum, Atom **atom_arr,
                         int coordAtomNum, int *coordAtomIds, int *stemAtomIds,
                         int *planeAtomIds)
{
    Molecule *newMol = malloc(sizeof(Molecule));
    newMol->name = strdup(name);
    newMol->atomNum = atomNum;
    newMol->atom_arr = calloc(atomNum, sizeof(Atom *));
    memcpy(newMol->atom_arr, atom_arr, sizeof(Atom *) * atomNum);
    newMol->coordAtomNum = coordAtomNum;
    newMol->coordAtomIds = malloc(coordAtomNum * sizeof(int));
    memcpy(newMol->coordAtomIds, coordAtomIds, sizeof(int) * coordAtomNum);
    newMol->stemAtomIds = malloc(2 * sizeof(int));
    memcpy(newMol->stemAtomIds, stemAtomIds, 2 * sizeof(int));
    newMol->planeAtomIds = malloc(3 * sizeof(int));
    memcpy(newMol->planeAtomIds, planeAtomIds, 3 * sizeof(int));
    return newMol;
}

void destroyMolecule(Molecule *molPtr)
{
    free(molPtr->name);
    for (int i = 0; i < molPtr->atomNum; ++i)
    {
        destroyAtom(molPtr->atom_arr[i]);
    }
    free(molPtr->atom_arr);
    free(molPtr->coordAtomIds);
    free(molPtr->stemAtomIds);
    free(molPtr->planeAtomIds);
    free(molPtr);
}

// Methods
Atom *Molecule_get_Atom_by_Id(Molecule *mPtr, int treeId)
{
    return mPtr->atom_arr[treeId];
}

Matrix **Molecule_get_coords(Molecule *mPtr)
{
    Matrix **MolCoords = calloc(mPtr->atomNum, sizeof(Matrix *));
    for (int i = 0; i < mPtr->atomNum; ++i)
    {
        MolCoords[i] = Atom_get_coord(mPtr->atom_arr[i]);
    }
    return MolCoords;
}

void Molecule_update_Atom_coords(Molecule *mPtr, Matrix **MolCoords)
{
    for (int i = 0; i < mPtr->atomNum; ++i)
    {
        Atom_update_coord(mPtr->atom_arr[i], MolCoords[i]);
    }
}

Matrix *Molecule_get_vector_ab(Molecule *mPtr, int aId, int bId)
{
    Atom *a = Molecule_get_Atom_by_Id(mPtr, aId);
    Atom *b = Molecule_get_Atom_by_Id(mPtr, bId);
    Matrix *a_coord = Atom_get_coord(a);
    Matrix *b_coord = Atom_get_coord(b);
    Matrix *minus_b = create_matrix(b_coord->lines, b_coord->columns);
    copy_matrix(b_coord, &minus_b);
    multiply_matrix_with_scalar(minus_b, -1.0);
    Matrix *res;
    add_matrices(a_coord, minus_b, &res);
    destroy_matrix(minus_b);
    res->value[0][3] = 1;
    return res;
}

Matrix *Molecule_get_stem_vector(Molecule *mPtr)
{
    return Molecule_get_vector_ab(mPtr, mPtr->stemAtomIds[0],
                                  mPtr->stemAtomIds[1]);
}

/*Matrix *Molecule_get_plane_normal(Molecule *mPtr)*/
/*{*/
/*Matrix *ba = Molecule_get_vector_ab(mPtr, mPtr->planeAtomIds[0],*/
/*mPtr->planeAtomIds[1]);*/
/*Matrix *ca = Molecule_get_vector_ab(mPtr, mPtr->planeAtomIds[0],*/
/*mPtr->planeAtomIds[2]);*/
/*}*/
