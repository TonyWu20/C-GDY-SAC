#include "molecule.h"
#include "atom.h"
#include "matrix.h"
#include "misc.h"
#include "my_maths.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#define NULLSITE 255
enum
{
    C1 = 41,
    C2 = 42,
    C3 = 54,
    C4 = 53,
    FR = 52,
    NR = 40,
    M = 73
};

int task_cd2_sym[][2] = {{C1, C2}, {C2, C3}, {C3, C4}, {C4, FR},
                         {NR, C1}, {C1, M},  {C2, M}};
int task_cd2_asym[][2] = {{C1, C2}, {C2, C1}, {C2, C3}, {C3, C2}, {C3, C4},
                          {C4, C3}, {C4, FR}, {FR, C4}, {NR, C1}, {C1, NR},
                          {C1, M},  {M, C1},  {C2, M},  {M, C2}};
int task_cd1[][2] = {{C1, NULLSITE}, {C2, NULLSITE}, {C3, NULLSITE},
                     {C4, NULLSITE}, {FR, NULLSITE}, {NR, NULLSITE},
                     {M, NULLSITE}};

// Implementation of Molecule struct

struct Molecule_vtable vtable = {Molecule_get_Atom_by_Id,
                                 Molecule_get_coords,
                                 Molecule_update_Atom_coords,
                                 Molecule_get_vector_ab,
                                 Molecule_get_centroid_ab,
                                 Molecule_apply_transformation,
                                 Molecule_textblock,
                                 Molecule_duplicate,
                                 destroyMolecule};
struct Adsorbate_vtable ads_vtable = {
    Adsorbate_get_stem_vector, Adsorbate_get_plane_normal,
    Adsorbate_make_upright,    Adsorbate_export_MSI,
    Adsorbate_duplicate,       destroyAdsorbate};

Molecule *createMolecule(char *name, int atomNum, Atom **atom_arr)
{
    Molecule *newMol = malloc(sizeof(Molecule));
    newMol->name = strdup(name);
    newMol->atomNum = atomNum;
    newMol->atom_arr = atom_arr;
    newMol->vtable = &vtable;
    return newMol;
}

Molecule *Molecule_duplicate(Molecule *self)
{
    Atom **dupAtoms = malloc(sizeof(Atom *) * self->atomNum);
    for (int i = 0; i < self->atomNum; ++i)
    {
        dupAtoms[i] = dupAtom(self->atom_arr[i]);
    }
    Molecule *dup = createMolecule(self->name, self->atomNum, dupAtoms);
    return dup;
}

Adsorbate *createAdsorbate(Molecule *newMol, int coordAtomNum,
                           int *coordAtomIds, int *stemAtomIds,
                           int *planeAtomIds, int bSym)
{
    Adsorbate *ads = malloc(sizeof(Adsorbate));
    ads->_mol = newMol;
    ads->coordAtomNum = coordAtomNum;
    ads->coordAtomIds = malloc(coordAtomNum * sizeof(int));
    memcpy(ads->coordAtomIds, coordAtomIds, sizeof(int) * coordAtomNum);
    memcpy(ads->stemAtomIds, stemAtomIds, 2 * sizeof(int));
    memcpy(ads->planeAtomIds, planeAtomIds, 3 * sizeof(int));
    ads->ads_vtable = &ads_vtable;
    ads->bSym = bSym;
    ads->taskLists = createTasks(ads);
    return ads;
}

Adsorbate *Adsorbate_duplicate(Adsorbate *self)
{
    Molecule *molCopy = self->_mol->vtable->duplicate(self->_mol);
    Adsorbate *dup =
        createAdsorbate(molCopy, self->coordAtomNum, self->coordAtomIds,
                        self->stemAtomIds, self->planeAtomIds, self->bSym);
    return dup;
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
    for (int i = 0; i < ads->taskLists->taskNum; ++i)
        free(ads->taskLists->tasks[i]);
    free(ads->taskLists->tasks);
    free(ads->taskLists);
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

char **Molecule_textblock(Molecule *self)
{
    char **atom_blocks = malloc(sizeof(char *) * self->atomNum);
    for (int i = 0; i < self->atomNum; ++i)
    {
        atom_blocks[i] = self->atom_arr[i]->vtable->export_text(
            self->atom_arr[i]); // malloc from func
    }
    return atom_blocks;
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
    double rot_angle;
    plane_normal->value[0][0] = 0;

    Matrix *stemVector = adsPtr->ads_vtable->get_stem_vector(adsPtr);
    if (stemVector->value[0][0] < 0)
        rot_angle = 2 * PI - vector_angle(plane_normal, y_base);
    else
        rot_angle = vector_angle(plane_normal, y_base);

    destroy_matrix(y_base);
    destroy_matrix(plane_normal);
    free(y_base);
    free(plane_normal);

    Matrix *rot_mat = rotate_angle_around_axis(stemVector, rot_angle);
    mPtr->vtable->apply_transformation(mPtr, rot_mat, rotate_around_origin);
    // Release rot_mat
    destroy_matrix(rot_mat);
    free(rot_mat);
    Matrix *after_cd = mPtr->vtable->get_mol_coords(mPtr);
    double *centroid = centroid_of_points(after_cd);
    destroy_matrix(after_cd);
    free(after_cd);
    double cd_to_centroid_z =
        mPtr->vtable->get_atom_by_Id(mPtr, adsPtr->coordAtomIds[0])
            ->coord->value[2][0] -
        centroid[2];
    free(centroid);
    if (cd_to_centroid_z > 0)
    {
        Matrix *invert = rotate_angle_around_axis(stemVector, PI);
        mPtr->vtable->apply_transformation(mPtr, invert, rotate_around_origin);
        destroy_matrix(invert);
        free(invert);
    }
    destroy_matrix(stemVector);
    free(stemVector);
}

void Adsorbate_export_MSI(Adsorbate *self, char *dest)
{
    char header_line[] = "# MSI CERIUS2 DataModel File Version 4 0\n";
    char model_start[] = "(1 Model\n";
    char model_end[] = ")\n";
    int lineSize = self->_mol->atomNum + 3;
    char **content_lines = malloc(sizeof(char *) * (lineSize));
    content_lines[0] = strdup(header_line);
    content_lines[1] = strdup(model_start);
    char **atoms = self->_mol->vtable->export_text(self->_mol);
    for (int i = 0; i < self->_mol->atomNum; ++i)
    {
        content_lines[i + 2] = atoms[i];
    }
    content_lines[lineSize - 1] = strdup(model_end);
    int destUndefined = 0;
    if (!dest)
    {
        dest = strdup("./C2_pathways_ads/exported/");
        destUndefined = 1;
    }
    createDirectory(dest);
    int dirLen = strlen(dest);
    int adsNameLen = strlen(self->_mol->name);
    char *exportName = malloc(dirLen + adsNameLen + 5);
    snprintf(exportName, dirLen + adsNameLen + 5, "%s%s.msi", dest,
             self->_mol->name);
    FILE *writeFile = fopen(exportName, "w");
    for (int i = 0; i < lineSize; ++i)
    {
        fputs(content_lines[i], writeFile);
        free(content_lines[i]);
    }
    fclose(writeFile);
    free(content_lines);
    free(exportName);
    free(atoms);
    if (destUndefined)
        free(dest);
}

struct taskTable *createTasks(Adsorbate *self)
{
    struct taskTable *tab = NULL;
    switch (self->coordAtomNum)
    {
    case 2:
        switch (self->bSym)
        {
        case 1:
        {
            tab = malloc(sizeof(*tab));
            tab->taskNum = sizeof(task_cd2_sym) / sizeof(task_cd2_sym[0]);
            tab->tasks = malloc(sizeof(int *) * tab->taskNum);
            for (int i = 0; i < tab->taskNum; ++i)
            {
                tab->tasks[i] = malloc(sizeof(int) * 2);
                for (int j = 0; j < 2; ++j)
                {
                    tab->tasks[i][j] = task_cd2_sym[i][j];
                }
            }
            break;
        }
        case 0:
        {
            tab = malloc(sizeof(*tab));
            tab->taskNum = sizeof(task_cd2_asym) / sizeof(task_cd2_asym[0]);
            tab->tasks = malloc(sizeof(int *) * tab->taskNum);
            for (int i = 0; i < tab->taskNum; ++i)
            {
                tab->tasks[i] = malloc(sizeof(int) * 2);
                for (int j = 0; j < 2; ++j)
                    tab->tasks[i][j] = task_cd2_asym[i][j];
            }
            break;
        }
        default:
            break;
        }
        break;
    case 1:
        tab = malloc(sizeof(*tab));
        tab->taskNum = sizeof(task_cd1) / sizeof(task_cd1[0]);
        tab->tasks = malloc(sizeof(int *) * 7);
        for (int i = 0; i < 7; ++i)
        {
            tab->tasks[i] = malloc(sizeof(int) * 2);
            for (int j = 0; j < 2; ++j)
                tab->tasks[i][j] = task_cd1[i][j];
        }
        break;
    default:
        break;
    }
    return tab;
}
