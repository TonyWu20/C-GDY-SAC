#include "molecule.h"
#include "atom.h"
#include "database/ads_database.h"
#include "misc.h"
#include "my_maths.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#define NULLSITE 255

static bool faceUp(Adsorbate *self);
enum
{
    C1 = 41,
    C2 = 42,
    C3 = 54,
    C4 = 53,
    FR = 52,
    NR = 40,
    M = 73
} site_codes;

// Implementation of Molecule struct

struct Molecule_vtable vtable = {
    Molecule_get_Atom_by_Id,    Molecule_get_vector_ab,
    Molecule_get_centroid_ab,   Molecule_apply_rotation,
    Molecule_apply_translation, Molecule_textblock,
    Molecule_duplicate,         destroyMolecule};
struct Adsorbate_vtable ads_vtable = {
    Adsorbate_get_stem_vector, Adsorbate_get_plane_normal,
    Adsorbate_make_upright,    Adsorbate_export_MSI,
    Adsorbate_duplicate,       destroyAdsorbate};

static bool faceUp(Adsorbate *self)
{
    Molecule *mol = self->mol;
    Atom *cdAtom = mol->vtable->get_atom_by_Id(mol, self->coordAtomIds[0]);
    Atom *upAtom = mol->vtable->get_atom_by_Id(mol, self->upperAtomId);
    if (cdAtom->coord.z < upAtom->coord.z)
        return true;
    else
        return false;
}

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
                           int *planeAtomIds, bool bVer, bool bSym,
                           int upperAtomId, char *pathName)
{
    Adsorbate *ads = malloc(sizeof(Adsorbate));
    ads->mol = newMol;
    ads->coordAtomNum = coordAtomNum;
    ads->coordAtomIds = malloc(coordAtomNum * sizeof(int));
    memcpy(ads->coordAtomIds, coordAtomIds, sizeof(int) * coordAtomNum);
    memcpy(ads->stemAtomIds, stemAtomIds, 2 * sizeof(int));
    memcpy(ads->planeAtomIds, planeAtomIds, 3 * sizeof(int));
    ads->vtable = &ads_vtable;
    ads->bVertical = bVer;
    ads->bSym = bSym;
    ads->taskLists = createTasks(ads);
    ads->upperAtomId = upperAtomId;
    ads->pathName = strdup(pathName);
    return ads;
}

Adsorbate *Adsorbate_duplicate(Adsorbate *self)
{
    Molecule *molCopy = self->mol->vtable->duplicate(self->mol);
    Adsorbate *dup =
        createAdsorbate(molCopy, self->coordAtomNum, self->coordAtomIds,
                        self->stemAtomIds, self->planeAtomIds, self->bVertical,
                        self->bSym, self->upperAtomId, self->pathName);
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
    ads->mol->vtable->destroy(ads->mol);
    free(ads->coordAtomIds);
    for (int i = 0; i < ads->taskLists->taskNum; ++i)
        free(ads->taskLists->tasks[i]);
    free(ads->taskLists->tasks);
    free(ads->taskLists);
    free(ads->pathName);
    free(ads);
}

// Methods
Atom *Molecule_get_Atom_by_Id(Molecule *mPtr, int atomId)
{
    return mPtr->atom_arr[atomId - 1];
}

simd_double3 Molecule_get_vector_ab(Molecule *mPtr, int aId, int bId)
{
    Atom *a = mPtr->vtable->get_atom_by_Id(mPtr, aId);
    Atom *b = mPtr->vtable->get_atom_by_Id(mPtr, bId);
    return b->coord - a->coord;
}

simd_double3 Molecule_get_centroid_ab(Molecule *mPtr, int aId, int bId)
{
    Atom *a = mPtr->vtable->get_atom_by_Id(mPtr, aId);
    Atom *b = mPtr->vtable->get_atom_by_Id(mPtr, bId);
    simd_double3 ab[2] = {a->coord, b->coord};
    simd_double3 c_ab = simd_centroid_of_points(ab, 2);
    return c_ab;
}

void Molecule_apply_rotation(Molecule *self, simd_quatd rotation)
{
    for (int i = 0; i < self->atomNum; ++i)
    {
        Atom *curAtom = self->atom_arr[i];
        curAtom->coord = simd_act(rotation, curAtom->coord);
    }
}

void Molecule_apply_translation(Molecule *self, simd_double4x4 *transMat)
{
    for (int i = 0; i < self->atomNum; ++i)
    {
        Atom *curAtom = self->atom_arr[i];
        simd_double4 tmp_4dCoord = simd_make_double4(curAtom->coord, 1);
        simd_double4 translated = simd_mul(*transMat, tmp_4dCoord);
        curAtom->coord = simd_make_double3(translated);
    }
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

simd_double3 Adsorbate_get_stem_vector(Adsorbate *adsPtr)
{
    Molecule *mPtr = adsPtr->mol;
    return mPtr->vtable->get_vector_ab(mPtr, adsPtr->stemAtomIds[0],
                                       adsPtr->stemAtomIds[1]);
}

simd_double3 Adsorbate_get_plane_normal(Adsorbate *adsPtr)
{
    Molecule *mPtr = adsPtr->mol;
    simd_double3 ba = mPtr->vtable->get_vector_ab(mPtr, adsPtr->planeAtomIds[0],
                                                  adsPtr->planeAtomIds[1]);
    simd_double3 ca = mPtr->vtable->get_vector_ab(mPtr, adsPtr->planeAtomIds[0],
                                                  adsPtr->planeAtomIds[2]);
    simd_double3 pNormal = simd_normalize(simd_cross(ba, ca));
    return pNormal;
}

void Adsorbate_make_upright(Adsorbate *adsPtr)
{
    Molecule *mPtr = adsPtr->mol;
    simd_double3 stemVector = adsPtr->vtable->get_stem_vector(adsPtr);
    if (adsPtr->bVertical == false)
    {
        simd_double3 zAxis = {0, 0, 1};
        double angle = simd_vector_angle(stemVector, zAxis);
        simd_double3 rotAxis = simd_normalize(simd_cross(stemVector, zAxis));
        simd_quatd rotation = simd_quaternion(angle, rotAxis);
        mPtr->vtable->rotateMol(mPtr, rotation);
    }
    else
    {
        simd_double3 plane_normal =
            adsPtr->vtable->get_plane_normal(adsPtr); // malloced
        printf("%f, %f, %f\n", plane_normal.x, plane_normal.y, plane_normal.z);
        simd_double3 plane_proj_XY = simd_normalize(
            simd_make_double3(plane_normal.x, plane_normal.y, 0));
        double rot_angle = simd_vector_angle(plane_normal, plane_proj_XY);
        simd_double3 rotAxis = simd_cross(plane_normal, plane_proj_XY);
        double rotAxis_stem_angle = simd_vector_angle(rotAxis, stemVector);
        simd_quatd rotation_upright_q;
        if (rotAxis_stem_angle < PI / 2)
        {
            rotation_upright_q =
                simd_quaternion(rot_angle, simd_normalize(stemVector));
        }
        else
        {
            rotation_upright_q =
                simd_quaternion(rot_angle, simd_normalize(-1 * stemVector));
        }
        mPtr->vtable->rotateMol(mPtr, rotation_upright_q);
        if (!faceUp(adsPtr))
        {
            stemVector = adsPtr->vtable->get_stem_vector(adsPtr);
            simd_quatd invert = simd_quaternion(PI, simd_normalize(stemVector));
            mPtr->vtable->rotateMol(mPtr, invert);
        }
    }
}

void Adsorbate_export_MSI(Adsorbate *self, char *dest)
{
    char header_line[] = "# MSI CERIUS2 DataModel File Version 4 0\n";
    char model_start[] = "(1 Model\n";
    char model_end[] = ")\n";
    int lineSize = self->mol->atomNum + 3;
    char **content_lines = malloc(sizeof(char *) * (lineSize));
    content_lines[0] = strdup(header_line);
    content_lines[1] = strdup(model_start);
    char **atoms = self->mol->vtable->export_text(self->mol);
    for (int i = 0; i < self->mol->atomNum; ++i)
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
    int exportLen = 1 + snprintf(NULL, 0, "%s%s.msi", dest, self->mol->name);
    char *exportName = malloc(exportLen);
    snprintf(exportName, exportLen, "%s%s.msi", dest, self->mol->name);
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
        if (self->bSym)
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
        }
        else
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
