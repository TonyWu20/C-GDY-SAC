#include "molecule.h"

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

int task_cd2[][2] = {{C1, C2}, {C2, C3}, {C3, C4}, {C4, FR},
                     {NR, C1}, {C1, M},  {C2, M}};
int task_cd1[][2] = {{C1, NULLSITE}, {C2, NULLSITE}, {C3, NULLSITE},
                     {C4, NULLSITE}, {FR, NULLSITE}, {NR, NULLSITE},
                     {M, NULLSITE}};

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
                           int *planeAtomIds, bool bSym, bool bVer,
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
    ads->bSym = bSym;
    ads->bVer = bVer;
    ads->upperAtomId = upperAtomId;
    ads->pathName = strdup(pathName);
    return ads;
}

Adsorbate *Adsorbate_duplicate(Adsorbate *self)
{
    Molecule *molCopy = self->mol->vtable->duplicate(self->mol);
    Adsorbate *dup =
        createAdsorbate(molCopy, self->coordAtomNum, self->coordAtomIds,
                        self->stemAtomIds, self->planeAtomIds, self->bSym,
                        self->bVer, self->upperAtomId, self->pathName);
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
    free(ads->pathName);
    free(ads);
}

// Methods
Atom *Molecule_get_Atom_by_Id(Molecule *mPtr, int atomId)
{
    return mPtr->atom_arr[atomId - 1];
}

vec_double3 Molecule_get_vector_ab(Molecule *mPtr, int aId, int bId)
{
    Atom *a = mPtr->vtable->get_atom_by_Id(mPtr, aId);
    Atom *b = mPtr->vtable->get_atom_by_Id(mPtr, bId);
    vec_double3 a_coord = a->vtable->get_coord(a);
    vec_double3 b_coord = b->vtable->get_coord(b);
    vec_double3 res = vec_sub(b_coord, a_coord);
    return res;
}

vec_double3 Molecule_get_centroid_ab(Molecule *mPtr, int aId, int bId)
{
    Atom *a = mPtr->vtable->get_atom_by_Id(mPtr, aId);
    Atom *b = mPtr->vtable->get_atom_by_Id(mPtr, bId);
    vec_double3 points[] = {a->coord, b->coord};
    vec_double3 c_ab = vec_centroid(points, 2);
    return c_ab;
}
void Molecule_apply_rotation(Molecule *self, vec_quatd q)
{
    for (int i = 0; i < self->atomNum; ++i)
    {
        Atom *cur = self->atom_arr[i];
        cur->coord = vec_act(q, cur->coord);
    }
    return;
}
void Molecule_apply_translation(Molecule *self, matrix_double4x4 T)
{
    for (int i = 0; i < self->atomNum; ++i)
    {
        Atom *cur = self->atom_arr[i];
        cur->coord =
            vec_make_double3(vec_mul(T, vec_make_double4(cur->coord, 1)));
    }
    return;
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

vec_double3 Adsorbate_get_stem_vector(Adsorbate *adsPtr)
{
    Molecule *mPtr = adsPtr->mol;
    return mPtr->vtable->get_vector_ab(mPtr, adsPtr->stemAtomIds[0],
                                       adsPtr->stemAtomIds[1]);
}

vec_double3 Adsorbate_get_plane_normal(Adsorbate *adsPtr)
{
    Molecule *mPtr = adsPtr->mol;
    vec_double3 ba = mPtr->vtable->get_vector_ab(mPtr, adsPtr->planeAtomIds[0],
                                                 adsPtr->planeAtomIds[1]);
    vec_double3 ca = mPtr->vtable->get_vector_ab(mPtr, adsPtr->planeAtomIds[0],
                                                 adsPtr->planeAtomIds[2]);
    vec_double3 normal = vec_cross(ba, ca);
    return normal;
}

void Adsorbate_make_upright(Adsorbate *adsPtr)
{
    Molecule *mPtr = adsPtr->mol;
    vec_double3 stemVector = adsPtr->vtable->get_stem_vector(adsPtr);
    vec_double3 planeNormal = adsPtr->vtable->get_plane_normal(adsPtr);
    vec_double3 planeNormalProj_XY =
        vec_make_double3(planeNormal.x, planeNormal.y, 0);
    double toRotate = vec_angle_uv(planeNormal, planeNormalProj_XY);
    vec_double3 rotAxis = vec_cross(planeNormal, planeNormalProj_XY);
    vec_quatd quatd_rot;
    if (vec_dot(stemVector, rotAxis) >
        0) // stemVector points to correct direction
    {
        quatd_rot = vec_make_quaternion(toRotate, vec_normalize(stemVector));
    }
    else // the reverse stemVector is the true rotate axis
    {
        quatd_rot = vec_make_quaternion(toRotate,
                                        vec_normalize(vec_mul(-1, stemVector)));
    }
    mPtr->vtable->rotateMol(mPtr, quatd_rot);
    if (!faceUp(adsPtr))
    {
        stemVector = adsPtr->vtable->get_stem_vector(adsPtr);
        quatd_rot = vec_make_quaternion(PI, stemVector);
        mPtr->vtable->rotateMol(mPtr, quatd_rot);
    }
    return;
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
    int dirLen = strlen(dest);
    int adsNameLen = strlen(self->mol->name);
    char *exportName = malloc(dirLen + adsNameLen + 5);
    snprintf(exportName, dirLen + adsNameLen + 5, "%s%s.msi", dest,
             self->mol->name);
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
