#include "lattice.h"
#include "misc.h"
#include "molecule.h"
#include "my_maths.h"
#include <matrix.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>

struct carbon_site siteDict[] = {{"c1", 41}, {"c2", 42}, {"c3", 54}, {"c4", 43},
                                 {"FR", 52}, {"NR", 40}, {"M", 73}};
struct Lattice_vtable lat_vtable = {lattice_get_carbon_chain_vector,
                                    lattice_get_carbon_metal_vector,
                                    lattice_attach_molecule,
                                    lattice_rotate_to_standard_orientation,
                                    lattice_modify_metal_element,
                                    lattice_export_dest,
                                    lattice_export_MSI,
                                    destroyLattice};
Lattice *createLattice(Molecule *mol, matrix_double3x3 lattice_vectors)
{
    Lattice *new = malloc(sizeof(Lattice));
    new->mol = mol;
    new->lattice_vectors = lattice_vectors;
    new->metal_site_id = 73;
    new->vtable = &lat_vtable;
    new->attached_adsName = NULL;
    new->pathName = NULL;
    return new;
}

void lattice_modify_metal_element(Lattice *self, const char *metal_symbol,
                                  int elementId)
{
    Atom *metalAtom =
        self->mol->vtable->get_atom_by_Id(self->mol, self->metal_site_id);
    strncpy(metalAtom->element, metal_symbol, strlen(metal_symbol));
    metalAtom->element[strlen(metal_symbol)] = '\0';
    metalAtom->elementId = elementId;
    int len = 1 + snprintf(NULL, 0, "SAC_GDY_%s", metal_symbol);
    snprintf(self->mol->name, len, "SAC_GDY_%s", metal_symbol);
    self->mol->name[len] = '\0';
}

void destroyLattice(Lattice *self)
{
    Molecule *mPtr = self->mol;
    mPtr->vtable->destroy(mPtr);
    free(self->attached_adsName);
    free(self->pathName);
    free(self);
}

/* Attach adsorbate to the lattice and update their atomId folloing the lattice
 * atoms */
Lattice *lattice_attach_molecule(Lattice *self, Adsorbate *ads, char *newName)
{
    Molecule *lat_mol = self->mol, *ads_mol = ads->mol;
    Atom **cur_arr = lat_mol->atom_arr;
    int new_atomNum = lat_mol->atomNum + ads_mol->atomNum;
    Atom *lastAtom = lat_mol->vtable->get_atom_by_Id(lat_mol, lat_mol->atomNum);
    int lastId = lastAtom->atomId;
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
    }
    Molecule *resMol = createMolecule(newName, new_atomNum, new_arr);
    Lattice *new = createLattice(resMol, self->lattice_vectors);
    new->attached_adsName = strdup(ads->mol->name);
    new->pathName = strdup(ads->pathName);
    return new;
}

void lattice_rotate_to_standard_orientation(Lattice *self)
{
    matrix_double3x3 lat_vectors = self->lattice_vectors;
    vec_double3 xAxis = vec_make_double3(1, 0, 0);
    double a_to_x = vec_angle_uv(lat_vectors.i, xAxis);
    if (a_to_x == 0)
    {
        return;
    }
    vec_double3 zAxis = vec_make_double3(0, 0, 1);
    vec_quatd rotateQ = vec_make_quaternion(a_to_x, zAxis);
    self->lattice_vectors = vec_act(rotateQ, lat_vectors);
    self->mol->vtable->rotateMol(self->mol, rotateQ);
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
    case 53:
        return siteDict[3].name;
    case 52:
        return siteDict[4].name;
    case 40:
        return siteDict[5].name;
    case 73:
        return siteDict[6].name;
    default:
        printf("Incorrect carbon site Id %d\n", siteId);
        return NULL;
    }
}

void lattice_export_MSI(Lattice *self)
{
    char header_line[] = "# MSI CERIUS2 DataModel File Version 4 0\n";
    char model_start[] = "(1 Model\n";
    char model_misc[] = "  (A I CRY/DISPLAY (192 256))\n  (A I PeriodicType "
                        "100)\n  (A C SpaceGroup \" 1 1\")\n  (A D A3 "
                        "(16.39518593025 -9.465765010246 0))\n  (A D B3 (0 "
                        "18.93153002049 0))\n  (A D C3 (0 0 9.999213039981))\n "
                        " (A D CRY/TOLERANCE 0.05)\n";
    char model_end[] = ")\n";
    Molecule *mol = self->mol;
    int lineSize = mol->atomNum + 4;
    char **atoms = mol->vtable->export_text(mol);
    char **contentLines = malloc(sizeof(char *) * lineSize);
    contentLines[0] = strdup(header_line);
    contentLines[1] = strdup(model_start);
    contentLines[2] = strdup(model_misc);
    for (int i = 0; i < self->mol->atomNum; ++i)
    {
        contentLines[i + 3] = atoms[i];
    }
    contentLines[lineSize - 1] = strdup(model_end);
    char *dest = lattice_export_dest(self);
    createDirectory(dest);
    int exportNameLen = snprintf(NULL, 0, "%s%s.msi", dest, mol->name) + 1;
    char *exportName = malloc(exportNameLen);
    snprintf(exportName, exportNameLen, "%s%s.msi", dest, mol->name);
    FILE *exportFile = fopen(exportName, "w");
    for (int i = 0; i < lineSize; ++i)
    {
        fputs(contentLines[i], exportFile);
        free(contentLines[i]);
    }
    fclose(exportFile);
    free(contentLines);
    free(exportName);
    free(atoms);
    free(dest);
}

static inline char *metal_family(int elementId)
{
    int metal_atomic_num = elementId;
    char *family = NULL;
    if (metal_atomic_num <= 30 && metal_atomic_num > 20)
        family = strdup("3d");
    else if (metal_atomic_num > 30 && metal_atomic_num <= 48)
        family = strdup("4d");
    else if (metal_atomic_num > 71 && metal_atomic_num <= 80)
        family = strdup("5d");
    else if (metal_atomic_num > 56 && metal_atomic_num <= 71)
        family = strdup("lm");
    else
    {
        printf("Not metal!\n");
    };
    return family;
}

char *lattice_export_dest(Lattice *self)
{
    Molecule *mPtr = self->mol;
    Atom *metal = mPtr->vtable->get_atom_by_Id(mPtr, self->metal_site_id);
    char *family = metal_family(metal->elementId);
    char *metal_symbol = metal->element;

    int destLen = 1 + snprintf(NULL, 0, "./C2_CO2RR_models/%s/%s/%s/%s/%s_opt/",
                               self->pathName, family, metal_symbol,
                               self->attached_adsName, self->mol->name);

    if (self->attached_adsName == NULL)
    {
        printf("No attached adsorbate!\n");
        return 0;
    }
    char *buffer = malloc(destLen);
    snprintf(buffer, destLen, "./C2_CO2RR_models/%s/%s/%s/%s/%s_opt/",
             self->pathName, family, metal_symbol, self->attached_adsName,
             self->mol->name);
    free(family);
    return buffer;
}
