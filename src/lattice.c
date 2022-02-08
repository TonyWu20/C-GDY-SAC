#include "lattice.h"
#include "misc.h"
#include "molecule.h"
#include "my_maths.h"
#include <matrix.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>

struct carbon_site siteDict[] = {
    {"c1", 41},       {"c2", 42},        {"c3", 54},   {"c4", 43},
    {"far_ring", 52}, {"near_ring", 40}, {"metal", 73}};
struct Lattice_vtable lat_vtable = {lattice_get_carbon_chain_vector,
                                    lattice_get_carbon_metal_vector,
                                    lattice_attach_molecule,
                                    lattice_rotate_to_standard_orientation,
                                    lattice_export_dest,
                                    lattice_export_MSI,
                                    destroyLattice};
Lattice *createLattice(Molecule *mol, Matrix *lattice_vectors)
{
    Lattice *new = malloc(sizeof(Lattice));
    new->_mol = mol;
    new->lattice_vectors = lattice_vectors;
    memcpy(new->carbon_sites, siteDict, sizeof(struct carbon_site) * 7);
    new->metal_site_id = 73;
    new->vtable = &lat_vtable;
    new->attached_adsName = NULL;
    new->pathName = NULL;
    lattice_metal_info(new);
    return new;
}

void destroyLattice(Lattice *self)
{
    Molecule *mPtr = self->_mol;
    mPtr->vtable->destroy(mPtr);
    destroy_matrix(self->lattice_vectors);
    free(self->lattice_vectors);
    free(self->metal_family);
    free(self->attached_adsName);
    free(self->pathName);
    free(self);
}

Matrix *lattice_get_carbon_chain_vector(Lattice *self)
{
    return self->_mol->vtable->get_vector_ab(self->_mol, 41, 42);
}

/* Return a vector pointing to metal from given carbon atom */
Matrix *lattice_get_carbon_metal_vector(Lattice *self, int cId)
{
    return self->_mol->vtable->get_vector_ab(self->_mol, cId, 73);
}

/* Attach adsorbate to the lattice and update their atomId folloing the lattice
 * atoms */
Lattice *lattice_attach_molecule(Lattice *self, Adsorbate *ads, char *newName,
                                 char *pathName)
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
    Matrix *lattice_vectors = create_matrix(4, 3);
    copy_matrix(self->lattice_vectors, &lattice_vectors);
    Lattice *new = createLattice(resMol, lattice_vectors);
    new->attached_adsName = strdup(ads->_mol->name);
    new->pathName = strdup(pathName);
    return new;
}

void lattice_rotate_to_standard_orientation(Lattice *self)
{
    Molecule *mol = self->_mol;
    Matrix *lat_vectors = self->lattice_vectors;
    Matrix *a = create_matrix(4, 1);
    Matrix *b = create_matrix(4, 1);
    for (int i = 0; i < 3; ++i)
    {
        a->value[i][0] = lat_vectors->value[i][0];
        b->value[i + 1][0] = 0;
    }
    a->value[3][0] = 1;
    b->value[0][0] = 1;
    b->value[3][0] = 1;
    double a_to_x = vector_angle(a, b);
    if (a_to_x == 0)
    {
        destroy_matrix(a);
        destroy_matrix(b);
        free(a);
        free(b);
        return;
    }
    Matrix *rot_mat = rotationMatrix(a_to_x, 'Z');
    Matrix *new_lat_vec;
    multiply_matrices(rot_mat, lat_vectors, &new_lat_vec);
    self->lattice_vectors = new_lat_vec;
    mol->vtable->apply_transformation(mol, rot_mat, rotate_around_origin);
    destroy_matrix(lat_vectors);
    destroy_matrix(a);
    destroy_matrix(b);
    destroy_matrix(rot_mat);
    free(lat_vectors);
    free(a);
    free(b);
    free(rot_mat);
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

void lattice_export_MSI(Lattice *self, char *pathName)
{
    char header_line[] = "# MSI CERIUS2 DataModel File Version 4 0\n";
    char model_start[] = "(1 Model\n";
    char model_misc[] = "  (A I CRY/DISPLAY (192 256))\n  (A I PeriodicType "
                        "100)\n  (A C SpaceGroup \" 1 1\")\n  (A D A3 "
                        "(16.39518593025 -9.465765010246 0))\n  (A D B3 (0 "
                        "18.93153002049 0))\n  (A D C3 (0 0 9.999213039981))\n "
                        " (A D CRY/TOLERANCE 0.05)\n";
    char model_end[] = ")\n";
    Molecule *mol = self->_mol;
    int lineSize = mol->atomNum + 4;
    char **atoms = mol->vtable->export_text(mol);
    char **contentLines = malloc(sizeof(char *) * lineSize);
    contentLines[0] = strdup(header_line);
    contentLines[1] = strdup(model_start);
    contentLines[2] = strdup(model_misc);
    for (int i = 0; i < self->_mol->atomNum; ++i)
    {
        contentLines[i + 3] = atoms[i];
    }
    contentLines[lineSize - 1] = strdup(model_end);
    char *dest = lattice_export_dest(self, pathName);
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

void lattice_metal_info(Lattice *self)
{
    Atom *metal =
        self->_mol->vtable->get_atom_by_Id(self->_mol, self->metal_site_id);
    self->metal_symbol = metal->element;
    char *token = NULL;
    char *buffer = strdup(metal->ACL_Label);
    char *rest;
    token = strtok_r(buffer, " ", &rest);
    int metal_atomic_num = atoi(token);
    if (metal_atomic_num <= 30)
        self->metal_family = strdup("3d");
    else if (metal_atomic_num > 30 && metal_atomic_num <= 48)
        self->metal_family = strdup("4d");
    else if (metal_atomic_num > 71 && metal_atomic_num <= 80)
        self->metal_family = strdup("5d");
    else if (metal_atomic_num > 56 && metal_atomic_num <= 71)
        self->metal_family = strdup("lm");
    self->metal_order = metal_atomic_num;
    free(buffer);
}

char *lattice_export_dest(Lattice *self, char *pathName)
{
    char *metal_family = self->metal_family;
    char *metal_symbol = self->metal_symbol;

    int destLen = 1 + snprintf(NULL, 0, "./C2_CO2RR_models/%s/%s/%s/%s/%s_opt/",
                               pathName, metal_family, metal_symbol,
                               self->attached_adsName, self->_mol->name);

    if (self->attached_adsName == NULL)
    {
        printf("No attached adsorbate!\n");
        return 0;
    }
    char *buffer = malloc(destLen);
    snprintf(buffer, destLen, "./C2_CO2RR_models/%s/%s/%s/%s/%s_opt/", pathName,
             metal_family, metal_symbol, self->attached_adsName,
             self->_mol->name);
    return buffer;
}
