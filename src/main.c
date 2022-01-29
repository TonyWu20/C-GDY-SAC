#include "assemble.h"
#include "atom.h"
#include "castep_output.h"
#include "matrix.h"
#include "misc.h"
#include "molecule.h"
#include "my_maths.h"
#include "parser.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define PI (atan(1) * 4)
Lattice *load_lat(char *fileName, char *name)
{
    Lattice *lat = parse_lattice_from_file(fileName, name);
    return lat;
}

void test_add(char *latFile, char *latName, char *adsFile, char *adsName,
              int c1, int c2)
{
    Lattice *lat = load_lat(latFile, latName);
    Adsorbate *ads = parse_molecule_from_file(adsFile, adsName);
    Lattice *result = Add_mol_to_lattice(lat, ads, ads->taskLists->tasks[0][0],
                                         ads->taskLists->tasks[0][1]);
    printf("%s\n", result->_mol->name);
    result->vtable->export_msi(result, "ethylene");
    result->vtable->destroy(result);
    lat->vtable->destroy(lat);
    ads->ads_vtable->destroy(ads);
}

void test_build()
{
    allocateTasks("ethylene");
}

void test_fracCoordMat()
{
    Lattice *t = load_lat("./SAC_GDY_V.msi", "SAC_GDY_V");
    Adsorbate *ads = parse_molecule_from_file("./C2_pathways_ads/ethylene_path/C2H4.msi", "C2H4");
    Lattice *result = Add_mol_to_lattice(t, ads, ads->taskLists->tasks[0][0], ads->taskLists->tasks[0][1]);
    Matrix *fracCoordMat = fractionalCoordMatrix(result->lattice_vectors);
    print_matrix(fracCoordMat);
    Matrix *t_coord = result->_mol->vtable->get_mol_coords(result->_mol);
    Matrix *res;
    multiply_matrices(fracCoordMat, t_coord, &res);
    for (int i = 0; i < result->_mol->atomNum; ++i)
    {
        Atom *cur = result->_mol->atom_arr[i];
        printf("%d %s: (%.15f, %.15f, %.15f)\n", cur->atomId, cur->element,
               res->value[0][i], res->value[1][i], res->value[2][i]);
    }
    destroy_matrix(t_coord);
    destroy_matrix(res);
    destroy_matrix(fracCoordMat);
    free(fracCoordMat);
    free(t_coord);
    free(res);
    t->vtable->destroy(t);
    ads->ads_vtable->destroy(ads);
    result->vtable->destroy(result);
}

void test_table()
{
    CastepInfo *table = initTable();
    CastepInfo *item;
    item = find_item(table, "Sc");
    printf("%s:\n\tLCAO:%d\n\tMass:%.12f\n\tPot:%s\n\tSpin:%d\n",
           item->info->elm, item->info->LCAO, item->info->mass,
           item->info->potential_file, item->info->spin);
    int needed = snprintf(NULL, 0, "head %s\n", item->info->potential_file);
    char buf[needed+1];
    snprintf(buf, needed+1, "head %s\n", item->info->potential_file);
    system(buf);
    delete_all(&table);
}

void test_cell()
{
    Lattice *t = load_lat("./SAC_GDY_Sc.msi", "SAC_GDY_Sc");
    Adsorbate *ads = parse_molecule_from_file("./C2_pathways_ads/ethylene_path/C2H4.msi", "C2H4");
    Lattice *result = Add_mol_to_lattice(t, ads, ads->taskLists->tasks[0][0], ads->taskLists->tasks[0][1]);
    t->vtable->destroy(t);
    ads->ads_vtable->destroy(ads);
    CastepInfo *table = initTable();
    Cell *cell = createCell(result, table);
    cell->vtable->sortAtoms(cell);
    for (int i = 0; i < cell->lattice->_mol->atomNum; ++i)
    {
        Atom *cur = cell->lattice->_mol->atom_arr[i];
        printf("%d %s %d: (%.15f, %.15f, %.15f)\n", cur->atomId, cur->element, cur->elementId,
               cur->coord->value[0][0], cur->coord->value[1][0], cur->coord->value[2][0]);
    }
    int elmNum = 0;
    char **elmList = cell->vtable->sortElmList(cell, &elmNum);
    for (int i = 0; i < elmNum; ++i)
    {
        CastepInfo *item = find_item(table, elmList[i]);
        printf("%8s%18.10f\n", elmList[i], item->info->mass);
    }
    for (int i = 0; i < elmNum; ++i)
    {
        CastepInfo *item = find_item(table, elmList[i]);
        char *pos = strrchr(item->info->potential_file, '/');
        printf("%8s  %s\n", elmList[i], pos+1);
        free(elmList[i]);
    }
    char *vec = cell->textTable->blockWriter(cell, "LATTICE_CART", cell_latticeVector_writer);
    printf("%s\n", vec);
    free(vec);
    free(elmList);
    cell->destroy(cell);
    delete_all(&table);
}

int main(int argc, char *argv[])
{
    test_cell();
    return 0;
}
