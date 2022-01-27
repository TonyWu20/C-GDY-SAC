#include "assemble.h"
#include "atom.h"
#include "matrix.h"
#include "misc.h"
#include "molecule.h"
#include "my_maths.h"
#include "parser.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

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
    Matrix *fracCoordMat = fractionalCoordMatrix(t->lattice_vectors);
    print_matrix(fracCoordMat);
    Matrix *t_coord = t->_mol->vtable->get_mol_coords(t->_mol);
    Matrix *res;
    multiply_matrices(fracCoordMat, t_coord, &res);
    for (int i = 0; i < t->_mol->atomNum; ++i)
    {
        Atom *cur = t->_mol->atom_arr[i];
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
}

int main(int argc, char *argv[])
{
    test_fracCoordMat();
    return 0;
}
