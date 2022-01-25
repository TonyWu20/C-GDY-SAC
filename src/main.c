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

int main(int argc, char *argv[])
{
    allocateTasks("ethylene");
    return 0;
}
