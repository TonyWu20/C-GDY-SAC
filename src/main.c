#include "assemble.h"
#include "atom.h"
#include "matrix.h"
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
    Lattice *result = Add_mol_to_lattice(lat, ads, c1, c2);
    result->vtable->export_msi(result, NULL);
    result->vtable->destroy(result);
    lat->vtable->destroy(lat);
    ads->ads_vtable->destroy(ads);
}

int main(int argc, char *argv[])
{
    test_add("SAC_GDY_V.msi", "SAC_GDY_V", "C2_pathways_ads/COCHO.msi", "COCHO",
             41, NULLSITE);
    test_add("SAC_GDY_V.msi", "SAC_GDY_V", "C2_pathways_ads/OCH2CH.msi",
             "OCH2CH", 41, 73);
    return 0;
}
