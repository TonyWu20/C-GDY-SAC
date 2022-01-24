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

int main(int argc, char *argv[])
{
    FILE *file = fopen("C2_pathways_ads/C2H4.msi", "rb");
    fseek(file, 0, SEEK_END);
    long fsize = ftell(file);
    rewind(file);
    char *content = malloc(fsize + 1);
    fread(content, fsize, 1, file);
    fclose(file);
    content[fsize] = 0;
    Adsorbate *ads =
        parse_molecule_from_file("C2_pathways_ads/C2H4.msi", "C2H4");
    free(content);
    Lattice *GDY_V = load_lat("SAC_GDY_V.msi", "SAC_GDY_V");
    Lattice *result = Add_cd2_mol_to_carbon_chain(GDY_V, ads, 42, 54);
    ads->ads_vtable->export_msi(ads, NULL);
    GDY_V->vtable->destroy(GDY_V);
    ads->ads_vtable->destroy(ads);
    printf("%s\n", result->_mol->name);
    result->vtable->export_msi(result, NULL);
    result->vtable->destroy(result);
    return 0;
}
