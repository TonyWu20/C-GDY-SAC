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
    FILE *f = fopen(fileName, "r");
    fseek(f, 0, SEEK_END);
    long fsize = ftell(f);
    rewind(f);
    char *body = malloc(fsize + 1);
    fread(body, fsize, 1, f);
    fclose(f);
    body[fsize] = 0;
    Lattice *lat = parse_lattice_from_file(fileName, name);
    for (int i = 0; i < 6; ++i)
    {
        printf("carbon site: %s, id: %d\n", lat->carbon_sites[i].name,
               lat->carbon_sites[i].id);
    }
    free(body);
    lat->vtable->destroy(lat);
    return lat;
}

int main(int argc, char *argv[])
{
    Matrix *rot_mat = rotationMatrix(PI, 'Z');
    print_matrix(rot_mat);
    destroy_matrix(rot_mat);
    free(rot_mat);
    FILE *file = fopen("C2_pathways_ads/C2H4.msi", "rb");
    fseek(file, 0, SEEK_END);
    long fsize = ftell(file);
    rewind(file);
    char *content = malloc(fsize + 1);
    fread(content, fsize, 1, file);
    fclose(file);
    content[fsize] = 0;
    printf("%s\n", content);
    Adsorbate *ads =
        parse_molecule_from_file("C2_pathways_ads/C2H4.msi", "C2H4");
    free(content);
    Lattice *GDY_V = load_lat("SAC_GDY_V.msi", "SAC_GDY_V");
    Add_mol_to_carbon_chain(GDY_V, ads, 41, 42);
    ads->ads_vtable->destroy(ads);
    GDY_V->vtable->destroy(GDY_V);
    return 0;
}
