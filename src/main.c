#include "atom.h"
#include "matrix.h"
#include "molecule.h"
#include "my_maths.h"
#include "parser.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define PI (atan(1) * 4)
int test_lat(char *fileName)
{
    FILE *f = fopen(fileName, "r");
    fseek(f, 0, SEEK_END);
    long fsize = ftell(f);
    rewind(f);
    char *body = malloc(fsize + 1);
    fread(body, fsize, 1, f);
    fclose(f);
    body[fsize] = 0;
    Matrix *res = NULL;
    get_lattice_vectors(body, &res);
    free(body);
    destroy_matrix(res);
    free(res);
    return 0;
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
    Molecule *mol = ads->_mol;
    Matrix *coord = mol->vtable->get_mol_coords(mol);
    print_matrix(coord);
    destroy_matrix(coord);
    free(coord);
    ads->ads_vtable->make_upright(ads);
    Matrix *after = mol->vtable->get_mol_coords(mol);
    print_matrix(after);
    destroy_matrix(after);
    free(after);
    free(content);
    ads->ads_vtable->destroy(ads);
    test_lat("SAC_GDY_V.msi");
    return 0;
}
