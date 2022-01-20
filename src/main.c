#include "atom.h"
#include "matrix.h"
#include "molecule.h"
#include "my_maths.h"
#include "parser.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define PI (atan(1) * 4)

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
    Molecule *mol =
        parse_molecule_from_file("C2_pathways_ads/C2H4.msi", "C2H4");
    Matrix *stem_vect = Molecule_get_stem_vector(mol);
    Matrix *a_vect = Molecule_get_vector_ab(mol, 1, 3);
    print_matrix(stem_vect);
    print_matrix(a_vect);
    printf("%f\n", vector_angle(stem_vect, a_vect) / PI * 180);
    cross_product(stem_vect, a_vect);
    free(content);
    destroy_matrix(stem_vect);
    destroy_matrix(a_vect);
    free(stem_vect);
    free(a_vect);
    destroyMolecule(mol);
    return 0;
}
