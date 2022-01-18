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
    FILE *file = fopen("COCHOH.msi", "rb");
    fseek(file, 0, SEEK_END);
    long fsize = ftell(file);
    rewind(file);
    char *content = malloc(fsize + 1);
    fread(content, fsize, 1, file);
    fclose(file);
    content[fsize] = 0;
    printf("%s\n", content);
    int atom_block_size = 0;
    Atom **atom_blocks = atom_block(content, &atom_block_size);
    for (int i = 0; i < atom_block_size; ++i)
    {
        print_matrix(Atom_get_coord(atom_blocks[i]));
    }
    for (int i = 0; i < atom_block_size; ++i)
    {
        destroyAtom(atom_blocks[i]);
    }
    free(content);
    free(atom_blocks);
    return 0;
}
