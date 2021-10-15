#include "main.h"
#include "maths/MyMaths.h"
#include "msiParser/msiParser.h"
#include "msiWriter/msiWriter.h"
#include <stdio.h>

void init_mol_direction(FILE *file);
int main(int argc, char *argv[])
{
    char *fileName = argv[1];
    FILE *file = fopen(fileName, "r");
    if (file == NULL)
    {
        printf("Cannot open file.\n");
        return 1;
    }
    BASE_LATTICE *lat = parseBase(file);
    free(lat);
    fclose(file);
}

void init_mol_direction(FILE *file)
{
    MOLECULE *mol;
    mol = parseMol(file);
    fclose(file);
    double *u, *v, *a;
    u = mol->molAtoms[0].coord;
    v = mol->molAtoms[1].coord;
    a = malloc(sizeof(double) * 3);
    initVector(u, v, a);
    double a_norm = NormVector(a, 3);
    double b[] = {a_norm, 0, 0};
    double theta = VecAngle(a, b);
    free(a);
    rotMol(mol, -theta, 'z');
    rotMol(mol, -90, 'x');
    MSI_FILE *rotated_mol;
    rotated_mol = build_MolMsi(mol);
    free(mol);
    for (int i = 0; i < rotated_mol->ItemNum; i++)
    {
        fprintf(stdout, "%s", rotated_mol->lines[i]);
    }
    free(rotated_mol);
}
