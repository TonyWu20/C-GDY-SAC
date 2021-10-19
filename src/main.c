#include "main.h"
#include "maths/MyMaths.h"
#include "msiParser/msiParser.h"
#include "msiWriter/msiWriter.h"
#include <stdio.h>

void init_mol_direction(FILE *file);
void test_latMsi_writer(FILE *file);
int main(int argc, char *argv[])
{
    char *fileName = argv[1];
    FILE *file = fopen(fileName, "r");
    if (file == NULL)
    {
        printf("Cannot open file.\n");
        return 1;
    }
    init_mol_direction(file);
}

void init_mol_direction(FILE *file)
{
    MOLECULE *mol;
    mol = parseMol(file);
    fclose(file);
    double *a;
    a = malloc(3 * sizeof(double));
    load_StemVector(mol, a);
    double a_norm = NormVector(a, 3);
    double b[] = {a_norm, 0, 0};
    double theta = VecAngle(a, b);
    free(a);
    rotMol(mol, -theta, 'z');
    /*rotMol(mol, 90, 'x');*/
    MSI_FILE *rotated_mol;
    rotated_mol = build_MolMsi(mol);
    free(mol);
    for (int i = 0; i < rotated_mol->ItemNum; i++)
    {
        printf("%s", rotated_mol->lines[i]);
    }
    free_MSI_MOL(rotated_mol);
}

void test_latMsi_writer(FILE *file)
{
    BASE_LATTICE *lat;
    lat = parseBase(file);
    fclose(file);
    MSI_FILE *lat_file;
    lat_file = build_LatMsi(lat);
    free(lat);
    free_MSI_LAT(lat_file);
}
