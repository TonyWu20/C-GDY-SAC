#include "main.h"
#include "msiParser/msiParser.h"
#include <stdio.h>

ATOM_BLOCK ads[10];

int main(int argc, char *argv[])
{
    char *fileName = argv[1];
    FILE *file = fopen(fileName, "r");
    if (file == NULL)
    {
        printf("Cannot open file.\n");
        return 1;
    }
    int atomCount;
    atomCount = countAtoms(file);
    MOLECULE *mol = init_mol(atomCount);
    scanAtom(file, mol->molAtoms);
    printf("Parse molecule file %s successfully! Has %d atoms.\n", fileName,
           mol->atomNum);
    resetXYZ(atomCount, mol->molAtoms);
    printf("Reset coordinates\n");
    for (int i = 0; i < atomCount; i++)
    {
        printf("Atom %d %s XYZ: %lf, %lf, %lf\n", mol->molAtoms[i].itemId,
               mol->molAtoms[i].elm, mol->molAtoms[i].coord[0],
               mol->molAtoms[i].coord[1], mol->molAtoms[i].coord[2]);
    }
    fclose(file);
    BASE_LATTICE model;
    fileName = "SAC_GDY_V.msi";
    file = fopen(fileName, "r");
    model = parseBase(file);
    for (int i = 0; i < 3; i++)
    {
        printf("Vector: %lf, %lf, %lf\n", model.latVector[i][0],
               model.latVector[i][1], model.latVector[i][2]);
    }
    fclose(file);
}
