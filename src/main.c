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
    ATOM_BLOCK atoms[10];
    int atomCount;
    atomCount = scanAtom(file, atoms);
    printf("Parse molecule file %s successfully!\n", fileName);
    resetXYZ(atomCount, atoms);
    printf("Reset coordinates\n");
    for (int i = 0; i < atomCount; i++)
    {
        printf("Atom %d %s XYZ: %lf, %lf, %lf\n", atoms[i].itemId, atoms[i].elm,
               atoms[i].coord[0], atoms[i].coord[1], atoms[i].coord[2]);
    }
}
