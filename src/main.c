#include "main.h"
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
    scanAtom(file);
    printf("Parse molecule file %s successfully!\n", fileName);
    resetXYZ();
    printf("Reset coordinates\n");
    for (int i = 0; i < atomCount; i++)
    {
        printf("Atom %d %s XYZ: %lf, %lf, %lf\n", ads[i].itemId, ads[i].elm,
               ads[i].x, ads[i].y, ads[i].z);
    }
}
