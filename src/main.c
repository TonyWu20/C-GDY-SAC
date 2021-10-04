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
    printf("%d Atom: %s, XYZ: %lf, %lf, %lf\n", ads[3].itemId, ads[3].elm,
           ads[3].x, ads[3].y, ads[3].z);
}
