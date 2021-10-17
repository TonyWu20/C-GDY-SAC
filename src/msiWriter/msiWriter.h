#include "../main.h"

typedef struct
{
    int ItemNum;
    char *lines[];
} MSI_FILE;

MSI_FILE *init_MSI_FILE(int ItemNum);
MSI_FILE *build_MolMsi(MOLECULE *mol);

void write_atomBlock(ATOM_BLOCK, char **);
void free_MSI_FILE(MSI_FILE *);
