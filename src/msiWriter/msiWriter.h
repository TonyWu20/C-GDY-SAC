#include "../main.h"

typedef struct
{
    int ItemNum;
    char *lines[];
} MSI_FILE;

MSI_FILE *init_MSI_FILE(int ItemNum);
MSI_FILE *build_MolMsi(MOLECULE *mol);
MSI_FILE *build_LatMsi(BASE_LATTICE *);

void write_atomBlock(ATOM_BLOCK, char **);
char *write_vector(double *, char VectorName);
void load_atoms(ATOM_BLOCK atoms[], int atomNum, MSI_FILE *, int startPos);
void load_headers(MSI_FILE *);
void load_vectors(MSI_FILE *, BASE_LATTICE *);

void free_MSI_MOL(MSI_FILE *);
void free_MSI_LAT(MSI_FILE *);
