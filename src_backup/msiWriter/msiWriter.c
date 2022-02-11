#include "msiWriter.h"
#include <stdio.h>

MSI_FILE *init_MSI_FILE(int ItemNum)
{
    MSI_FILE *model;
    model = malloc(sizeof(MSI_FILE) + ItemNum * 150 * sizeof(char));
    model->ItemNum = ItemNum;
    return model;
}

MSI_FILE *build_MolMsi(MOLECULE *mol)
{
    int ItemNum = mol->atomNum + 3;
    MSI_FILE *mol_file = init_MSI_FILE(ItemNum);
    load_headers(mol_file);
    load_atoms(mol->totalAtoms, mol->atomNum, mol_file, 2);
    return mol_file;
}

MSI_FILE *build_LatMsi(BASE_LATTICE *lat)
{
    int ItemNum =
        lat->atomNum + 1 + 2 + 7; /* declaration line; (1 Model + ); crystal
                                     settings and lattice parameters */
    MSI_FILE *lat_file = init_MSI_FILE(ItemNum);
    load_headers(lat_file);
    load_vectors(lat_file, lat);
    load_atoms(lat->totalAtoms, lat->atomNum, lat_file, 9);
    return lat_file;
}

void write_atomBlock(ATOM_BLOCK atom, char **line)
{
    asprintf(line,
             "  (%d Atom\n    (A C ACL \"%d %s\")\n    (A C Label "
             "\"%s\")\n    (A D XYZ (%f %f %f))\n    (A I Id %d)\n  )\n",
             atom.itemId, atom.elmId, atom.elm, atom.elm, atom.coord[0],
             atom.coord[1], atom.coord[2], atom.itemId - 1);
}

void load_atoms(ATOM_BLOCK *atoms, int atomNum, MSI_FILE *mfile, int startPos)
{
    for (int i = 0; i < atomNum; ++i)
    {
        write_atomBlock(atoms[i], &mfile->lines[i + startPos]);
    }
}

void load_headers(MSI_FILE *mfile)
{
    mfile->lines[0] = "# MSI CERIUS2 DataModel File Version 4.0\n";
    mfile->lines[1] = "(1 Model\n";
    mfile->lines[mfile->ItemNum - 1] = ")";
}

char *write_vector(double *coord, char VectorName)
{
    char *vector;
    asprintf(&vector, "  (A D %c3 (%13.11lf %13.11lf %13.11lf))\n", VectorName,
             coord[0], coord[1], coord[2]);
    return vector;
}
void load_vectors(MSI_FILE *lat_file, BASE_LATTICE *lat)
{
    char *Avector, *Bvector, *Cvector;
    Avector = write_vector(lat->latVector[0], 'A');
    Bvector = write_vector(lat->latVector[1], 'B');
    Cvector = write_vector(lat->latVector[2], 'C');
    char *parameters[] = {"  (A I CRY/DISPLAY (192 256))\n",
                          "  (A I PeriodicType 100)\n",
                          "  (A C SpaceGroup \"1 1\")\n",
                          Avector,
                          Bvector,
                          Cvector,
                          "  (A D CRY/TOLERANCE 0.05)\n"};
    ;
    for (int i = 0; i < 7; ++i)
    {
        lat_file->lines[i + 2] = parameters[i];
    }
}

void free_MSI_MOL(MSI_FILE *mfile)
{
    for (int i = 0; i < mfile->ItemNum - 3; ++i)
    {
        free(mfile->lines[i + 2]);
    }
    free(mfile);
}
void free_MSI_LAT(MSI_FILE *mfile)
{
    for (int i = 0; i < 3; ++i)
    {
        free(mfile->lines[i + 5]);
    }
    for (int i = 0; i < mfile->ItemNum - 10; ++i)
    {
        free(mfile->lines[i + 9]);
    }
    free(mfile);
}
