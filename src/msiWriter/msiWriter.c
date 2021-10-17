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
    mol_file->lines[0] = "# MSI CERIUS2 DataModel File Version 4.0\n";
    char *begin = "(1 Model\n";
    mol_file->lines[1] = begin;
    load_atoms(mol->totalAtoms, mol->atomNum, mol_file, 2);
    mol_file->lines[mol->atomNum + 2] = ")";
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

void load_headers(MSI_FILE *lat_file)
{
    ;
}
void load_vectors(MSI_FILE *lat_file, BASE_LATTICE *lat)
{
    ;
}

void free_MSI_FILE(MSI_FILE *mfile)
{
    for (int i = 0; i < mfile->ItemNum - 3; ++i)
    {
        free(mfile->lines[i + 2]);
    }
    free(mfile);
}
