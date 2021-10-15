#include "msiWriter.h"

MSI_FILE *init_MSI_FILE(int ItemNum)
{
    MSI_FILE *model;
    model = malloc(sizeof(MSI_FILE) + ItemNum * 2 * sizeof(char));
    model->ItemNum = ItemNum;
    return model;
}

MSI_FILE *build_MolMsi(MOLECULE *mol)
{
    int ItemNum = mol->atomNum + 3;
    MSI_FILE *mol_file = init_MSI_FILE(ItemNum);
    mol_file->lines[0] = strdup("# MSI CERIUS2 DataModel File Version 4.0\n");
    char *begin = "(1 Model\n";
    mol_file->lines[1] = strdup(begin);
    for (int i = 0; i < mol->atomNum; i++)
    {
        mol_file->lines[i + 2] = strdup(write_atomBlock(mol->molAtoms[i]));
    }
    mol_file->lines[mol->atomNum + 2] = strdup(")");
    return mol_file;
}

char *write_atomBlock(ATOM_BLOCK atom)
{
    char *buffer;
    asprintf(&buffer,
             "  (%d Atom\n    (A C ACL \"%d %s\")\n    (A C Label "
             "\"%s\")\n    (A D XYZ (%f %f %f))\n    (A I Id %d)\n  )\n",
             atom.itemId, atom.elmId, atom.elm, atom.elm, atom.coord[0],
             atom.coord[1], atom.coord[2], atom.itemId - 1);
    return buffer;
}
