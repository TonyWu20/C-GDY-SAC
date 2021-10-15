#include "msiWriter.h"

MSI_FILE *init_MSI_FILE(int ItemNum)
{
    MSI_FILE *model;
    model = malloc(sizeof(MSI_FILE) + ItemNum * 80 * sizeof(char));
    model->ItemNum = ItemNum;
    return model;
}

MSI_FILE *build_MolMsi(MOLECULE *mol)
{
    int ItemNum = mol->atomNum + 3;
    MSI_FILE *mol_file = init_MSI_FILE(ItemNum);
    mol_file->lines[0] = "# MSI CERIUS2 DataModel File Version 4.0\n";
    char *begin = strdup("(1 Model\n");
    mol_file->lines[1] = begin;
    for (int i = 0; i < mol->atomNum; i++)
    {
        asprintf(&mol_file->lines[i + 2],
                 "  (%d Atom\n    (A C ACL \"%d %s\")\n    (A C Label "
                 "\"%s\")\n    (A D XYZ (%f %f %f))\n    (A I Id %d)\n  )\n",
                 mol->molAtoms[i].itemId, mol->molAtoms[i].elmId,
                 mol->molAtoms[i].elm, mol->molAtoms[i].elm,
                 mol->molAtoms[i].coord[0], mol->molAtoms[i].coord[1],
                 mol->molAtoms[i].coord[2], mol->molAtoms[i].itemId - 1);
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
