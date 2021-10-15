#include "msiParser.h"
MOLECULE *init_mol(int atomNum)
{
    MOLECULE *mol;
    mol = malloc(sizeof(MOLECULE) + atomNum * sizeof(ATOM_BLOCK));
    mol->atomNum = atomNum;
    return mol;
}

int find_CdAtom(MOLECULE *mol)
{
    int i;
    for (i = 0; i < mol->atomNum && mol->molAtoms[i].bCdSite == 0; i++)
    {
        ;
    }
    return i;
}

int find_Stem(MOLECULE *mol)
{
    int i;
    for (i = 0; i < mol->atomNum && mol->molAtoms[i].bStem == 0; i++)
    {
        ;
    }
    return i;
}

/** parse .msi to MOLECULE pointer.
 *  Attention! malloc is used
 */
MOLECULE *parseMol(FILE *file)
{
    int atomNum = 0;
    atomNum = countAtoms(file);
    MOLECULE *mol = init_mol(atomNum);
    scanAtom(file, mol->molAtoms);
    mol->CdSiteId = find_CdAtom(mol);
    mol->StemId = find_Stem(mol);
    resetXYZ(atomNum, mol->molAtoms);
    return mol;
}
