#include "../maths/MyMaths.h"
#include "msiParser.h"

MOLECULE *init_mol(int atomNum)
{
    MOLECULE *mol;
    mol = malloc(sizeof(MOLECULE) + atomNum * sizeof(ATOM_BLOCK));
    mol->atomNum = atomNum;
    mol->StemAtomId[0] = 0;
    mol->StemAtomId[1] = 0;
    return mol;
}

int find_CdAtom(MOLECULE *mol)
{
    int i;
    for (i = 0; i < mol->atomNum && mol->totalAtoms[i].bCdSite == 0; i++)
    {
        ;
    }
    return i;
}

int find_Stem(MOLECULE *mol)
{
    int StemNum = 0;
    for (int i = 0; i < mol->atomNum; i++)
    {
        switch (mol->totalAtoms[i].StemInfo)
        {
        case 1:
            mol->StemAtomId[0] = i;
            StemNum++;
            break;
        case 2:
            mol->StemAtomId[1] = i;
            StemNum++;
            break;
        default:
            break;
        }
    }
    return StemNum;
}

/** parse .msi to MOLECULE pointer.
 *  Attention! malloc is used
 */
MOLECULE *parseMol(FILE *file)
{
    int atomNum = 0;
    int StemNum = 0;
    atomNum = countAtoms(file);
    MOLECULE *mol = init_mol(atomNum);
    scanAtom(file, mol->totalAtoms);
    mol->CdSiteId = find_CdAtom(mol);
    StemNum = find_Stem(mol);
    switch (StemNum)
    {
    case 0:
        printf("No atom marked as stem is found!\n");
    case 1:
        printf("Warning: only one atom is marked as stem!\n");
    case 2:
        break;
    default:
        printf("Error: more than two atoms are marked as stem atoms, please "
               "check model .msi file.\n");
    }
    resetXYZ(atomNum, mol->totalAtoms);
    return mol;
}

void load_StemVector(MOLECULE *mol, double *stemRes)
{
    int Stem_a, Stem_b;
    Stem_a = mol->StemAtomId[0];
    Stem_b = mol->StemAtomId[1];
    double *a, *b;
    a = mol->totalAtoms[Stem_a].coord;
    b = mol->totalAtoms[Stem_b].coord;
    initVector(a, b, stemRes);
}
