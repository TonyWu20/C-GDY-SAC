#include "msiParser.h"

/** Append passed ATOM_BLOCK to the source BASE_LATTICE, make changes to the
 * passed target BASE_LATTICE pointer
 * The *target should be initialized before passed into the function.
 * Args:
 * 		BASE_LATTICE source; source lattice
 * 		MOLECULE add_mol; appended molecule
 * 		BASE_LATTICE *target; target lattice pointer
 */
BASE_LATTICE *init_adsorbed_lat(BASE_LATTICE *source, MOLECULE *mol)
{
    int src_atomNum = source->atomNum;
    int added_atomNum = mol->atomNum;
    BASE_LATTICE *target;
    target = init_lattice(src_atomNum + added_atomNum);
    target->atomNum = src_atomNum + added_atomNum;
    /* copy latVector */
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            target->latVector[i][j] = source->latVector[i][j];
        }
    }
    return target;
}
void appendMolAtoms(BASE_LATTICE *source, MOLECULE *add_mol,
                    BASE_LATTICE *target)
{
    ATOM_BLOCK *s, *t;
    s = source->totalAtoms;
    t = target->totalAtoms;
    for (int i = 0; i < source->atomNum; i++)
    {
        memcpy(t++, s++, sizeof(ATOM_BLOCK));
    }
    s = add_mol->molAtoms;
    for (int i = 0; i < add_mol->atomNum; i++)
    {
        memcpy(t, s++, sizeof(ATOM_BLOCK));
        t++->itemId += source->atomNum;
    }
}
