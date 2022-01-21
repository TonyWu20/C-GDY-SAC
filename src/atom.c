#include "atom.h"
#include "matrix.h"
#include <stdlib.h>
#include <string.h>

// Implementation of Atom structure
struct _Atom
{
    char *element;
    char *ACL_Label;
    Matrix *coord;
    int atomId;
    int treeId;
};

Atom *createAtom(char *element, char *label, Matrix *coord, int atomId,
                 int treeId)
{
    Atom *newAtom = malloc(sizeof(Atom));
    newAtom->element = element;
    newAtom->ACL_Label = label;
    newAtom->coord = coord;
    newAtom->atomId = atomId;
    newAtom->treeId = treeId;
    return newAtom;
}
void destroyAtom(Atom *atomPtr)
{
    free(atomPtr->element);
    free(atomPtr->ACL_Label);
    destroy_matrix(atomPtr->coord);
    free(atomPtr->coord);
    free(atomPtr);
}

Matrix *Atom_get_coord(Atom *aPtr)
{
    return aPtr->coord;
}

void Atom_update_coord(Atom *aPtr, double x, double y, double z)
{
    aPtr->coord->value[0][0] = x;
    aPtr->coord->value[1][0] = y;
    aPtr->coord->value[2][0] = z;
}

int Atom_get_atomId(Atom *aPtr)
{
    return aPtr->atomId;
}

int Atom_get_treeId(Atom *aPtr)
{
    return aPtr->treeId;
}

// End of Atom
