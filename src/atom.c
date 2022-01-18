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
    newAtom->element = strdup(element);
    newAtom->ACL_Label = strdup(label);
    newAtom->coord = create_matrix(4, 4);
    copy_matrix(coord, &newAtom->coord);
    newAtom->atomId = atomId;
    newAtom->treeId = treeId;
    return newAtom;
}
void destroyAtom(Atom *atomPtr)
{
    free(atomPtr->element);
    free(atomPtr->ACL_Label);
    destroy_matrix(atomPtr->coord);
    free(atomPtr);
}

Matrix *Atom_get_coord(Atom *aPtr)
{
    Matrix *coord = aPtr->coord;
    return coord;
}

void Atom_update_coord(Atom *aPtr, Matrix *newCoord)
{
    copy_matrix(newCoord, &aPtr->coord);
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
