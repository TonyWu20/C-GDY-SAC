#include "atom.h"
#include "matrix.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// Implementation of Atom methods

struct Atom_vtable atom_vtable = {
    Atom_get_coord,  Atom_update_coord, Atom_get_atomId,
    Atom_set_atomId, Atom_get_treeId,   Atom_set_treeId,
    dupAtom,         Atom_textblock,    destroyAtom};
Atom *createAtom(char *element, char *label, Matrix *coord, int atomId,
                 int treeId)
{
    Atom *newAtom = malloc(sizeof(Atom));
    newAtom->element = element;
    newAtom->ACL_Label = label;
    newAtom->coord = coord;
    newAtom->atomId = atomId;
    newAtom->treeId = treeId;
    newAtom->vtable = &atom_vtable;
    return newAtom;
}

Atom *dupAtom(Atom *self)
{
    Atom *dup = malloc(sizeof(Atom));
    dup->element = strdup(self->element);
    dup->ACL_Label = strdup(self->ACL_Label);
    dup->coord = create_matrix(4, 1);
    copy_matrix(self->coord, &dup->coord);
    dup->atomId = self->atomId;
    dup->treeId = self->treeId;
    dup->vtable = &atom_vtable;
    return dup;
}
void destroyAtom(Atom *atomPtr)
{
    free(atomPtr->element);
    free(atomPtr->ACL_Label);
    destroy_matrix(atomPtr->coord);
    free(atomPtr->coord);
    free(atomPtr);
}

Matrix *Atom_get_coord(Atom *self)
{
    return self->coord;
}

void Atom_update_coord(Atom *self, double x, double y, double z)
{
    self->coord->value[0][0] = x;
    self->coord->value[1][0] = y;
    self->coord->value[2][0] = z;
}

int Atom_get_atomId(Atom *self)
{
    return self->atomId;
}

void Atom_set_atomId(Atom *self, int newId)
{
    self->atomId = newId;
}

int Atom_get_treeId(Atom *self)
{
    return self->treeId;
}

void Atom_set_treeId(Atom *self, int newId)
{
    self->atomId = newId;
}

char *Atom_textblock(Atom *self)
{
    char *buffer = malloc(144);
    snprintf(buffer, 144,
             "  (%d Atom\n    (A C ACL \"%s\")\n    (A C Label \"%s\")\n    (A "
             "D XYZ (%.12f %.12f %.12f))\n    (A I Id %d)\n  )",
             self->treeId, self->ACL_Label, self->element,
             self->coord->value[0][0], self->coord->value[1][0],
             self->coord->value[2][0], self->atomId);
    buffer = realloc(buffer, strlen(buffer) + 1);
    return buffer;
}
// End of Atom
