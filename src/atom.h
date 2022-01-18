#pragma once
#include <matrix.h>
struct _Atom;
typedef struct _Atom Atom;

// Memory management
Atom *createAtom(char *element, char *label, Matrix *coord, int atomId,
                 int treeId);

void destroyAtom(Atom *);

// Methods
Matrix *Atom_get_coord(Atom *);
void Atom_update_coord(Atom *, Matrix *);
int Atom_get_atomId(Atom *);
int Atom_get_treeId(Atom *);
