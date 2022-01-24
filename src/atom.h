#pragma once
#include <matrix.h>
struct _Atom
{
    char *element;
    char *ACL_Label;
    Matrix *coord;
    int atomId;
    int treeId;
    struct Atom_vtable *vtable;
};
typedef struct _Atom Atom;
struct Atom_vtable
{
    Matrix *(*get_coord)(Atom *);
    void (*update_coord)(Atom *, double, double, double);
    int (*get_atomId)(Atom *);
    void (*set_atomId)(Atom *, int);
    int (*get_treeId)(Atom *);
    void (*set_treeId)(Atom *, int);
    Atom *(*dupAtom)(Atom *self);
    char *(*export_text)(Atom *self);
    void (*destroy)(Atom *);
};

// Memory management
Atom *createAtom(char *element, char *label, Matrix *coord, int atomId,
                 int treeId);

Atom *dupAtom(Atom *self);
void destroyAtom(Atom *);

// Methods
Matrix *Atom_get_coord(Atom *);
void Atom_update_coord(Atom *, double x, double y, double z);
int Atom_get_atomId(Atom *);
void Atom_set_atomId(Atom *, int);
int Atom_get_treeId(Atom *);
void Atom_set_treeId(Atom *, int);
char *Atom_textblock(Atom *self);
