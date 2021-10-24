#include "../maths/MyMaths.h"
#include "msiParser.h"
#include <math.h>
#include <stdio.h>

static void rotationMatrix(double degree, char axis, double (*mat)[3]);

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
    s = add_mol->totalAtoms;
    for (int i = 0; i < add_mol->atomNum; i++)
    {
        memcpy(t, s++, sizeof(ATOM_BLOCK));
        t++->itemId += source->atomNum;
    }
}

void get_CoordMat(MOLECULE *s, double matrix[][3])
{
    for (int i = 0; i < s->atomNum; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            memcpy(matrix[i], &s->totalAtoms[i].coord, 3 * sizeof(double));
        }
    }
}

void assignCoordtoMol(double coords[][3], MOLECULE *mol)
{
    for (int i = 0; i < mol->atomNum; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            memcpy(&mol->totalAtoms[i].coord, coords[i], 3 * sizeof(double));
        }
    }
}

static void rotationMatrix(double degree, char axis, double (*mat)[3])
{
    double theta;
    theta = (degree * M_PI) / 180;
    double cos_theta;
    double sin_theta;
    cos_theta = (((int)(degree * 1000000) % 90000000 == 0) &&
                 ((int)(degree * 1000000) % 180000000 != 0))
                    ? 0
                    : cos(theta);
    sin_theta = ((int)(degree * 1000000) % 180000000 == 0) ? 0 : sin(theta);
    double rot_x[3][3] = {{1.0, 0, 0},
                          {0, cos_theta, -1 * sin_theta},
                          {0.0, sin_theta, cos_theta}};
    double rot_y[3][3] = {
        {cos_theta, 0, sin_theta}, {0, 1, 0}, {-1 * sin_theta, 0, cos_theta}};
    double rot_z[3][3] = {{cos_theta, -1 * sin_theta, 0}, // line 1
                          {sin_theta, cos_theta, 0},      // line 2
                          {0, 0, 1}};
    switch (axis)
    {
    case 'x':
        for (int i = 0; i < 3; i++)
        {
            memcpy(mat[i], rot_x[i], sizeof(double) * 3);
        }
        break;
    case 'y':
        for (int i = 0; i < 3; i++)
        {
            memcpy(mat[i], rot_y[i], sizeof(double) * 3);
        }
        break;
    case 'z':
        for (int i = 0; i < 3; i++)
        {
            memcpy(mat[i], rot_z[i], sizeof(double) * 3);
        }
        break;
    default:
        perror("Wrong axis input.\n");
        break;
    }
}

void rotMol(MOLECULE *mol, double degree, char axis)
{
    double matrix[mol->atomNum][3];
    double rotated_coord[mol->atomNum][3];
    get_CoordMat(mol, matrix);
    double rot_mat[3][3];
    rotationMatrix(degree, axis, rot_mat);
    double *center = mol->totalAtoms[mol->StemAtomId[1]].coord;
    for (int i = 0; i < 3; ++i)
    {
        center[i] /= -2;
    }
    moveMatrix(matrix, center, mol->atomNum);
    for (int i = 0; i < 3; ++i)
    {
        center[i] *= -1;
    }
    multiplyMatrices_mx3(matrix, rot_mat, rotated_coord, mol->atomNum);
    moveMatrix(matrix, center, mol->atomNum);
    assignCoordtoMol(rotated_coord, mol);
}

void placeMol(MOLECULE *mol, BASE_LATTICE *lat, int destId,
              BASE_LATTICE *target)
{
    double *u = malloc(3 * sizeof(double));
    memcpy(u, lat->totalAtoms[destId].coord, 3 * sizeof(double));
    double coord[mol->atomNum][3];
    get_CoordMat(mol, coord);
    u[2] += 1.54221;
    moveMatrix(coord, u, mol->atomNum);
    assignCoordtoMol(coord, mol);
    appendMolAtoms(lat, mol, target);
    free(u);
}

void align_carbon_chain(MOLECULE *mol, double *chain_vec)
{
    double *mol_stem; // cc_vec is carbon chain vector
    double theta;
    mol_stem = malloc(3 * sizeof(double));
    load_StemVector(mol, mol_stem);
    theta = VecAngle(mol_stem, chain_vec);
    free(mol_stem);
    rotMol(mol, 180 - theta, 'z');
    /*resetXYZ(mol->atomNum, mol->totalAtoms);*/
}

void attach_carbon_chain(MOLECULE *mol, BASE_LATTICE *lat, int carbonId)
{
    double coord[mol->atomNum][3];
    align_carbon_chain(mol, lat->carbon_chain_vec);
    get_CoordMat(mol, coord);
    moveMatrix(coord, lat->totalAtoms[carbonId - 1].coord, mol->atomNum);
    assignCoordtoMol(coord, mol);
}
