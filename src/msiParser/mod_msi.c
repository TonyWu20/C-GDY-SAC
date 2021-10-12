#include "../maths/MyMaths.h"
#include "msiParser.h"
#include <math.h>
#include <stdio.h>
#define PI 3.141592654

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

void get_CoordMat(MOLECULE *s, double matrix[][3])
{
    for (int i = 0; i < s->atomNum; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            memcpy(matrix[i], &s->molAtoms[i].coord, 3 * sizeof(double));
        }
    }
}

void assignCoordtoMol(double coords[][3], MOLECULE *mol)
{
    for (int i = 0; i < mol->atomNum; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            memcpy(&mol->molAtoms[i].coord, coords[i], 3 * sizeof(double));
        }
    }
}

void rotMol(MOLECULE *mol, double angle, char axis)
{
    double matrix[mol->atomNum][3];
    get_CoordMat(mol, matrix);
    double theta;
    theta = (angle * M_PI) / 180;
    printf("theta = %lf\n", theta);
    double cos_theta;
    double sin_theta;
    cos_theta = (((int)(angle * 1000000) % 90000000 == 0) &&
                 ((int)(angle * 1000000) % 180000000 != 0))
                    ? 0
                    : cos(theta);
    sin_theta = ((int)(angle * 1000000) % 180000000 == 0) ? 0 : sin(theta);
    printf("cos(theta): %lf, sin(theta): %lf\n", cos_theta, sin_theta);
    double rotated_coord[mol->atomNum][3];
    double rot_x[3][3] = {{1.0, 0, 0},
                          {0, cos_theta, -1.0 * sin_theta},
                          {0.0, sin_theta, cos_theta}};
    double rot_y[3][3] = {
        {cos_theta, 0, sin_theta}, {0, 1, 0}, {-1 * sin_theta, 0, cos_theta}};
    double rot_z[3][3] = {
        {cos_theta, -1 * sin_theta, 0}, {sin_theta, cos_theta, 0}, {0, 0, 1}};
    switch (axis)
    {
    case 'x':
        multiplyMatrices_mx3(matrix, rot_x, rotated_coord, mol->atomNum);
        assignCoordtoMol(rotated_coord, mol);
        break;
    case 'y':
        multiplyMatrices_mx3(matrix, rot_y, rotated_coord, mol->atomNum);
        assignCoordtoMol(rotated_coord, mol);
        break;
    case 'z':
        multiplyMatrices_mx3(matrix, rot_z, rotated_coord, mol->atomNum);
        assignCoordtoMol(rotated_coord, mol);
        break;
    default:
        perror("Wrong axis.\n");
        break;
    }
}
