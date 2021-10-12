#include "main.h"
#include "maths/MyMaths.h"
#include "msiParser/mod_msi.h"
#include <stdio.h>

ATOM_BLOCK ads[10];

int main(int argc, char *argv[])
{
    char *fileName = argv[1];
    FILE *file = fopen(fileName, "r");
    if (file == NULL)
    {
        printf("Cannot open file.\n");
        return 1;
    }
    int atomCount;
    atomCount = countAtoms(file);
    MOLECULE *mol = init_mol(atomCount);
    scanAtom(file, mol->molAtoms);
    printf("Parse molecule file %s successfully! Has %d atoms.\n", fileName,
           mol->atomNum);
    resetXYZ(atomCount, mol->molAtoms);
    printf("Reset coordinates\n");
    fclose(file);
    BASE_LATTICE *model;
    fileName = argv[2];
    file = fopen(fileName, "r");
    if (file == NULL)
    {
        printf("Cannot open file.\n");
        return 1;
    }
    model = parseBase(file);
    fclose(file);
    BASE_LATTICE *ads_model;
    ads_model = init_adsorbed_lat(model, mol);
    for (int i = 0; i < 3; i++)
    {
        printf("Vector: %lf, %lf, %lf\n", ads_model->latVector[i][0],
               ads_model->latVector[i][1], ads_model->latVector[i][2]);
    }
    double *u, *v, *a, *b;
    u = mol->molAtoms[0].coord;
    v = mol->molAtoms[1].coord;
    a = malloc(3 * sizeof(double));
    b = malloc(3 * sizeof(double));
    initVector(u, v, a);
    u = model->totalAtoms[41].coord;
    v = model->totalAtoms[42].coord;
    initVector(u, v, b);
    rotMol(mol, 149.9999, 'z');
    rotMol(mol, 270.0, 'x');
    double theta = VecAngle(a, b);
    printf("Angle deviated %lf\n", theta);
    for (int i = 0; i < mol->atomNum; i++)
    {
        printf("%lf, %lf, %lf\n", mol->molAtoms[i].coord[0],
               mol->molAtoms[i].coord[1], mol->molAtoms[i].coord[2]);
    }
    rotMol(mol, theta, 'z');
    for (int i = 0; i < mol->atomNum; i++)
    {
        printf("%lf, %lf, %lf\n", mol->molAtoms[i].coord[0],
               mol->molAtoms[i].coord[1], mol->molAtoms[i].coord[2]);
    }
    appendMolAtoms(model, mol, ads_model);
    for (int i = model->atomNum; i < ads_model->atomNum; i++)
    {
        printf("  (%d Atom\n    (A C ACL \"%d %s\")\n    (A C Label \"%s\")\n  "
               "  (A D XYZ (%15.13lf %15.13lf %15.13lf))\n    (A I Id %d)\n)\n",
               ads_model->totalAtoms[i].itemId, ads_model->totalAtoms[i].elmId,
               ads_model->totalAtoms[i].elm, ads_model->totalAtoms[i].elm,
               ads_model->totalAtoms[i].coord[0],
               ads_model->totalAtoms[i].coord[1],
               ads_model->totalAtoms[i].coord[2],
               ads_model->totalAtoms[i].itemId - 1);
    }
    free(mol);
    free(ads_model);
    free(model);
    free(a);
}
