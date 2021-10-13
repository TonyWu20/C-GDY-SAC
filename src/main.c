#include "main.h"
#include "maths/MyMaths.h"
#include "msiParser/msiParser.h"
#include <stdio.h>

int main(int argc, char *argv[])
{
    char *fileName = argv[1];
    FILE *file = fopen(fileName, "r");
    if (file == NULL)
    {
        printf("Cannot open file.\n");
        return 1;
    }
    MOLECULE *mol;
    mol = parseMol(file);
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
    rotMol(mol, 149.9999, 'z');
    rotMol(mol, 270.0, 'x');
    double *u, *v, *a, *b;
    u = mol->molAtoms[4].coord;
    v = mol->molAtoms[5].coord;
    a = malloc(3 * sizeof(double));
    initVector(u, v, a);

    u = model->totalAtoms[40].coord;
    v = model->totalAtoms[41].coord;
    b = malloc(3 * sizeof(double));
    initVector(u, v, b);
    double theta = VecAngle(a, b);
    rotMol(mol, -1 * theta, 'z');
    placeMol(mol, model, 40, ads_model);
    for (int i = model->atomNum; i < ads_model->atomNum; i++)
    {
        printf(
            "  (%d Atom\n    (A C ACL \"%d %s\")\n    (A C Label \"%s\")\n  "
            "  (A D XYZ (%15.13lf %15.13lf %15.13lf))\n    (A I Id %d)\n  )\n",
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
