#include "main.h"
#include "msiParser/msiParser.h"
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
    fileName = "SAC_GDY_V.msi";
    file = fopen(fileName, "r");
    model = parseBase(file);
    fclose(file);
    BASE_LATTICE *ads_model;
    ads_model = init_adsorbed_lat(model, mol);
    for (int i = 0; i < 3; i++)
    {
        printf("Vector: %lf, %lf, %lf\n", ads_model->latVector[i][0],
               ads_model->latVector[i][1], ads_model->latVector[i][2]);
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
}
