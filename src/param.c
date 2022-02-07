#include "param.h"

int getFinalCutoffEnergy(Cell *cell)
{
    CastepInfo *table = cell->infoTab;
    int energy = 0;
    for (int i = 0; i < cell->elmNums; ++i)
    {
        CastepInfo *item = find_item(table, cell->elmLists[i]);
        int fineEnergy = parse_fineCutoffEnergy(item->info->potential_file);
        int ultraFineEnergy = roundupBiggerTenth(fineEnergy) * 1.1;
        energy = (ultraFineEnergy > energy) ? ultraFineEnergy : energy;
    }
    return energy;
}

int parse_fineCutoffEnergy(const char *fileName)
{
    FILE *f = fopen(fileName, "r");
    int n = 4;
    char str[16];
    for (int i = 0; i < n; ++i)
    {
        fgets(str, 16, f);
    }
    char quality[4];
    int energy;
    sscanf(str, "%d %s", &energy, quality);
    fclose(f);
    return energy;
}

char *load_geomParamFile(void)
{
    const char geomParamPath[] = "./resources/geom.param";
    char *ret = readWholeFile(geomParamPath);
    return ret;
}

char *load_dosParamFile(void)
{
    const char dosParamPath[] = "./resources/dos.param";
    char *ret = readWholeFile(dosParamPath);
    return ret;
}
