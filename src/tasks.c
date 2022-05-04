#include "tasks.h"
#include "misc.h"
#include "param.h"
#include "parser.h"
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <unistd.h>
/* static char *findBaseByElementId(int i); */

#define TOTAL_ELEMENT_NUM 1
#define TOTAL_MODELS (231 * TOTAL_ELEMENT_NUM)

enum
{
    C1 = 41,
    C2 = 42,
    C3 = 54,
    C4 = 53,
    FR = 52,
    NR = 40,
    M = 73
} site_codes;

int task_cd[][2] = {{C1, C2}, {C2, C3}, {C3, C4}, {C4, FR},
                    {NR, C1}, {C1, M},  {C2, M}};
int task_cd1[][2] = {{C1, NULLSITE}, {C2, NULLSITE}, {C3, NULLSITE},
                     {C4, NULLSITE}, {NR, NULLSITE}, {FR, NULLSITE},
                     {M, NULLSITE}};

Adsorbate *getCurrentAds(AdsTableYAML *adsTableYAML, int adsId)
{
    AdsInfo curAdsInfo = adsTableYAML->adsInfoItem[adsId];
    char *path;
    asprintf(&path, "./C2_pathways_ads/%s_path/%s.msi", curAdsInfo.pathName,
             curAdsInfo.name);
    Adsorbate *ads =
        parse_adsorbate_from_file(path, curAdsInfo.name, &curAdsInfo);
    free(path);
    return ads;
}

void exportAll(Lattice *base, Adsorbate *ads, int c1, int c2,
               HashNode *elmTable, int *progress)
{
    Lattice *result = Add_mol_to_lattice(base, ads, c1, c2, 1.4);
    result->vtable->export_msi(result);
    Cell *cell = createCell(result, elmTable);
    cell->vtable->exportCell(cell);
    cell->vtable->exportSeeds(cell);
    (*progress)++;
    double percentage = ((double)*progress) / (double)(TOTAL_MODELS);
    printProgress(*progress, TOTAL_MODELS, percentage,
                  cell->lattice->mol->name);
    cell->destroy(cell);
}

void generator(ElmTableYAML *elmTableYAML, char **elements,
               AdsTableYAML *adsTableYAML, int *progress)
{
    HashNode *elmTable = init_ElmInfoTable(elmTableYAML);
    int adsNums = adsTableYAML->adsInfoItem_count;
    int curElmId = 0;
// clang-format off
    #pragma omp parallel shared(adsNums, curElmId, elmTable, adsTableYAML)
    // clang-format on
    {
        Lattice *base = parse_lattice_from_file(
            "./msi_models/3d/SAC_GDY_Sc.msi", "SAC_GDY_Sc");
        // clang-format off
        #pragma omp for collapse(3)
        // clang-format on 
        for (int j = 0; j < adsNums; ++j)
        {
            for (int i = 0; i < TOTAL_ELEMENT_NUM; ++i)
            {
                for (int k = 0; k < 7; ++k)
                {
                    HashNode *nextElm =
                        find_item_by_str(elmTable, elements[i]);
                    int nextElmId = ((ElmInfo *)(nextElm->val))->atomicNum;
                    base->vtable->modify_metal(base, nextElm->key,
                                                nextElmId);
                    Adsorbate *ads = getCurrentAds(adsTableYAML, j);
                    switch (ads->coordAtomNum)
                    {
                    case 2:
                    {
                        exportAll(base, ads, task_cd[k][0], task_cd[k][1],
                                  elmTable, progress);
                        if (ads->bSym==false)
                        {
                            exportAll(base, ads, task_cd[k][1], task_cd[k][0],
                                      elmTable, progress);
                        }
                        break;
                    }
                    case 1:
                        exportAll(base, ads, task_cd1[k][0], task_cd1[k][1], elmTable,
                                  progress);
                        break;
                    default:
                        break;
                    }
                    ads->vtable->destroy(ads);
                }
            }
        }
        base->vtable->destroy(base);
    }
    delete_all(&elmTable);
}
