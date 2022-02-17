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

#define TOTAL_ELEMENT_NUM 44
#define TOTAL_MODELS (231 * TOTAL_ELEMENT_NUM)
enum
{
    T3D,
    T4D,
    T5D,
    LM
};

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
char *pathways[] = {"ethylene", "acetic_acid", "ethanol", "ethanol_other_ads"};
char *ethylene_ads[] = {
    "CO.msi",      "CHO.msi",      "COCHO.msi",  "COCHOH.msi",  "OCH2CO.msi",
    "OCH2COH.msi", "OCH2CHOH.msi", "OCH2CH.msi", "OCH2CH2.msi", "C2H4.msi"};
char *acetic_acid_ads[] = {"OCH2C_cyc_OH.msi", "CH3COOH.msi"};
char *ethanol_ads[] = {
    "Glyoxal.msi",        "HOCCHO.msi",   "HOHCCHO.msi",
    "Glycolaldehyde.msi", "CH2CHO.msi",   "CH2CHOH.msi",
    "CH2CH2OH.msi",       "CH3CH2OH.msi", "acetaldehyde.msi"};
char *ethanol_other_ads[] = {"CH2OHCH2O.msi", "ethylene_glycol.msi"};

double heightChoice[] = {1.4, 1.6, 1.4, 1.4};

Adsorbate *getAds_from_i(int i, char **adsList, int adsListLen,
                         HashNode *adsTable)
{
    int currAds = i % adsListLen;
    char *ads_name = extractStemName(adsList[currAds]);
    AdsInfo *info = (AdsInfo *)find_item_by_str(adsTable, ads_name)->val;
    Adsorbate *ads =
        parse_adsorbate_from_file(adsList[currAds], ads_name, info);
    free(ads_name);
    return ads;
}

char **pathway_adsLists(int pathNameCode, int *adsListLen)
{
    *adsListLen = 0;
    char **list = NULL;
    switch (pathNameCode)
    {
    case ETHYLENE:
        *adsListLen = sizeof(ethylene_ads) / sizeof(char *);
        list = adsListBuild("ethylene", ethylene_ads, *adsListLen);
        break;
    case ACETIC_ACID:
        *adsListLen = sizeof(acetic_acid_ads) / sizeof(char *);
        list = adsListBuild("acetic_acid", acetic_acid_ads, *adsListLen);
        break;
    case ETHANOL:
        *adsListLen = sizeof(ethanol_ads) / sizeof(char *);
        list = adsListBuild("ethanol", ethanol_ads, *adsListLen);
        break;
    case ETHANOL_OTHER:
        *adsListLen = sizeof(ethanol_other_ads) / sizeof(char *);
        list = adsListBuild("ethanol_other", ethanol_other_ads, *adsListLen);
        break;
    default:
        break;
    }
    return list;
}

char **adsListBuild(char *pathName, char *adsList[], int listLen)
{
    char **list = malloc(sizeof(char *) * listLen);
    for (int i = 0; i < listLen; ++i)
    {
        int filePathLen = 1 + snprintf(NULL, 0, "./C2_pathways_ads/%s_path/%s",
                                       pathName, adsList[i]);
        list[i] = malloc(filePathLen);
        snprintf(list[i], filePathLen, "./C2_pathways_ads/%s_path/%s", pathName,
                 adsList[i]);
    }
    return list;
}

Adsorbate *getCurrendAds(AdsTableYAML *adsTableYAML, int adsId)
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
                    Adsorbate *ads = getCurrendAds(adsTableYAML, j);
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
