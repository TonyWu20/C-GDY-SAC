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

#define TOTAL_ELEMENT_NUM 3
#define TOTAL_MODELS (231 * TOTAL_ELEMENT_NUM)
enum
{
    T3D,
    T4D,
    T5D,
    LM
};
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

char *elements[] = {"Sc", "Ti", "V",  "Cr", "Mn", "Fe", "Co", "Ni", "Cu",
                    "Zn", "Y",  "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd",
                    "Ag", "Cd", "Hf", "Ta", "W",  "Re", "Os", "Ir", "Pt",
                    "Au", "Hg", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu",
                    "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu"};

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

void generateModels(int total_tasks, int pathNameCode, int adsListLen,
                    char **adsList, char **elements, char **pathways,
                    HashNode *elmTable, HashNode *adsTable, int *progress)
{
    /* int i, k; */
    // clang-format off
    Lattice *lat = parse_lattice_from_file("./msi_models/3d/SAC_GDY_Sc.msi", "SAC_GDY_Sc");
    /* #pragma omp parallel private(lat) shared (adsList, elements, elmTable, adsTable, progress) */
    // clang-format on
    {
        // clang-format off
        /* #pragma omp for */
        // clang-format on
        for (int i = 0; i < total_tasks; ++i)
        {
            char *curElement = elements[i / adsListLen];
            ElmInfo *curElmInfo = find_item_by_str(elmTable, curElement)->val;
            lat->vtable->modify_metal(lat, curElement, curElmInfo->atomicNum);
            Adsorbate *ads = getAds_from_i(i, adsList, adsListLen, adsTable);
            for (int k = 0; k < ads->taskLists->taskNum; ++k)
            {
                Adsorbate *ads_copy = ads->vtable->duplicate(ads);
                Lattice *result = Add_mol_to_lattice(
                    lat, ads_copy, ads_copy->taskLists->tasks[k][0],
                    ads_copy->taskLists->tasks[k][1],
                    heightChoice[pathNameCode]);
                result->vtable->export_msi(result);
                Cell *cell = createCell(result, elmTable);
                cell->vtable->exportCell(cell);
                write_param(cell);
                write_kptaux(cell);
                write_trjaux(cell);
                write_pbsScript(cell);
                write_SMCastepExtension(cell);
                /* copy_potentials(cell, potTable); */
                (*progress)++;
                double percentage =
                    (double)(*progress) / (double)(TOTAL_MODELS);
                printProgress(*progress, TOTAL_MODELS, percentage,
                              result->mol->name);
                cell->destroy(cell);
                ads_copy->vtable->destroy(ads_copy);
            }
            // clang-format off
            ads->vtable->destroy(ads);
            // clang-format on
        }
    }
    lat->vtable->destroy(lat);
}

void allocateTasks(int pathNameCode, int *progress, HashNode *elmTable,
                   HashNode *adsTable)
{
    int adsListLen = 0;
    char **adsList = pathway_adsLists(pathNameCode, &adsListLen);
    int total_tasks = TOTAL_ELEMENT_NUM * adsListLen;
    generateModels(total_tasks, pathNameCode, adsListLen, adsList, elements,
                   pathways, elmTable, adsTable, progress);
    for (int i = 0; i < adsListLen; ++i)
    {
        free(adsList[i]);
    }
    free(adsList);
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
