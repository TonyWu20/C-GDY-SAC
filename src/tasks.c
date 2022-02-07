#include "tasks.h"
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
static char *findBaseByElementId(int i);

#define TOTAL_ELEMENT_NUM 44
enum
{
    T3D,
    T4D,
    T5D,
    LM
};
char *pathways[] = {"ethylene", "acetic_acid", "ethanol", "ethanol_other_ads"};
char *ethylene_ads[] = {"COCHO.msi",   "COCHOH.msi",   "OCH2CO.msi",
                        "OCH2COH.msi", "OCH2CHOH.msi", "OCH2CH.msi",
                        "OCH2CH2.msi", "C2H4.msi"};
char *acetic_acid_ads[] = {"OCH2C_cyc_OH.msi", "CH3COOH.msi"};
char *ethanol_ads[] = {"Glyoxal.msi",  "HOCCHO.msi",   "HOHCCHO.msi",
                       "CH2OHCHO.msi", "CH2CHO.msi",   "CH2CHOH.msi",
                       "CH2CH2OH.msi", "CH3CH2OH.msi", "acetaldehyde.msi"};
char *ethanol_other_ads[] = {"Glycolaldehyde.msi", "CH2OHCH2OH.msi",
                             "ethlylene_glycol.msi"};

char *elements[] = {"Sc", "Ti", "V",  "Cr", "Mn", "Fe", "Co", "Ni", "Cu",
                    "Zn", "Y",  "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd",
                    "Ag", "Cd", "Hf", "Ta", "W",  "Re", "Os", "Ir", "Pt",
                    "Au", "Hg", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu",
                    "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu"};

void allocateTasks(int pathNameCode)
{
    int adsListLen = 0;
    char **adsList = pathway_adsLists(pathNameCode, &adsListLen);
    int total_tasks = TOTAL_ELEMENT_NUM * adsListLen;
    CastepInfo *table = initTable();
    int i, k;
    int pg = 0;
    // clang-format off
    #pragma omp parallel private(k) shared(table, adsList)
    // clang-format on
    {
        // clang-format off
        #pragma omp for
        // clang-format on
        for (i = 0; i < total_tasks; ++i)
        {
            int currElement = i / adsListLen;
            int currAds = i % adsListLen;
            int baseNameLen = 8 + strlen(elements[currElement]);
            char *baseName = malloc(baseNameLen + 1);
            snprintf(baseName, baseNameLen + 1, "SAC_GDY_%s",
                     elements[currElement]);
            char *basePath = findBaseByElementId(currElement);
            Lattice *lat = parse_lattice_from_file(basePath, baseName);
            char *ads_name = extractStemName(adsList[currAds]);
            Adsorbate *ads =
                parse_molecule_from_file(adsList[currAds], ads_name);
            for (k = 0; k < ads->taskLists->taskNum; ++k)
            {
                Lattice *result = Add_mol_to_lattice(
                    lat, ads, ads->taskLists->tasks[k][0],
                    ads->taskLists->tasks[k][1], pathways[pathNameCode]);
                result->vtable->export_msi(result, pathways[pathNameCode]);
                Cell *cell = createCell(result, table);
                cell->vtable->exportCell(cell);
                cell->destroy(cell);
            }
            pg++;
// clang-format off
            #pragma omp critical
            // clang-format on
            {
                double percentage = (double)(pg) / (double)total_tasks;
                printProgress(pg, total_tasks, percentage, lat->_mol->name);
                ads->ads_vtable->destroy(ads);
                free(ads_name);
                free(basePath);
                free(baseName);
                lat->vtable->destroy(lat);
            }
        }
    }
    for (int i = 0; i < adsListLen; ++i)
    {
        free(adsList[i]);
    }
    free(adsList);
    delete_all(&table);
}

/* Malloc file path of the base model */
static char *findBaseByElementId(int i)
{
    char *family;
    if (i < 10)
        family = strdup("3d");
    else if (i >= 10 && i < 20)
        family = strdup("4d");
    else if (i >= 20 && i < 29)
        family = strdup("5d");
    else
        family = strdup("lm");

    int filePathLen = 1 + snprintf(NULL, 0, "./msi_models/%s/SAC_GDY_%s.msi",
                                   family, elements[i]);
    char *filePath = malloc(filePathLen);
    snprintf(filePath, filePathLen, "./msi_models/%s/SAC_GDY_%s.msi", family,
             elements[i]);
    free(family);
    return filePath;
}

char **pathway_adsLists(int pathNameCode, int *adsListLen)
{
    *adsListLen = 0;
    char **list = NULL;
    switch (pathNameCode)
    {
    case ETHYLENE:
        *adsListLen = 8;
        list = adsListBuild("ethylene", ethylene_ads, *adsListLen);
        break;
    case ACETIC_ACID:
        *adsListLen = 2;
        list = adsListBuild("acetic_acid", acetic_acid_ads, *adsListLen);
        break;
    case ETHANOL:
        *adsListLen = 9;
        list = adsListBuild("ethanol", ethanol_ads, *adsListLen);
        break;
    case ETHANOL_OTHER:
        *adsListLen = 3;
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
