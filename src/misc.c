#include "misc.h"
#include "assemble.h"
#include "lattice.h"
#include "molecule.h"
#include "parser.h"
#include <stdio.h>
#include <stdlib.h>
#define PBSTR "||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||"
#define PBWIDTH 60
static void printProgress(int cur, int total, double percentage, char *name)
{
    int val = (int)(percentage * 100);
    int lpad = (int)(percentage * PBWIDTH);
    int rpad = PBWIDTH - lpad;
    int strPad = 32 - strlen(name);
    printf("\e[32m \rNow %.*s%*s %d/%d %3d%% [%.*s%*s]\e[m", 32, name, strPad,
           "", cur, total, val, lpad, PBSTR, rpad, "");
    fflush(stdout);
}
enum
{
    T3D,
    T4D,
    T5D,
    LM
};
char *pathways[] = {"ethylene"};
char *ethylene_ads[] = {"COCHO.msi",   "COCHOH.msi",   "OCH2CO.msi",
                        "OCH2COH.msi", "OCH2CHOH.msi", "OCH2CH.msi",
                        "OCH2CH2.msi", "C2H4.msi"};
char *elements[] = {"Sc", "Ti", "V",  "Cr", "Mn", "Fe", "Co", "Ni", "Cu",
                    "Zn", "Y",  "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd",
                    "Ag", "Cd", "Hf", "Ta", "W",  "Re", "Os", "Ir", "Pt",
                    "Au", "Hg", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu",
                    "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu"};

void allocateTasks(char *pathName)
{
    int adsListLen = 0;
    char **adsList = pathway_adsLists(pathName, &adsListLen);
    int total_tasks = TOTAL_ELEMENT_NUM * adsListLen;
    for (int i = 0; i < total_tasks; ++i)
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
        Adsorbate *ads = parse_molecule_from_file(adsList[currAds], ads_name);
        for (int k = 0; k < ads->taskLists->taskNum; ++k)
        {
            Lattice *result =
                Add_mol_to_lattice(lat, ads, ads->taskLists->tasks[k][0],
                                   ads->taskLists->tasks[k][1], pathName);
            result->vtable->export_msi(result, pathName);
            double percentage = (double)(i + 1) / (double)total_tasks;
            printProgress(i + 1, total_tasks, percentage, result->_mol->name);
            result->vtable->destroy(result);
        }
        ads->ads_vtable->destroy(ads);
        free(ads_name);
        free(basePath);
        free(baseName);
        lat->vtable->destroy(lat);
    }
    for (int i = 0; i < adsListLen; ++i)
    {
        free(adsList[i]);
    }
    free(adsList);
}

void createDirectory(char *dest)
{
    struct stat s;
    int err = stat(dest, &s);
    if (err == -1)
    {
        int destStringLen = strlen(dest);
        int commandLen = destStringLen + 9;
        char mkdir_command[commandLen + 1];
        snprintf(mkdir_command, commandLen + 1, "mkdir -p %s", dest);
        system(mkdir_command);
    }
}

char *extractStemName(char *filepath)
{
    char *buffer = strdup(filepath);
    char *last = strrchr(buffer, '/');
    char *stempath = NULL;
    if (last != NULL)
        stempath = last + 1;
    char *token, *rest;
    token = strtok_r(stempath, ".", &rest);
    char *ret = strdup(token);
    free(buffer);
    return ret;
}

char *findBaseByElementId(int i)
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

    int filePathLen = 14 + 3 + 13 + strlen(elements[i]);
    char *filePath = malloc(filePathLen + 1);
    snprintf(filePath, filePathLen, "./msi_models/%s/SAC_GDY_%s.msi", family,
             elements[i]);
    free(family);
    return filePath;
}

char **pathway_adsLists(char *pathName, int *adsListLen)
{
    *adsListLen = 0;
    char **list = NULL;
    int baseLen = 25 + strlen(pathName);
    if (!strcmp(pathName, "ethylene"))
    {
        *adsListLen = 8;
        list = malloc(sizeof(char *) * (*adsListLen));
        for (int i = 0; i < 8; ++i)
        {
            int filePathLen = baseLen + strlen(ethylene_ads[i]);
            list[i] = malloc(filePathLen + 1);
            snprintf(list[i], filePathLen + 1, "./C2_pathways_ads/%s_path/%s",
                     pathName, ethylene_ads[i]);
        }
    }
    return list;
}
