#pragma once
#include "database.h"
#include "uthash.h"
#include <stdbool.h>

typedef struct
{
    char *name;
    int *coordAtoms;
    unsigned coordAtoms_count;
    int *stemAtoms;
    unsigned stemAtoms_count;
    int *planeAtoms;
    unsigned planeAtoms_count;
    bool vertical;
    bool bSym;
    int upperAtomId;
    char *pathName;
} AdsInfo;

typedef struct AdsTableYAML AdsTableYAML;
struct AdsTableYAML
{
    AdsInfo *adsInfoItem;
    unsigned adsInfoItem_count;
    void (*destroy)(AdsTableYAML **self);
};
AdsTableYAML *load_adsTableYAML(void);
void destroy_adsTableYAML(AdsTableYAML **adsTableYAML);
HashNode *init_adsInfoTable(AdsTableYAML *adsTableYAML);
