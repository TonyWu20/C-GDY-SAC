#pragma once
#include "database.h"
#include "uthash.h"
#include <stdbool.h>
extern int task_cd2_sym[7][2];
extern int task_cd2_asym[14][2];
extern int task_cd1[7][2];

typedef struct {
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
    int atomNums;
    char *pathName;
} AdsInfo;

typedef struct AdsTableYAML AdsTableYAML;
struct AdsTableYAML {
    AdsInfo *adsInfoItem;
    unsigned adsInfoItem_count;
    void (*destroy)(AdsTableYAML **self);
};
AdsTableYAML *load_adsTableYAML(void);
void destroy_adsTableYAML(AdsTableYAML **adsTableYAML);
HashNode *init_adsInfoTable(AdsTableYAML *adsTableYAML);
