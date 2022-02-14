#pragma once
#include "database.h"
#include "uthash.h"
#include <stdbool.h>
extern int task_cd2_sym[7][2];
extern int task_cd2_asym[14][2];
extern int task_cd1[7][2];

typedef struct
{
    char *name;
    int coordAtomNum;
    int coordAtomIds[2];
    int stemAtomIds[2];
    int planeAtomIds[3];
    bool vertical;
    bool bSym;
    int upperAtomId;
} AdsorbateInfo;
HashNode *init_adsInfoTable();
