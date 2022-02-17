#pragma once
#include "database.h"

typedef struct
{
    char *name;
    int atomicNum;
    int LCAO;
    double mass;
    char *potFile;
    int spin;
} ElmInfo;

typedef struct element_table_yaml ElmTableYAML;
struct element_table_yaml
{
    ElmInfo *infoItems;
    int infoItems_count;
    void (*destroy)(ElmTableYAML **self);
};
void test_elm_table(void);
ElmTableYAML *load_elmTableYAML(void);
void destroy_ElmTableYAML(ElmTableYAML **elmTableYAML);
HashNode *init_ElmInfoTable(ElmTableYAML *elmTableYAML);
HashNode *init_ElmInfoTable_internal();
