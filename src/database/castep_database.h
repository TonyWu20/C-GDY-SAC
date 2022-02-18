#pragma once
#include "uthash.h"
typedef struct
{
    const char *name;
    struct ElmItem *info;
    UT_hash_handle hh;
} CastepInfo;

typedef struct
{
    const char *elm;
    const char *potential_file;
    char *fileContent;
    UT_hash_handle hh;
} PotentialFile;

struct ElmItem
{
    const char *elm;
    int LCAO;
    double mass;
    const char *potential_file;
    int spin;
};
void add_item(CastepInfo **table, struct ElmItem *item);

CastepInfo *find_item(CastepInfo *table, const char *elm);

void delete_all(CastepInfo **table);

CastepInfo *initTable(void);

PotentialFile *initPotTable(void);

void add_PotItem(PotentialFile **table, struct ElmItem *item);
PotentialFile *find_PotItem(PotentialFile *table, const char *elm);
void delete_PotAll(PotentialFile **table);
