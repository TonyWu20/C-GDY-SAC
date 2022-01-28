#include "castep_output.h"

void add_item(CastepInfo **table, struct ElmItem *item)
{
    CastepInfo *tableItem = malloc(sizeof(CastepInfo));
    tableItem->name = item->elm;
    tableItem->info = malloc(sizeof(struct ElmItem));
    memcpy(tableItem->info, item, sizeof(struct ElmItem));
    HASH_ADD_KEYPTR(hh, *table, tableItem->name, strlen(tableItem->name),
                    tableItem);
}

CastepInfo *find_item(CastepInfo *table, const char *elm)
{
    CastepInfo *ret;
    HASH_FIND_STR(table, elm, ret);
    return ret;
}

void delete_all(CastepInfo **table)
{
    CastepInfo *currItem, *tmp;
    HASH_ITER(hh, *table, currItem, tmp)
    {
        HASH_DEL(*table, currItem);
        free(currItem->info);
        free(currItem);
    }
}

CastepInfo *initTable()
{
    CastepInfo *table = NULL;
    for (int i = 0; i < 3; ++i)
    {
        struct ElmItem cur = ads[i];
        add_item(&table, &cur);
    }
    for (int i = 0; i < 10; ++i)
    {
        struct ElmItem cur = metal3D[i];
        add_item(&table, &cur);
    }
    for (int i = 0; i < 10; ++i)
    {
        struct ElmItem cur = metal4D[i];
        add_item(&table, &cur);
    }
    for (int i = 0; i < 9; ++i)
    {
        struct ElmItem cur = metal5D[i];
        add_item(&table, &cur);
    }
    for (int i = 0; i < 15; ++i)
    {
        struct ElmItem cur = metalLM[i];
        add_item(&table, &cur);
    }
    return table;
}
