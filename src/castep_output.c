#include "castep_output.h"

void add_item(CastepInfo **table, const char *elm, int lcaoVal, double massVal,
              const char *potential_file, int spinVal)
{
    CastepInfo *item = malloc(sizeof(CastepInfo));
    item->elm = elm;
    item->LCAO = lcaoVal;
    item->mass = massVal;
    item->potential_file = potential_file;
    item->spin = spinVal;
    HASH_ADD_KEYPTR(hh, *table, item->elm, strlen(item->elm), item);
}

CastepInfo *find_item(CastepInfo *table, char *elm)
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
        free(currItem);
    }
}

CastepInfo *initTable()
{
    CastepInfo *table = NULL;
    for (int i = 0; i < 3; ++i)
    {
        struct ElmItem cur = ads[i];
        add_item(&table, cur.elm, cur.LCAO, cur.mass, cur.potential_file,
                 cur.spin);
    }
    for (int i = 0; i < 10; ++i)
    {
        struct ElmItem cur = metal3D[i];
        add_item(&table, cur.elm, cur.LCAO, cur.mass, cur.potential_file,
                 cur.spin);
    }
    for (int i = 0; i < 10; ++i)
    {
        struct ElmItem cur = metal4D[i];
        add_item(&table, cur.elm, cur.LCAO, cur.mass, cur.potential_file,
                 cur.spin);
    }
    for (int i = 0; i < 9; ++i)
    {
        struct ElmItem cur = metal5D[i];
        add_item(&table, cur.elm, cur.LCAO, cur.mass, cur.potential_file,
                 cur.spin);
    }
    for (int i = 0; i < 15; ++i)
    {
        struct ElmItem cur = metalLM[i];
        add_item(&table, cur.elm, cur.LCAO, cur.mass, cur.potential_file,
                 cur.spin);
    }
    return table;
}
