#pragma once
#include "uthash.h"

typedef struct
{
    const void *key;
    void *val;
    UT_hash_handle hh;
} HashNode;

void add_str_keyptr_item(HashNode **hashTab, char *key, void *hashValItem);
HashNode *find_item_by_str(HashNode *hashTab, const char *key);
void delete_all(HashNode **hashTab);
