#include "database.h"
#include <stdio.h>

void add_str_keyptr_item(HashNode **hashTab, const char *key, void *hashValItem)

{
    HashNode *np = malloc(sizeof(HashNode));
    if (!np)
    {
        printf("Malloc memory for HashNode failed!\n");
        return;
    }
    np->key = key;
    np->val = hashValItem;
    HASH_ADD_KEYPTR(hh, *hashTab, np->key, strlen((const char *)np->key), np);
}

HashNode *find_item_by_str(HashNode *hashTab, char *key)
{
    HashNode *np;
    printf("%s\n", key);
    HASH_FIND_STR(hashTab, key, np);
    if (!np)
    {
        printf("Cannot find hash item %s\n", key);
    }
    return np;
}

void delete_all(HashNode **hashTab)
{
    HashNode *curr, *tmp;
    HASH_ITER(hh, *hashTab, curr, tmp)
    {
        HASH_DEL(*hashTab, curr);
        free(curr);
    }
}
