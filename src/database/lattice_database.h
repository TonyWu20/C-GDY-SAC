#include "database.h"

typedef struct
{
    const char *name;
    int atomicNum;
    int LCAO;
    double mass;
    const char *potPath;
    int spin;
} ElmInfo;

struct element_table_yaml
{
    ElmInfo *infoItems;
    int infoItems_count;
};
void test_elm_table(void);
HashNode *init_ElmInfoTable();
