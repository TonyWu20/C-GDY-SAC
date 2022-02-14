#include "database.h"

typedef struct
{
    char *name;
    int atomicNum;
    int LCAO;
    double mass;
    char *potPath;
    int spin;
} ElmInfo;

struct element_table_yaml
{
    ElmInfo *infoItems;
    int infoItems_count;
};
void test_elm_table(void);
struct element_table_yaml *load_elmTableYAML(void);
void destroy_element_table_yaml(struct element_table_yaml **elmTableYAML);
HashNode *init_ElmInfoTable(struct element_table_yaml *elmTableYAML);
