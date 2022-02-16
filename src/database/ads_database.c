#include "ads_database.h"
#include "cyaml/cyaml.h"
#include <stdio.h>
#define NULLSITE 255
enum
{
    C1 = 41,
    C2 = 42,
    C3 = 54,
    C4 = 53,
    FR = 52,
    NR = 40,
    M = 73
} site_codes;

int task_cd2_sym[][2] = {{C1, C2}, {C2, C3}, {C3, C4}, {C4, FR},
                         {NR, C1}, {C1, M},  {C2, M}};
int task_cd2_asym[][2] = {{C1, C2}, {C2, C1}, {C2, C3}, {C3, C2}, {C3, C4},
                          {C4, C3}, {C4, FR}, {FR, C4}, {NR, C1}, {C1, NR},
                          {C1, M},  {M, C1},  {C2, M},  {M, C2}};
int task_cd1[][2] = {{C1, NULLSITE}, {C2, NULLSITE}, {C3, NULLSITE},
                     {C4, NULLSITE}, {FR, NULLSITE}, {NR, NULLSITE},
                     {M, NULLSITE}};

static const cyaml_schema_value_t id_entry = {
    CYAML_VALUE_INT(CYAML_FLAG_DEFAULT, int),
};

static const cyaml_schema_field_t ads_info_field_schema[] = {
    CYAML_FIELD_STRING_PTR("name", CYAML_FLAG_POINTER, AdsInfo, name, 0,
                           CYAML_UNLIMITED),
    CYAML_FIELD_SEQUENCE("coordAtomIds", CYAML_FLAG_POINTER, AdsInfo,
                         coordAtoms, &id_entry, 0, CYAML_UNLIMITED),
    CYAML_FIELD_SEQUENCE("stemAtomIds", CYAML_FLAG_POINTER, AdsInfo, stemAtoms,
                         &id_entry, 0, CYAML_UNLIMITED),
    CYAML_FIELD_SEQUENCE("planeAtomIds", CYAML_FLAG_POINTER, AdsInfo,
                         planeAtoms, &id_entry, 0, CYAML_UNLIMITED),
    CYAML_FIELD_BOOL("vertical", CYAML_FLAG_DEFAULT, AdsInfo, vertical),
    CYAML_FIELD_BOOL("bSym", CYAML_FLAG_DEFAULT, AdsInfo, bSym),
    CYAML_FIELD_UINT("upperAtomId", CYAML_FLAG_DEFAULT, AdsInfo, upperAtomId),
    CYAML_FIELD_STRING_PTR("pathName", CYAML_FLAG_POINTER, AdsInfo, pathName, 0,
                           CYAML_UNLIMITED),
    CYAML_FIELD_END,
};

static const cyaml_schema_value_t ads_info_schema = {

    CYAML_VALUE_MAPPING(CYAML_FLAG_DEFAULT, AdsInfo, ads_info_field_schema)};

static const cyaml_schema_field_t ads_table_fields_schema[] = {
    CYAML_FIELD_SEQUENCE("Adsorbates", CYAML_FLAG_POINTER, AdsTableYAML,
                         adsInfoItem, &ads_info_schema, 0, CYAML_UNLIMITED),
    CYAML_FIELD_END,
};

static const cyaml_schema_value_t ads_table_schema = {
    CYAML_VALUE_MAPPING(CYAML_FLAG_POINTER, AdsTableYAML,
                        ads_table_fields_schema),
};

static const cyaml_config_t config = {
    .log_fn = cyaml_log,            /* Use the default logging function. */
    .mem_fn = cyaml_mem,            /* Use the default memory allocator. */
    .log_level = CYAML_LOG_WARNING, /* Logging errors and warnings only. */
};

AdsTableYAML *load_adsTableYAML(void)
{
    cyaml_err_t err;
    AdsTableYAML *adsTableYAML;
    err = cyaml_load_file("./src/database/ads_table.yaml", &config,
                          &ads_table_schema, (void **)&adsTableYAML, NULL);
    if (err != CYAML_OK)
    {
        printf("ERROR: %s\n", cyaml_strerror(err));
        return NULL;
    }
    adsTableYAML->destroy = destroy_adsTableYAML;
    return adsTableYAML;
}
void destroy_adsTableYAML(AdsTableYAML **adsTableYAML)
{
    cyaml_free(&config, &ads_table_schema, *adsTableYAML, 0);
}

HashNode *init_adsInfoTable(AdsTableYAML *adsTableYAML)
{
    HashNode *hashTab = NULL;
    for (int i = 0; i < adsTableYAML->adsInfoItem_count; ++i)
    {
        AdsInfo *curPtr = &adsTableYAML->adsInfoItem[i];
        add_str_keyptr_item(&hashTab, curPtr->name, (void *)curPtr);
    }
    return hashTab;
}
