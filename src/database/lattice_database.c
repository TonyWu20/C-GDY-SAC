#include "lattice_database.h"
#include "cyaml/cyaml.h"
#include <stdio.h>

static const cyaml_schema_field_t element_info_field_schema[] = {
    CYAML_FIELD_STRING_PTR("element", CYAML_FLAG_POINTER, ElmInfo, name, 0,
                           CYAML_UNLIMITED),
    CYAML_FIELD_UINT("atomic_num", CYAML_FLAG_DEFAULT, ElmInfo, atomicNum),
    CYAML_FIELD_UINT("LCAO", CYAML_FLAG_DEFAULT, ElmInfo, LCAO),
    CYAML_FIELD_FLOAT("mass", CYAML_FLAG_DEFAULT, ElmInfo, mass),
    CYAML_FIELD_STRING_PTR("pot", CYAML_FLAG_POINTER, ElmInfo, potPath, 0,
                           CYAML_UNLIMITED),
    CYAML_FIELD_UINT("spin", CYAML_FLAG_DEFAULT, ElmInfo, spin),
    CYAML_FIELD_END};

static const cyaml_schema_value_t element_info_schema = {CYAML_VALUE_MAPPING(
    CYAML_FLAG_DEFAULT, ElmInfo, element_info_field_schema)};

static const cyaml_schema_field_t elm_table_fields_schema[] = {

    CYAML_FIELD_SEQUENCE("Element_info", CYAML_FLAG_POINTER,

                         struct element_table_yaml, infoItems,
                         &element_info_schema,

                         0, CYAML_UNLIMITED),

    CYAML_FIELD_END};

static const cyaml_schema_value_t elm_table_schema = {CYAML_VALUE_MAPPING(

    CYAML_FLAG_POINTER, struct element_table_yaml, elm_table_fields_schema)};

static const cyaml_config_t config = {
    .log_fn = cyaml_log,            /* Use the default logging function. */
    .mem_fn = cyaml_mem,            /* Use the default memory allocator. */
    .log_level = CYAML_LOG_WARNING, /* Logging errors and warnings only. */
};

struct element_table_yaml *load_elmTableYAML(void)
{
    cyaml_err_t err;
    struct element_table_yaml *elmTableYAML;
    err = cyaml_load_file("./src/database/element_table.yaml", &config,
                          &elm_table_schema, (void **)&elmTableYAML, NULL);
    if (err != CYAML_OK)
    {
        printf("ERROR: %s\n", cyaml_strerror(err));
        return NULL;
    }
    return elmTableYAML;
}

void destroy_element_table_yaml(struct element_table_yaml **elmTableYAML)
{
    cyaml_free(&config, &elm_table_schema, *elmTableYAML, 0);
}

HashNode *init_ElmInfoTable(struct element_table_yaml *elmTableYAML)
{
    HashNode *hashTab = NULL;
    for (int i = 0; i < elmTableYAML->infoItems_count; ++i)
    {
        ElmInfo *curItem = &elmTableYAML->infoItems[i];
        add_str_keyptr_item(&hashTab, curItem->name, (void *)curItem);
    }
    return hashTab;
}
