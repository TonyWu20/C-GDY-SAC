#include "lattice_database.h"
#include "cyaml/cyaml.h"
#include <stdio.h>

ElmInfo elements[] = {
    {"C", 6, 2, 12.0109996796, "C_00PBE.usp", 0},
    {"H", 0, 1, 1.0080000162, "H_00PBE.usp", 0},
    {"O", 8, 2, 15.9989995956, "O_00PBE.usp", 0},
    {"Co", 27, 3, 58.9329986572, "Co_00PBE.uspcc", 3},
    {"Cr", 24, 5, 51.9959983826, "Cr_00PBE.usp", 3},
    {"Cu", 29, 3, 63.5460014343, "Cu_00PBE.usp", 1},
    {"Fe", 26, 3, 55.8470001221, "Fe_00PBE.uspcc", 4},
    {"Mn", 25, 3, 54.9379997253, "Mn_00PBE.uspcc", 5},
    {"Ni", 28, 3, 58.7099990845, "Ni_00PBE.uspcc", 2},
    {"Sc", 21, 5, 44.9560012817, "Sc_00PBE.usp", 0},
    {"Ti", 22, 5, 47.9000015259, "Ti_00PBE.usp", 0},
    {"V", 23, 5, 50.9410018921, "V_00PBE.usp", 2},
    {"Zn", 30, 4, 65.3799972534, "Zn_00PBE.usp", 0},
    {"Ag", 47, 3, 107.8679962158, "Ag_00PBE.usp", 0},
    {"Cd", 48, 4, 112.4100036621, "Cd_00PBE.usp", 0},
    {"Mo", 42, 5, 95.9400024414, "Mo_00PBE.usp", 2},
    {"Nb", 41, 5, 92.90599823, "Nb_00PBE.usp", 1},
    {"Pd", 46, 3, 106.4000015259, "Pd_00PBE.usp", 2},
    {"Rh", 45, 3, 102.90599823, "Rh_00PBE.usp", 4},
    {"Ru", 44, 5, 101.0699996948, "Ru_00PBE.usp", 5},
    {"Tc", 43, 5, 98.90599823, "Tc_00PBE.usp", 3},
    {"Y", 39, 3, 88.90599823, "Y_00PBE.uspcc", 0},
    {"Zr", 40, 5, 91.2200012207, "Zr_00PBE.usp", 0},
    {"Au", 79, 3, 196.966003418, "Au_00PBE.usp", 2},
    {"Hf", 72, 3, 178.4900054932, "Hf_00PBE.uspcc", 0},
    {"Hg", 80, 4, 200.5899963379, "Hg_00PBE.usp", 0},
    {"Ir", 77, 3, 192.2200012207, "Ir_00PBE.usp", 4},
    {"Os", 76, 5, 190.1999969482, "Os_00PBE.usp", 5},
    {"Pt", 78, 3, 195.0899963379, "Pt_00PBE.usp", 4},
    {"Re", 75, 5, 186.2070007324, "Re_00PBE.usp", 3},
    {"Ta", 73, 3, 180.9479980469, "Ta_00PBE.usp", 1},
    {"W", 74, 5, 183.8500061035, "W_00PBE.usp", 2},
    {"Ce", 58, 6, 140.1199951172, "Ce_00PBE.usp", 1},
    {"Dy", 66, 6, 162.5, "Dy_00.usp", 5},
    {"Er", 68, 6, 167.2599945068, "Er_00.usp", 3},
    {"Eu", 63, 6, 151.9600067139, "Eu_00.usp", 6},
    {"Gd", 64, 6, 157.25, "Gd_00.usp", 7},
    {"Ho", 67, 6, 164.9299926758, "Ho_00PBE.usp", 4},
    {"La", 57, 6, 138.9049987793, "La_00PBE.usp", 0},
    {"Lu", 71, 4, 174.9700012207, "Lu_00.usp", 0},
    {"Nd", 60, 6, 144.2400054932, "Nd_00.usp", 3},
    {"Pm", 61, 6, 147.0, "Pm_00.usp", 4},
    {"Pr", 59, 6, 140.9080047607, "Pr_00.usp", 2},
    {"Sm", 62, 6, 150.3999938965, "Sm_00.usp", 5},
    {"Tb", 65, 6, 158.9250030518, "Tb_00.usp", 6},
    {"Tm", 69, 6, 168.9340057373, "Tm_00.usp", 2},
    {"Yb", 70, 6, 173.0399932861, "Yb_00PBE.usp", 1},
};

static const cyaml_schema_field_t element_info_field_schema[] = {
    CYAML_FIELD_STRING_PTR("element", CYAML_FLAG_POINTER, ElmInfo, name, 0,
                           CYAML_UNLIMITED),
    CYAML_FIELD_UINT("atomic_num", CYAML_FLAG_DEFAULT, ElmInfo, atomicNum),
    CYAML_FIELD_UINT("LCAO", CYAML_FLAG_DEFAULT, ElmInfo, LCAO),
    CYAML_FIELD_FLOAT("mass", CYAML_FLAG_DEFAULT, ElmInfo, mass),
    CYAML_FIELD_STRING_PTR("pot", CYAML_FLAG_POINTER, ElmInfo, potFile, 0,
                           CYAML_UNLIMITED),
    CYAML_FIELD_UINT("spin", CYAML_FLAG_DEFAULT, ElmInfo, spin),
    CYAML_FIELD_END};

static const cyaml_schema_value_t element_info_schema = {CYAML_VALUE_MAPPING(
    CYAML_FLAG_DEFAULT, ElmInfo, element_info_field_schema)};

static const cyaml_schema_field_t elm_table_fields_schema[] = {

    CYAML_FIELD_SEQUENCE("Element_info", CYAML_FLAG_POINTER,

                         ElmTableYAML, infoItems, &element_info_schema,

                         0, CYAML_UNLIMITED),

    CYAML_FIELD_END};

static const cyaml_schema_value_t elm_table_schema = {CYAML_VALUE_MAPPING(

    CYAML_FLAG_POINTER, ElmTableYAML, elm_table_fields_schema)};

static const cyaml_config_t config = {
    .log_fn = cyaml_log,            /* Use the default logging function. */
    .mem_fn = cyaml_mem,            /* Use the default memory allocator. */
    .log_level = CYAML_LOG_WARNING, /* Logging errors and warnings only. */
};

ElmTableYAML *load_elmTableYAML(void)
{
    cyaml_err_t err;
    ElmTableYAML *elmTableYAML;
    err = cyaml_load_file("./src/database/element_table.yaml", &config,
                          &elm_table_schema, (void **)&elmTableYAML, NULL);
    if (err != CYAML_OK)
    {
        printf("ERROR: %s\n", cyaml_strerror(err));
        return NULL;
    }
    elmTableYAML->destroy = destroy_ElmTableYAML;
    return elmTableYAML;
}

void destroy_ElmTableYAML(ElmTableYAML **elmTableYAML)
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

HashNode *init_ElmInfoTable_internal()
{
    HashNode *hashTab = NULL;
    for (int i = 0; i < 47; ++i)
    {
        add_str_keyptr_item(&hashTab, elements[i].name, (void *)&elements[i]);
    }
    return hashTab;
}
