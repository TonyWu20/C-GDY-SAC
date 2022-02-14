#include "ads_database.h"
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
AdsorbateInfo ethylene[] = {
    {"CO", 1, {1, -1}, {1, 2}, {1, 2, 2}, true, false, 2},
    {"CHO", 1, {1, -1}, {2, 3}, {1, 2, 3}, true, false, 2},
    {"COCHO", 1, {3, -1}, {1, 3}, {1, 2, 3}, true, false, 5},
    {"COCHOH", 1, {3, -1}, {1, 3}, {1, 2, 3}, true, false, 6},
    {"OCH2CO", 2, {1, 3}, {1, 3}, {1, 2, 3}, true, false, 6},
    {"OCH2COH", 2, {1, 3}, {1, 3}, {1, 2, 3}, true, false, 6},
    {"OCH2CHOH", 2, {1, 3}, {1, 3}, {1, 2, 3}, true, false, 4},
    {"OCH2CH", 2, {1, 3}, {1, 3}, {1, 2, 3}, true, false, 4},
    {"OCH2CH2", 2, {1, 3}, {1, 3}, {1, 2, 3}, true, false, 4},
    {"C2H4", 2, {1, 2}, {1, 2}, {1, 2, 3}, false, true, 0}};

AdsorbateInfo acetic_acid[] = {
    {"OCH2C_cyc_OH", 1, {1, -1}, {2, 4}, {1, 2, 3}, true, false, 6},
    {"CH3COOH", 1, {2, -1}, {2, 3}, {2, 3, 4}, true, false, 7}};

AdsorbateInfo ethanol[] = {
    {"Glyoxal", 2, {1, 6}, {1, 6}, {1, 2, 3}, true, false, 5},
    {"HOCCHO", 2, {1, 6}, {1, 6}, {1, 2, 3}, true, true, 4},
    {"HOHCCHO", 2, {1, 6}, {1, 6}, {1, 2, 3}, true, false, 4},
    {"Glycolaldehyde", 2, {1, 6}, {1, 6}, {1, 2, 3}, true, false, 4},
    {"CH2CHO", 1, {5, -1}, {1, 2}, {1, 2, 5}, true, false, 3},
    {"CH2CHOH", 2, {1, 2}, {1, 2}, {1, 2, 3}, true, false, 4},
    {"CH2CH2OH", 1, {1, -1}, {1, 2}, {1, 2, 3}, true, false, 4},
    {"CH3CH2OH", 1, {1, -1}, {1, 2}, {1, 2, 3}, true, false, 9},
    {"acetaldehyde", 1, {1, -1}, {1, 2}, {1, 2, 6}, true, false, 6}};

AdsorbateInfo ethanol_other[] = {
    {"CH2OHCH2O", 2, {1, 6}, {1, 6}, {1, 2, 3}, true, false, 4},
    {"ethylene_glycol", 2, {1, 6}, {1, 6}, {1, 2, 3}, true, true, 4}};

HashNode *init_adsInfoTable()
{
    HashNode *hashTab = NULL;
    for (int i = 0; i < sizeof(ethylene) / sizeof(ethylene[0]); ++i)
    {
        add_str_keyptr_item(&hashTab, ethylene[i].name, (void *)&ethylene[i]);
    }
    for (int i = 0; i < sizeof(ethanol) / sizeof(ethanol[0]); ++i)
    {
        add_str_keyptr_item(&hashTab, ethanol[i].name, (void *)&ethanol[i]);
    }
    for (int i = 0; i < sizeof(ethanol_other) / sizeof(ethanol_other[0]); ++i)
    {
        add_str_keyptr_item(&hashTab, ethanol_other[i].name,
                            (void *)&ethanol_other[i]);
    }
    for (int i = 0; i < sizeof(acetic_acid) / sizeof(acetic_acid[0]); ++i)
    {
        add_str_keyptr_item(&hashTab, acetic_acid[i].name,
                            (void *)&acetic_acid[i]);
    }
    return hashTab;
}
