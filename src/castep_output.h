#include "atom.h"
#include "lattice.h"
#include "molecule.h"
#include "uthash.h"
#include <stdio.h>

typedef struct Cell Cell;

typedef struct
{
    const char *name;
    struct ElmItem *info;
    UT_hash_handle hh;
} CastepInfo;

struct Cell
{
    Lattice *lattice;
    struct ElmItem *info;
    char *(*write_block)(char *block_name, char *(*block_text)(Cell *self));
};

struct ElmItem
{
    const char *elm;
    int LCAO;
    double mass;
    const char *potential_file;
    int spin;
};

struct ElmItem ads[] = {{"C", 2, 12.0109996796, "./Potentials/C_00PBE.usp", 0},
                        {"H", 1, 1.0080000162, "./Potentials/H_00PBE.usp", 0},
                        {"O", 2, 15.9989995956, "./Potentials/O_00PBE.usp", 0}};

struct ElmItem metal3D[] = {
    {"Co", 3, 58.9329986572, "./Potentials/Co_00PBE.uspcc", 3},
    {"Cr", 5, 51.9959983826, "./Potentials/Cr_00PBE.usp", 3},
    {"Cu", 3, 63.5460014343, "./Potentials/Cu_00PBE.usp", 1},
    {"Fe", 3, 55.8470001221, "./Potentials/Fe_00PBE.uspcc", 4},
    {"Mn", 3, 54.9379997253, "./Potentials/Mn_00PBE.uspcc", 5},
    {"Ni", 3, 58.7099990845, "./Potentials/Ni_00PBE.uspcc", 2},
    {"Sc", 5, 44.9560012817, "./Potentials/Sc_00PBE.usp", 0},
    {"Ti", 5, 47.9000015259, "./Potentials/Ti_00PBE.usp", 0},
    {"V", 5, 50.9410018921, "./Potentials/V_00PBE.usp", 2},
    {"Zn", 4, 65.3799972534, "./Potentials/Zn_00PBE.usp", 0}};

struct ElmItem metal4D[] = {
    {"Ag", 3, 107.8679962158, "./Potentials/Ag_00PBE.usp", 0},
    {"Cd", 4, 112.4100036621, "./Potentials/Cd_00PBE.usp", 0},
    {"Mo", 5, 95.9400024414, "./Potentials/Mo_00PBE.usp", 2},
    {"Nb", 5, 92.90599823, "./Potentials/Nb_00PBE.usp", 1},
    {"Pd", 3, 106.4000015259, "./Potentials/Pd_00PBE.usp", 2},
    {"Rh", 3, 102.90599823, "./Potentials/Rh_00PBE.usp", 4},
    {"Ru", 5, 101.0699996948, "./Potentials/Ru_00PBE.usp", 5},
    {"Tc", 5, 98.90599823, "./Potentials/Tc_00PBE.usp", 3},
    {"Y", 3, 88.90599823, "./Potentials/Y_00PBE.uspcc", 0},
    {"Zr", 5, 91.2200012207, "./Potentials/Zr_00PBE.usp", 0}};

struct ElmItem metal5D[] = {
    {"Au", 3, 196.966003418, "./Potentials/Au_00PBE.usp", 2},
    {"Hf", 3, 178.4900054932, "./Potentials/Hf_00PBE.uspcc", 0},
    {"Hg", 4, 200.5899963379, "./Potentials/Hg_00PBE.usp", 0},
    {"Ir", 3, 192.2200012207, "./Potentials/Ir_00PBE.usp", 4},
    {"Os", 5, 190.1999969482, "./Potentials/Os_00PBE.usp", 5},
    {"Pt", 3, 195.0899963379, "./Potentials/Pt_00PBE.usp", 4},
    {"Re", 5, 186.2070007324, "./Potentials/Re_00PBE.usp", 3},
    {"Ta", 3, 180.9479980469, "./Potentials/Ta_00PBE.usp", 1},
    {"W", 5, 183.8500061035, "./Potentials/W_00PBE.usp", 2},
};

struct ElmItem metalLM[] = {
    {"Ce", 6, 140.1199951172, "./Potentials/Ce_00PBE.usp", 1},
    {"Dy", 6, 162.5, "./Potentials/Dy_00.usp", 5},
    {"Er", 6, 167.2599945068, "./Potentials/Er_00.usp", 3},
    {"Eu", 6, 151.9600067139, "./Potentials/Eu_00.usp", 6},
    {"Gd", 6, 157.25, "./Potentials/Gd_00.usp", 7},
    {"Ho", 6, 164.9299926758, "./Potentials/Ho_00PBE.usp", 4},
    {"La", 6, 138.9049987793, "./Potentials/La_00PBE.usp", 0},
    {"Lu", 4, 174.9700012207, "./Potentials/Lu_00.usp", 0},
    {"Nd", 6, 144.2400054932, "./Potentials/Nd_00.usp", 3},
    {"Pm", 6, 147.0, "./Potentials/Pm_00.usp", 4},
    {"Pr", 6, 140.9080047607, "./Potentials/Pr_00.usp", 2},
    {"Sm", 6, 150.3999938965, "./Potentials/Sm_00.usp", 5},
    {"Tb", 6, 158.9250030518, "./Potentials/Tb_00.usp", 6},
    {"Tm", 6, 168.9340057373, "./Potentials/Tm_00.usp", 2},
    {"Yb", 6, 173.0399932861, "./Potentials/Yb_00PBE.usp", 1}};

/* Add CastepInfo item to hash table
 */
void add_item(CastepInfo **table, struct ElmItem *item);

CastepInfo *find_item(CastepInfo *table, const char *elm);

void delete_all(CastepInfo **table);

CastepInfo *initTable();
