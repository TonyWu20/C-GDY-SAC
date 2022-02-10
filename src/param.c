#include "param.h"
#include "parser.h"

int getFinalCutoffEnergy(Cell *cell)
{
    CastepInfo *table = cell->infoTab;
    int energy = 0;
    for (int i = 0; i < cell->elmNums; ++i)
    {
        CastepInfo *item = find_item(table, cell->elmLists[i]);
        int fineEnergy = parse_fineCutoffEnergy(item->info->potential_file);
        int ultraFineEnergy = roundupBiggerTenth(fineEnergy * 1.1);
        energy = (ultraFineEnergy > energy) ? ultraFineEnergy : energy;
    }
    return energy;
}

int parse_fineCutoffEnergy(const char *fileName)
{
    FILE *f = fopen(fileName, "r");
    int n = 4;
    char str[16];
    for (int i = 0; i < n; ++i)
    {
        fgets(str, 16, f);
    }
    char quality[4];
    int energy;
    sscanf(str, "%d %s", &energy, quality);
    fclose(f);
    return energy;
}

char *formatParam(char *template, int spin, double cutoff_energy)
{
    int len = 1 + snprintf(NULL, 0, template, spin, cutoff_energy);
    char *formatted = malloc(len);
    snprintf(formatted, len, template, spin, cutoff_energy);
    return formatted;
}

void write_param(Cell *self)
{
    double cutoff_energy = getFinalCutoffEnergy(self);
    CastepInfo *metal = find_item(self->infoTab, self->lattice->metal_symbol);
    int spin = metal->info->spin;
    char template_geom[] =
        "task : GeometryOptimization\n"
        "comment : CASTEP calculation from Materials Studio\n"
        "xc_functional : PBE\n"
        "spin_polarized : true\n"
        "spin : %8d\n"
        "opt_strategy : Speed\n"
        "page_wvfns : 0\n"
        "cut_off_energy : %18.15f\n"
        "grid_scale : 1.500000000000000\n"
        "fine_grid_scale : 1.500000000000000\n"
        "finite_basis_corr : 0\n"
        "elec_energy_tol : 1.000000000000000e-005\n"
        "max_scf_cycles : 6000\n"
        "fix_occupancy : false\n"
        "metals_method : dm\n"
        "mixing_scheme : Pulay\n"
        "mix_charge_amp : 0.500000000000000\n"
        "mix_spin_amp : 2.000000000000000\n"
        "mix_charge_gmax : 1.500000000000000\n"
        "mix_spin_gmax : 1.500000000000000\n"
        "mix_history_length : 20\n"
        "perc_extra_bands : 72\n"
        "smearing_width : 0.100000000000000\n"
        "spin_fix : 6\n"
        "num_dump_cycles : 0\n"
        "geom_energy_tol : 5.000000000000000e-005\n"
        "geom_force_tol : 0.100000000000000\n"
        "geom_stress_tol : 0.200000000000000\n"
        "geom_disp_tol : 0.005000000000000\n"
        "geom_max_iter : 6000\n"
        "geom_method : BFGS\n"
        "fixed_npw : false\n"
        "calculate_ELF : false\n"
        "calculate_stress : false\n"
        "popn_calculate : true\n"
        "calculate_hirshfeld : true\n"
        "calculate_densdiff : false\n"
        "popn_bond_cutoff : 3.000000000000000\n"
        "pdos_calculate_weights : true";
    char *newGeomParam = formatParam(template_geom, spin, cutoff_energy);
    char *exportDir = self->lattice->vtable->exportDir(self->lattice,
                                                       self->lattice->pathName);
    char *stemName = self->lattice->_mol->name;
    char *geomParamName;
    int fileNameLen = 1 + snprintf(NULL, 0, "%s%s.param", exportDir, stemName);
    geomParamName = malloc(fileNameLen);
    snprintf(geomParamName, fileNameLen, "%s%s.param", exportDir, stemName);
    FILE *geomParam = fopen(geomParamName, "w");
    fputs(newGeomParam, geomParam);
    fclose(geomParam);
    free(geomParamName);
    free(newGeomParam);
    char template_dos[] = "task : BandStructure\n"
                          "continuation : default\n"
                          "comment : CASTEP calculation from Materials Studio\n"
                          "xc_functional : PBE\n"
                          "spin_polarized : true\n"
                          "spin :        %d\n"
                          "opt_strategy : Speed\n"
                          "page_wvfns :        0\n"
                          "cut_off_energy :      %18.15f\n"
                          "grid_scale :        1.500000000000000\n"
                          "fine_grid_scale :        1.500000000000000\n"
                          "finite_basis_corr :        0\n"
                          "elec_energy_tol :   1.000000000000000e-005\n"
                          "max_scf_cycles :     6000\n"
                          "fix_occupancy : false\n"
                          "metals_method : dm\n"
                          "mixing_scheme : Pulay\n"
                          "mix_charge_amp :        0.500000000000000\n"
                          "mix_spin_amp :        2.000000000000000\n"
                          "mix_charge_gmax :        1.500000000000000\n"
                          "mix_spin_gmax :        1.500000000000000\n"
                          "mix_history_length :       20\n"
                          "perc_extra_bands :      72\n"
                          "smearing_width :        0.100000000000000\n"
                          "spin_fix :        6\n"
                          "num_dump_cycles : 0\n"
                          "bs_nextra_bands :       72\n"
                          "bs_xc_functional : PBE\n"
                          "bs_eigenvalue_tol :   1.000000000000000e-005\n"
                          "calculate_stress : false\n"
                          "calculate_ELF : false\n"
                          "popn_calculate : false\n"
                          "calculate_hirshfeld : false\n"
                          "calculate_densdiff : false\n"
                          "pdos_calculate_weights : true\n"
                          "bs_write_eigenvalues : true\n";
    char *newDosParam = formatParam(template_dos, spin, cutoff_energy);
    char *dosParamName;
    fileNameLen = 1 + snprintf(NULL, 0, "%s%s_DOS.param", exportDir, stemName);
    dosParamName = malloc(fileNameLen);
    snprintf(dosParamName, fileNameLen, "%s%s_DOS.param", exportDir, stemName);
    free(exportDir);
    FILE *dosParamFile = fopen(dosParamName, "w");
    fputs(newDosParam, dosParamFile);
    fclose(dosParamFile);
    free(dosParamName);
    free(newDosParam);
}

void write_kptaux(Cell *self)
{
    char kptauxText[] = "MP_GRID :        1       1       1\n"
                        "MP_OFFSET :   0.000000000000000e+000  "
                        "0.000000000000000e+000  0.000000000000000e+000\n"
                        "%BLOCK KPOINT_IMAGES\n"
                        "   1   1\n"
                        "%ENDBLOCK KPOINT_IMAGES";
    char *exportDir = self->lattice->vtable->exportDir(self->lattice,
                                                       self->lattice->pathName);
    char *stemName = self->lattice->_mol->name;
    char *kptauxName;
    int len = 1 + snprintf(NULL, 0, "%s%s.kptaux", exportDir, stemName);
    kptauxName = malloc(len);
    snprintf(kptauxName, len, "%s%s.kptaux", exportDir, stemName);
    FILE *kptauxFile = fopen(kptauxName, "w");
    free(kptauxName);
    fputs(kptauxText, kptauxFile);
    fclose(kptauxFile);
    len = 1 + snprintf(NULL, 0, "%s%s_DOS.kptaux", exportDir, stemName);
    char *kptauxDosName;
    kptauxDosName = malloc(len);
    snprintf(kptauxDosName, len, "%s%s_DOS.kptaux", exportDir, stemName);
    free(exportDir);
    FILE *kptauxDosFile = fopen(kptauxDosName, "w");
    free(kptauxDosName);
    fputs(kptauxText, kptauxDosFile);
    fclose(kptauxDosFile);
}

void write_trjaux(Cell *self)
{
    char *exportDir = self->lattice->vtable->exportDir(self->lattice,
                                                       self->lattice->pathName);
    char *stemName = self->lattice->_mol->name;
    char *trjauxName;
    int len = 1 + snprintf(NULL, 0, "%s%s.trjaux", exportDir, stemName);
    trjauxName = malloc(len);
    snprintf(trjauxName, len, "%s%s.trjaux", exportDir, stemName);
    free(exportDir);
    FILE *trjauxFile = fopen(trjauxName, "w");
    char trjauxHeader[] =
        "# Atom IDs to appear in any .trj file to be generated.\n"
        "# Correspond to atom IDs which will be used in exported .msi file\n"
        "# required for animation/analysis of trajectory within Cerius2.\n";
    fputs(trjauxHeader, trjauxFile);
    int atomNum = self->lattice->_mol->atomNum;
    Molecule *mPtr = self->lattice->_mol;
    for (int i = 0; i < atomNum; ++i)
    {
        Atom *curAtom = mPtr->atom_arr[i];
        int lineLen = 1 + snprintf(NULL, 0, "%d\n", curAtom->atomId);
        char *line = malloc(lineLen);
        snprintf(line, lineLen, "%d\n", curAtom->atomId);
        fputs(line, trjauxFile);
        free(line);
    }
    char trjauxEnding[] = "#Origin  0.000000000000000e+000  "
                          "0.000000000000000e+000  0.000000000000000e+000";
    fputs(trjauxEnding, trjauxFile);
    free(trjauxName);
    fclose(trjauxFile);
}

void copy_potentials(Cell *self, PotentialFile *table)
{
    char *exportDir = self->lattice->vtable->exportDir(self->lattice,
                                                       self->lattice->pathName);
    for (int i = 0; i < self->elmNums; ++i)
    {
        PotentialFile *potItem = find_PotItem(table, self->elmLists[i]);
        char *potential_stem = strrchr(potItem->potential_file, '/') + 1;
        int pathLen = 1 + snprintf(NULL, 0, "%s%s", exportDir, potential_stem);
        char *potPath = malloc(pathLen);
        snprintf(potPath, pathLen, "%s%s", exportDir, potential_stem);
        FILE *copiedPotFile = fopen(potPath, "w");
        free(potPath);
        fputs(potItem->fileContent, copiedPotFile);
        fclose(copiedPotFile);
    }
    free(exportDir);
}
