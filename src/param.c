#include "param.h"
#include "parser.h"

static char *exportFileName(Cell *self, const char *format,
                            const char *desiredName)
{
    char *exportDir = self->lattice->vtable->exportDir(
        self->lattice); /* Malloced String Created */
    int fileNameLen = 1 + snprintf(NULL, 0, format, exportDir, desiredName);
    char *ret = malloc(fileNameLen);
    snprintf(ret, fileNameLen, format, exportDir, desiredName);
    free(exportDir);
    return ret;
}

int getFinalCutoffEnergy(Cell *cell)
{
    HashNode *elmTable = cell->lookupTable;
    int energy = 0;
    for (int i = 0; i < cell->elmNums; ++i)
    {
        ElmInfo *item =
            find_item_by_str(elmTable, (const char *)cell->elmLists[i])->val;
        char *potPath = NULL;
        asprintf(&potPath, "./Potentials/%s", item->potFile);
        int fineEnergy = parse_fineCutoffEnergy(potPath);
        int ultraFineEnergy = roundupBiggerTenth(fineEnergy * 1.1);
        energy = (ultraFineEnergy > energy) ? ultraFineEnergy : energy;
        free(potPath);
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
    Atom *metalAtom = self->lattice->mol->vtable->get_atom_by_Id(
        self->lattice->mol, self->lattice->metal_site_id);
    ElmInfo *metal =
        find_item_by_str(self->lookupTable, (const char *)metalAtom->element)
            ->val;
    int spin = metal->spin;
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
    char *exportDir = self->lattice->vtable->exportDir(self->lattice);
    char *stemName = self->lattice->mol->name;
    char *geomParamName = exportFileName(self, "%s%s.param", stemName);
    FILE *geomParam = fopen(geomParamName, "w");
    if (!geomParam)
    {
        printf("Open file %s failed!\n", geomParamName);
    }
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
    char *dosParamName = exportFileName(self, "%s%s_DOS.param", stemName);
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
    char *exportDir = self->lattice->vtable->exportDir(self->lattice);
    char *stemName = self->lattice->mol->name;
    char *kptauxName = exportFileName(self, "%s%s.kptaux", stemName);
    FILE *kptauxFile = fopen(kptauxName, "w");
    free(kptauxName);
    fputs(kptauxText, kptauxFile);
    fclose(kptauxFile);
    char *kptauxDosName = exportFileName(self, "%s%s_DOS.kptaux", stemName);
    free(exportDir);
    FILE *kptauxDosFile = fopen(kptauxDosName, "w");
    free(kptauxDosName);
    fputs(kptauxText, kptauxDosFile);
    fclose(kptauxDosFile);
}

void write_trjaux(Cell *self)
{
    char *exportDir = self->lattice->vtable->exportDir(self->lattice);
    char *stemName = self->lattice->mol->name;
    char *trjauxName = exportFileName(self, "%s%s.trjaux", stemName);
    free(exportDir);
    FILE *trjauxFile = fopen(trjauxName, "w");
    free(trjauxName);
    char trjauxHeader[] =
        "# Atom IDs to appear in any .trj file to be generated.\n"
        "# Correspond to atom IDs which will be used in exported .msi file\n"
        "# required for animation/analysis of trajectory within Cerius2.\n";
    fputs(trjauxHeader, trjauxFile);
    int atomNum = self->lattice->mol->atomNum;
    Molecule *mPtr = self->lattice->mol;
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
    fclose(trjauxFile);
}

void copy_potentials(Cell *self)
{
    char *exportDir = self->lattice->vtable->exportDir(self->lattice);
    for (int i = 0; i < self->elmNums; ++i)
    {
        ElmInfo *elmInfo =
            find_item_by_str(self->lookupTable, (const char *)self->elmLists[i])
                ->val;
        char *potPath = exportFileName(self, "%s%s", elmInfo->potFile);
        struct stat s;
        if (stat(potPath, &s) == 0)
        {
            free(potPath);
            continue;
        }
        char *potSrc;
        asprintf(&potSrc, "./Potentials/%s", elmInfo->potFile);
        char *potContent = readWholeFile(potSrc);
        free(potSrc);
        FILE *copiedPotFile = fopen(potPath, "w");
        free(potPath);
        fputs(potContent, copiedPotFile);
        fclose(copiedPotFile);
        free(potContent);
    }
    free(exportDir);
}

void write_pbsScript(Cell *self)
{
    char *exportDir = self->lattice->vtable->exportDir(
        self->lattice); /* Malloced String Created */
    char *pbsScriptTemplate = readWholeFile(
        "./resource/pbs_template.sh"); /* Malloced String Created */
    char *pathName = exportFileName(self, "%s%s",
                                    "hpc.pbs.sh"); /* Malloced String Created */
    free(exportDir);                               /* Malloced String Freed */
    FILE *script = fopen(pathName, "w");           /* FILE opened */
    free(pathName);                                /* Malloced String Freed */
    fputs(pbsScriptTemplate, script);
    free(pbsScriptTemplate); /* Malloced String Freed */
    char cmdFormat[] = "mpirun --mca btl ^tcp --hostfile hostfile "
                       "/home/bhuang/castep.mpi %s\n";
    int cmdLen = 1 + snprintf(NULL, 0, cmdFormat, self->lattice->mol->name);
    char *cmdLine = malloc(cmdLen); /* Malloced String Created */
    snprintf(cmdLine, cmdLen, cmdFormat, self->lattice->mol->name);
    fputs(cmdLine, script);
    free(cmdLine); /* Malloced String Freed */
    fputs("rm ./hostfile", script);
    fclose(script); /* FILE closed */
}

void write_SMCastepExtension(Cell *self)
{
    char *fileName = exportFileName(self, "%sSMCastep_Extension_%s.xms",
                                    self->lattice->mol->name);
    char *content = readWholeFile("./resource/SMCastep_Extension.xms");
    FILE *ext = fopen(fileName, "w");
    free(fileName);
    fputs(content, ext);
    free(content);
    fclose(ext);
}
