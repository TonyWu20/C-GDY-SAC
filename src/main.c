#include "assemble.h"
#include "atom.h"
#include "cell.h"
#include "misc.h"
#include "molecule.h"
#include "my_maths.h"
#include "parser.h"
#include "tasks.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define PI (atan(1) * 4)
Lattice *load_lat(char *fileName, char *name)
{
    Lattice *lat = parse_lattice_from_file(fileName, name);
    return lat;
}

void test_tasks()
{
    char *elements[] = {"Sc", "Ti", "V",  "Cr", "Mn", "Fe", "Co", "Ni", "Cu",
                        "Zn", "Y",  "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd",
                        "Ag", "Cd", "Hf", "Ta", "W",  "Re", "Os", "Ir", "Pt",
                        "Au", "Hg", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu",
                        "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu"};
    AdsTableYAML *adsTableYAML = load_adsTableYAML();
    ElmTableYAML *elmTableYAML = load_elmTableYAML();
    int prog = 0;
    generator(elmTableYAML, elements, adsTableYAML, &prog);
    elmTableYAML->destroy(&elmTableYAML);
    adsTableYAML->destroy(&adsTableYAML);
}

int main(int argc, char *argv[])
{
    test_tasks();
    return 0;
}
