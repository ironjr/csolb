/**
 * solb-single.cpp
 *
 * Standalone version of SolB solenoid magnetic flux density calculator.
 *
 * Version 1.0 @ 10/26/2018
 *
 * Jaerin Lee
 * Applied Superconductivity Laboratory
 * Dept. of Electrical and Computer Engineering
 * Seoul National University
 */

#include <stdlib.h>
#include <argp.h>

#include "../solb.h"

const char *argp_program_version = "csolb 1.0";
const char *argp_program_bug_address = "<jarin.lee@gmail.com>";

static char doc[] =
    "cSolB -- Standalone calculator for magnetic flux density in solenoidal magnets.";

static char args_doc[] = "ARG1 ARG2";

static struct argp_option options[] = {
    { "verbose",        'v', 0,         0, "Produce verbose output" },
    { "interactive",    't', 0,         0, "Run in command line mode" },
    { "coil",           'c', "FILE",    0, "Coil data input" },
    { "radial",         'r', 0,         0, "Radial coordinate of the probe" },
    { "azimuthal",      'z', 0,         0, "Azimuthal coordinate of the probe" },
    { "probe",          'p', "FILE",    0, "File of list of probes" },
    { "output",         'o', "FILE",    0, "Output file of B field strength" },
    { 0 }
};

struct arguments
{
    char *args[2];
    int verbose;
    int interactive;
    double r;
    double z;
    char *coil_file;
    char *probe_file;
    char *output_file;
}

static error_t
parse_opt (int key, char *arg, struct argp_state *state)
{
    struct arguments *arguments = state->input;
    switch (key)
    {
        case 'v':
            arguments->verbose = 1;
            break;
        case 't':
            arguments->interactive = 1;
            break;
            case 

int main(int argc, char** argv)
{
    struct arguments arguements;

    arguments.verbose = 0;
    arguments.interactive = 0;
    arguments.r = 0;
    arguments.z = 0;
    arguments.coil_file = "-";
    arguments.probe_file = "-";
    arguments.output_file = "-";

    argp_parse(&argp, argc, argv, 0, 0, &arguments);



    return 0;
}
