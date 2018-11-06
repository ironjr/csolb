/**
 * solb_app.cpp
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

#include <vector>
#include <cstdlib>
#include <cstring>
#include <argp.h>
#include <omp.h>

#include "../core/solb.h"

#define BUF_SIZE 200 // Buffer size in bytes for parsing of input files
#define CHUNK_SIZE (1<<8) // Number of point probes to be fully evaluated before written on the disk


const char *argp_program_version = "csolb 1.0";
const char *argp_program_bug_address = "<jarin.lee@gmail.com>";

static char doc[] =
    "cSolB -- Standalone calculator for magnetic flux density in solenoidal magnets.";

static struct argp_option options[] = {
    { "verbose",        'v', 0,         0, "Produce verbose output" },
    { "interactive",    't', 0,         0, "Run in command line mode" },
    { "coil",           'c', "FILE",    0, "Coil data input" },
    { "probe",          'p', "FILE",    0, "File of list of probes" },
    { "output",         'o', "FILE",    0, "Output file of B field strength" },
    { 0 }
};

struct arguments
{
    int verbose;
    int interactive;
    char *coil_file;
    char *probe_file;
    char *output_file;
};

static error_t
parse_opt(int key, char *arg, struct argp_state *state)
{
    struct arguments *arguments = (struct arguments *)state->input;
    switch (key)
    {
        case 'v':
            arguments->verbose = 1;
            break;
        case 't':
            arguments->interactive = 1;
            break;
        case 'c':
            arguments->coil_file = arg;
            break;
        case 'p':
            arguments->probe_file = arg;
            break;
        case 'o':
            arguments->output_file = arg;
            break;
        default:
            return ARGP_ERR_UNKNOWN;
    }
    return 0;
}

/* Parser */
static struct argp argp = {
    options,
    parse_opt,
    0,
    doc
};

int parse_coil(FILE *, std::vector<top_solenoid_t *> *);
int parse_single_coil(char *, top_solenoid_t *);
int parse_probe(FILE *, std::vector<vec2d_t> *);
void run_interactive(struct arguments *);

int
main(int argc, char** argv)
{
    struct arguments arguments;

    /* Default option values */
    arguments.verbose = 0;
    arguments.interactive = 0;
    arguments.coil_file = NULL;
    arguments.probe_file = NULL;
    arguments.output_file = NULL;

    argp_parse(&argp, argc, argv, 0, 0, &arguments);

    /* Change to interactive mode if input files are insufficient */
    if (arguments.coil_file == NULL || arguments.probe_file == NULL)
    {
        fprintf(stderr, "No input files are detected\nChange to interactive mode");
        arguments.interactive = 1;
    }

    /* Interactive mode */
    if (arguments.interactive)
    {
        run_interactive(&arguments);
        return 0;
    }

    /* Open the necessary file for auto-calculation */
    FILE *c_fp = NULL;
    if ((c_fp = fopen(arguments.coil_file, "rt")) == NULL)
    {
        fprintf(stderr, "%s: No such file or directory", arguments.coil_file);
        exit(0);
    }

    FILE *p_fp = NULL;
    if ((p_fp = fopen(arguments.probe_file, "rt")) == NULL)
    {
        fprintf(stderr, "%s: No such file or directory", arguments.probe_file);
        exit(0);
    }

    /* Output file is optional, it is stdout by default */
    FILE *o_fp = stdout;
    if (arguments.output_file != NULL)
    {
        if ((o_fp = fopen(arguments.output_file, "wt")) == NULL)
        {
            fprintf(stderr, "%s: No such file or directory", arguments.output_file);
            exit(0);
        }
    }

    /* Get coil configuration from the file */
    std::vector<top_solenoid_t *> coils;
    parse_coil(c_fp, &coils);
    size_t ncoil = coils.size();
    if (ncoil == 0)
    {
        fprintf(stderr, "%s: No coils are properly specified", arguments.coil_file);
        exit(0);
    }
    fclose(c_fp);

    std::vector<vec2d_t> probes;
    parse_probe(p_fp, &probes);
    size_t nprobe = probes.size();
    if (nprobe == 0)
    {
        fprintf(stderr, "%s: No probes are properly specified", arguments.probe_file);
        exit(0);
    }
    fclose(p_fp);

    /* Run the main program
     * We assume that typically there are more probes than coils
     * RAM may handle 100s of millions of results, though we need to consult with the size of
     * cache memory to optimize the performance
     */
    mag_field_2d_t res_chunk[CHUNK_SIZE]; 
    size_t rem = nprobe % CHUNK_SIZE;
    size_t full_chunks = nprobe / CHUNK_SIZE;
    fprintf(o_fp, "%16s%16s%16s%16s\n", "Coord_R", "Coord_Z", "Br", "Bz");
    if (full_chunks == 0)
    {
        for (size_t i = 0; i < nprobe; ++i)
        {
            res_chunk[i].Br = 0;
            res_chunk[i].Bz = 0;
        }
#pragma omp parallel for collapse(2)
        for (size_t i = 0; i < nprobe; ++i)
        {
            for (size_t j = 0; j < ncoil; ++j)
            {
                mag_field_2d_t res = solb_single(coils[j], probes[i].r, probes[i].z);
                res_chunk[i].Br += res.Br;
                res_chunk[i].Bz += res.Bz;
            }
        }
        //fprintf(stderr, "INFO: Evaluation finished; begin writing to the result in file");
        for (size_t i = 0; i < nprobe; ++i)
        {
            fprintf(o_fp, "%16lf%16lf%16lf%16lf\n", probes[i].r, probes[i].z,
                    res_chunk[i].Br, res_chunk[i].Bz);
        }
    }
    else
    {
        for (size_t chunk = 0; chunk < full_chunks; ++chunk)
        {
            for (size_t i = 0; i < CHUNK_SIZE; ++i)
            {
                res_chunk[i].Br = 0;
                res_chunk[i].Bz = 0;
            }
#pragma omp parallel for collapse(2)
            for (size_t i = 0; i < CHUNK_SIZE; ++i)
            {
                for (size_t j = 0; j < ncoil; ++j)
                {
                    mag_field_2d_t res = solb_single(coils[j],
                            probes[i + CHUNK_SIZE * chunk].r,
                            probes[i + CHUNK_SIZE * chunk].z);
                    res_chunk[i].Br += res.Br;
                    res_chunk[i].Bz += res.Bz;
                }
            }
            for (size_t i = 0; i < CHUNK_SIZE; ++i)
            {
                fprintf(o_fp, "%16lf%16lf%16lf%16lf\n",
                        probes[i + CHUNK_SIZE * chunk].r,
                        probes[i + CHUNK_SIZE * chunk].z,
                        res_chunk[i].Br, res_chunk[i].Bz);
            }
        }
        for (size_t i = 0; i < rem; ++i)
        {
            res_chunk[i].Br = 0;
            res_chunk[i].Bz = 0;
        }
#pragma omp parallel for collapse(2)
        for (size_t i = 0; i < rem; ++i)
        {
            for (size_t j = 0; j < ncoil; ++j)
            {
                mag_field_2d_t res = solb_single(coils[j],
                        probes[i + CHUNK_SIZE * full_chunks].r,
                        probes[i + CHUNK_SIZE * full_chunks].z);
                res_chunk[i].Br += res.Br;
                res_chunk[i].Bz += res.Bz;
            }
        }
        for (size_t i = 0; i < rem; ++i)
        {
            fprintf(o_fp, "%16lf%16lf%16lf%16lf\n",
                    probes[i + CHUNK_SIZE * full_chunks].r,
                    probes[i + CHUNK_SIZE * full_chunks].z,
                    res_chunk[i].Br, res_chunk[i].Bz);
        }
    }


    return 0;
}

/*
 * parse_coil
 * Parse coils from the coil configuration input file.
 * returns 1 on success, 0 on failure
 */
int
parse_coil(FILE *coil_fp, std::vector<top_solenoid_t *> *coils)
{
    top_solenoid_t sol;
    char buf[BUF_SIZE];
    while (fgets(buf, BUF_SIZE, coil_fp) != NULL)
    {
        if (!parse_single_coil(buf, &sol))
            return 0;
        top_solenoid_t *item = new top_solenoid_t(&sol);
        std::vector<top_solenoid_t *>::iterator it = coils->end();
        coils->insert(it, item);
    }
    return 1;
}

/*
 * parse_single_coil
 * Parse a single line of coil configuration file into solenoid data.
 * returns 1 on success, 0 on failure
 */
int
parse_single_coil(char *str, top_solenoid_t *sol)
{
    static const char* label = "parse_single_coil";

    char *tok;
    double base;
    char de;
    int exponent;

    int idx = 0;
    tok = strtok(str, " \t\r\n");
    while (tok != NULL)
    {
        /* Parse fortran/matlab like floats as well as c-type floats */
        int res = sscanf(tok, "%lf%c%d", &base, &de, &exponent);
        switch (res)
        {
            case 0:
                fprintf(stderr, "%s: Wrong floating point expression in the file, try %%lf[d|D|e|E]%%d", label);
                return 0;
            case 1:
                exponent = 3;
                break;
            case 2:
                exponent = 3;
            case 3:
                if (de != 'e' && de != 'E' && de != 'd' && de != 'D')
                {
                    fprintf(stderr, "%s: Wrong floating point expression in the file, try %%lf[d|D|e|E]%%d", label);
                    return 0;
                }
                break;
            default:
                fprintf(stderr, "%s: Unknown error from sscanf", label);
                return 0;
        }
        switch (idx)
        {
            case 0:
            case 1:
            case 2:
            case 3:
                /* Coil dimension in mm; should be translated into m */
                if (exponent != 3)
                    base *= pow(10, exponent - 3);
                break;
            case 4:
                /* Current density is given in A/mm^2; should be translated into SI units */
                if (exponent != -6)
                    base *= pow(10, exponent + 6);
                break;
            default:
                fprintf(stderr, "%s: Wrong number of coil parameters", label);
                return 0;
        }
        switch (idx)
        {
            case 0:
                sol->a1 = base;
                break;
            case 1:
                sol->a2 = base;
                break;
            case 2:
                sol->b1 = base;
                break;
            case 3:
                sol->b2 = base;
                break;
            case 4:
                sol->j = base;
                break;
        }
        tok = strtok(NULL, " \t\r\n");
        ++idx;
    }
    return 1;
}

/*
 * parse_probe
 * Parse probes from the probe configuration input file.
 * returns 1 on success, 0 on failure
 */
int
parse_probe(FILE *probe_fp, std::vector<vec2d_t> *probes)
{
    vec2d_t probe;
    char buf[BUF_SIZE];
    while (fgets(buf, BUF_SIZE, probe_fp) != NULL)
    {
        if (sscanf(buf, "%lf%lf", &probe.r, &probe.z) != 2)
            return 0;
        vec2d_t *item = new vec2d_t(&probe);
        std::vector<vec2d_t>::iterator it = probes->end();
        probes->insert(it, item);
    }
    return 1;
}

/*
 * run_interactive
 * Run in interactive mode.
 * TODO
 */
void
run_interactive(struct arguments *arguments)
{
    /* Initial setup */
    FILE *c_fp = NULL;
    FILE *p_fp = NULL;
    if (arguments->coil_file != NULL)
    {
        if ((c_fp = fopen(arguments->coil_file, "rt")) == NULL)
        {
            fprintf(stderr, "%s: No such file or directory", arguments->coil_file);
            return;
        }
    }
    if (arguments->probe_file != NULL)
    {
        if ((p_fp = fopen(arguments->probe_file, "rt")) == NULL)
        {
            fprintf(stderr, "%s: No such file or directory", arguments->probe_file);
            return;
        }
    }

    /* In the interactive mode, output file saves the entire log, though it is
     * optional */
    FILE *o_fp = stdout;
    if (arguments->output_file != NULL)
    {
        if ((o_fp = fopen(arguments->output_file, "wt")) == NULL)
        {
            fprintf(stderr, "%s: No such file or directory", arguments->output_file);
            return;
        }
    }

    /* Buffers and other stuffs for interactive calculation */
    //char buf[BUF_SIZE];

    /* Interactive loop */
    // TODO Use bison/yacc
    while (1)
    {
        break;

    }
    return;
}
