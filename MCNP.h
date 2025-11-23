/*
 * Program Name: MCNPHeavy
 * Program Release Year: 2025
 * Program Author: Steven S.
 * Program Link: https://github.com/BitEU/MCNPHeavy
 * Purpose: A Monte Carlo Neutron Transport program for k-effective criticality calculation with Windows Console/UNIVAC 1219 support
 */

#ifndef MCNP_H
#define MCNP_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

/* Platform-specific includes */
#ifndef UNIVAC
#include <windows.h>
#endif

/* Constants */
#define MAX_MATERIALS 10
#define MAX_REGIONS 20
#define MAX_GENERATIONS 100
#define MAX_PARTICLES 10000
#define PI 3.14159265358979323846

/* Multi-group energy structure */
#define NUM_ENERGY_GROUPS 3
#define ENERGY_GROUP_FAST 0      /* 0.1 - 20 MeV */
#define ENERGY_GROUP_EPITHERMAL 1 /* 1 eV - 0.1 MeV */
#define ENERGY_GROUP_THERMAL 2   /* 0.025 eV - 1 eV */

/* Energy boundaries (eV) */
#define ENERGY_FAST_MAX 20.0e6
#define ENERGY_FAST_MIN 0.1e6
#define ENERGY_EPI_MAX 0.1e6
#define ENERGY_EPI_MIN 1.0
#define ENERGY_THERMAL_MAX 1.0
#define ENERGY_THERMAL_MIN 0.025

/* Neutron speeds by group (cm/s) */
#define SPEED_FAST 1.38e9        /* ~14 MeV neutron */
#define SPEED_EPITHERMAL 4.4e7   /* ~1 keV neutron */
#define SPEED_THERMAL 2.2e5      /* ~0.025 eV neutron */

/* Material types */
#define MAT_VOID 0
#define MAT_U235 1               /* Fissile material (Uranium-235) */
#define MAT_LEAD 2               /* Lead shielding */
#define MAT_CONCRETE 3           /* Concrete shielding */
#define MAT_WATER 4              /* Water moderator */
#define MAT_GRAPHITE 5           /* Graphite moderator */

/* Interaction types */
#define INTERACT_ABSORB 0
#define INTERACT_SCATTER 1
#define INTERACT_FISSION 2
#define INTERACT_LEAK 3

/* Geometry types */
#define GEOM_SPHERE 0
#define GEOM_SLAB 1
#define GEOM_CYLINDER 2

/* Energy-dependent cross-section data (in barns, 1 barn = 1e-24 cm^2) */
typedef struct {
    double sigma_absorption;     /* Absorption cross-section (barns) */
    double sigma_scatter;        /* Scattering cross-section (barns) */
    double sigma_fission;        /* Fission cross-section (barns) */
    double sigma_total;          /* Total cross-section (barns) */
    double nu;                   /* Average neutrons per fission */
    double alpha;                /* Scattering anisotropy parameter (avg cosine) */
} CrossSectionData;

/* Scattering matrix for down-scattering between groups */
typedef struct {
    double prob[NUM_ENERGY_GROUPS][NUM_ENERGY_GROUPS]; /* Probability of scatter from group i to j */
} ScatterMatrix;

/* Material with energy-dependent nuclear data */
typedef struct {
    int type;                    /* Material type ID */
    char name[32];               /* Material name */
    double density;              /* g/cm^3 */
    double atomic_mass;          /* Atomic mass (amu) */
    double temperature;          /* Temperature (Kelvin) */
    CrossSectionData xs[NUM_ENERGY_GROUPS]; /* Cross-sections by energy group */
    ScatterMatrix scatter_matrix; /* Down-scattering probabilities */
} Material;

/* 3D Vector */
typedef struct {
    double x;
    double y;
    double z;
} Vector3;

/* Neutron particle */
typedef struct {
    Vector3 position;            /* Current position (cm) */
    Vector3 direction;           /* Direction unit vector */
    double energy;               /* Energy in eV */
    int energy_group;            /* Current energy group (0=fast, 1=epi, 2=thermal) */
    double weight;               /* Statistical weight (for variance reduction) */
    int alive;                   /* 1 = alive, 0 = absorbed/leaked */
    int generation;              /* Generation number for k-eff */
    int collisions;              /* Number of collisions */
} Neutron;

/* Geometric region */
typedef struct {
    int geometry_type;           /* GEOM_SPHERE, GEOM_SLAB, etc. */
    int material_id;             /* Material filling this region */
    Vector3 center;              /* Center point (cm) */
    double radius;               /* For sphere (cm) */
    double thickness;            /* For slab (cm) */
    double height;               /* For cylinder (cm) */
} Region;

/* Simulation configuration */
typedef struct {
    int num_particles;           /* Particles per generation */
    int num_generations;         /* Number of generations */
    int num_regions;             /* Number of geometric regions */
    Region regions[MAX_REGIONS]; /* Geometric regions */
    Material materials[MAX_MATERIALS]; /* Material database */
    int num_materials;           /* Number of materials defined */
    double critical_radius;      /* Critical sphere radius estimate */
} Config;

/* Simulation tallies/results */
typedef struct {
    int absorptions;             /* Number of absorptions */
    int scatters;                /* Number of scatters */
    int fissions;                /* Number of fissions */
    int leakages;                /* Number of leakages */
    int scatters_by_group[NUM_ENERGY_GROUPS]; /* Scatters per energy group */
    int fissions_by_group[NUM_ENERGY_GROUPS]; /* Fissions per energy group */
    double total_distance;       /* Total distance traveled */
    int total_collisions;        /* Total collisions */
    double k_eff;                /* k-effective (criticality) */
    double k_eff_array[MAX_GENERATIONS]; /* k-eff per generation */
    double shannon_entropy[MAX_GENERATIONS]; /* Spatial entropy per generation */
} Tallies;

/* Main simulation state */
typedef struct {
    Config config;
    Tallies tallies;
    unsigned long rng_seed;      /* Random number generator seed */
    int current_generation;
    double neutrons_this_gen;
    double neutrons_next_gen;
} SimState;

/* Function declarations */

/* Initialization */
void init_simulation(SimState* state);
void init_materials(SimState* state);
void init_geometry(SimState* state);
void init_tallies(Tallies* tallies);
void setup_test_problem(SimState* state);

/* Random number generation (Linear Congruential Generator for portability) */
void seed_rng(unsigned long seed);
double rand_uniform(void);
double rand_exponential(double lambda);
void rand_isotropic_direction(Vector3* dir);

/* Geometry and transport */
int find_region(const Vector3* pos, const Config* config);
double distance_to_boundary(const Neutron* n, const Config* config, int* next_region);
void move_neutron(Neutron* n, double distance);

/* Physics interactions */
int sample_interaction(const Material* mat, int energy_group);
void handle_absorption(Neutron* n, Tallies* tallies);
void handle_scatter(Neutron* n, const Material* mat, Tallies* tallies);
void handle_fission(Neutron* n, Tallies* tallies, SimState* state);

/* Energy physics */
int get_energy_group(double energy);
double sample_fission_energy(void);
void anisotropic_scatter(Vector3* dir, double alpha, double A);
int sample_scatter_group(int current_group, const ScatterMatrix* matrix);
double get_neutron_speed(int energy_group);

/* Monte Carlo tracking */
void track_neutron(Neutron* n, SimState* state);
void run_generation(SimState* state);
void run_simulation(SimState* state);

/* Analysis and output */
void calculate_k_effective(SimState* state);
void print_results(const SimState* state);
void print_material_info(const Material* mat);
void print_header(void);

/* Vector operations */
double vector_length(const Vector3* v);
void vector_normalize(Vector3* v);
double vector_dot(const Vector3* a, const Vector3* b);
void vector_add(Vector3* result, const Vector3* a, const Vector3* b);
void vector_scale(Vector3* result, const Vector3* v, double scalar);

/* Utility functions */
double cm_to_mfp(double distance, double sigma_total, double density);
void print_separator(void);

/* Platform-specific functions */
#ifndef UNIVAC
void console_setup(void);
#endif

/* Platform-specific string copy */
#ifdef UNIVAC
#define SAFE_STRCPY(dest, src, size) do { strncpy(dest, src, (size)-1); (dest)[(size)-1] = '\0'; } while(0)
#else
#define SAFE_STRCPY(dest, src, size) strcpy_s(dest, size, src)
#endif

#endif /* MCNP_H */
