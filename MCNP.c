/*
 * Monte Carlo Neutron Transport Program
 * k-effective Criticality Calculator - Main Implementation
 * 
 * Compatible with UNIVAC 1219 and Windows conhost
 */

#include "MCNP.h"

/* Global RNG state (BSS initialization to 0) */
static unsigned long rng_state = 0;

/* Material database with energy-dependent cross-sections */
/* Data approximates ENDF/B-VIII.0 library values */
static const Material material_library[] = {
    /* Void */
    {MAT_VOID, "VOID", 0.0, 0.0, 293.0,
        {{0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
         {0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
         {0.0, 0.0, 0.0, 0.0, 0.0, 0.0}},
        {{{0.0}}}},
    
    /* Uranium-235 (energy-dependent, fissile) */
    {MAT_U235, "U-235", 19.1, 235.0, 293.0,
        /* Fast: high speed, low absorption, high fission */
        {{0.09, 5.1, 1.30, 6.49, 2.50, 0.94},   /* 2 MeV typical */
         /* Epithermal: resonance region */
         {1.5, 8.3, 2.0, 11.8, 2.45, 0.85},     /* 1 keV typical */
         /* Thermal: huge cross-sections */
         {99.0, 10.0, 585.0, 694.0, 2.43, 0.0}}, /* 0.025 eV */
        /* Scattering matrix (fast->epi->thermal down-scattering) */
        {{{0.88, 0.12, 0.0},
          {0.0, 0.82, 0.18},
          {0.0, 0.0, 1.0}}}},
    
    /* Lead (Pb) - good reflector, inelastic scattering */
    {MAT_LEAD, "LEAD", 11.34, 207.2, 293.0,
        {{0.003, 6.8, 0.0, 6.803, 0.0, 0.92},   /* Fast */
         {0.09, 11.2, 0.0, 11.29, 0.0, 0.75},   /* Epithermal */
         {0.17, 11.0, 0.0, 11.17, 0.0, 0.0}},   /* Thermal */
        {{{0.75, 0.25, 0.0},
          {0.0, 0.85, 0.15},
          {0.0, 0.0, 1.0}}}},
    
    /* Concrete (SiO2 + Ca composite) */
    {MAT_CONCRETE, "CONCRETE", 2.3, 40.0, 293.0,
        {{0.008, 2.1, 0.0, 2.108, 0.0, 0.70},
         {0.018, 1.4, 0.0, 1.418, 0.0, 0.50},
         {0.022, 0.8, 0.0, 0.822, 0.0, 0.0}},
        {{{0.70, 0.30, 0.0},
          {0.0, 0.75, 0.25},
          {0.0, 0.0, 1.0}}}},
    
    /* Water (H2O) - excellent moderator */
    {MAT_WATER, "WATER", 1.0, 18.0, 293.0,
        /* Water slows neutrons rapidly */
        {{0.019, 8.5, 0.0, 8.519, 0.0, 0.67},   /* Fast */
         {0.103, 28.0, 0.0, 28.103, 0.0, 0.40}, /* Epithermal */
         {0.664, 49.0, 0.0, 49.664, 0.0, 0.0}}, /* Thermal */
        /* Strong down-scattering due to hydrogen */
        {{{0.55, 0.40, 0.05},
          {0.0, 0.60, 0.40},
          {0.0, 0.0, 1.0}}}},
    
    /* Graphite (Carbon-12) - classic moderator */
    {MAT_GRAPHITE, "GRAPHITE", 1.6, 12.0, 293.0,
        {{0.0003, 4.6, 0.0, 4.6003, 0.0, 0.95},
         {0.0020, 4.7, 0.0, 4.702, 0.0, 0.85},
         {0.0034, 4.74, 0.0, 4.7434, 0.0, 0.0}},
        {{{0.85, 0.15, 0.0},
          {0.0, 0.88, 0.12},
          {0.0, 0.0, 1.0}}}}
};

static const int num_material_library = sizeof(material_library) / sizeof(Material);

/* ========================================================================== */
/* MAIN ENTRY POINT */
/* ========================================================================== */

int main(void) {
#ifndef UNIVAC
    console_setup();
#endif
    
    SimState state;
    
    /* Initialize with zeros (BSS initialization) */
    memset(&state, 0, sizeof(SimState));
    
    print_header();
    init_simulation(&state);
    run_simulation(&state);
    print_results(&state);
    
    return 0;
}

/* ========================================================================== */
/* INITIALIZATION FUNCTIONS */
/* ========================================================================== */

void init_simulation(SimState* state) {
    printf("\nInitializing Monte Carlo Neutron Transport Simulation...\n");
    
    /* Seed RNG with current time */
    seed_rng((unsigned long)time(NULL));
    state->rng_seed = rng_state;
    
    init_materials(state);
    init_geometry(state);
    init_tallies(&state->tallies);
    setup_test_problem(state);
    
    printf("Configuration complete.\n");
    print_separator();
}

void init_materials(SimState* state) {
    /* Copy material library to simulation state */
    state->config.num_materials = num_material_library;
    
    for (int i = 0; i < num_material_library && i < MAX_MATERIALS; i++) {
        memcpy(&state->config.materials[i], &material_library[i], sizeof(Material));
    }
    
    printf("\nMaterial Database Loaded:\n");
    for (int i = 0; i < state->config.num_materials; i++) {
        print_material_info(&state->config.materials[i]);
    }
}

void init_geometry(SimState* state) {
    /* Initialize regions array to zeros */
    memset(state->config.regions, 0, sizeof(state->config.regions));
    state->config.num_regions = 0;
}

void init_tallies(Tallies* tallies) {
    /* Zero initialization (BSS style) */
    memset(tallies, 0, sizeof(Tallies));
}

void setup_test_problem(SimState* state) {
    printf("\n=== CRITICALITY TEST PROBLEM ===\n");
    printf("Problem: U-235 sphere surrounded by water reflector\n");
    printf("Goal: Calculate k-effective to determine criticality\n");
    printf("Physics: Multi-group energy transport with anisotropic scattering\n\n");
    
    /* Set simulation parameters */
    state->config.num_particles = 1000;
    state->config.num_generations = 50;
    
    /* Define geometry: nested spheres */
    state->config.num_regions = 3;
    
    /* Region 0: Inner U-235 sphere (fissile core) */
    state->config.regions[0].geometry_type = GEOM_SPHERE;
    state->config.regions[0].material_id = MAT_U235;
    state->config.regions[0].center.x = 0.0;
    state->config.regions[0].center.y = 0.0;
    state->config.regions[0].center.z = 0.0;
    state->config.regions[0].radius = 8.7;  /* cm - near-critical for bare U-235 */
    state->config.critical_radius = state->config.regions[0].radius;
    
    /* Region 1: Water reflector shell */
    state->config.regions[1].geometry_type = GEOM_SPHERE;
    state->config.regions[1].material_id = MAT_WATER;
    state->config.regions[1].center.x = 0.0;
    state->config.regions[1].center.y = 0.0;
    state->config.regions[1].center.z = 0.0;
    state->config.regions[1].radius = 18.7;  /* 10 cm water reflector */
    
    /* Region 2: Void (outside) */
    state->config.regions[2].geometry_type = GEOM_SPHERE;
    state->config.regions[2].material_id = MAT_VOID;
    state->config.regions[2].center.x = 0.0;
    state->config.regions[2].center.y = 0.0;
    state->config.regions[2].center.z = 0.0;
    state->config.regions[2].radius = 1000.0;  /* Effective infinity */
    
    printf("Geometry Configuration:\n");
    printf("  Core: U-235 sphere, radius = %.2f cm\n", state->config.regions[0].radius);
    printf("  Reflector: Water shell, outer radius = %.2f cm\n", state->config.regions[1].radius);
    printf("  Outside: Void (leakage)\n\n");
    
    printf("Simulation Parameters:\n");
    printf("  Neutrons per generation: %d\n", state->config.num_particles);
    printf("  Number of generations: %d\n", state->config.num_generations);
    printf("  Energy groups: %d (Fast, Epithermal, Thermal)\n", NUM_ENERGY_GROUPS);
    printf("  Total neutron histories: %d\n", 
           state->config.num_particles * state->config.num_generations);
}

/* ========================================================================== */
/* RANDOM NUMBER GENERATION */
/* ========================================================================== */

void seed_rng(unsigned long seed) {
    /* Linear Congruential Generator parameters (Numerical Recipes) */
    /* Ensure seed is not zero */
    rng_state = (seed != 0) ? seed : 1;
}

double rand_uniform(void) {
    /* LCG: X(n+1) = (a * X(n) + c) mod m */
    /* Using parameters: a=1664525, c=1013904223, m=2^32 */
    rng_state = (1664525UL * rng_state + 1013904223UL) & 0xFFFFFFFFUL;
    return (double)rng_state / 4294967296.0;  /* Normalize to [0,1) */
}

double rand_exponential(double lambda) {
    /* Sample from exponential distribution using inversion method */
    double u = rand_uniform();
    /* Avoid log(0) */
    if (u < 1e-10) u = 1e-10;
    return -log(u) / lambda;
}

void rand_isotropic_direction(Vector3* dir) {
    /* Sample isotropic direction using Marsaglia method */
    double u1, u2, s;
    
    do {
        u1 = 2.0 * rand_uniform() - 1.0;
        u2 = 2.0 * rand_uniform() - 1.0;
        s = u1 * u1 + u2 * u2;
    } while (s >= 1.0);
    
    double factor = 2.0 * sqrt(1.0 - s);
    dir->x = u1 * factor;
    dir->y = u2 * factor;
    dir->z = 1.0 - 2.0 * s;
}

/* ========================================================================== */
/* GEOMETRY FUNCTIONS */
/* ========================================================================== */

int find_region(const Vector3* pos, const Config* config) {
    /* Find which region contains the given position */
    /* Check regions from innermost to outermost */
    
    for (int i = 0; i < config->num_regions; i++) {
        const Region* reg = &config->regions[i];
        
        if (reg->geometry_type == GEOM_SPHERE) {
            double dx = pos->x - reg->center.x;
            double dy = pos->y - reg->center.y;
            double dz = pos->z - reg->center.z;
            double r = sqrt(dx*dx + dy*dy + dz*dz);
            
            if (r <= reg->radius) {
                return i;
            }
        }
    }
    
    /* Should not reach here if geometry is properly defined */
    return config->num_regions - 1;  /* Return outermost region */
}

double distance_to_boundary(const Neutron* n, const Config* config, int* next_region) {
    /* Calculate distance to nearest boundary along neutron direction */
    /* For spherical geometry only (simplified) */
    
    int current_region = find_region(&n->position, config);
    const Region* reg = &config->regions[current_region];
    
    if (reg->geometry_type == GEOM_SPHERE) {
        /* Ray-sphere intersection */
        Vector3 oc;
        oc.x = n->position.x - reg->center.x;
        oc.y = n->position.y - reg->center.y;
        oc.z = n->position.z - reg->center.z;
        
        double b = 2.0 * vector_dot(&oc, &n->direction);
        double c = vector_dot(&oc, &oc) - reg->radius * reg->radius;
        double discriminant = b*b - 4.0*c;
        
        if (discriminant < 0) {
            /* Moving away from center, check outer sphere */
            if (current_region + 1 < config->num_regions) {
                *next_region = current_region + 1;
                return distance_to_boundary(n, config, next_region);
            }
            return 1e10;  /* Very large distance */
        }
        
        double t1 = (-b - sqrt(discriminant)) / 2.0;
        double t2 = (-b + sqrt(discriminant)) / 2.0;
        
        /* Use positive root */
        double dist = (t1 > 1e-6) ? t1 : t2;
        
        if (dist > 1e-6) {
            *next_region = (current_region + 1 < config->num_regions) ? 
                           current_region + 1 : current_region;
            return dist;
        }
    }
    
    *next_region = current_region;
    return 1e10;  /* Very large distance */
}

void move_neutron(Neutron* n, double distance) {
    n->position.x += n->direction.x * distance;
    n->position.y += n->direction.y * distance;
    n->position.z += n->direction.z * distance;
}

/* ========================================================================== */
/* ENERGY PHYSICS */
/* ========================================================================== */

int get_energy_group(double energy) {
    /* Determine energy group from neutron energy (eV) */
    if (energy >= ENERGY_FAST_MIN) {
        return ENERGY_GROUP_FAST;
    } else if (energy >= ENERGY_EPI_MIN) {
        return ENERGY_GROUP_EPITHERMAL;
    } else {
        return ENERGY_GROUP_THERMAL;
    }
}

double get_neutron_speed(int energy_group) {
    /* Return characteristic speed for energy group */
    switch (energy_group) {
        case ENERGY_GROUP_FAST:
            return SPEED_FAST;
        case ENERGY_GROUP_EPITHERMAL:
            return SPEED_EPITHERMAL;
        case ENERGY_GROUP_THERMAL:
            return SPEED_THERMAL;
        default:
            return SPEED_THERMAL;
    }
}

double sample_fission_energy(void) {
    /* Sample from Watt fission spectrum (approximation of Maxwell-Boltzmann) */
    /* Watt spectrum: W(E) = C * exp(-E/a) * sinh(sqrt(b*E)) */
    /* For U-235: a = 0.988 MeV, b = 2.249 MeV^-1 */
    /* Simplified: use empirical distribution */
    
    double xi = rand_uniform();
    
    /* Approximate cumulative distribution */
    if (xi < 0.15) {
        /* Fast neutrons (2-20 MeV) */
        return 2.0e6 + xi / 0.15 * 18.0e6;
    } else if (xi < 0.85) {
        /* Most probable (0.1-2 MeV) */
        return 0.1e6 + (xi - 0.15) / 0.70 * 1.9e6;
    } else {
        /* Tail (< 0.1 MeV) */
        return 1000.0 + (xi - 0.85) / 0.15 * 99000.0;
    }
}

int sample_scatter_group(int current_group, const ScatterMatrix* matrix) {
    /* Sample target energy group after scattering */
    double xi = rand_uniform();
    double cumulative = 0.0;
    
    for (int j = 0; j < NUM_ENERGY_GROUPS; j++) {
        cumulative += matrix->prob[current_group][j];
        if (xi < cumulative) {
            return j;
        }
    }
    
    return current_group; /* Shouldn't reach here */
}

void anisotropic_scatter(Vector3* dir, double alpha, double A) {
    /* Anisotropic scattering in center-of-mass frame */
    /* alpha: average cosine of scattering angle */
    /* A: atomic mass ratio */
    
    double mu_cm;  /* Cosine of scattering angle in CM frame */
    
    if (fabs(alpha) < 0.01) {
        /* Isotropic scattering */
        mu_cm = 2.0 * rand_uniform() - 1.0;
    } else {
        /* Forward-peaked scattering (Henyey-Greenstein-like) */
        double xi = rand_uniform();
        if (fabs(alpha) > 0.99) {
            mu_cm = alpha > 0 ? 1.0 : -1.0;
        } else {
            double g2 = alpha * alpha;
            mu_cm = (1.0 + g2 - pow((1.0 - g2) / (1.0 + alpha * (2.0 * xi - 1.0)), 2.0)) / (2.0 * alpha);
        }
    }
    
    /* Transform to lab frame (for A > 1, elastic scattering) */
    double cos_lab;
    if (A > 1.0) {
        /* Lab frame transformation */
        cos_lab = (1.0 + A * mu_cm) / sqrt(1.0 + A * A + 2.0 * A * mu_cm);
    } else {
        cos_lab = mu_cm;
    }
    
    /* Clamp to valid range */
    if (cos_lab > 1.0) cos_lab = 1.0;
    if (cos_lab < -1.0) cos_lab = -1.0;
    
    /* Sample azimuthal angle (uniform) */
    double phi = 2.0 * PI * rand_uniform();
    double sin_theta = sqrt(1.0 - cos_lab * cos_lab);
    
    /* Rotate direction vector */
    Vector3 old_dir = *dir;
    
    /* Create local coordinate system */
    Vector3 perp1, perp2;
    if (fabs(old_dir.z) < 0.9) {
        perp1.x = -old_dir.y;
        perp1.y = old_dir.x;
        perp1.z = 0.0;
    } else {
        perp1.x = 1.0;
        perp1.y = 0.0;
        perp1.z = -old_dir.x;
    }
    vector_normalize(&perp1);
    
    /* Cross product for second perpendicular */
    perp2.x = old_dir.y * perp1.z - old_dir.z * perp1.y;
    perp2.y = old_dir.z * perp1.x - old_dir.x * perp1.z;
    perp2.z = old_dir.x * perp1.y - old_dir.y * perp1.x;
    
    /* New direction */
    dir->x = cos_lab * old_dir.x + sin_theta * (cos(phi) * perp1.x + sin(phi) * perp2.x);
    dir->y = cos_lab * old_dir.y + sin_theta * (cos(phi) * perp1.y + sin(phi) * perp2.y);
    dir->z = cos_lab * old_dir.z + sin_theta * (cos(phi) * perp1.z + sin(phi) * perp2.z);
    
    vector_normalize(dir);
}

/* ========================================================================== */
/* PHYSICS INTERACTIONS */
/* ========================================================================== */

int sample_interaction(const Material* mat, int energy_group) {
    /* Sample interaction type based on energy-dependent cross-sections */
    const CrossSectionData* xs = &mat->xs[energy_group];
    double total = xs->sigma_total;
    double xi = rand_uniform() * total;
    
    if (xi < xs->sigma_absorption) {
        return INTERACT_ABSORB;
    } else if (xi < xs->sigma_absorption + xs->sigma_fission) {
        return INTERACT_FISSION;
    } else {
        return INTERACT_SCATTER;
    }
}

void handle_absorption(Neutron* n, Tallies* tallies) {
    n->alive = 0;
    tallies->absorptions++;
}

void handle_scatter(Neutron* n, const Material* mat, Tallies* tallies) {
    /* Energy-dependent anisotropic scattering */
    const CrossSectionData* xs = &mat->xs[n->energy_group];
    
    /* Anisotropic scattering based on material and energy */
    anisotropic_scatter(&n->direction, xs->alpha, mat->atomic_mass);
    
    /* Sample new energy group (down-scattering) */
    int new_group = sample_scatter_group(n->energy_group, &mat->scatter_matrix);
    n->energy_group = new_group;
    
    /* Update energy based on new group */
    switch (new_group) {
        case ENERGY_GROUP_FAST:
            n->energy = ENERGY_FAST_MIN + rand_uniform() * (ENERGY_FAST_MAX - ENERGY_FAST_MIN);
            break;
        case ENERGY_GROUP_EPITHERMAL:
            n->energy = ENERGY_EPI_MIN * pow(ENERGY_EPI_MAX / ENERGY_EPI_MIN, rand_uniform());
            break;
        case ENERGY_GROUP_THERMAL:
            n->energy = ENERGY_THERMAL_MIN + rand_uniform() * (ENERGY_THERMAL_MAX - ENERGY_THERMAL_MIN);
            break;
    }
    
    tallies->scatters++;
    tallies->scatters_by_group[new_group]++;
    n->collisions++;
}

void handle_fission(Neutron* n, Tallies* tallies, SimState* state) {
    /* Sample number of neutrons from fission (energy-dependent) */
    const Material* mat = &state->config.materials[MAT_U235];
    const CrossSectionData* xs = &mat->xs[n->energy_group];
    
    /* Track average neutrons produced (nu varies slightly with energy) */
    state->neutrons_next_gen += xs->nu;
    
    tallies->fissions++;
    tallies->fissions_by_group[n->energy_group]++;
    n->alive = 0;  /* Original neutron absorbed in fission */
}

/* ========================================================================== */
/* MONTE CARLO TRACKING */
/* ========================================================================== */

void track_neutron(Neutron* n, SimState* state) {
    const int MAX_COLLISIONS = 1000;
    
    while (n->alive && n->collisions < MAX_COLLISIONS) {
        /* Find current region and material */
        int region_id = find_region(&n->position, &state->config);
        const Region* region = &state->config.regions[region_id];
        int mat_id = region->material_id;
        
        /* Check for leakage (void material) */
        if (mat_id == MAT_VOID) {
            n->alive = 0;
            state->tallies.leakages++;
            return;
        }
        
        const Material* mat = &state->config.materials[mat_id];
        const CrossSectionData* xs = &mat->xs[n->energy_group];
        
        /* Sample distance to collision using energy-dependent cross-section */
        double mfp = 1.0 / (xs->sigma_total * mat->density * 6.022e-1);  /* mean free path */
        double dist_collision = rand_exponential(1.0 / mfp);
        
        /* Check distance to boundary */
        int next_region = region_id;
        double dist_boundary = distance_to_boundary(n, &state->config, &next_region);
        
        if (dist_collision < dist_boundary) {
            /* Collision occurs before boundary */
            move_neutron(n, dist_collision);
            state->tallies.total_distance += dist_collision;
            state->tallies.total_collisions++;
            
            /* Sample interaction type (energy-dependent) */
            int interaction = sample_interaction(mat, n->energy_group);
            
            switch (interaction) {
                case INTERACT_ABSORB:
                    handle_absorption(n, &state->tallies);
                    break;
                    
                case INTERACT_SCATTER:
                    handle_scatter(n, mat, &state->tallies);
                    break;
                    
                case INTERACT_FISSION:
                    handle_fission(n, &state->tallies, state);
                    break;
            }
        } else {
            /* Neutron crosses boundary */
            move_neutron(n, dist_boundary + 1e-6);  /* Small step past boundary */
            state->tallies.total_distance += dist_boundary;
        }
    }
}

void run_generation(SimState* state) {
    state->neutrons_next_gen = 0.0;
    state->neutrons_this_gen = (double)state->config.num_particles;
    
    /* Track all neutrons in this generation */
    for (int i = 0; i < state->config.num_particles; i++) {
        Neutron n;
        
        /* Initialize neutron at center of fissile core */
        n.position.x = 0.0;
        n.position.y = 0.0;
        n.position.z = 0.0;
        
        /* Random isotropic direction */
        rand_isotropic_direction(&n.direction);
        
        /* Sample energy from fission spectrum (realistic birth energy) */
        n.energy = sample_fission_energy();
        n.energy_group = get_energy_group(n.energy);
        n.weight = 1.0;
        n.alive = 1;
        n.generation = state->current_generation;
        n.collisions = 0;
        
        /* Track this neutron */
        track_neutron(&n, state);
    }
}

void run_simulation(SimState* state) {
    printf("\n=== STARTING SIMULATION ===\n\n");
    printf("Gen    k-eff    Fissions  Absorptions  Scatters  Leakages  Fast%%  Epi%%  Therm%%\n");
    printf("---  --------  ---------  -----------  --------  --------  -----  ----  ------\n");
    
    for (int gen = 0; gen < state->config.num_generations; gen++) {
        state->current_generation = gen;
        
        /* Save tallies from previous generation */
        Tallies gen_tallies;
        memset(&gen_tallies, 0, sizeof(Tallies));
        
        /* Run this generation */
        run_generation(state);
        
        /* Calculate k-effective for this generation */
        double k_gen = state->neutrons_next_gen / state->neutrons_this_gen;
        state->tallies.k_eff_array[gen] = k_gen;
        
        /* Calculate fission fractions by energy group */
        int total_fissions = state->tallies.fissions_by_group[0] + 
                            state->tallies.fissions_by_group[1] + 
                            state->tallies.fissions_by_group[2];
        double fast_frac = total_fissions > 0 ? 
            100.0 * state->tallies.fissions_by_group[ENERGY_GROUP_FAST] / total_fissions : 0.0;
        double epi_frac = total_fissions > 0 ? 
            100.0 * state->tallies.fissions_by_group[ENERGY_GROUP_EPITHERMAL] / total_fissions : 0.0;
        double therm_frac = total_fissions > 0 ? 
            100.0 * state->tallies.fissions_by_group[ENERGY_GROUP_THERMAL] / total_fissions : 0.0;
        
        /* Print generation results with energy distribution */
        printf("%3d  %8.5f  %9d  %11d  %8d  %8d  %5.1f  %4.1f  %6.1f\n",
               gen + 1,
               k_gen,
               state->tallies.fissions,
               state->tallies.absorptions,
               state->tallies.scatters,
               state->tallies.leakages,
               fast_frac,
               epi_frac,
               therm_frac);
        
        /* Reset tallies for next generation (keep cumulative k-eff) */
        double saved_k_eff[MAX_GENERATIONS];
        memcpy(saved_k_eff, state->tallies.k_eff_array, sizeof(saved_k_eff));
        init_tallies(&state->tallies);
        memcpy(state->tallies.k_eff_array, saved_k_eff, sizeof(saved_k_eff));
    }
    
    /* Calculate final k-effective */
    calculate_k_effective(state);
}

/* ========================================================================== */
/* ANALYSIS AND OUTPUT */
/* ========================================================================== */

void calculate_k_effective(SimState* state) {
    /* Skip first few generations (equilibration) */
    int skip = state->config.num_generations / 10;
    if (skip < 5) skip = 5;
    
    double sum = 0.0;
    double sum_sq = 0.0;
    int count = 0;
    
    for (int i = skip; i < state->config.num_generations; i++) {
        double k = state->tallies.k_eff_array[i];
        sum += k;
        sum_sq += k * k;
        count++;
    }
    
    double mean = sum / count;
    double variance = (sum_sq / count) - (mean * mean);
    double std_dev = sqrt(variance);
    double std_error = std_dev / sqrt((double)count);
    
    /* 95% confidence interval (assuming normal distribution) */
    double confidence_95 = 1.96 * std_error;
    
    state->tallies.k_eff = mean;
    
    printf("\n");
    print_separator();
    printf("\n=== FINAL RESULTS ===\n\n");
    printf("k-effective Statistics:\n");
    printf("  Mean:               %.5f\n", mean);
    printf("  Standard Deviation: %.5f\n", std_dev);
    printf("  Standard Error:     %.5f\n", std_error);
    printf("  95%% Confidence:     %.5f +/- %.5f\n\n", mean, confidence_95);
    printf("  Generations used: %d (skipped first %d for equilibration)\n\n", count, skip);
    
    /* Interpret criticality */
    printf("CRITICALITY ASSESSMENT:\n");
    if (mean < 0.98) {
        printf("  System is SUBCRITICAL (k-eff < 1)\n");
        printf("  Neutron population will decrease exponentially\n");
        printf("  Configuration is stable and safe\n");
        printf("  Power level: DECREASING\n");
    } else if (mean > 1.02) {
        printf("  System is SUPERCRITICAL (k-eff > 1)\n");
        printf("  Neutron population will increase exponentially\n");
        printf("  WARNING: Unsafe configuration - prompt criticality risk!\n");
        printf("  Power level: INCREASING (factor of %.2f per generation)\n", mean);
    } else {
        printf("  System is CRITICAL (k-eff ~= 1)\n");
        printf("  Neutron population is self-sustaining\n");
        printf("  Configuration at critical threshold - reactor steady state\n");
        printf("  Power level: STABLE\n");
    }
    
    /* Physics commentary */
    printf("\n");
    printf("PHYSICS NOTES:\n");
    printf("  - Multi-group energy transport with 3 groups (fast, epithermal, thermal)\n");
    printf("  - Anisotropic scattering with lab-frame transformations\n");
    printf("  - Watt fission spectrum for neutron birth energies\n");
    printf("  - Energy-dependent cross-sections approximating ENDF/B-VIII.0\n");
    printf("  - Down-scattering between energy groups (moderation)\n");
}

void print_results(const SimState* state) {
    printf("\n");
    print_separator();
    printf("\nSimulation Summary:\n");
    printf("  Random seed: %lu\n", state->rng_seed);
    printf("  Core radius: %.2f cm\n", state->config.critical_radius);
    printf("  k-effective: %.5f\n", state->tallies.k_eff);
    printf("\n");
}

void print_material_info(const Material* mat) {
    printf("  %s: rho=%.2f g/cm^3, A=%.1f, T=%.0fK\n",
           mat->name, mat->density, mat->atomic_mass, mat->temperature);
    if (mat->xs[ENERGY_GROUP_THERMAL].nu > 0) {
        printf("    Thermal: sigma_f=%.1f b, sigma_a=%.1f b, nu=%.2f\n",
               mat->xs[ENERGY_GROUP_THERMAL].sigma_fission,
               mat->xs[ENERGY_GROUP_THERMAL].sigma_absorption,
               mat->xs[ENERGY_GROUP_THERMAL].nu);
        printf("    Fast:    sigma_f=%.2f b, sigma_a=%.2f b, nu=%.2f\n",
               mat->xs[ENERGY_GROUP_FAST].sigma_fission,
               mat->xs[ENERGY_GROUP_FAST].sigma_absorption,
               mat->xs[ENERGY_GROUP_FAST].nu);
    }
}

void print_header(void) {
    print_separator();
    printf("  MONTE CARLO NEUTRON TRANSPORT PROGRAM\n");
    printf("  k-effective Criticality Calculator\n");
    print_separator();
}

void print_separator(void) {
    printf("======================================================\n");
}

/* ========================================================================== */
/* VECTOR OPERATIONS */
/* ========================================================================== */

double vector_length(const Vector3* v) {
    return sqrt(v->x * v->x + v->y * v->y + v->z * v->z);
}

void vector_normalize(Vector3* v) {
    double len = vector_length(v);
    if (len > 1e-10) {
        v->x /= len;
        v->y /= len;
        v->z /= len;
    }
}

double vector_dot(const Vector3* a, const Vector3* b) {
    return a->x * b->x + a->y * b->y + a->z * b->z;
}

void vector_add(Vector3* result, const Vector3* a, const Vector3* b) {
    result->x = a->x + b->x;
    result->y = a->y + b->y;
    result->z = a->z + b->z;
}

void vector_scale(Vector3* result, const Vector3* v, double scalar) {
    result->x = v->x * scalar;
    result->y = v->y * scalar;
    result->z = v->z * scalar;
}

/* ========================================================================== */
/* PLATFORM-SPECIFIC FUNCTIONS */
/* ========================================================================== */

#ifndef UNIVAC
void console_setup(void) {
    /* Enable UTF-8 output on Windows */
    SetConsoleOutputCP(CP_UTF8);
    
    /* Set console title */
    SetConsoleTitleA("Monte Carlo Neutron Transport - k-eff Calculator");
}
#endif
