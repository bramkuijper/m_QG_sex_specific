// Coevolution of maternal and paternal effects
// in silver spoon contexts
//
#include <ctime>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <string>
#include <cmath>
#include <cassert>
#include <random>



// set up the random number generator 
std::random_device rd;
unsigned seed = rd();
std::mt19937 rng_r(seed);

// generate a uniform distribution between 0 and 1
std::uniform_real_distribution<> uniform(0.0,1.0);

// generate a bernouilli distribution
// for meiotic segregation that spits out 
// either 0 or 1
std::bernoulli_distribution random_allele(0.5);


// population size
const int N = 5000;

// maximum number of generations
int max_generations = 50000;
int skip = 10;

// the baseline environment
double eps = 0.0;

// the selective optimum
double theta = 0.0;

double sigma_e = 1.0;

// environmental variables
double sinusoidal_frequency = 0.0;
double gauss_autocorr = 0.0;
double gauss_sd_eps = 0.0;

// effect of environment on theta
//
// see Lande (2009) J Evol Biol
double A = 0.0;
double B = 0.0;

// selection coefficients
double omega_y = 0.0;
double omega_hm = 0.0;
double omega_hp = 0.0;
double omega_m0 = 0.0;
double omega_m1 = 0.0;

// minimal survival probability
double s0_init = 0.0;

// which we put to 0 after a while
double s0 = 0.0;

// mutation rates
double mu_hm = 0.0;
double mu_hp = 0.0;
double mu_m0 = 0.0;
double mu_m1 = 0.0;

// mutational variances
double sdmu_hm = 0.0;
double sdmu_hp = 0.0;
double sdmu_m0 = 0.0;
double sdmu_m1 = 0.0;

// initial values
double init_hm = 0;
double init_hp = 0;
double init_m0 = 0;
double init_m1 = 0;

// maternal error in environmental maternal effect
double sd_mat_error = 0.0;

// stats
int Nf = 0;
int Nm = 0;
int Nfsurv = 0;
int Nmsurv = 0;

// empty file name
std::string file_name{};

// diploid individuals
struct Individual
{
    double hm[2]; // maternal phenotypic effect
    double hp[2]; // paternal phenotypic effect
    double m0[2]; // baseline phenotype
    double m1[2]; // maternal environmental effect
    double y; // an individual's phenotype

    double eps; // an individual's environment
};

typedef Individual Pop[N];

// Make subpopulations
Pop Females, Males, FemaleSurvivors, MaleSurvivors;

// initializes parameters using command line arguments
void init_arguments(int argc, char *argv[])
{
    init_hm = std::stod(argv[1]);
    init_hp = std::stod(argv[2]);
    init_m0 = std::stod(argv[3]);
    init_m1 = std::stod(argv[4]);

    mu_hm = std::stod(argv[5]);
    mu_hp = std::stod(argv[6]);
    mu_m0 = std::stod(argv[7]);
    mu_m1 = std::stod(argv[8]);
    
    sdmu_hm = std::stod(argv[9]);
    sdmu_hp = std::stod(argv[10]);
    sdmu_m0 = std::stod(argv[11]);
    sdmu_m1 = std::stod(argv[12]);

    omega_y = std::stod(argv[13]);
    omega_hm = std::stod(argv[14]);
    omega_hp = std::stod(argv[15]);
    omega_m0 = std::stod(argv[16]);
    omega_m1 = std::stod(argv[17]);

    sd_mat_error = std::stod(argv[18]); 

    sinusoidal_frequency = std::stod(argv[19]);
    gauss_autocorr = std::stod(argv[20]);
    gauss_sd_eps = std::stod(argv[21]);

    A = std::stod(argv[22]);
    B = std::stod(argv[23]);

    s0_init = s0 = std::stod(argv[24]);
    sigma_e = std::stod(argv[25]);

    file_name = argv[26];
}

// init_ialization function, runs at start of sim
// gives all the individuals some starting
// value
void init_pop()
{ 
    std::normal_distribution<double> developmental_noise_dist(0.0, sigma_e);

    // init_ialize the population
    for (size_t i = 0; i < N/2; ++i)
    {
        for (size_t allele_i = 0; allele_i < 2; ++allele_i)
        {
            Females[i].hm[allele_i] = init_hm;
            Females[i].hp[allele_i] = init_hp;
            Females[i].m0[allele_i] = init_m0;
            Females[i].m1[allele_i] = init_m1;
            
            Males[i].hm[allele_i] = init_hm;
            Males[i].hp[allele_i] = init_hp;
            Males[i].m0[allele_i] = init_m0;
            Males[i].m1[allele_i] = init_m1;
        }

        Females[i].y = developmental_noise_dist(rng_r);
        Males[i].y = developmental_noise_dist(rng_r);
        Females[i].eps = 0;
        Males[i].eps = 0;
    }

    Nf = Nm = N/2;
}


// mutate allelic values when trait value is bounded between 0 and 1
void mut01(double &val, double const mu, double const sdmu)
{
    // see whether we want t
    if (uniform(rng_r) < mu)
    {
        // initialize distribution of mutational effect sizes
        std::normal_distribution<double> mutational_effect_size(0.0, sdmu);

        val += mutational_effect_size(rng_r);

        val = val < 0 ? 0 : val > 1.0 ? 1.0 : val;
    } 
}

// mutate allelic values for traits that are not bounded
void mut(double &val, double const mu, double const sdmu)
{
    if (uniform(rng_r) < mu)
    {
        std::normal_distribution<double> mutational_effect_size(0.0, sdmu);
        val += mutational_effect_size(rng_r);
    } 
}

// make a kid out of parent asexually
// add it either to local offspring pool
// or to pool of dispersers
void create_offspring(
        Individual &mother, 
        Individual &father, 
        Individual &Kid)
{
    // inherit maternal effect loci
    Kid.hm[0] = mother.hm[random_allele(rng_r)];
	mut(Kid.hm[0], mu_hm, sdmu_hm);

    Kid.hm[1] = father.hm[random_allele(rng_r)];
	mut(Kid.hm[1], mu_hm, sdmu_hm);

    // inherit paternal effect loci
    Kid.hp[0] = mother.hp[random_allele(rng_r)];
	mut(Kid.hp[0], mu_hp, sdmu_hp);

    Kid.hp[1] = father.hp[random_allele(rng_r)];
	mut(Kid.hp[1], mu_hp, sdmu_hp);

    // elevation
    Kid.m0[0] = mother.m0[random_allele(rng_r)];
	mut(Kid.m0[0], mu_m0, sdmu_m0);

    Kid.m0[1] = father.m0[random_allele(rng_r)];
	mut(Kid.m0[1], mu_m0, sdmu_m0);
    
    // environmental maternal effect
    Kid.m1[0] = mother.m1[random_allele(rng_r)];
	mut(Kid.m1[0], mu_m1, sdmu_m1);

    Kid.m1[1] = father.m1[random_allele(rng_r)];
	mut(Kid.m1[1], mu_m1, sdmu_m1);


    // distributions to sample noise from,
    // either developmental or maternal
    std::normal_distribution<double> developmental_noise(0, sigma_e);
    std::normal_distribution<double> maternal_environmental_noise(0, sd_mat_error);

    // build child's phenotype
    Kid.y = 0.5 * (Kid.hm[0] + Kid.hm[1]) * mother.y // maternal phenotypic effect
            + 0.5 * (Kid.hp[0] + Kid.hp[1]) * father.y // paternal phenotypic effect
            + 0.5 * (Kid.m0[0] + Kid.m0[1]) // baseline phenotype
            + developmental_noise(rng_r) // developmental error
            + 0.5 * (Kid.m1[0] + Kid.m1[1]) * // maternal enviromental effect
            (mother.eps + maternal_environmental_noise(rng_r));

    assert(std::isnan(Kid.y) == 0);

    Kid.eps = eps;
}

// do the stats on the data
void write_data(int const generation, std::ofstream &data_file)
{
    double mean_y = 0.0;
    double mean_hm = 0.0;
    double mean_hp = 0.0;
    double mean_m0 = 0.0;
    double mean_m1 = 0.0;
    
    double ss_y = 0.0;
    double ss_hm = 0.0;
    double ss_hp = 0.0;
    double ss_m0 = 0.0;
    double ss_m1 = 0.0;

    // loop through patches and individuals
    // and get stats on genotypes and patch freqs.
    //
    double val;

    // average over females
    for (int i = 0; i < Nf; ++i)
    {
        val = Females[i].y;
        mean_y += val;
        ss_y += val * val;

        val = 0.5 * (Females[i].hm[0] + Females[i].hm[1]);
        mean_hm += val;
        ss_hm += val * val;

        val = 0.5 * (Females[i].hp[0] + Females[i].hp[1]);
        mean_hp += val;
        ss_hp += val * val;

        val = 0.5 * (Females[i].m0[0] + Females[i].m0[1]);
        mean_m0 += val;
        ss_m0 += val * val;
        
        val = 0.5 * (Females[i].m1[0] + Females[i].m1[1]);
        mean_m1 += val;
        ss_m1 += val * val;
    }
    
    // average over males
    for (int i = 0; i < Nm; ++i)
    {
        val = Males[i].y;
        mean_y += val;
        ss_y += val * val;

        val = 0.5 * (Males[i].hm[0] + Males[i].hm[1]);
        mean_hm += val;
        ss_hm += val * val;

        val = 0.5 * (Males[i].hp[0] + Males[i].hp[1]);
        mean_hp += val;
        ss_hp += val * val;

        val = 0.5 * (Males[i].m0[0] + Males[i].m0[1]);
        mean_m0 += val;
        ss_m0 += val * val;
        
        val = 0.5 * (Males[i].m1[0] + Males[i].m1[1]);
        mean_m1 += val;
        ss_m1 += val * val;
    }

    data_file << generation << ";"; 

    mean_y /= Nf + Nm;
    mean_hm /= Nf + Nm;
    mean_hp /= Nf + Nm;
    mean_m0 /= Nf + Nm;
    mean_m1 /= Nf + Nm;

    double var_y = ss_y / (Nf + Nm) - mean_y * mean_y;
    double var_hm = ss_hm / (Nf + Nm) - mean_hm * mean_hm;
    double var_hp = ss_hp / (Nf + Nm) - mean_hp * mean_hp;
    double var_m0 = ss_m0 / (Nf + Nm) - mean_m0 * mean_m0;
    double var_m1 = ss_m1 / (Nf + Nm) - mean_m1 * mean_m1;

    data_file 
        << eps << ";"
        << theta << ";"
        << mean_y << ";"
        << mean_hm << ";"
        << mean_hp << ";"
        << mean_m0 << ";"
        << mean_m1 << ";"
        << var_y << ";"
        << var_hm << ";"
        << var_hp << ";"
        << var_m0 << ";"
        << var_m1 << ";"
        << Nf << ";"
        << Nm << ";"
        << Nfsurv << ";"
        << Nmsurv << ";" << std::endl;
}

// headers of the datafile
void write_data_headers(std::ofstream &data_file)
{
    data_file 
        << "generation;"
        << "eps;"
        << "theta;"
        << "meanphen;"
        << "meanhm;"
        << "meanhp;"
        << "meanm0;"
        << "meanm1;"
        << "varphen;"
        << "varhm;"
        << "varhp;"
        << "varm0;"
        << "varm1;"
        << "Nf;"
        << "Nm;"
        << "Nfsurv;"
        << "Nmsurv;" << std::endl;
}

// parameters at the end of the sim
void write_parameters(std::ofstream &DataFile)
{
    DataFile << std::endl << std::endl 
        << "init_hm;" << init_hm << std::endl
        << "init_hp;" << init_hp << std::endl
        << "init_m0;" << init_m0 << std::endl
        << "init_m1;" << init_m1 << std::endl
        << "mu_hm;" << mu_hm << std::endl
        << "mu_hp;" << mu_hp << std::endl
        << "mu_m0;" << mu_m0 << std::endl
        << "mu_m1;" << mu_m1 << std::endl
        << "sdmu_hm;" << sdmu_hm << std::endl
        << "sdmu_hp;" << sdmu_hp << std::endl
        << "sdmu_m0;" << sdmu_m0 << std::endl
        << "sdmu_m1;" << sdmu_m1 << std::endl
        << "omega_y;" << omega_y << std::endl
        << "omega_hm;" << omega_hm << std::endl
        << "omega_hp;" << omega_hp << std::endl
        << "omega_m0;" << omega_m0 << std::endl
        << "omega_m1;" << omega_m1 << std::endl
        << "s0_init;" << s0_init << std::endl
        << "s0;" << s0 << std::endl
        << "sigma_e;" << sigma_e << std::endl
        << "A;" << A << std::endl
        << "B;" << B << std::endl
        << "sinusoidal_frequency;" << sinusoidal_frequency << std::endl
        << "gauss_autocorr;" << gauss_autocorr << std::endl
        << "gauss_sd_eps;" << gauss_sd_eps << std::endl
        << "seed;" << seed << std::endl;
}

// survival probability
double w(Individual &ind)
{
    double retval = s0 + (1.0 - s0) * 
            exp(-0.5 * ( 
                pow((ind.y - theta)/omega_y,2.0)
                + pow(0.5 *(ind.hm[0] + ind.hm[1])/omega_hm,2.0)
                + pow(0.5 *(ind.hp[0] + ind.hp[1])/omega_hp,2.0)
                + pow(0.5 *(ind.m1[0] + ind.m1[1])/omega_m1,2.0)));

    // bound values - when exponentials get really large all sorts of numerical misery happens
    if (std::isnormal(retval) == 0)
    {
        retval = 0.0;
    }

    return(retval);
}

void survival()
{
    Nfsurv = 0;
    Nmsurv = 0;
    // female survival
    for (int i = 0; i < Nf; ++i)
    {
        if (uniform(rng_r) < w(Females[i]))
        {
            FemaleSurvivors[Nfsurv++] = Females[i];
        }
    }

    // male survival
    for (int i = 0; i < Nm; ++i)
    {
        if (uniform(rng_r) < w(Males[i]))
        {
            MaleSurvivors[Nmsurv++] = Males[i];
        }
    }

    Nf = 0;
    Nm = 0;
}

void reproduction(int const generation, std::ofstream &data_file)
{
    int mother_id, father_id;

    if (Nfsurv < 1 || Nmsurv < 1)
    {
        write_parameters(data_file);
        exit(1);
    }

    // set up samplers to sample from the population of survivors
    std::uniform_int_distribution<int> female_survivor_sampler(0, Nfsurv - 1);
    std::uniform_int_distribution<int> male_survivor_sampler(0, Nmsurv - 1);

    // sample survivors and have them reproduce
    for (int i = 0; i < N; ++i)
    {
        mother_id = female_survivor_sampler(rng_r);
        father_id = male_survivor_sampler(rng_r);

        if (uniform(rng_r) < 0.5)
        {
            create_offspring(FemaleSurvivors[mother_id]
                                ,MaleSurvivors[father_id]
                                ,Females[Nf++]);
        }
        else
        {
            create_offspring(FemaleSurvivors[mother_id]
                                ,MaleSurvivors[father_id]
                                ,Males[Nm++]);
        }
    }

    assert(Nm + Nf <= N);
}

void change_environment(int const generation)
{
    std::normal_distribution<double> noise_in_envtal_signal(
            0.0, 
            sqrt(1.0 - gauss_autocorr * gauss_autocorr) * gauss_sd_eps);

    eps = sin(sinusoidal_frequency * generation) 
        + noise_in_envtal_signal(rng_r)
        + gauss_autocorr * eps;

    theta = A + B * eps;
}

// main function body
int main(int argc, char **argv)
{
    // initialize the arguments
    init_arguments(argc, argv);

    std::ofstream data_file(file_name.c_str());

    init_pop();

    write_data_headers(data_file);

    for (int generation = 0; generation < max_generations; ++generation)
    {
        change_environment(generation);

        survival();

        reproduction(generation, data_file);

        if (generation % skip == 0)
        {
            write_data(generation, data_file);
        }

//        if (generation == 0.5 * max_generations)
//        {
//            s0 = s0_final;
//        }
    }

    write_parameters(data_file);
}

