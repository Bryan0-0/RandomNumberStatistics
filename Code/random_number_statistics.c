 /******************************************************************************
 File: random_number_statistics.c
 Purpose:
 =========
This program generates random numbers using uniform and normal 
distributions,writes the generated sequences to text files, and calculates the 
sample mean and standard deviation for each sequence. It also generates 
histograms for the generated sequences.
 Features:
 =========
 1.**Uniform Integer Distribution**: Generates random integers within a specified 
range [m, M], where all integers have equal probability.
 2.**Uniform Real Distribution**: Generates random real numbers uniformly 
distributed in [m, M].
 3.**Normal Distribution (Real)**: Generates random real numbers from a normal 
distribution with specified mean and standard deviation.
 4.**Normal Distribution (Integer)**: Generates random integer numbers from a 
normal distribution with specified mean and standard deviation.
 5.**Truncated Normal Distribution (Real)**: Generates random real numbers from a 
normal distribution truncated to lie within [m, M].
 6.**Truncated Normal Distribution (Integer)**: Generates random integers from a 
normal distribution truncated to lie within [m, M].
 7.To compute the sample mean and standard deviation of the generated sequences.
 8.Generates histograms for all the 6 distributions in Scenario 3.
 Scenarios:
 ==========
 1.Scenario 1: mu=5, sigma=1, m=1, M=8, N=20
 2.Scenario 2: mu=2^10, sigma=2^8, m=1, M=2000, N=200,000
 3.Scenario 3: mu=2^12, sigma=1.3*(2^10), m=1, M=8100, N=2,000,000
 Output:
 ==========- Separate .txt files are created for each generator in the respective subfolders 
(Scenario1, Scenario2, Scenario3) under the DATA directory.- The program calculates the sample mean and standard deviation for all the six 
generators of the three given scenarios.- Histograms are generated for all the generators in Scenario 3 and written to 
the HISTOGRAM folder.
 ********************************************************************************
 */
 #include <stdio.h>    // For standard I/O operations
 #include <math.h>     // For mathematical functions
 #include <stdlib.h>   // For memory allocation and exit() function
 #include <time.h>     // For seeding the random number generator  
 #include <string.h>   // For string manipulation
 #include <sys/stat.h> // For mkdir()
 #include <errno.h>    // For error handling
 

 // Macros for generating random numbers
 #define frand() (rand() / (double)RAND_MAX) // Uniform random number in [0, 1)
 #define nrand() (sqrt(-2 * log(frand())) * cos(2 * M_PI * frand())) // Normal random number
 #ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif
 // Number of bins for histograms
 #define HISTOGRAM_BINS 50
// Function prototypes
 void generate_random_numbers_to_file(const char* filename, int type, double m, 
double M, double mu, double sigma, int N);
 void calculate_statistics_from_file(const char* filename, int N);
 void generate_histogram_from_file(const char* input_filename, const char* 
output_filename, int bins, double min, double max);
 int generate_uniform_integer(double m, double M);
 double generate_uniform_real(double m, double M);
 double generate_normal_real(double mu, double sigma);
 int generate_normal_integer(double mu, double sigma);
 double generate_truncated_normal_real(double m, double M, double mu, double 
sigma);
 int generate_truncated_normal_integer(double m, double M, double mu, double 
sigma);
 double calculate_mean(double* data, int n);
 double calculate_std_dev(double* data, int n, double mean);
 int create_directory(const char* path);
 int main() {
    srand(time(NULL)); // Seed the random number generator
    // Define scenarios: {mu, sigma, m, M, N}
    double scenarios[3][5] = {
        {5, 1, 1, 8, 20},
        {pow(2, 10), pow(2, 8), 1, 2000, 200000},
        {pow(2, 12), 1.3 * pow(2, 10), 1, 8100, 2000000}
    };
    // Paths for directories and files
    char* subfolders[3] = {"DATA/Scenario1", "DATA/Scenario2", "DATA/Scenario3"};
    char* histogram_folder = "Histogram";
    // Directory creation (e.g., DATA, HISTOGRAM, subfolders)
    // Ensure to use `create_directory` function and validate each directory creation.
    // Loop through each scenario to process random numbers
     if (create_directory("DATA") != 0) {
        fprintf(stderr, "Failed to create base DATA directory. Exiting.\n");
        return -1;
    }

    if (create_directory(histogram_folder) != 0) {
        fprintf(stderr, "Failed to create HISTOGRAM directory. Exiting.\n");
        return -1;
    }

    // File names for random number types
    char* file_names[6] = {
        "uniform_integers.txt",
        "uniform_reals.txt",
        "normal_integers.txt",
        "normal_reals.txt",
        "truncated_normal_integers.txt",
        "truncated_normal_reals.txt"
    };

    // Process each scenario
    for (int i = 0; i < 3; i++) {
        // Create subfolder for the scenario
        if (create_directory(subfolders[i]) != 0) {
            fprintf(stderr, "Failed to create subfolder %s. Exiting.\n", subfolders[i]);
            return -1;
        }
    
    for (int o = 0; o < 3; o++) {
        double mu = scenarios[o][0];
        double sigma = scenarios[o][1];
        double m = scenarios[o][2];
        double M = scenarios[o][3];
        int N = (int)scenarios[o][4];
        // Process random number generation, statistics, and histograms
        // Utilize appropriate function calls (detailed below) for each step
        for (int type = 1; type <= 6; type++) {
            char filename[256];
            snprintf(filename, sizeof(filename), "%s/%s", subfolders[i], file_names[type - 1]);

            // Generate random numbers and save to file
            printf("Generating random numbers for %s...\n", filename);
            generate_random_numbers_to_file(filename, type, m, M, mu, sigma, N);

            // Calculate statistics for the generated file
            printf("Calculating statistics for %s...\n", filename);
            calculate_statistics_from_file(filename, N);

            // For Scenario 3, generate histograms
            if (i == 2) {
                char histogram_filename[256];
                snprintf(histogram_filename, sizeof(histogram_filename), "%s/hist_%s", histogram_folder, file_names[type - 1]);

                printf("Generating histogram for %s...\n", filename);
                generate_histogram_from_file(filename, histogram_filename, HISTOGRAM_BINS, m, M);
            }
    }
    }
        }
    

    printf("All random number sequences generated, statistics calculated, and histograms created.\n");
    return 0;
    
 }
 /*
 Function: create_directory
 ==========================
 Creates a directory if it does not already exist.
 Parameters:
    path - Path of the directory to create
 Returns:
    0 if the directory is successfully created or already exists,
    -1 if an error occurs.
 */
 int create_directory(const char* path) 
{ 
    struct stat st;
    int mkdir(const char *pathname, mode_t mode);

    // Check if the directory already exists
    if (stat(path, &st) == 0 && S_ISDIR(st.st_mode)) {
        return 0; // Directory already exists
    }

    // Try to create the directory
    if (mkdir(path, 0755) == 0) {
        return 0; // Successfully created the directory
    } else {
        // An error occurred while creating the directory
        perror("Error creating directory");
        return -1;
    }
}
   
/*
 Function: generate_random_numbers_to_file
 =========================================
 Generates random numbers based on the 6 generator types (e.g., uniform, normal) 
and writes them to a file.
 Parameters:
    filename - Name of the file to write the random numbers
    type - Type of random number generator (1=Uniform Int, 2=Uniform Real, etc.)
    m, M - Range for uniform generators or truncation
    mu, sigma - Mean and standard deviation for normal generators
    N - Number of random numbers to generate
 */
 void generate_random_numbers_to_file(const char* filename, int type, double m, 
double M, double mu, double sigma, int N) 
{
     FILE* file = fopen(filename, "w");
    if (!file) {
        perror("Error opening file");
        return;
    }

    for (int i = 0; i < N; i++) {
        double value;
        switch (type) {
            case 1:
                value = generate_uniform_integer(m, M);
                break;
            case 2:
                value = generate_uniform_real(m, M);
                break;
            case 3:
                value = generate_normal_integer(mu, sigma);
                break;
            case 4:
                value = generate_normal_real(mu, sigma);
                break;
            case 5:
                value = generate_truncated_normal_integer(m, M, mu, sigma);
                break;
            case 6:
                value = generate_truncated_normal_real(m, M, mu, sigma);
                break;
            default:
                fprintf(stderr, "Unknown generator type\n");
                fclose(file);
                return;
        }
        fprintf(file, "%f\n", value);
    }

    fclose(file);
}
 /*
 Function: calculate_statistics_from_file
 =========================================
 Reads random numbers from a file, computes mean and standard deviation, and 
prints results.
  
Parameters:
    filename - Name of the file to read random numbers from
    N - Number of random numbers in the file
*/
 void calculate_statistics_from_file(const char* filename, int N) 
{
   FILE* file = fopen(filename, "r");
    if (!file) {
        perror("Error opening file");
        return;
    }

    double* data = malloc(N * sizeof(double));
    if (!data) {
        perror("Memory allocation error");
        fclose(file);
        return;
    }

    for (int i = 0; i < N; i++) {
        fscanf(file, "%lf", &data[i]);
    }

    double mean = calculate_mean(data, N);
    double std_dev = calculate_std_dev(data, N, mean);

    printf("File: %s\nMean: %f\nStandard Deviation: %f\n", filename, mean, std_dev);

    free(data);
    fclose(file);
}
 /*
 Function: generate_histogram_from_file
 =======================================
 Reads numbers from a file, generates a histogram, and saves the histogram to a 
new file.
 Parameters:
    input_filename  - File containing the numbers to process
    output_filename - File to write the histogram to
    bins            - Number of bins in the histogram
    min, max        - Range of values to be included in the histogram
 */
 void generate_histogram_from_file(const char* input_filename, const char* 
output_filename, int bins, double min, double max) 
{
    FILE* input = fopen(input_filename, "r");
    FILE* output = fopen(output_filename, "w");

    if (!input || !output) {
        perror("Error opening file");
        if (input) fclose(input);
        if (output) fclose(output);
        return;
    }

    int* histogram = calloc(bins, sizeof(int));
    if (!histogram) {
        perror("Memory allocation error");
        fclose(input);
        fclose(output);
        return;
    }

    double range = max - min;
    double bin_width = range / bins;
    double value;

    while (fscanf(input, "%lf", &value) == 1) {
        if (value >= min && value < max) {
            int bin = (int)((value - min) / bin_width);
            if (bin >= bins) bin = bins - 1;
            histogram[bin]++;
        }
    }

    for (int i = 0; i < bins; i++) {
        fprintf(output,"Bin:[%d] ----> Bin Count:[%d]\n", i, histogram[i]);
    }

    free(histogram);
    fclose(input);
    fclose(output);
}
 /*
 Function: generate_uniform_integer
 ==================================
 Generates a random integer uniformly distributed in [m, M].
 Parameters:
    m - Lower bound of the range
    M - Upper bound of the range
 Returns:
    A random integer in the range [m, M].
 */
 int generate_uniform_integer(double m, double M)
 {
   return (int)m + rand() % ((int)M - (int)m + 1);
}
 /*
 Function: generate_uniform_real
 ================================
 Generates a random real number uniformly distributed in [m, M].
 Parameters:
    m - Lower bound of the range
    M - Upper bound of the range
 
Returns:
    A random real number in the range [m, M].
    
*/
 double generate_uniform_real(double m, double M) 
{
   return m + frand() * (M - m);
}
 /*
 Function: generate_normal_integer
 =================================
 Generates a random integer from a normal distribution with the specified mean and 
standard deviation.
 
Parameters:
    mu - Mean of the normal distribution
    sigma - Standard deviation of the normal distribution
 Returns:
    A random integer from the normal distribution.
 */
 int generate_normal_integer(double mu, double sigma) 
{
      return (int)round(generate_normal_real(mu, sigma));
}
 /*
 Function: generate_normal_real
 ===============================
 Generates a random real number from a normal distribution with the specified mean 
and standard deviation.
 Parameters:
    mu - Mean of the normal distribution
    sigma - Standard deviation of the normal distribution
 
Returns:
    A random real number from the normal distribution.
 */
 double generate_normal_real(double mu, double sigma)
 {
    return mu + sigma * nrand();
}
 /*
 Function: generate_truncated_normal_integer
 ============================================
 Generates a random integer from a normal distribution, truncated to lie within 
the range [m, M].
 
Parameters:
    m - Lower bound of the range
    M - Upper bound of the range
    mu - Mean of the normal distribution
    sigma - Standard deviation of the normal distribution
 
Returns:
    A random integer from the truncated normal distribution.
 */
 int generate_truncated_normal_integer(double m, double M, double mu, double 
sigma) {
{
    return (int)round(generate_truncated_normal_real(m, M, mu, sigma));
}
}
 /*
 Function: generate_truncated_normal_real
 =========================================
 Generates a random real number from a normal distribution, truncated to lie 
within the range [m, M].
 Parameters:
    m - Lower bound of the range
    M - Upper bound of the range
    mu - Mean of the normal distribution
    sigma - Standard deviation of the normal distribution
 Returns:
    A random real number from the truncated normal distribution.
 */
 double generate_truncated_normal_real(double m, double M, double mu, double 
sigma) 
{
    double value;
    do {
        value = generate_normal_real(mu, sigma);
    } while (value < m || value > M);
    return value;
}
 /*
 Function: calculate_mean
 ========================
 Calculates the sample mean of a sequence of numbers.
 Parameters:
    data - Pointer to an array of numbers
    n - Number of elements in the array
 Returns:
    The sample mean of the sequence.
 */
 double calculate_mean(double* data, int n) 
{

  double sum = 0.0;
    for (int i = 0; i < n; i++) {
        sum += data[i];
    }
    return sum / n;
 
}
 /*
 Function: calculate_std_dev
 ===========================
 Calculates the sample standard deviation of a sequence of numbers.
 
Parameters:
    data - Pointer to an array of numbers
    n - Number of elements in the array
    mean - The sample mean of the sequence
 
Returns:
    The sample standard deviation of the sequence.
 */
 double calculate_std_dev(double* data, int n, double mean)
 {
   double sum_sq_diff = 0.0;
    for (int i = 0; i < n; i++) {
        double diff = data[i] - mean;
        sum_sq_diff += diff * diff;
    }
   return sqrt(sum_sq_diff / n);
 }