/******************************************************************************
 * File: histogram_generator.c
 * 
 * Purpose:
 * This program generates a histogram representing the distribution of random 
 * numbers drawn from a normal distribution. It calculates how many numbers 
 * fall into each bin within a specified range and outputs the histogram.
 * 
 * Features:
 * 1. Generates random numbers from a normal distribution using the Box-Muller 
 *    transform.
 * 2. Defines a range of values for the distribution based on the mean (`mu`) 
 *    and standard deviation (`sigma`).
 * 3. Truncates values that fall outside the range to ensure the histogram is 
 *    accurate and meaningful.
 * 4. Divides the range into a fixed number of bins (`N_BINS`) and counts how 
 *    many numbers fall into each bin.
 * 5. Outputs the histogram data as a list of counts, where each count 
 *    corresponds to a bin.
 * 
 * Output:
 * - Prints the histogram data to the console as a series of integers, one per line.
 * 
 * Usage:
 * - Compile: gcc histogram_generator.c -lm -o histogram
 * - Execute: ./histogram
 * 
 * Experimentation Instructions:
 * 1. Modify `mu` (mean), `sigma` (standard deviation), and `n` (number of random 
 *    numbers) to experiment with different distributions.
 * 2. Change the `min` and `max` values to adjust the truncation range.
 * 3. Adjust `N_BINS` to increase or decrease the granularity of the histogram.
 * 
 *****************************************************************************/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>

// Macros for generating random numbers
#define frand() (rand() / (double)RAND_MAX) // Uniform random number in [0, 1)
#define nrand() (sqrt(-2 * log(frand())) * cos(2 * M_PI * frand())) // Normal random number

// Number of bins for the histogram
#define N_BINS 50 

int main() {
    // Variables
    int i, n = 20000;           // Number of random numbers to generate
    double temp;                // Temporary variable for each random number
    double mu = 100, sigma = 10; // Mean and standard deviation for normal distribution
    int histogram[N_BINS];      // Array to store the counts for each bin

    // Initialize random seed
    srand(time(NULL));

    // Truncation range: four standard deviations from the mean
    double min = mu - 4 * sigma; 
    double max = mu + 4 * sigma;

    // Initialize the histogram bins to 0
    for (i = 0; i < N_BINS; i++) {
        histogram[i] = 0;
    }

    // Generate `n` normally distributed random numbers and populate the histogram
    for (i = 0; i < n; i++) {
        temp = nrand() * sigma + mu; // Generate a normal random number

        // Skip numbers outside the truncation range
        if (temp < min || temp > max) {
            i--; // Retry for this iteration
            continue;
        }

        // Determine the bin index for the random number
        int bin_index = (int)floor((temp - min) * N_BINS / (max - min));

        // Increment the count for the corresponding bin
        histogram[bin_index]++;
    }

    // Print the histogram data
    for (i = 0; i < N_BINS; i++) {
        printf("Bin:[%d] ----> Bin Count:[%d]\n", i, histogram[i]); // Print the count for each bin
    }

    return 0;
}
