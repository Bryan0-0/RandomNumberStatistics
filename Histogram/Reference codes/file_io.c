/******************************************************************************
 * File: file_io_random_numbers.c
 * 
 * Purpose:
 * This program demonstrates file I/O operations in C by writing random numbers
 * to a text file and then reading the data back to display it on the console.
 * 
 * Features:
 * 1. Generates 100 random numbers using the C library's `rand()` function.
 * 2. Writes each random number to a text file in the format: "Random number X: [value]".
 * 3. Reads the random numbers back from the file line by line.
 * 4. Displays the contents of the file on the console.
 * 
 * Output:
 * - A text file named `random_numbers.txt` containing 100 random numbers.
 * - Console output displaying the file contents.
 * 
 * Usage:
 * - Compile: gcc file_io_random_numbers.c -o file_io_random_numbers
 * - Execute: ./file_io_random_numbers
 * 
 * Experimentation Instructions:
 * 1. Modify the number of random numbers (`NUM_RANDOMS`) to generate a different number.
 * 2. Use different random number generation logic if needed (e.g., normal distribution).
 * 3. Change the file name for saving output to test file handling.
 * 
 *****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

// Number of random numbers to generate
#define NUM_RANDOMS 100

// Function prototypes
void write_random_numbers_to_file(const char *filename, int count);
void read_and_display_file(const char *filename);

int main() {
    const char *filename = "random_numbers.txt"; // Output file name

    // Seed the random number generator with the current time
    srand(time(NULL));

    // Write random numbers to the file
    printf("Writing %d random numbers to the file: %s\n", NUM_RANDOMS, filename);
    write_random_numbers_to_file(filename, NUM_RANDOMS);

    // Read and display the file contents
    printf("\nReading and displaying the contents of the file:\n");
    read_and_display_file(filename);

    return 0;
}

/**
 * Writes `count` random numbers to the specified file.
 * 
 * Each number is written in the format: "Random number X: [random value]"
 * 
 * Parameters:
 *   filename - Name of the file to write to.
 *   count - Number of random numbers to generate and write.
 */
void write_random_numbers_to_file(const char *filename, int count) {
    FILE *file = fopen(filename, "w"); // Open the file in write mode
    if (file == NULL) {
        perror("Error opening file for writing");
        exit(EXIT_FAILURE);
    }

    // Write random numbers to the file
    for (int i = 1; i <= count; i++) {
        int random_number = rand(); // Generate a random number
        fprintf(file, "Random number [%d] --> %d\n", i, random_number); // Write to file
    }

    fclose(file); // Close the file after writing
    printf("Successfully wrote %d random numbers to the file: %s\n", count, filename);
}

/**
 * Reads the contents of the specified file and displays them on the console.
 * 
 * Each line of the file is printed as it is read.
 * 
 * Parameters:
 *   filename - Name of the file to read from.
 */
void read_and_display_file(const char *filename) {
    FILE *file = fopen(filename, "r"); // Open the file in read mode
    if (file == NULL) {
        perror("Error opening file for reading");
        exit(EXIT_FAILURE);
    }

    char line[256]; // Buffer to store each line read from the file

    // Read and print each line of the file
    while (fgets(line, sizeof(line), file) != NULL) {
        printf("%s", line); // Print the line to the console
    }

    fclose(file); // Close the file after reading
}
