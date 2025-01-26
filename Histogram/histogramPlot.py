import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import os

# Folder containing histogram files
FOLDER = "Histogram"  # The directory where .txt histogram files are stored
OUTPUT_PDF = "histograms.pdf"  # The name of the output PDF file containing all histograms
MAX_FILES = 6  # The maximum number of .txt files to process

def read_histogram_data(filename):
    """
    Reads histogram data from a .txt file. The file is expected to have lines formatted as:
    "Bin [0] ----> Count: 476"

    Parameters:
        filename (str): Path to the .txt file containing histogram data.

    Returns:
        bins (list of int): A list of bin indices (e.g., [0, 1, 2, ...]).
        counts (list of int): A list of corresponding bin counts.
    """
    bins = []  # List to store bin indices
    counts = []  # List to store bin counts

    # Open the file and read its content line by line
    with open(filename, "r") as file:
        for line in file:
            try:
                # Parse the line into bin index and count values
                parts = line.split("---->")
                bin_index = int(parts[0].split("[")[1].split("]")[0].strip())  # Extract the bin index
                count = int(parts[1].split(":")[1].strip())  # Extract the count
                bins.append(bin_index)  # Append the bin index to the list
                counts.append(count)  # Append the count to the list
            except (IndexError, ValueError):
                # Handle lines that don't match the expected format
                print(f"Skipping invalid line in {filename}: {line}")

    return bins, counts  # Return the parsed bin indices and counts

def plot_histogram(bins, counts, title, file_label):
    """
    Plots a single histogram based on the provided bin indices and counts.

    Parameters:
        bins (list of int): A list of bin indices.
        counts (list of int): A list of bin counts.
        title (str): Title of the histogram plot.
        file_label (str): File name to display below the plot as additional information.
    """
    # Create a bar chart with bins on the x-axis and counts on the y-axis
    plt.bar(bins, counts, color='blue', alpha=0.7, width=0.8)  # Bar chart appearance
    plt.xlabel("Bins")  # Label for the x-axis
    plt.ylabel("Counts")  # Label for the y-axis
    plt.title(title)  # Set the title of the plot
    plt.xticks(bins, fontsize=8)  # Customize x-axis tick labels with smaller font
    plt.tight_layout()  # Adjust layout to avoid clipping of labels

    # Add the file name below the plot as a figure text
    plt.figtext(0.5, -0.02, f"Source File: {file_label}", wrap=True, horizontalalignment='center', fontsize=10)

def get_histogram_files(folder, max_files):
    """
    Retrieves a list of up to `max_files` .txt files from the specified folder.

    Parameters:
        folder (str): Path to the folder containing histogram files.
        max_files (int): Maximum number of files to retrieve.

    Returns:
        list of str: List of full file paths to the .txt files.
    """
    files = []  # List to store paths to .txt files

    # Iterate over all files in the folder
    for file in os.listdir(folder):
        if file.endswith(".txt"):  # Check if the file has a .txt extension
            files.append(os.path.join(folder, file))  # Add the full path to the list

    files.sort()  # Sort the files alphabetically for consistent ordering
    return files[:max_files]  # Return only the first `max_files` files

# Main logic for processing and plotting histograms
with PdfPages(OUTPUT_PDF) as pdf:  # Open a PDF file to save the histograms
    histogram_files = get_histogram_files(FOLDER, MAX_FILES)  # Get the list of .txt files to process

    # Check if there are no .txt files in the folder
    if not histogram_files:
        print(f"No .txt files found in folder: {FOLDER}")
        exit()  # Exit the program if no files are found

    # Loop through each .txt file and process it
    for filename in histogram_files:
        print(f"Processing {filename}...")  # Notify which file is being processed

        # Read histogram data from the file
        bins, counts = read_histogram_data(filename)

        # Extract only the file name (excluding the folder path)
        file_label = os.path.basename(filename)

        # Create a new figure for the histogram plot
        plt.figure(figsize=(10, 5))
        plot_histogram(bins, counts, title=f"Histogram for {file_label}", file_label=file_label)

        # Save the current figure into the PDF
        pdf.savefig(bbox_inches='tight')  # Save the figure to the PDF file
        plt.close()  # Close the figure to free memory

# Notify the user that all histograms have been saved
print(f"All histograms saved to {OUTPUT_PDF}")
