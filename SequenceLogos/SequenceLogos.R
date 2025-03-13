# Install and load required packages
if (!requireNamespace("ggseqlogo", quietly = TRUE)) {
  install.packages("ggseqlogo")
}
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2")
}
if (!requireNamespace("reshape2", quietly = TRUE)) {
  install.packages("reshape2")
}
library(ggseqlogo)
library(ggplot2)
library(reshape2)

# Function to get the script's directory
get_script_dir <- function() {
  # Try to get the script path from RStudio
  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    return(dirname(rstudioapi::getActiveDocumentContext()$path))
  }

  # Try to get the script path from the command line
  cmd_args <- commandArgs(trailingOnly = FALSE)
  script_path <- grep("--file=", cmd_args, value = TRUE)
  if (length(script_path) > 0) {
    return(dirname(sub("--file=", "", script_path)))
  }

  # Fallback: Use the current working directory
  return(getwd())
}

# Get the directory where the script is located
script_dir <- get_script_dir()

# Set the working directory to the script's directory
setwd(script_dir)

# List all CSV files in the directory
csv_files <- list.files(path = script_dir, pattern = "\\.csv$", full.names = TRUE)

# Check if there are any CSV files
if (length(csv_files) == 0) {
  stop("No CSV files found in the directory. Exiting script.")
}

# Define the order of amino acids explicitly
amino_acid_order <- c("A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y")

# Function to process a single CSV file and create a sequence logo
process_csv_file <- function(input_file) {
  # Load the data from the CSV file
  data <- read.csv(input_file)

  # Ensure the data has the correct column names
  if (!all(c("Amino.Acid", "Position", "Frequency") %in% colnames(data))) {
    stop("The CSV file must contain columns: Amino.Acid, Position, Frequency")
  }

  # Sort the data by Position and Amino Acid
  data <- data[order(data$Position, data$Amino.Acid), ]

  # Debug: Print the sorted data
  print("Sorted Data:")
  print(head(data))

  # Pivot the data to wide format
  # Rows: Amino Acids
  # Columns: Positions
  # Values: Frequencies
  wide_data <- reshape2::dcast(data, Amino.Acid ~ Position, value.var = "Frequency", fill = 0)

  # Debug: Print the reshaped data
  print("Reshaped Data (Wide Format):")
  print(wide_data)

  # Ensure the amino acids are in the correct order
  wide_data <- wide_data[match(amino_acid_order, wide_data$Amino.Acid), ]

  # Convert to a matrix (required by ggseqlogo)
  wide_data <- as.matrix(wide_data[, -1])  # Remove the Amino Acid column
  rownames(wide_data) <- amino_acid_order  # Set row names to amino acids in the correct order

  # Debug: Print the matrix
  print("Matrix for ggseqlogo:")
  print(wide_data)

  # Normalize the frequencies to probabilities (ensure they sum to 1 at each position)
  wide_data <- sweep(wide_data, 2, colSums(wide_data), `/`)

  # Debug: Print the normalized matrix
  print("Normalized Matrix:")
  print(wide_data)

  # Debug: Verify column sums
  print("Column Sums (Should be 1):")
  print(colSums(wide_data))

  # Create the sequence logo using ggseqlogo
  p <- ggseqlogo(wide_data, method = "prob", seq_type = "aa", col_scheme = "chemistry") +
    theme_classic() +
    theme(
      axis.title.x = element_text(size = 14, face = "bold"),  # Customize x-axis title
      axis.title.y = element_text(size = 14, face = "bold"),  # Customize y-axis title
      axis.text.x = element_text(size = 12, face = "bold"),   # Customize x-axis labels
      axis.text.y = element_text(size = 12, face = "bold"),   # Customize y-axis labels
      axis.ticks = element_blank(),                          # Remove ticks
      panel.border = element_rect(color = "black", fill = NA, size = 1)  # Add a black frame
    ) +
    labs(
      x = "Position",  # X-axis label
      y = "Probability"  # Y-axis label
    )

  # Customize x-axis labels (P1, P2, P3, etc.)
  positions <- colnames(wide_data)
  p <- p + scale_x_continuous(breaks = 1:length(positions), labels = paste0("P", positions))

  # Remove the logo title
  p <- p + ggtitle("")

  # Get the base name of the input file (without extension)
  input_basename <- tools::file_path_sans_ext(basename(input_file))

  # Define the output file path
  output_file <- file.path(script_dir, paste0(input_basename, "_probability_logo.png"))

  # Save the plot
  ggsave(output_file, p, width = 10, height = 5, dpi = 300)

  # Print the path where the file was saved
  print(paste("Probability logo saved to:", output_file))
}

# Process all CSV files
for (csv_file in csv_files) {
  process_csv_file(csv_file)
}
