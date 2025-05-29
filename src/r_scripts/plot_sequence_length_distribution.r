suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(ggplot2))



# Define the options
option_list <- list(
    make_option(c("-s", "--stage"),
        type = "character", default = "",
        help = "Processing stage", metavar = "STAGE"
    ),
    make_option(c("-a", "--amplicon"),
        type = "character", default = "",
        help = "Amplicon", metavar = "AMPLICON"
    ),
    make_option(c("-i", "--input"),
        type = "character", default = "",
        help = "Input file", metavar = "INPUT_FILE"
    ),
    make_option(c("-o", "--output"),
        type = "character", default = "",
        help = "Output file", metavar = "OUTPUT_FILE"
    )
)

# Create parser and parse the arguments
opt_parser <- OptionParser(option_list = option_list)
options <- parse_args(opt_parser)

# Get parameters from options
stage <- options$stage
amplicon <- options$amplicon
input_file <- options$input
output_file <- options$output


amplicon_max_lengths <- list(
    demultiplexed = list("16s" = 2000, "its" = 4000, "16s_leaf" = 1600),
    noprimers = list("16s" = 1800, "its" = 1500, "16s_leaf" = 1500),
    size_filt = list("16s" = 1800, "its" = 1200, "16s_leaf" = 1200),
    chopper = list("16s" = 1800, "its" = 1200, "16s_leaf" = 1200),
    nochimera = list("16s" = 1800, "its" = 1200, "16s_leaf" = 1200),
    nohost = list("16s" = 1800, "its" = 1200, "16s_leaf" = 1200),
    itsxpress = list("16s" = 1800, "its" = 1200, "16s_leaf" = 1200)
)

if (!stage %in% names(amplicon_max_lengths)) {
    stop("Unrecognized stage.")
}

if (!amplicon %in% names(amplicon_max_lengths[[stage]])) {
    stop("Unrecognized amplicon type.")
}

max_length <- amplicon_max_lengths[[stage]][[amplicon]]

# Open a PDF device
pdf(output_file)

# Extract sample name from the file path
sample_name <- sub(".*/(.*?)_length_distribution\\.txt$", "\\1", input_file)

print(paste(
    "Graphing read length distribution of", sample_name, amplicon,
    stage, "sample",
    sep = " "
))

# Read in length distribution of sample as data-frame
read_len_distro <- as.data.frame(read.table(input_file, header = T, sep = "\t"))

# Set colnames to count and length
colnames(read_len_distro) <- c("length", "count")

total_length <- sum(read_len_distro$length * read_len_distro$count)
total_frequency <- sum(read_len_distro$count)

mean_read_length_round <- round(total_length / total_frequency, digits = 0)

# mean_read_length calculated previously
mean_read_length <- total_length / total_frequency

# Calculating the weighted variance
part1 <- sum(read_len_distro$count * (read_len_distro$length - mean_read_length)^2) / total_frequency
part2 <- sum(read_len_distro$count^2 * (read_len_distro$length - mean_read_length)^2) / total_frequency^2
weighted_variance <- part1 - part2

# Calculating the weighted standard deviation
weighted_std_dev <- round(sqrt(weighted_variance), digits = 0)

# Filter data to remove data over the max_length
filtered_data <- subset(read_len_distro, length <= max_length - 1)

seq_len_distro_plot <- ggplot(filtered_data, aes(x = length, y = count)) +
    geom_bar(stat = "identity", position = "dodge", fill = "blue") +
    labs(
        title = paste0("Sequence Length Distribution of\n", stage, " sample ", sample_name),
        subtitle = paste("Mean read length =", mean_read_length_round, "\nweighted standard deviation = ", weighted_std_dev, sep = " "),
        x = "Sequence length", y = "Frequency"
    ) +
    scale_x_continuous(breaks = seq(0, max_length, by = max_length / 10), labels = seq(0, max_length, by = max_length / 10), limits = c(0, max_length)) +
    theme(
        panel.background = element_rect(fill = "white"),
        plot.title = element_text(size = 12, hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(face = "italic", size = 12, hjust = 0.5, family = "Helvetica"),
        axis.line.y = element_line(linetype = 1, linewidth = 0.5, colour = "black"),
        axis.line.x = element_line(linetype = 1, linewidth = 0.5, colour = "black"),
        axis.text.x = element_text(vjust = 0.5, hjust = 0.5, face = "bold", size = 10, color = "black", angle = 45),
        axis.text.y = element_text(vjust = 0.5, face = "bold", size = 10, color = "black"),
        axis.title.x = element_text(face = "bold", size = 10),
        axis.title.y = element_text(face = "bold", size = 10),
        legend.position = "bottom",
        legend.text = element_text(colour = "black", size = 10, face = "italic")
    )

print(seq_len_distro_plot)

# Turn off device
dev.off()

print("Finished graphing")
