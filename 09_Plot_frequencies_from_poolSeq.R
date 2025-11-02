# Load necessary libraries
library(dplyr)
library(ggplot2)

# List all files matching the pattern
file_list <- list.files(pattern = "Helene_clone_magnapulex_singletons_mQ100_GC153502_ALL_freqs.txt")
file_list <- list.files(pattern = "Helene_clone_magnapulex_singletons_mQ100_GC159229_ALL_freqs.txt")
file_list <- list.files(pattern = "Helene_clone_magnapulex_singletons_mQ100_GC159230_ALL_freqs.txt")
file_list <- list.files(pattern = "Helene_clone_magnapulex_singletons_mQ100_GC159231_ALL_freqs.txt")
file_list <- list.files(pattern = "Helene_clone_magnapulex_singletons_mQ100_GC159232_ALL_freqs.txt")
file_list <- list.files(pattern = "Helene_clone_magnapulex_singletons_mQ100_GC159233_ALL_freqs.txt")
file_list <- list.files(pattern = "Helene_clone_magnapulex_singletons_mQ100_GC159234_ALL_freqs.txt")

# Read and combine all files into one data frame
combined_data <- lapply(file_list, read.table, sep = "\t", header = TRUE) %>%
  bind_rows()

combined_data$Alt_Allele_Frequency_corrected <- combined_data$Alt_Allele_Frequency * combined_data$Homozygous_or_Heterozygous




# Plot the density distribution of Alt_Allele_Frequency_corrected for each sample in the file
ggplot(combined_data, aes(x = Alt_Allele_Frequency_corrected, color = Sample, fill = Sample)) +
  geom_histogram(alpha = 0.5, bins = 100) +
  # scale_x_log10() +
  labs(title = "Density Distribution of Alt_Allele_Frequency_corrected by Sample (Log Scale)",
       x = "Alt Allele Frequency Corrected",
       y = "Density (Log Scale)") +
  theme_minimal() +
  theme(legend.position = "right") + facet_grid(Sample~., scales = "free_y") + xlim(-0.01,1)  #+ scale_y_log10(minor_breaks = NULL) #+ # Apply log scale to y-axis
  #scale_y_continuous(trans = 'log1p')  # Log transform adding a small constant + coord_cartesian(ylim = c(0, 100)) 
# 
# # Plot the density distribution of Alt_Allele_Frequency_corrected for each sample
# ggplot(combined_data, aes(x = Alt_Allele_Frequency_corrected, color = Sample, fill = Sample)) +
#   geom_histogram(alpha = 0.5, bins = 100) +
#   # scale_x_log10() +
#   labs(title = "Density Distribution of Alt_Allele_Frequency_corrected by Sample (Log Scale)",
#        x = "Alt Allele Frequency Corrected",
#        y = "Density (Log Scale)") +
#   theme_minimal() +
#   theme(legend.position = "right") + facet_grid(Sample~.) #+ xlim(0,2) #+ scale_y_log10(minor_breaks = NULL) #+ # Apply log scale to y-axis
# #scale_y_continuous(trans = 'log1p')  # Log transform adding a small constant + coord_cartesian(ylim = c(0, 100)) 


summary_stats <- combined_data %>%
  group_by(Sample) %>%
  summarise(
    Mean_Alt_Allele_Frequency = mean(Alt_Allele_Frequency_corrected, na.rm = TRUE),
    Median_Alt_Allele_Frequency = median(Alt_Allele_Frequency_corrected, na.rm = TRUE))

sum(summary_stats$Mean_Alt_Allele_Frequency)

sum(summary_stats$Median_Alt_Allele_Frequency)

write.table(summary_stats[,c(4,6,7)], file = "summary_stats_GC159234.txt", sep = '\t', quote = F, row.names = F)



summary_stats <- combined_data %>%
  group_by(Sample) %>%
  summarise(
    Occurrences = n(),  # Count the number of rows (occurrences) for each sample
    Non_Zero_Alt_Allele_Frequency = sum(Alt_Allele_Frequency_corrected != 0, na.rm = TRUE),  # Count how many times Alt_Allele_Frequency_corrected is not zero
    Mean_Alt_Allele_Frequency = if_else(
      Non_Zero_Alt_Allele_Frequency / Occurrences > 0.5, 
      mean(Alt_Allele_Frequency_corrected[Alt_Allele_Frequency_corrected != 0], na.rm = TRUE),  # Exclude zero if more than 40% of occurrences are non-zero
      mean(Alt_Allele_Frequency_corrected, na.rm = TRUE)  # Otherwise, calculate the mean normally
    ),
    Median_Alt_Allele_Frequency = if_else(
      Non_Zero_Alt_Allele_Frequency / Occurrences > 0.5, 
      median(Alt_Allele_Frequency_corrected[Alt_Allele_Frequency_corrected != 0], na.rm = TRUE),  # Exclude zero if more than 40% of occurrences are non-zero
      median(Alt_Allele_Frequency_corrected, na.rm = TRUE)  # Otherwise, calculate the mean normally
    ),
    CI_Lower_95 = if_else(
      Non_Zero_Alt_Allele_Frequency / Occurrences > 0.5, 
      mean(Alt_Allele_Frequency_corrected[Alt_Allele_Frequency_corrected != 0], na.rm = TRUE) - 1.96 * sd(Alt_Allele_Frequency_corrected[Alt_Allele_Frequency_corrected != 0], na.rm = TRUE) / sqrt(Occurrences),
      mean(Alt_Allele_Frequency_corrected, na.rm = TRUE) - 1.96 * sd(Alt_Allele_Frequency_corrected, na.rm = TRUE) / sqrt(Occurrences)
    ),
    CI_Upper_95 = if_else(
      Non_Zero_Alt_Allele_Frequency / Occurrences > 0.5, 
      mean(Alt_Allele_Frequency_corrected[Alt_Allele_Frequency_corrected != 0], na.rm = TRUE) + 1.96 * sd(Alt_Allele_Frequency_corrected[Alt_Allele_Frequency_corrected != 0], na.rm = TRUE) / sqrt(Occurrences),
      mean(Alt_Allele_Frequency_corrected, na.rm = TRUE) + 1.96 * sd(Alt_Allele_Frequency_corrected, na.rm = TRUE) / sqrt(Occurrences)
    )
  )


# Print the summary statistics for each input file (=sample)
print(summary_stats, n=30)

write.table(summary_stats[,c(4,6,7)], file = "summary_stats_GC159234.txt", sep = '\t', quote = F, row.names = F)

sum(summary_stats$Mean_Alt_Allele_Frequency)

sum(summary_stats$Median_Alt_Allele_Frequency)
