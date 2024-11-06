#===========Plot GONE results
# Set the working directory
setwd("/faststorage/project/Coregonus/Aja/Dragonfly/GONE/")

# Create the output PDF file
pdf(file="plots/GONE_7.37cM_hc0.05.NEZ+SJ.jackknife_CI.pdf", width = 100/25.4, height = 45/25.4)

# Create the layout for side-by-side plots
nf <- layout(matrix(c(1, 2), nrow=1, ncol=2, byrow=TRUE))

# Define populations and corresponding labels/colors
pops <- c("NEZALL", "SJALL")
labels <- c("NEZ", "SJ")
col.palette <- c("#efdf48", "#8bd346")

# Loop through each population
for (i in seq_along(pops)) {
  pop <- pops[i]
  col <- col.palette[i]  # Assign color
  label <- labels[i]     # Assign label
  
  # Plot main population results
  data <- read.table(paste(pop, "/GONE_hc0.05_7.37cM/Output_Ne_", pop, "_filtered.minGQ20.minDP6.maxDP45.HWE6.noSexChr.100longest.miss1.mac1", sep=""), header=TRUE, skip=1)
  
  # Read the jackknife file for the population to get individual names
  jackknife_file <- paste0(pop, ".txt")  # Replace with actual file name if different
  jackknife_data <- read.table(jackknife_file, header=F)
  individuals <- jackknife_data$V1  # Adjust column name if the file has headers
  
  # Initialize a list to store jackknife data for each individual
  jackknife_data <- list()
  
  # Loop through each individual to gather jackknife data
  for (indiv in individuals) {
    indiv_data_file <- paste(pop, "/", indiv, "/GONE_hc0.05_7.37cM/Output_Ne_", pop,"_", indiv, "_filtered.minGQ20.minDP6.maxDP45.HWE6.noSexChr.100longest.miss1.mac1", sep="")
    indiv_data <- read.table(indiv_data_file, header=TRUE, skip=1)
    jackknife_data[[indiv]] <- indiv_data$Geometric_mean  # Store only the Geometric_mean for each generation
  }
  
  # Combine jackknife data into a matrix
  jackknife_matrix <- do.call(cbind, jackknife_data)
  
  # Calculate mean, lower, and upper bounds (95% CI) for each generation
  jackknife_mean <- rowMeans(jackknife_matrix, na.rm=TRUE)
  jackknife_lower <- apply(jackknife_matrix, 1, function(x) quantile(x, 0.025, na.rm=TRUE))
  jackknife_upper <- apply(jackknife_matrix, 1, function(x) quantile(x, 0.975, na.rm=TRUE))
  
  # Plot the main population line
  par(mgp=c(1.5, 0.5, 0), mar=c(2.5, 3, 1, 0.5)+0.1, cex=0.6)
  plot(data$Generation, data$Geometric_mean,
       type = "l",
       xlab = "Years ago",
       ylab = "Effective population size",
       ylim = c(0, 4000),
       xlim = c(0, 400),
       col = col,
       lwd = 2)  # Main line width for emphasis
  
# Add dashed lines for the 95% CI
  lines(data$Generation, jackknife_lower, col=adjustcolor(col, alpha.f=0.5), lty=2)
  lines(data$Generation, jackknife_upper, col=adjustcolor(col, alpha.f=0.5), lty=2)
  
  # Add legend
  legend("topleft",
         legend = c(label, "95% CI"),
         col = c(col, adjustcolor(col, alpha.f=0.5)),
         lty = c(1, 2), 
         lwd = c(2, 1),
         cex = 0.75)
}

# Close the PDF device
dev.off()