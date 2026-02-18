library(data.table)
library(ggplot2)
library(dplyr)
library(scales)

#######
# MOI range
MOI_range  <- seq(0.05, 1.5, by = 0.05)

# Poisson-derived probabilities
df <- data.frame(
  MOI = MOI_range,
  Infected = 1 - exp(-MOI_range),
  Single = MOI_range * exp(-MOI_range)
)

# Multiple integrations (>1)
df$Multiple <- df$Infected - df$Single

# Convert to long format for ggplot
df_long <- df %>%
  tidyr::pivot_longer(cols = c("Infected", "Single", "Multiple"),
               names_to = "Category",
               values_to = "Probability")

pdf("out1.pdf")

ggplot(df_long, aes(x = MOI, y = Probability, color = Category)) +
  geom_line(size = 1.2) +
  labs(
    title = "Poisson Model of Viral Infection vs MOI",
    x = "MOI",
    y = "Fraction of Cells"
  ) +
  theme_minimal(base_size = 14)

dev.off()

################
# Parameters
gRNAs     <- 20e3
MOI_range <- seq(0.2, 1, by = 0.05)
coverages <- c(300, 500, 1000)

# Create full grid of MOI × coverage
dt <- CJ(MOI = MOI_range, Coverage = coverages)

# Fraction infected (Poisson model)
dt[, FractionInfected := 1 - exp(-MOI)]

# Required infected cells
dt[, InfectedCellsNeeded := gRNAs * Coverage]

# Total starting cells required
dt[, TotalCellsRequired := InfectedCellsNeeded / FractionInfected]

# Convert coverage to factor for plotting
dt[, Coverage := factor(Coverage)]

pdf("out2.pdf")

ggplot(dt, aes(x = MOI, y = TotalCellsRequired, color = Coverage)) +
  geom_line(size = 1.2) +
  geom_point() +
  scale_y_continuous(labels = comma) +
  labs(
    title = "Total Cells Required vs MOI",
    x = "MOI",
    y = "Total Starting Cells Required",
    linetype = "Coverage"
  ) +
  theme_minimal(base_size = 14)

dev.off()

##############
gRNAs     <- 20e3
MOI_range <- seq(0.2, 1, by = 0.05)
coverages <- c(300, 500, 1000)

# Create full grid of MOI × coverage
dt <- CJ(MOI = MOI_range, Coverage = coverages)

# Probability of exactly one integration
dt[, FractionSingle := MOI * exp(-MOI)]

# Required singly infected cells
dt[, SingleCellsNeeded := gRNAs * Coverage]

# Total starting cells required to achieve that many singly infected cells
dt[, TotalCellsRequired := SingleCellsNeeded / FractionSingle]

# Convert coverage to factor for plotting
dt[, Coverage := factor(Coverage)]

pdf("out3.pdf")

	ggplot(dt, aes(x = MOI, y = TotalCellsRequired, color = Coverage)) +
	  geom_line(size = 1.2) +
	  geom_point() +
	  scale_y_continuous(labels = comma) +
	  labs(
	    title = "Total Cells Required vs MOI",
	    subtitle = "Coverage Based on Singly Infected Cells",
	    x = "MOI",
	    y = "Total Starting Cells Required",
	    color = "Coverage"
	  ) +
	  theme_minimal(base_size = 14)

	ggplot(dt[MOI > 0.2 & MOI < 0.65], aes(x = MOI, y = TotalCellsRequired, color = Coverage)) +
	  geom_line(size = 1.2) +
	  geom_point() +
	  scale_y_continuous(
		  labels = comma,
		  breaks = scales::pretty_breaks(n = 12)
		) +
	  labs(
	    title = "Total Cells Required vs MOI",
	    subtitle = "Coverage Based on Singly Infected Cells",
	    x = "MOI",
	    y = "Total Starting Cells Required",
	    color = "Coverage"
	  ) +
	  theme_minimal(base_size = 14)

dev.off()








