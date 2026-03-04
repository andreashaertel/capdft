# SPDX-FileCopyrightText: 2022 Moritz Bültmann <moritz.bueltmann@physik.uni-freiburg.de>
# SPDX-License-Identifier: CC0-1.0
# ______________________________________________________________________________
# This script creates a heat map of two slices of the 3D density profile created
# by the example program.
# ______________________________________________________________________________
rm(list = ls())
library(tidyverse)
library(viridis)
data <- read_delim("3d_profile.dat", delim = " ", col_names = FALSE, comment = "#", trim_ws = TRUE)
names(data) <- c("x", "y", "z", "density")
# First plot
data_plot <- data %>% filter(y == .25) %>% filter(density != 0.)
my_plot1 <- ggplot(data_plot, aes(x = z, y = x, fill = density)) +
    geom_tile() +
    scale_fill_viridis(option = "turbo")
dev.new()
print(my_plot1)
ggsave("3d_profile1.pdf", my_plot1, width = 3, height = 3)
# Second plot
data_plot <- data %>% filter(y == 0.) %>% filter(density != 0.)
my_plot2 <- ggplot(data_plot, aes(x = z, y = x, fill = density)) +
    geom_tile() +
    scale_fill_viridis(option = "turbo")
dev.new()
print(my_plot2)
ggsave("3d_profile2.pdf", my_plot2, width = 3, height = 3)
