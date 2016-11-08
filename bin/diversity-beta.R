# diversity-beta.R
# Geoffrey Hannigan
# Schloss Lab
# University of Michigan

###################
# Set Environment #
###################

write("Calculating alpha diversity", stderr())

library("optparse")
library("ggplot2")
library("plotROC")
library("reshape2")
library("wesanderson")
library("vegan")
library("plyr")

option_list <- list(
  make_option(c("-i", "--input"),
    type = "character",
    default = NULL,
    help = "OTU abundance table straight from Mothur.",
    metavar = "character"),
  make_option(c("-m", "--metadata"),
    type = "character",
    default = NULL,
    help = "Metadata table with disease states.",
    metavar = "character"),
  make_option(c("-s", "--subsample"),
    type = "integer",
    default = 10000,
    help = "Metadata table with disease states.",
    metavar = "character"),
  make_option(c("-o", "--out"),
    type = "character",
    default = NULL,
    help = "Output pdf to draw results.",
    metavar = "character")
)
