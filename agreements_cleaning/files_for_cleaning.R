
# 1. select all rds files in 4.MB_indv_amt_1h that make it into the database 
# 2. move to folder for cleaning

library(tidyverse);library(bit64); library(geodata);library(dggridR);library(terra); library(sf); library(mapview);library(RColorBrewer);library(dplyr);library(stringr);library(purrr);library(ggspatial)

mt_ind <- readRDS("DATA/MoveTraitsData/8.MoveTraits_db/MoveTrait.v0.1_individual.sum_20260807.rds")

mt_ind <- mt_ind |> 
  filter(source %in% "movebank.august2026")

lookup <- paste(mt_ind$study_id,mt_ind$individual_id,sep="_")
lookup <- paste(lookup,"rds",sep=".")
length(lookup)

# Copy files for cleaning
src_dir  <- "DATA/MoveTraitsData/4.MB_indv_amt_1h"
dest_dir <- "DATA/MoveTraitsData/5.MB_indv_amt_1h_cleaning"
dir.create(dest_dir, showWarnings = FALSE, recursive = TRUE)

src_paths  <- file.path(src_dir, lookup)
dest_paths <- file.path(dest_dir, lookup)
exists     <- file.exists(src_paths)

cat(sprintf("%d of %d files found in source directory\n", sum(exists), length(lookup)))
file.copy(src_paths[exists], dest_paths[exists])
cat(sprintf("%d files copied to %s\n", sum(exists), dest_dir))

# 13360 of 13360 files found in source directory

# when were files downloaded? Need to also clean the data that were shared by dataowners (> 09.2026)

# Overview of file creation dates in 1.MB_indv_mv2
mv2_files <- list.files("DATA/MoveTraitsData/1.MB_indv_mv2", full.names = TRUE)
mv2_info  <- file.info(mv2_files)

mv2_info |>
  tibble::rownames_to_column("file") |>
  mutate(date = as.Date(mtime)) |>
  count(date, name = "n_files") |>
  arrange(date) 

# 3. move all individuals from one study into seperate folders 

# Check dest_dir and preview file names
files <- list.files(dest_dir, pattern = "\\.rds$", full.names = FALSE)
length(files)
head(files, 10)

# Create one subfolder per study and move files into it
study_ids <- sub("_.*$", "", files)

# Create all directories first
unique_studies <- unique(study_ids)
dirs_to_create <- file.path(dest_dir, unique_studies)
invisible(sapply(dirs_to_create, dir.create, showWarnings = FALSE))

# Move files into their respective study folder
src_paths  <- file.path(dest_dir, files)
dest_paths <- file.path(dest_dir, study_ids, files)

results <- file.rename(src_paths, dest_paths)

# 4. divide into 4 folders for sofia, diego, kevin and anne 

# Get file count per study folder
study_folders <- list.dirs(dest_dir, recursive = FALSE, full.names = TRUE)
study_names   <- basename(study_folders)

file_counts <- sapply(study_folders, function(d) {
  length(list.files(d, pattern = "\\.rds$"))
})

total_files <- sum(file_counts)

# Greedy bin-packing assignment
# Sort studies descending by file count (helps balance large studies)
order_idx    <- order(file_counts, decreasing = TRUE)
sorted_names <- study_names[order_idx]
sorted_counts <- file_counts[order_idx]

targets <- c(sofia = 0.30, diego = 0.30, kevin = 0.25, anne = 0.15) * total_files
assigned_counts <- c(sofia = 0, diego = 0, kevin = 0, anne = 0)
assignments <- character(length(sorted_names))

for (i in seq_along(sorted_names)) {
  # Assign to person with most remaining quota
  remaining <- targets - assigned_counts
  scorer <- names(which.max(remaining))
  assignments[i] <- scorer
  assigned_counts[scorer] <- assigned_counts[scorer] + sorted_counts[i]
}

# Summary
result_df <- data.frame(
  study   = sorted_names,
  files   = sorted_counts,
  scorer  = assignments
)

result_df |>
  dplyr::summarise(
    .by = scorer,
    n_studies = dplyr::n(),
    n_files   = sum(files),
    pct_files = round(sum(files) / total_files * 100, 1)
  ) |>
  dplyr::arrange(dplyr::desc(n_files))


# Create the four scorer folders
scorers <- c("sofia", "diego", "kevin", "anne")
scorer_dirs <- file.path(dest_dir, scorers)
invisible(sapply(scorer_dirs, dir.create, showWarnings = FALSE))

# Move each study folder into the assigned scorer folder
move_results <- mapply(function(study, scorer) {
  src  <- file.path(dest_dir, study)
  dst  <- file.path(dest_dir, scorer, study)
  file.rename(src, dst)
}, result_df$study, result_df$scorer)

cat("Study folders moved successfully:", sum(move_results), "\n")
cat("Failed:", sum(!move_results), "\n")

# 5. move to dropbox
# 6. write code to extract from each study folder the rds files in the "scrubdata/cleaned_files" 
#    folder into one overall location "cleaned amt"
# 

