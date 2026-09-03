library(tidyverse)

## ----Load data sharing agreements -------------------------------------------------------------

# Unzip each agreement into its own subfolder
zip_files <- list.files("./DATA/MoveTraitsData/data_agreements", pattern = "\\.zip$", full.names = TRUE)

for (f in zip_files) {
  out_dir <- tools::file_path_sans_ext(f)
  unzip(f, exdir = out_dir)
}

# Find all CSVs in the extracted subfolders (exclude output subfolders)
csv_files <- list.files(
  "./DATA/MoveTraitsData/data_agreements",
  pattern = "\\.csv$",
  full.names = TRUE,
  recursive = TRUE
) |>
  Filter(\(f) !grepl("data_agreements_pdf", f), x = _)

# Load and bind
agreements_bind <- csv_files |>
  set_names(csv_files) |>
  map(read_csv, show_col_types = FALSE) |>
  list_rbind(names_to = "source_file") |>
  mutate(source_file = basename(source_file))

## ----wrangle spatial resolutions to be shared and remove duplicates -------------------------------------------------------------

agreements <- agreements_bind |> 
  dplyr::select(1:10)

# Classify duplicated Study_ids as identical vs conflicting
dup_ids <- agreements |>
  group_by(Study_id) |>
  filter(n() > 1) |>
  pull(Study_id) |>
  unique()

id_classification <- agreements |>
  filter(Study_id %in% dup_ids) |>
  group_by(Study_id) |>
  summarise(all_identical = nrow(distinct(pick(-source_file))) == 1)

identical_ids <- id_classification |> filter(all_identical)  |> pull(Study_id)
conflict_ids  <- id_classification |> filter(!all_identical) |> pull(Study_id)

# Conflicts: Study_ids with differing rows — extract for manual review
agreements_conflicts <- agreements |> filter(Study_id %in% conflict_ids)

# Clean: keep one row for true duplicates, exclude conflicts pending review
agreements_clean <- agreements |>
  filter(!Study_id %in% conflict_ids) |>
  distinct(across(-source_file), .keep_all = TRUE)

# Expand peculiarity into binary yes/no columns (one per survey option)
agreements_clean <- agreements_clean |>
  separate_rows(peculiarity, sep = ";\\s*") |>
  mutate(present = "Yes") |>
  pivot_wider(
    names_from = peculiarity,
    values_from = present,
    values_fill = "No"
  ) |>
  select(-any_of("None")) |>
  rename(
    peculiarity_no_cleaning       = `Has not undergone thorough data cleaning`,
    peculiarity_experimental      = `Is recorded under experimentally altered conditions`,
    peculiarity_semi_domesticated = `Stems from (semi-)domesticated animals`
  )

agreements_clean2 <- agreements_clean |>
  mutate(
    species            = str_extract(species, "\\d+"),
    individual         = str_extract(individual, "\\d+"),
    individual_monthly = str_extract(individual_monthly, "\\d+"),
    within_individual  = if_else(
      str_detect(within_individual, "observation"),
      "observation",
      str_extract(within_individual, "\\d+")
    )
  )

saveRDS(agreements_clean2,"./DATA/MoveTraitsData/data_agreements/data_agreements.rds")

## ----Save pdfs for Kevin to send -------------------------------------------------------------

# Copy PDFs from each unzipped subfolder into a single pdf folder
pdf_out <- "./DATA/MoveTraitsData/data_agreements/data_agreements_pdf"
dir.create(pdf_out, showWarnings = FALSE)

pdf_files <- list.files(
  "./DATA/MoveTraitsData/data_agreements",
  pattern = "\\.pdf$",
  full.names = TRUE,
  recursive = TRUE
) |>
  Filter(\(f) !grepl("data_agreements_pdf", f), x = _)

file.copy(pdf_files, pdf_out)

# data owners
send_agreement_pdf <- agreements_bind |> 
  dplyr::select(Study_id,Study_Name,User_Surname,User_Firstname,User_Email,PI_Surname,PI_Firstname,PI_Email,PI_Alternate_Email,coauthor_Surname,coauthor_Firstname,coauthor_Email)

write_csv(send_agreement_pdf,"./DATA/MoveTraitsData/data_agreements/data_agreements_pdf/send_agreement_pdf.csv")
