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

# Only these columns define whether two agreements are substantively the same;
# differences in PI, submission date, etc. are ignored
key_cols <- c("Study_id", "Study_Name", "species", "individual",
              "individual_monthly", "within_individual")

# Collapse duplicates on key columns, keeping the first record
agreements_dedup <- agreements |>
  distinct(across(all_of(key_cols)), .keep_all = TRUE)

# Any Study_id still appearing more than once differs in at least one key column
conflict_ids <- agreements_dedup |>
  count(Study_id) |>
  filter(n > 1) |>
  pull(Study_id)

# Clean: one row per Study_id, conflicts excluded pending review
agreements_clean <- agreements_dedup |>
  filter(!Study_id %in% conflict_ids)

# Conflicts: one row per distinct version of each conflicting Study_id, for manual review
agreements_conflicts <- agreements_dedup |>
  filter(Study_id %in% conflict_ids) |>
  arrange(Study_id)
View(agreements_conflicts)

# Resolve conflicts by keeping the most conservative agreement per Study_id.
# Leniency rank per field (0 = most conservative):
#   "Do not include" = 0, 100 km = 1, 10 km = 2, 1 km = 3, observation = 4
res_cols <- c("species", "individual", "individual_monthly", "within_individual")

leniency_rank <- function(x) {
  case_when(
    str_detect(x, regex("do not include", ignore_case = TRUE)) ~ 0L,
    str_detect(x, regex("observation", ignore_case = TRUE))    ~ 4L,
    str_extract(x, "\\d+") == "100"                            ~ 1L,
    str_extract(x, "\\d+") == "10"                             ~ 2L,
    str_extract(x, "\\d+") == "1"                              ~ 3L
  )
}

# Called once per Study_id (see group_by below), so rows are only ever compared
# with other agreements for the same study. Flags the first row that is at least
# as strict as every other row in all four fields. If no such row exists (strict
# on one field, lenient on another), or a value could not be ranked (NA), no row
# is flagged.
flag_strictest <- function(ranks) {
  r <- as.matrix(ranks)
  out <- rep(FALSE, nrow(r))
  if (anyNA(r)) return(out)
  dominant <- map_lgl(seq_len(nrow(r)), \(i) all(t(r) - r[i, ] >= 0))
  if (any(dominant)) out[which(dominant)[1]] <- TRUE
  out
}

agreements_conflicts_ranked <- agreements_conflicts |>
  mutate(across(all_of(res_cols), leniency_rank, .names = "{.col}_rank")) |>
  group_by(Study_id) |>
  mutate(strictest = flag_strictest(pick(ends_with("_rank")))) |>
  ungroup()

# Conflicts with a clear stricter agreement: keep that row
agreements_conflicts_res <- agreements_conflicts_ranked |>
  filter(strictest) |>
  select(-ends_with("_rank"), -strictest)

# No agreement is stricter on every field: needs manual review
manual_review_ids <- setdiff(conflict_ids, agreements_conflicts_res$Study_id)
agreements_manual_review <- agreements_conflicts |>
  filter(Study_id %in% manual_review_ids)
if (nrow(agreements_manual_review) > 0) View(agreements_manual_review)

agreements_clean <- rbind(agreements_clean, agreements_conflicts_res)

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

nrow(agreements_clean2)

## ----Save pdfs for Kevin to send -------------------------------------------------------------

# # Copy PDFs from each unzipped subfolder into a single pdf folder
# pdf_out <- "./DATA/MoveTraitsData/data_agreements/data_agreements_pdf"
# dir.create(pdf_out, showWarnings = FALSE)
# 
# pdf_files <- list.files(
#   "./DATA/MoveTraitsData/data_agreements",
#   pattern = "\\.pdf$",
#   full.names = TRUE,
#   recursive = TRUE
# ) |>
#   Filter(\(f) !grepl("data_agreements_pdf", f), x = _)
# 
# file.copy(pdf_files, pdf_out)
# 
# # data owners
# send_agreement_pdf <- agreements_bind |> 
#   dplyr::select(Study_id,Study_Name,User_Surname,User_Firstname,User_Email,
#                 PI_Surname,PI_Firstname,PI_Email,PI_Alternate_Email,
#                 coauthor_Surname,coauthor_Firstname,coauthor_Email,
#                 Alternate_Name,Alternate_Email,Alternate_Contact)
# 
# write_csv(send_agreement_pdf,"./DATA/MoveTraitsData/data_agreements/data_agreements_pdf/send_agreement_pdf.csv")
