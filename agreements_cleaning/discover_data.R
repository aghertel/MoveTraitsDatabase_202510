library(move2)
library(bit64)
library(units)
library(R.utils)
library(dplyr)

# specify account to use in the R session
#keyring::key_list()
options("move2_movebank_key_name" = "Anne.Hertel")
open_studies <- movebank_download_study_info()

options("move2_movebank_key_name" = "MoveTraits")
all_studies <- movebank_download_study_info()

table(all_studies$study_permission)
table(all_studies$license_type)


# collab_studies <- studies %>%
#   filter(study_permission == "collaborator" | study_permission == "data_manager" | i_am_collaborator == TRUE)

collab_studies <- all_studies %>%
  filter(study_permission == "collaborator" | study_permission == "data_manager" | i_am_collaborator == TRUE)

overview <- collab_studies %>%
  select(id, name, principal_investigator_name, 
         number_of_individuals,
         sensor_type_ids, taxon_ids) %>%
  arrange(name)

View(overview)
nrow(overview)

studies$sensor_type_ids |> strsplit(",") |> unlist() |> table()

overview_GPS <- overview |> 
  filter(grepl("\\bGPS\\b", sensor_type_ids))

sum(as.numeric(overview_GPS$number_of_individuals),na.rm=T)
nrow(overview_GPS)

#4012435577_4023828432
class(Ids_toDo) <- "integer64"
Ids_toDo %in% 4012435577

length(unique(overview$principal_investigator_name))


agreements_clean <- readRDS("./DATA/MoveTraitsData/data_agreements/data_agreements.rds")

## ---- Compare access (overview_GPS) with data sharing agreements ----

# Movebank ids are integer64, agreement ids are numeric; ids are < 2^53 so
# converting to numeric is lossless
access <- overview_GPS |>
  mutate(id = as.numeric(id))

# One row per study, in case an id still has multiple agreement rows
agreement_studies <- agreements_clean |>
  distinct(Study_id, .keep_all = TRUE)

# Access and a completed agreement
agreements_complete <- access |>
  inner_join(
    select(agreement_studies, Study_id, agreement_study_name = Study_Name),
    by = c("id" = "Study_id")
  )

# Access but no agreement
agreements_missing <- access |>
  anti_join(agreement_studies, by = c("id" = "Study_id"))

# Agreement but no access (expected to be empty); reason column says why
no_access <- agreement_studies |>
  select(Study_id, Study_Name) |>
  anti_join(access, by = c("Study_id" = "id")) |>
  left_join(
    studies |>
      transmute(
        Study_id = as.numeric(id), movebank_name = name,
        study_permission, i_am_collaborator, sensor_type_ids
      ),
    by = "Study_id"
  ) |>
  mutate(
    reason = case_when(
      is.na(movebank_name)                       ~ "id not found in Movebank study list",
      !grepl("\\bGPS\\b", sensor_type_ids)       ~ "no GPS sensor listed",
      TRUE                                       ~ "GPS study, but no collaborator/data manager access"
    )
  )

tibble(
  group = c("access + agreement", "access, no agreement", "agreement, no access"),
  n_studies = c(nrow(agreements_complete), nrow(agreements_missing), nrow(no_access)),
  n_individuals = c(
    sum(as.numeric(agreements_complete$number_of_individuals), na.rm = TRUE),
    sum(as.numeric(agreements_missing$number_of_individuals), na.rm = TRUE),
    NA
  )
)

count(no_access, reason)

# Individuals in GPS studies with access and a sharing agreement
sum(as.numeric(agreements_complete$number_of_individuals), na.rm = TRUE)

