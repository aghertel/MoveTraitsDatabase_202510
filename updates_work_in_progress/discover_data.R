library(move2)
library(bit64)
library(units)
library(R.utils)
library(dplyr)

# specify account to use in the R session
keyring::key_list()
options("move2_movebank_key_name" = "MoveTraits")

studies <- movebank_download_study_info()
names(studies)

table(studies$study_permission)

collab_studies <- studies %>%
  filter(study_permission == "collaborator" | i_am_collaborator == TRUE)

overview <- collab_studies %>%
  select(id, name, principal_investigator_name, 
         number_of_individuals,
         sensor_type_ids, taxon_ids) %>%
  arrange(name)

View(overview)

sum(as.numeric(overview$number_of_individuals),na.rm=T)
