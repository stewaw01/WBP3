library(tidyverse)
library(janitor)
library(purrr)
library(lubridate)

# Make sure to include the correct cohort_start and cohort_end dates.
cohort_start <- as.Date("2024-10-01")
cohort_end <- as.Date("2025-09-30")

# Load in synthetic data approximating an extract from the National SACT data
# Only contains fields required by this script
raw_pts <- read_rds("SyntheticData/Step_01_raw_pts.rds")

n1 <- n_distinct(raw_pts$patient_upi) # 137 patients


## Flag regimen intervals based on gap (if regimen_name differs)
regimen_intervals <- raw_pts |> 
  distinct(patient_upi, regimen_name, appointment_date) |> 
  arrange(patient_upi, regimen_name, appointment_date) |> 
  group_by(patient_upi, regimen_name) |> 
  mutate(
    gap_days = as.numeric(appointment_date - lag(appointment_date, default = first(appointment_date))),
    is_new_episode = if_else(gap_days > as.numeric(168), 1, 0),
    base_ep_id = cumsum(is_new_episode)
  ) |> 
  group_by(patient_upi, regimen_name, base_ep_id) |> 
  summarise(
    reg_start = min(appointment_date),
    .groups = "drop"
  )

## Define "eras" per patient based on regimen starts and re-number
regimen_eras <- regimen_intervals |> 
  distinct(patient_upi, reg_start) |> 
  arrange(patient_upi, reg_start) |> 
  group_by(patient_upi) |> 
  mutate(
    new_regimen_number = row_number(),
    era_start = reg_start,
    ### Era ends the day before a new regimen starts
    era_end = lead(reg_start) - days(1),
    era_end = coalesce(era_end, as.Date("2099-12-31"))
  ) |> 
  ungroup()

## Join these back to the raw_pts dataframe
mapped_pts <- raw_pts |> 
  left_join(regimen_eras, 
            by = join_by(patient_upi, 
                         between(appointment_date, era_start, era_end)))


## Generate new_regimen_name/start_date/latest_date
era_summaries <- mapped_pts |> 
  group_by(patient_upi, new_regimen_number) |> 
  summarise(
    new_regimen_name = paste(sort(unique(regimen_name)), collapse = "_"),
    new_tumour_reporting_group = paste(sort(unique(tumour_reporting_group)), collapse = "_"),
    new_regimen_start_date = as.Date(min(appointment_date)),
    new_regimen_latest_date = as.Date(max(appointment_date)), 
    .groups = "drop"
  )

## Create the final regimen dataframe by joining back
## Re-number the cycles in each new regimen block
regimens <- mapped_pts |> 
  left_join(era_summaries,
            by = c("patient_upi", "new_regimen_number")) |> 
  mutate(local_cycle_regimen_name = paste(local_cycle_number, local_regimen_name, sep = "_")) |> 
  arrange(patient_upi, new_regimen_number, appointment_date, local_cycle_number, local_regimen_name, drug_name) |> 
  group_by(patient_upi, new_regimen_number) |> 
  mutate(new_cycle_number = rle(local_cycle_regimen_name)$lengths %>%
           {rep(seq(length(.)), .)}) |> 
  ungroup() |> 
  mutate(local_cycle_length = as.numeric(local_cycle_length),
         ecog_ps = if_else(is.na(performance_status_cycle), 
                           performance_status_care_episode,
                           performance_status_cycle) |> as.numeric()
  ) |> 
  group_by(patient_upi) |> 
  mutate(cancer_network = first(cancer_network), 
         source = first(source)) |> 
  ungroup() |> 
  select(-era_start, -era_end)


## Remove all these temporary dataframes
rm(regimen_intervals, regimen_eras, mapped_pts, era_summaries)



# Filter to patients who first recieved our comparator within our date range ----
cohort_date_range_upis <- regimens |>
  group_by(patient_upi, new_regimen_name) |> 
  arrange(patient_upi, appointment_date) |> 
  slice(1) |> 
  ungroup() |> 
  filter(appointment_date >= cohort_start 
         & appointment_date <= cohort_end
         & str_detect(drug_name, "ELRANATAMAB")
         & regimen_trial_flag == "N") |> 
  select(patient_upi)

cohort_in_date_range <- regimens |> 
  filter(patient_upi %in% cohort_date_range_upis$patient_upi)

n2 <- n_distinct(cohort_in_date_range$patient_upi) # 69


# Add flags for comparator drugs ----
cohort_med_flag <- cohort_in_date_range |> 
  group_by(patient_upi) |>
  mutate(
    # This was a much nicer way of flagging when we were still including Talquetamab
    # As of now it's a bit overengineered - but useful in future somewhere?
    # Easily extensible by just adding a new "had_x" and the drug name to the vector
    comparator_med = {
      had_elranatamab = any(grepl("ELRANATAMAB", drug_name))
      had_teclistamab = any(grepl("TECLISTAMAB", drug_name))
      
      found_comparator <- c("ELRANATAMAB", "TECLISTAMAB")[c(had_elranatamab, had_teclistamab)]
      
      if (length(found_comparator) > 0){
        paste(found_comparator, collapse = "/")
      } else {
        NA_character_
      }
    }) |> 
  ungroup()


## Get first regimen number of comparator drug for each patient
## This is a workaround for NA's which were appearing when the mutate was included in the above chunk
comparator_first_regimens <- regimens |> 
  group_by(patient_upi, new_regimen_name) |>
  filter(appointment_date >= cohort_start 
         & appointment_date <= cohort_end) |> 
  filter(drug_name == "ELRANATAMAB") |> 
  mutate(comparator_first_regimen = first(new_regimen_number[drug_name == "ELRANATAMAB"])) |> 
  ungroup() |> 
  group_by(patient_upi) |> 
  slice_min(comparator_first_regimen) |> 
  select(patient_upi, comparator_first_regimen) |> 
  distinct()

## Join back on to cohort_med_flag dataframe
cohort_med_flag <- cohort_med_flag |> 
  left_join(comparator_first_regimens)


# Line of Treatment pipeline ----
## Create summarised dataframe with a 'drug_list' created from the new_regimen_name
## components (this allows for neater merging after the lines have been identified)
df_pre_rules <- cohort_med_flag |> 
  group_by(patient_upi, new_regimen_number, new_regimen_name, new_regimen_start_date, new_regimen_latest_date, patient_date_of_death) |>
  summarise() |> 
  mutate(drug_list = strsplit(trimws(as.character(new_regimen_name)), "\\s*[\\+\\_]\\s*")) |> 
  arrange(patient_upi, new_regimen_start_date) |> 
  ungroup()


## List of rules which determine line of treatment groupings as agreed with Rhona
line_of_treatment_rules <- list(
  
  rule_28_days_from_start <- function(current_line, next_regimen) {
    
    if (current_line$line_number > 1) return(FALSE)
    
    days_from_start <- as.numeric(difftime(next_regimen$new_regimen_start_date, current_line$line_start, units = "days"))
    return(days_from_start <= 28)
  },
  
  rule_100_days_overlap <- function(current_line, next_regimen) {
    days_from_start <- as.numeric(difftime(next_regimen$new_regimen_start_date, current_line$line_start, units = "days"))
    has_overlap <- next_regimen$new_regimen_start_date <= current_line$line_end
    return(days_from_start <= 100 && has_overlap)
  },
  
  rule_no_new_drugs_under_180 <- function(current_line, next_regimen) {
    
    next_drugs <- next_regimen$drug_list[[1]]
    if (length(next_drugs) == 0 || all(is.na(next_drugs))) return(FALSE)
    
    gap <- as.numeric(difftime(next_regimen$new_regimen_start_date, current_line$line_end, units = "days"))
    no_new_drugs <- all(next_drugs %in% current_line$drug_list)
    
    return(no_new_drugs && gap < 180)
  }
  
)


## Evaluation function to apply the rules to each patient
evaluate_lines_of_treatment <- function(patient_data, rules) {
  patient_data <- patient_data |> arrange(new_regimen_start_date)
  
  ### If there is no data, return early (this silenced a warning I was getting
  ### but basically doesn't do anything)
  if(nrow(patient_data) == 0) return (patient_data)
  
  ### Split the patient data into individual rows
  regimens_list <- split(patient_data, seq_len(nrow(patient_data)))
  
  init_regimen <- regimens_list[[1]]
  
  ### Create the "state" of the lines of treatment from the first regimen to
  ### give us a starting point
  init_state <- list(
    line_number = 1,
    line_start = init_regimen$new_regimen_start_date,
    line_end = init_regimen$new_regimen_latest_date,
    drug_list = init_regimen$drug_list[[1]]
  )
  
  ### Function to process each regimen, applying the rules and updating the state
  ### If any rule returns TRUE, the regimen is deemed to be part of the same line
  process_regimen <- function(current_state, next_regimen_df) {
    combine_regimens <- any(map_lgl(rules, ~ .x(current_state, next_regimen_df)),
                            na.rm = TRUE)
    
    ### If any rule returns true, update the "state"
    if (combine_regimens) {
      if (is.na(current_state$line_end) && is.na(next_regimen_df$new_regimen_latest_date)) {
        new_end <- NA_Date_
      } else {
        new_end <- max(current_state$line_end, next_regimen_df$new_regimen_latest_date, na.rm = TRUE)
      }
      
      list(
        line_number = current_state$line_number,
        line_start = current_state$line_start,
        line_end = new_end,
        drug_list = unique(c(current_state$drug_list, next_regimen_df$drug_list[[1]]))
      )
    } else {
      list(
        line_number = current_state$line_number + 1,
        line_start = next_regimen_df$new_regimen_start_date,
        line_end = next_regimen_df$new_regimen_latest_date,
        drug_list = next_regimen_df$drug_list[[1]]
      )
    }
  }
  
  if (length(regimens_list) > 1) {
    ### Accumulate is the correct choice, since it uses the prior result as
    ### the argument for the next iteration, which is exactly what we need
    state_history <- accumulate(
      regimens_list[-1],
      process_regimen,
      .init = init_state
    )
  } else {
    state_history <- list(init_state)
  }
  
  patient_data$line_of_treatment <- map_dbl(state_history, "line_number")
  
  return(patient_data)
}


# Apply the "rules engine" function to the prepared dataframe
df_rules_applied <- df_pre_rules |>
  group_by(patient_upi) |> 
  group_modify(~ evaluate_lines_of_treatment(.x, line_of_treatment_rules)) |> 
  ungroup()


# Reformat the data frame to required specification. In the real use case, the
# line_of_treatment and line_of_treatment_name would be joined onto back onto
# the main data for later stages of analysis. These steps are omitted.
final_lines <- df_rules_applied |> 
  group_by(patient_upi, line_of_treatment) |> 
  reframe(
    new_regimen_number = new_regimen_number,
    new_regimen_name = new_regimen_name,
    line_start_date = min(new_regimen_start_date, na.rm = TRUE),
    last_sact_appt = max(new_regimen_latest_date, na.rm = TRUE),
    date_of_death = first(patient_date_of_death),
    line_of_treatment_name = paste(unique(unlist(drug_list)), collapse = "+")
  ) |> 
  group_by(patient_upi) |> 
  mutate(next_line_start = lead(line_start_date),
         day_before_next = next_line_start - days(1),
         
         line_end_date = pmin(
           last_sact_appt,
           day_before_next,
           date_of_death,
           na.rm = TRUE)
  ) |> 
  ungroup() |> 
  select(patient_upi, new_regimen_number, new_regimen_name, line_of_treatment, line_of_treatment_name, line_start_date, line_end_date)

