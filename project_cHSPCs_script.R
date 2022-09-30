#############################    PROJECT cHSPCs    #############################
# Latest update: 08/07/2022
# @Author: Giulia Pais
################################################################################
library(ISAnalytics)
library(dplyr)
library(data.table)
library(openxlsx)

date <- "20220708"
date_month <- "202207"
project_name <- "cHSPCs"
file_prefix <- paste(date, project_name, sep = "_")

# FOLDER STRUCTURE ----
project_folder <- fs::path(fs::path_wd(), date_month)
report_folder <- fs::path(project_folder, "report")
sharing_fold <- fs::path(project_folder, "sharing")
plots_fold <- fs::path(project_folder, "plots")

for (fold in c(project_folder, report_folder, sharing_fold, plots_fold)) {
  fs::dir_create(fold) # does nothing if folder already exists
}

# ADDITIONAL FILES ----
selected_samples_file <- "./20210628_WAS_cHSCPs_project_association_file_complete_aggiornato.txt"
selected_samples <- readr::read_tsv(file = selected_samples_file)

blood_lineages_file <- "./blood_lineages_update_FINAL.xlsx"
blood_lineages_tab <- openxlsx::read.xlsx(blood_lineages_file)

# AF IMPORT ----
af_path <- fs::path("./WAS.assofile.extract.July2022.txt")
af <- import_association_file(af_path, root = "/mnt/OracleResults/",
                              import_iss = TRUE, 
                              report_path = report_folder)

# OUTLIER DETECTION ----
outliers_removed <- outlier_filter(
  af, report_path = report_folder
)

### Removing also sample WASHE02 - TP 90 - PB - CD3 as per request
outliers_removed <- outliers_removed %>%
  filter(!(SubjectID == "WASHE02" & TimePoint == "0090" & Tissue == "PB" & 
             CellMarker == "CD3"))

# MATRIX IMPORT ----
matrices <- import_parallel_Vispa2Matrices(
  outliers_removed, workers = 5, report_path = report_folder, 
  mode = "AUTO", quantification_type = c("seqCount", "fragmentEstimate",
                                         "ShsCount")
)
### * Keeping only samples with associated metadata
matrices <- matrices %>%
  semi_join(outliers_removed, by = "CompleteAmplificationID")

# RECALIBRATION ----
recalibr <- compute_near_integrations(matrices,
                                      value_columns = c("seqCount", "fragmentEstimate",
                                                        "ShsCount"),
                                      max_workers = 15,
                                      map_as_file = TRUE,
                                      file_path = report_folder)

# COLLISIONS ----
coll <- remove_collisions(recalibr, outliers_removed, 
                          quant_cols = c(seqCount = "seqCount", 
                                         fragmentEstimate = "fragmentEstimate",
                                         ShsCount = "ShsCount"),
                          max_workers = 15,
                          report_path = report_folder)


# AGGREGATION ----
## to detect duplicates between PCRmethods
agg_methods <- aggregate_values_by_key(coll,
                                       outliers_removed,
                                       value_cols = c("seqCount", 
                                                      "fragmentEstimate", 
                                                      "ShsCount"), 
                                       key = c("SubjectID", "CellMarker",
                                               "Tissue", "TimePoint",
                                               "PCRMethod"))
agg_methods <- agg_methods %>%
  semi_join(selected_samples, by = c("SubjectID", "CellMarker",
                                     "Tissue", "TimePoint"))

detected_dupl <- agg_methods %>%
  group_by(SubjectID, CellMarker, Tissue, TimePoint) %>%
  summarise(n_methods = n_distinct(PCRMethod), .groups = "drop") %>%
  filter(n_methods > 1)

is_counts_dupl <- agg_methods %>%
  semi_join(detected_dupl, by = c("SubjectID", "CellMarker",
                                  "Tissue", "TimePoint")) %>%
  group_by(SubjectID, CellMarker,
             Tissue, TimePoint,
             PCRMethod) %>%
  summarise(nIS = n_distinct(chr, integration_locus, strand),
            seqCount = sum(seqCount_sum),
            fragmentEstimate = sum(fragmentEstimate_sum, na.rm = TRUE),
            ShsCount = sum(ShsCount_sum, na.rm = TRUE), .groups = "drop")

is_counts_dupl %>%
  readr::write_tsv(file = fs::path(project_folder, paste(file_prefix, 
                                "_is-counts-for-duplicate-removal.tsv"))
  )

# REMOVAL OF PCR METHODS DUPLICATES ----
# * As per request keep LAM samples
final_af <- outliers_removed %>%
  group_by(across(all_of(c("SubjectID", "CellMarker",
             "Tissue", "TimePoint")))) %>%
  group_modify(.f = ~ {
    if (length(unique(.x$PCRMethod)) == 2) {
      return(
        .x %>%
          filter(PCRMethod == "LAM-PCR")
      )
    } else {
      return(.x)
    }
  }) %>%
  ungroup()

## * Complete AF with blood lineages info and time point type info
final_af <- final_af %>%
  left_join(blood_lineages_tab, by = "CellMarker") %>%
  left_join(selected_samples %>% 
              distinct(CompleteAmplificationID, TimePointType),
            by = "CompleteAmplificationID")

final_af %>%
  readr::write_tsv(file = fs::path(project_folder, paste(file_prefix, "final_af.tsv", sep = "_")))
final_af <- final_af %>%
  rename(Cell.type = "CellType")

## * Matrix with duplicates removed
final_matrix <- coll %>%
  semi_join(final_af, by = "CompleteAmplificationID")


# ----- SHARING ANALYSES -------------------------------------------------------
## General utils for xlsx format
perc_format <- createStyle(numFmt = "0.00\\%")
write_section <- function(tbl, tables, wb, sheet, start_row) {
  current_tab <- do.call(`$`, args = list(tables, tbl))
  if (!is.null(current_tab) && !purrr::is_empty(current_tab)) {
    writeData(wb, sheet, startRow = start_row, x = paste("TABLE", tbl))
    mergeCells(wb, sheet, cols = c(1:ncol(current_tab)), rows = start_row)
    writeDataTable(wb, sheet, startRow = start_row + 1,
                   x = current_tab)
    addStyle(wb = wb, sheet = sheet, style = perc_format,
             cols = which(stringr::str_starts(
               colnames(current_tab), "on_")),
             rows = start_row + 1: nrow(current_tab) + 1, 
             gridExpand = TRUE
    )
    end_row <- start_row + nrow(current_tab) + 3
    return(end_row)
  } else {
    return(start_row)
  }
}

## SHARING PT1 ----
pt1_patients <- c("WAS1009", "WASHE02", "WASCUP04")
data_selection_1 <- final_matrix %>%
  left_join(final_af, by = "CompleteAmplificationID") %>%
  filter(SubjectID %in% pt1_patients)
setDT(data_selection_1)

## * Excluding some samples as per request
## Excluding BM-CD34-WASHE02-0090-DX/SX fram early tps
data_selection_1[Tissue == "BM" & CellMarker == "CD34" & 
                   SubjectID == "WASHE02" & 
                   TimePoint == "0090" & 
                   (AddedField2 == "DX" | AddedField2 == "SX"), 
                 TimePointType := NA_character_]
## Excluding BM-CD34-WASCUP04-0030 fram early tps
data_selection_1[Tissue == "BM" & 
                   CellMarker == "CD34" & SubjectID == "WASCUP04" & 
                   TimePoint == "0030", TimePointType := NA_character_]
## Excluding BM-CFCBULK-WAS1009-0090 fram early tps
data_selection_1[Tissue == "BM" & 
                   TimePoint == "0090" & CellMarker == "CFCBULK", 
                 TimePointType := NA_character_]
## Excluding WAS1009 CD34 amplified
data_selection_1[Tissue == "BM" & TimePoint == "0030" &
                   CellMarker == "CD34" & AddedField3 == "AMP", 
                 TimePointType := NA_character_]

### Define a function for analysis table 1 - process and save results
pt_1_sharing <- function(df, tps = c("early", "late")) {
  tbl_pt1 <- purrr::map_df(tps, function(tp) {
    left_tbl_a <- df[Tissue == "PB" & HematoLineage == "CFCBULK" &
                       TimePointType == tp]
    right_tbl_a_1 <- df[Tissue == "BM" & AddedField2 %in% c("DX", "SX") &
                          HematoLineage %in% c("CFCBULK", "CD34", "PrimHSPC", 
                                               "MyeloHSPC", 
                                               "LymphoHSPC", "ErythroHSPC") &
                          TimePointType == tp]
    tbl_a_1 <- if (nrow(left_tbl_a) == 0 | nrow(right_tbl_a_1) == 0) {
      NULL
    } else {
      is_sharing(left_tbl_a, right_tbl_a_1, 
                 group_keys = list(g1 = c("Tissue", "HematoLineage", "TimePointType"),
                                   g2 = c("Tissue", "AddedField2", "HematoLineage", 
                                          "TimePointType")), 
                 minimal = TRUE)
    }
    
    right_tbl_a_2 <- df[Tissue == "BM" & 
                          HematoLineage %in% c("CFCBULK", "CD34", "PrimHSPC", 
                                               "MyeloHSPC", "LymphoHSPC", 
                                               "ErythroHSPC") &
                          Cell.type %in% c("Other", "CD34", "HSPC") &
                          TimePointType == tp]
    tbl_a_2 <- if (nrow(left_tbl_a) == 0 | nrow(right_tbl_a_2) == 0) {
      NULL
    } else {
      is_sharing(left_tbl_a, right_tbl_a_2, 
                 group_keys = list(g1 = c("Tissue", "HematoLineage", "TimePointType"),
                                   g2 = c("Tissue", "HematoLineage",  "Cell.type",
                                          "TimePointType")), 
                 minimal = TRUE)
    }
    
    right_tbl_a_3 <- df[Tissue == "BM" & Cell.type == "HSPC" & 
                          TimePointType == tp]
    tbl_a_3 <- if (nrow(left_tbl_a) == 0 | nrow(right_tbl_a_3) == 0) {
      NULL
    } else {
      is_sharing(left_tbl_a, right_tbl_a_3, 
                 group_keys = list(g1 = c("Tissue", "HematoLineage", "TimePointType"),
                                   g2 = c("Tissue", "Cell.type",
                                          "TimePointType")), 
                 minimal = TRUE)
    }
    
    right_tbl_a_4 <- df[Tissue == "BM" & AddedField2 %in% c("DX", "SX") &
                          Cell.type == "HSPC" & 
                          TimePointType == tp]
    
    tbl_a_4 <- if (nrow(left_tbl_a) == 0 | nrow(right_tbl_a_4) == 0) {
      NULL
    } else {
      is_sharing(left_tbl_a, right_tbl_a_4, 
                 group_keys = list(g1 = c("Tissue", "HematoLineage", "TimePointType"),
                                   g2 = c("Tissue", "AddedField2", "Cell.type",
                                          "TimePointType")), 
                 minimal = TRUE)
    }
    rbindlist(list(tbl_a_1, tbl_a_2, tbl_a_3, tbl_a_4))
  })
  
  tbl_pt2 <- purrr::map_df(tps, function(tp) {
    left_tbl_b <- df[Tissue == "BM" & HematoLineage == "CFCBULK" & 
                       TimePointType == tp]
    right_tbl_b_1 <- df[Tissue == "BM" & HematoLineage == "CD34" & 
                          TimePointType == tp]
    
    right_tbl_b_2 <- df[Tissue == "BM" & Cell.type == "HSPC" &
                          TimePointType == tp]
    
    tbl_b_1 <- if (nrow(left_tbl_b) == 0 | nrow(right_tbl_b_1) == 0) {
      NULL
    } else {
      is_sharing(left_tbl_b, right_tbl_b_1,
                 group_keys = list(g1 = c("Tissue", "HematoLineage", "TimePointType"),
                                   g2 = c("Tissue", "HematoLineage",
                                          "TimePointType")), 
                 minimal = TRUE)
    }
    
    tbl_b_2 <- if (nrow(left_tbl_b) == 0 | nrow(right_tbl_b_2) == 0) {
      NULL
    } else {
      is_sharing(left_tbl_b, right_tbl_b_1,
                 group_keys = list(g1 = c("Tissue", "HematoLineage", "TimePointType"),
                                   g2 = c("Tissue", "Cell.type",
                                          "TimePointType")), 
                 minimal = TRUE)
    }
    
    left_tbl_c <- df[Tissue == "BM" & AddedField2 %in% c("DX", "SX") &
                       HematoLineage == "CFCBULK" & 
                       TimePointType == tp]
    
    right_tbl_c_1 <- df[Tissue == "BM" & HematoLineage == "CD34" & 
                          AddedField2 %in% c("DX", "SX") &
                          TimePointType == tp]
    
    right_tbl_c_2 <- df[Tissue == "BM" & AddedField2 %in% c("DX", "SX") &
                          Cell.type == "HSPC" & TimePointType == tp]
    
    tbl_c_1 <- if (nrow(left_tbl_c) == 0 | nrow(right_tbl_c_1) == 0) {
      NULL
    } else {
      is_sharing(left_tbl_c, right_tbl_c_1,
                 group_keys = list(g1 = c("Tissue", "AddedField2",
                                          "HematoLineage", "TimePointType"),
                                   g2 = c("Tissue", "HematoLineage", "AddedField2",
                                          "TimePointType")), 
                 minimal = TRUE)
    }
    
    tbl_c_2 <- if (nrow(left_tbl_c) == 0 | nrow(right_tbl_c_2) == 0) {
      NULL
    } else {
      is_sharing(left_tbl_c, right_tbl_c_2,
                 group_keys = list(g1 = c("Tissue", "AddedField2",
                                          "HematoLineage", "TimePointType"),
                                   g2 = c("Tissue", "AddedField2", "Cell.type",
                                          "TimePointType")), 
                 minimal = TRUE)
    }
    
    left_tbl_d <- df[Tissue == "BM" & Cell.type == "HSPC" & 
                       TimePointType == tp]
    
    tbl_d_1 <- if (nrow(left_tbl_d) == 0 | nrow(right_tbl_b_1) == 0) {
      NULL
    } else {
      is_sharing(left_tbl_d, right_tbl_b_1,
                 group_keys = list(g1 = c("Tissue", "Cell.type", "TimePointType"),
                                   g2 = c("Tissue", "HematoLineage",
                                          "TimePointType")), 
                 minimal = TRUE)
    }
    
    left_tbl_e <- df[Tissue == "BM" & AddedField2 %in% c("DX", "SX") &
                       Cell.type == "HSPC" & 
                       TimePointType == tp]
    tbl_e_1 <- if (nrow(left_tbl_e) == 0 | nrow(right_tbl_c_1) == 0) {
      NULL
    } else {
      is_sharing(left_tbl_e, right_tbl_c_1,
                 group_keys = list(g1 = c("Tissue", "Cell.type", "TimePointType"),
                                   g2 = c("Tissue", "AddedField2", "HematoLineage",
                                          "TimePointType")), 
                 minimal = TRUE)
    }
    rbindlist(list(
      tbl_b_1, tbl_b_2, tbl_c_1, tbl_c_2, tbl_d_1, tbl_e_1
    ))
  })
  
  rbindlist(list(tbl_pt1, tbl_pt2))
}

pt_1_sharing_results <- pt_1_sharing(data_selection_1)
pt_1_sharing_results_per_patient <- purrr::map(pt1_patients, ~ {
  pt_1_sharing(data_selection_1[SubjectID == .x])
}) %>%
  purrr::set_names(pt1_patients)

wb1 <- createWorkbook()
addWorksheet(wb1, "Overall")
write_section("A", tables = list(A = pt_1_sharing_results), wb = wb1, 
              sheet = "Overall", start_row = 1)

for (pt in names(pt_1_sharing_results_per_patient)) {
  addWorksheet(wb1, pt)
  write_section("A", tables = list(A = pt_1_sharing_results[[pt]]), wb = wb1, 
                sheet = pt, start_row = 1)
}
saveWorkbook(wb1, file = fs::path(sharing_fold, paste(file_prefix, "sharing_analyses_pt1.xlsx", 
                                                      sep = "_")),
             overwrite = TRUE)

## SHARING PT2 ----
### Define functions for the various analyses tables
pt_2A_sharing <- function(df) {
  
  lineages_comp_early <- function(lineage, side) {
    left_tbl <- df[Tissue == "BM" & AddedField2 == side & 
                      HematoLineage == lineage & TimePointType == "early"]
    right_tbl <- if (side == "DX") {
      df[(Tissue == "BM" & AddedField2 == "SX" & HematoLineage == lineage &
            TimePointType == "early") | 
           (Tissue == "BM" & AddedField2 %in% c("DX", "SX") & 
              HematoLineage == lineage &
              TimePointType == "late")]
    } else {
      df[Tissue == "BM" & AddedField2 %in% c("DX", "SX") & 
           HematoLineage == lineage & TimePointType == "late"]
    }
    
    tbl <- if (nrow(left_tbl) == 0 | nrow(right_tbl) == 0) {
      NULL
    } else {
      is_sharing(left_tbl, right_tbl,
                 group_key = c("Tissue", "AddedField2",
                               "HematoLineage", "TimePointType"))
    }
    return(tbl)
  }
  
  tbl_pt1 <- purrr::map2_df(rep(c("CFCBULK", "CD34", "PrimHSPC", 
                                  "MyeloHSPC", 
                                  "LymphoHSPC", "ErythroHSPC"), each = 2),
                            rep(c("DX", "SX"), 6),
                            ~ lineages_comp_early(.x, .y))
  
  lineages_comp_late <- function(lineage) {
    left_tbl <- df[Tissue == "BM" & AddedField2 == "DX" & 
                     HematoLineage == lineage & TimePointType == "late"]
    right_tbl <- df[Tissue == "BM" & AddedField2 == "SX" & 
           HematoLineage == lineage & TimePointType == "late"]
    
    tbl <- if (nrow(left_tbl) == 0 | nrow(right_tbl) == 0) {
      NULL
    } else {
      is_sharing(left_tbl, right_tbl,
                 group_key = c("Tissue", "AddedField2",
                               "HematoLineage", "TimePointType"))
    }
    return(tbl)
  }
  
  tbl_pt2 <- purrr::map_df(c("CFCBULK", "CD34", "PrimHSPC", 
                             "MyeloHSPC", 
                             "LymphoHSPC", "ErythroHSPC"), 
                           ~ lineages_comp_late(.x))
  
  left_tbl_pt3_early_1 <- df[Tissue == "BM" & AddedField2 == "DX" &
                               Cell.type == "HSPC" & TimePointType == "early"]
  left_tbl_pt3_early_2 <- df[Tissue == "BM" & AddedField2 == "SX" &
                               Cell.type == "HSPC" & TimePointType == "early"]
  
  right_tbl_pt3_early_1 <- df[(Tissue == "BM" & 
                                 AddedField2 == "SX" & Cell.type == "HSPC" &
                                 TimePointType == "early") | 
                                (Tissue == "BM" & 
                                   AddedField2 %in% c("DX", "SX") & 
                                   Cell.type == "HSPC" &
                                   TimePointType == "late")]
  right_tbl_pt3_early_2 <- df[Tissue == "BM" & AddedField2 %in% c("DX", "SX") &
                                Cell.type == "HSPC" & TimePointType == "late"]
  tbl_pt3_1 <- if (nrow(left_tbl_pt3_early_1) == 0 | 
                 nrow(right_tbl_pt3_early_1) == 0) {
    NULL
  } else {
    is_sharing(left_tbl_pt3_early_1, right_tbl_pt3_early_1,
               group_key = c("Tissue", "AddedField2", "Cell.type", 
                             "TimePointType"),
               minimal = TRUE)
  }
  tbl_pt3_2 <- if (nrow(left_tbl_pt3_early_2) == 0 | 
                   nrow(right_tbl_pt3_early_2) == 0) {
    NULL
  } else {
    is_sharing(left_tbl_pt3_early_2, right_tbl_pt3_early_2,
               group_key = c("Tissue", "AddedField2", "Cell.type", 
                             "TimePointType"),
               minimal = TRUE)
  }
  left_tbl_pt3_late <- df[Tissue == "BM" & AddedField2 == "DX" & 
                            Cell.type == "HSPC" & TimePointType == "late"]
  right_tbl_pt3_late <- df[Tissue == "BM" & AddedField2 == "SX" & 
                            Cell.type == "HSPC" & TimePointType == "late"]
  tbl_pt3_3 <- if (nrow(left_tbl_pt3_late) == 0 | 
                   nrow(right_tbl_pt3_late) == 0) {
    NULL
  } else {
    is_sharing(left_tbl_pt3_late, right_tbl_pt3_late,
               group_key = c("Tissue", "AddedField2", "Cell.type", 
                             "TimePointType"),
               minimal = TRUE)
  }
  
  left_tbl_pt4 <- df[Tissue == "PB" & HematoLineage == "CFCBULK" & 
                       TimePointType == "early"]
  right_tbl_pt4 <- df[Tissue == "PB" & HematoLineage == "CFCBULK" & 
                       TimePointType == "late"]
  tbl_pt4 <- if (nrow(left_tbl_pt4) == 0 | 
                   nrow(right_tbl_pt4) == 0) {
    NULL
  } else {
    is_sharing(left_tbl_pt4, right_tbl_pt4,
               group_key = c("Tissue", "AddedField2", "Cell.type", 
                             "TimePointType"),
               minimal = TRUE)
  }
  
  rbindlist(list(tbl_pt1, tbl_pt2, tbl_pt3_1, tbl_pt3_2, tbl_pt3_3, tbl_pt4))
}

pt_2B_sharing <- function(df) {
  right_tbl_a <- df[Tissue == "PB" & 
                      HematoLineage %in% c("Myeloid", "Lymphoid", 
                                           "Erythroid") &
                      TimePointType == "early"]
  right_tbl_b <- df[Tissue == "PB" & 
                      HematoLineage %in% c("Myeloid", "Lymphoid", 
                                           "Erythroid") &
                      TimePointType == "late"]
  left_tbl_a <- df[Tissue %in% c("PB", "BM") & HematoLineage == "CFCBULK" &
                     TimePointType == "early"]
  left_tbl_b <- df[Tissue %in% c("PB", "BM") & HematoLineage == "CFCBULK" &
                     TimePointType == "late"]
  
  tbl_a <- if (nrow(right_tbl_a) == 0 | nrow(left_tbl_a) == 0) {
    NULL
  } else {
    is_sharing(left_tbl_a, right_tbl_a,
               group_key = c("Tissue", "HematoLineage", "TimePointType"),
               minimal = TRUE)
  }
  tbl_b <- if (nrow(right_tbl_b) == 0 | nrow(left_tbl_b) == 0) {
    NULL
  } else {
    is_sharing(left_tbl_b, right_tbl_b,
               group_key = c("Tissue", "HematoLineage", "TimePointType"),
               minimal = TRUE)
  }
  
  left_tbl_c <- df[Tissue == "BM" & 
                     AddedField2 %in% c("DX", "SX") &
                      HematoLineage == "CFCBULK" &
                      TimePointType == "early"]
  left_tbl_d <- df[Tissue == "BM" & 
                     AddedField2 %in% c("DX", "SX") &
                   HematoLineage == "CFCBULK" &
                     TimePointType == "late"]
  tbl_c <- if (nrow(right_tbl_a) == 0 | nrow(left_tbl_c) == 0) {
    NULL
  } else {
    is_sharing(left_tbl_c, right_tbl_a,
               group_keys = list(
                 g1 = c("Tissue", "AddedField2", "HematoLineage", "TimePointType"),
                 g2 = c("Tissue", "HematoLineage", "TimePointType")
               ),
               minimal = TRUE)
  }
  tbl_d <- if (nrow(right_tbl_b) == 0 | nrow(left_tbl_d) == 0) {
    NULL
  } else {
    is_sharing(left_tbl_d, right_tbl_b,
               group_keys = list(
                 g1 = c("Tissue", "AddedField2", "HematoLineage", "TimePointType"),
                 g2 = c("Tissue", "HematoLineage", "TimePointType")
               ),
               minimal = TRUE)
  }
  
  right_tbl_e <- df[Tissue == "PB" & HematoLineage %in% c("Myeloid", "Lymphoid") &
                      Cell.type %in% c("Granulo", "Mono", "B", "NK", "T") &
                      TimePointType == "early"]
  right_tbl_f <- df[Tissue == "PB" & HematoLineage %in% c("Myeloid", "Lymphoid") &
                      Cell.type %in% c("Granulo", "Mono", "B", "NK", "T") &
                      TimePointType == "late"]
  tbl_e_1 <- if (nrow(right_tbl_e) == 0 | nrow(left_tbl_a) == 0) {
    NULL
  } else {
    is_sharing(left_tbl_a, right_tbl_e,
               group_keys = list(
                 g1 = c("Tissue", "HematoLineage", "TimePointType"),
                 g2 = c("Tissue", "HematoLineage", "Cell.type","TimePointType")
               ),
               minimal = TRUE)
  }
  tbl_e_2 <- if (nrow(right_tbl_e) == 0 | nrow(left_tbl_c) == 0) {
    NULL
  } else {
    is_sharing(left_tbl_c, right_tbl_e,
               group_keys = list(
                 g1 = c("Tissue", "AddedField2", "HematoLineage", "TimePointType"),
                 g2 = c("Tissue", "HematoLineage", "Cell.type","TimePointType")
               ),
               minimal = TRUE)
  }
  tbl_f_1 <- if (nrow(right_tbl_f) == 0 | nrow(left_tbl_b) == 0) {
    NULL
  } else {
    is_sharing(left_tbl_b, right_tbl_f,
               group_keys = list(
                 g1 = c("Tissue", "HematoLineage", "TimePointType"),
                 g2 = c("Tissue", "HematoLineage", "Cell.type","TimePointType")
               ),
               minimal = TRUE)
  }
  tbl_f_2 <- if (nrow(right_tbl_f) == 0 | nrow(left_tbl_d) == 0) {
    NULL
  } else {
    is_sharing(left_tbl_d, right_tbl_f,
               group_keys = list(
                 g1 = c("Tissue", "AddedField2", "HematoLineage", "TimePointType"),
                 g2 = c("Tissue", "HematoLineage", "Cell.type","TimePointType")
               ),
               minimal = TRUE)
  }
  
  rbindlist(list(tbl_a, tbl_b, tbl_c, tbl_d, tbl_e_1, tbl_e_2, tbl_f_1, tbl_f_2))
}

pt_2C_sharing <- function(df) {
  # 3-way comparisons
  g1_1 <- df[Tissue == "PB" & HematoLineage == "CFCBULK" & TimePointType == "early"]
  g1_2 <- df[Tissue == "PB" & HematoLineage == "CFCBULK" & TimePointType == "late"]
  
  g2_1 <- df[Tissue == "BM" & AddedField2 == "DX" & 
               HematoLineage %in% c("CFCBULK", "CD34") & TimePointType == "early"]
  g2_2 <- df[Tissue == "BM" & AddedField2 == "DX" & 
               HematoLineage %in% c("CFCBULK", "CD34") & TimePointType == "late"]
  g2_3 <- df[Tissue == "BM" & AddedField2 == "DX" & 
               Cell.type == "HSPC" & TimePointType == "early"]
  g2_4 <- df[Tissue == "BM" & AddedField2 == "DX" & 
               Cell.type == "HSPC" & TimePointType == "late"]
  
  g3_1 <- df[Tissue == "BM" & AddedField2 == "SX" & 
               HematoLineage %in% c("CFCBULK", "CD34") & TimePointType == "early"]
  g3_2 <- df[Tissue == "BM" & AddedField2 == "SX" & 
               HematoLineage %in% c("CFCBULK", "CD34") & TimePointType == "late"]
  
  g3_3 <- df[Tissue == "BM" & AddedField2 == "SX" & 
               Cell.type == "HSPC" & TimePointType == "early"]
  g3_4 <- df[Tissue == "BM" & AddedField2 == "SX" & 
               Cell.type == "HSPC" & TimePointType == "late"]
  
  tbl_1 <- if (nrow(g1_1) == 0 | nrow(g2_1) == 0 | nrow(g3_1) == 0) {
    NULL
  } else {
    is_sharing(g1_1, g2_1, g3_1,
               group_keys = list(
                 g1 = c("Tissue", "HematoLineage", "TimePointType"),
                 g2 = c("Tissue", "AddedField2", "HematoLineage", "TimePointType"),
                 g3 = c("Tissue", "AddedField2", "HematoLineage", "TimePointType")
               ))
  }
  tbl_2 <- if (nrow(g1_2) == 0 | nrow(g2_2) == 0 | nrow(g3_2) == 0) {
    NULL
  } else {
    is_sharing(g1_2, g2_2, g3_2,
               group_keys = list(
                 g1 = c("Tissue", "HematoLineage", "TimePointType"),
                 g2 = c("Tissue", "AddedField2", "HematoLineage", "TimePointType"),
                 g3 = c("Tissue", "AddedField2", "HematoLineage", "TimePointType")
               ))
  }
  tbl_3 <- if (nrow(g1_1) == 0 | nrow(g2_3) == 0 | nrow(g3_3) == 0) {
    NULL
  } else {
    is_sharing(g1_1, g2_3, g3_3,
               group_keys = list(
                 g1 = c("Tissue", "HematoLineage", "TimePointType"),
                 g2 = c("Tissue", "AddedField2", "Cell.type", "TimePointType"),
                 g3 = c("Tissue", "AddedField2", "Cell.type", "TimePointType")
               ))
  }
  tbl_4 <- if (nrow(g1_2) == 0 | nrow(g2_4) == 0 | nrow(g3_4) == 0) {
    NULL
  } else {
    is_sharing(g1_1, g2_3, g3_3,
               group_keys = list(
                 g1 = c("Tissue", "HematoLineage", "TimePointType"),
                 g2 = c("Tissue", "AddedField2", "Cell.type", "TimePointType"),
                 g3 = c("Tissue", "AddedField2", "Cell.type", "TimePointType")
               ))
  }
  
  rbindlist(list(tbl_1, tbl_2, tbl_3, tbl_4))
}

pt_2D_sharing <- function(df) {
  # 4 way comparisons
  g1 <- df[Tissue == "PB" & HematoLineage == "CFCBULK" & TimePointType == "early"]
  g2 <- df[Tissue == "PB" & HematoLineage == "CFCBULK" & TimePointType == "late"]
  g3_1 <- df[Tissue == "BM" & HematoLineage %in% c("CFCBULK", "CD34") &
               TimePointType == "early"]
  g3_2 <- df[Tissue == "BM" & Cell.type == "HSPC" &
               TimePointType == "early"]
  g4_1 <- df[Tissue == "BM" & HematoLineage %in% c("CFCBULK", "CD34") &
               TimePointType == "late"]
  g4_2 <- df[Tissue == "BM" & Cell.type == "HSPC" &
               TimePointType == "late"]
  
  tbl_1 <- if (nrow(g1) == 0 | nrow(g2) == 0 | nrow(g3_1) == 0 | nrow(g4_1) == 0) {
    NULL
  } else {
    is_sharing(g1, g2, g3_1, g4_1,
               group_key = c("Tissue", "HematoLineage", "TimePointType"),
               minimal = TRUE
               )
  }
  tbl_2 <- if (nrow(g1) == 0 | nrow(g2) == 0 | nrow(g3_1) == 0 | nrow(g4_1) == 0) {
    NULL
  } else {
    is_sharing(g1, g2, g3_2, g4_2,
               group_keys = list(
                 g1 = c("Tissue", "HematoLineage", "TimePointType"),
                 g2 = c("Tissue", "HematoLineage", "TimePointType"),
                 g3 = c("Tissue", "Cell.type", "TimePointType"),
                 g4 = c("Tissue", "Cell.type", "TimePointType")
                 ),
               minimal = TRUE
    )
  }
  rbindlist(list(tbl_1, tbl_2))
}

pt_2E_sharing <- function(df) {
  # 6 way comparison
  g1 <- df[Tissue == "PB" & HematoLineage == "CFCBULK" & TimePointType == "early"]
  g2 <- df[Tissue == "PB" & HematoLineage == "CFCBULK" & TimePointType == "late"]
  g3 <- df[Tissue == "BM" & AddedField2 == "DX" &
             HematoLineage == "CFCBULK" & TimePointType == "early"]
  g4_1 <- df[Tissue == "BM" & AddedField2 == "SX" &
             HematoLineage %in% c("CFCBULK", "CD34") & TimePointType == "early"]
  g4_2 <- df[Tissue == "BM" & AddedField2 == "SX" &
               Cell.type == "HSPC" & TimePointType == "early"]
  g5_1 <- df[Tissue == "BM" & AddedField2 == "DX" &
               HematoLineage %in% c("CFCBULK", "CD34") & TimePointType == "late"]
  g5_2 <- df[Tissue == "BM" & AddedField2 == "DX" &
               Cell.type == "HSPC" & TimePointType == "late"]
  g6_1 <- df[Tissue == "BM" & AddedField2 == "SX" &
               HematoLineage %in% c("CFCBULK", "CD34") & TimePointType == "late"]
  g6_2 <- df[Tissue == "BM" & AddedField2 == "SX" &
               Cell.type == "HSPC" & TimePointType == "late"]
  
  tbl_1 <- if (nrow(g1) == 0 | nrow(g2) == 0 | nrow(g3) == 0 |
               nrow(g4_1) == 0 | nrow(g5_1) == 0 | nrow(g6_1) == 0) {
    NULL
  } else {
    is_sharing(g1, g2, g3, g4_1, g5_1, g6_1,
               group_keys = list(
                 g1 = c("Tissue", "HematoLineage", "TimePointType"),
                 g2 = c("Tissue", "HematoLineage", "TimePointType"),
                 g3 = c("Tissue", "AddedField2", "HematoLineage", "TimePointType"),
                 g4 = c("Tissue", "AddedField2", "HematoLineage", "TimePointType"),
                 g5 = c("Tissue", "AddedField2", "HematoLineage", "TimePointType"),
                 g6 = c("Tissue", "AddedField2", "HematoLineage", "TimePointType")
               ))
  }
  tbl_2 <- if (nrow(g1) == 0 | nrow(g2) == 0 | nrow(g3) == 0 |
               nrow(g4_2) == 0 | nrow(g5_2) == 0 | nrow(g6_2) == 0) {
    NULL
  } else {
    is_sharing(g1, g2, g3, g4_1, g5_1, g6_1,
               group_keys = list(
                 g1 = c("Tissue", "HematoLineage", "TimePointType"),
                 g2 = c("Tissue", "HematoLineage", "TimePointType"),
                 g3 = c("Tissue", "AddedField2", "HematoLineage", "TimePointType"),
                 g4 = c("Tissue", "AddedField2", "Cell.type", "TimePointType"),
                 g5 = c("Tissue", "AddedField2", "Cell.type", "TimePointType"),
                 g6 = c("Tissue", "AddedField2", "Cell.type", "TimePointType")
               ))
  }
  rbindlist(list(tbl_1, tbl_2))
}

pt_2_sharing <- function(df) {
  list(
    A = pt_2A_sharing(df),
    B = pt_2B_sharing(df),
    C = pt_2C_sharing(df),
    D = pt_2D_sharing(df),
    E = pt_2E_sharing(df)
  )
}

write_sheet_2 <- function(wb, sheet, tables) {
  curr_row <- 1
  curr_row <- write_section("A", tables, wb, sheet, curr_row)
  curr_row <- write_section("B", tables, wb, sheet, curr_row)
  curr_row <- write_section("C", tables, wb, sheet, curr_row)
  curr_row <- write_section("D", tables, wb, sheet, curr_row)
  curr_row <- write_section("E", tables, wb, sheet, curr_row)
  return(NULL)
}

pt_2_sharing_results <- pt_2_sharing(data_selection_1)

pt_2_sharing_results_per_patient <- purrr::map(pt1_patients, ~ {
  pt_2_sharing(data_selection_1[SubjectID == .x])
}) %>%
  purrr::set_names(pt1_patients)

wb2 <- createWorkbook()
addWorksheet(wb2, "Overall")
write_sheet_2(wb2, "Overall", pt_2_sharing_results)

for (pt in names(pt_2_sharing_results_per_patient)) {
  addWorksheet(wb2, pt)
  write_sheet_2(wb2, pt, pt_2_sharing_results_per_patient[[pt]])
}

saveWorkbook(wb2, file = fs::path(sharing_fold, 
                                  paste(file_prefix, "sharing_analyses_pt2.xlsx", 
                                                      sep = "_")),
             overwrite = TRUE) 


## SHARING PT3 ----
### Define functions for the different tables
pt_3C_sharing <- function(df) {
  print("TABLE C")
  g1 <- df[Tissue == "PB" & HematoLineage == "CFCBULK" & TimePointType == "very late"]
  g2_1 <- df[Tissue == "BM" & AddedField2 == "DX" & HematoLineage %in% c("CFCBULK", "CD34") &
               TimePointType == "very late"]
  g2_2 <- df[Tissue == "BM" & AddedField2 == "DX" & Cell.type == "HSPC" &
               TimePointType == "very late"]
  g3_1 <- df[Tissue == "BM" & AddedField2 == "SX" & HematoLineage %in% c("CFCBULK", "CD34") &
               TimePointType == "very late"]
  g3_2 <- df[Tissue == "BM" & AddedField2 == "SX" & Cell.type == "HSPC" &
               TimePointType == "very late"]
  
  tbl_1 <- NULL
  try({
    tbl_1 <- is_sharing(g1, g2_1, g3_1,
               group_keys = list(
                 g1 = c("Tissue", "HematoLineage", "TimePointType"),
                 g2 = c("Tissue", "AddedField2", "HematoLineage", "TimePointType"),
                 g3 = c("Tissue", "AddedField2", "HematoLineage", "TimePointType")
               ), minimal = TRUE)
  })
  
  tbl_2 <- NULL
  try({
    tbl_2 <- is_sharing(g1, g2_2, g3_2,
               group_keys = list(
                 g1 = c("Tissue", "HematoLineage", "TimePointType"),
                 g2 = c("Tissue", "AddedField2", "Cell.type", "TimePointType"),
                 g3 = c("Tissue", "AddedField2", "Cell.type", "TimePointType")
               ), minimal = TRUE)
  })
  
  way3 <- rbindlist(list(tbl_1, tbl_2))
  
  combo1 <- NULL
  try({
    combo1 <- is_sharing(g1, rbindlist(list(g2_1, g3_1)),
               group_keys = list(
                 g1 = c("Tissue", "HematoLineage", "TimePointType"),
                 g2 = c("Tissue", "AddedField2", "HematoLineage", "TimePointType")
               ))
  })
  combo2 <- NULL
  try({
    combo2 <- is_sharing(g1, rbindlist(list(g2_2, g3_2)),
               group_keys = list(
                 g1 = c("Tissue", "HematoLineage", "TimePointType"),
                 g2 = c("Tissue", "AddedField2", "Cell.type", "TimePointType")
               ))
  })
  combo3 <- NULL
  try({
    combo3 <- is_sharing(rbindlist(list(g2_1, g3_1)),
                         group_key = c("Tissue", "AddedField2", "HematoLineage",
                                       "TimePointType"))
  })
  combo4 <- NULL
  try({
    combo4 <- is_sharing(
      rbindlist(list(g2_1, g3_1)),
      rbindlist(list(g2_2, g3_2)),
      group_keys = list(
        g1 = c("Tissue", "AddedField2", "HematoLineage",
               "TimePointType"),
        g2 = c("Tissue", "AddedField2", "Cell.type",
               "TimePointType")
      )
    )
  })
  rbindlist(list(way3, combo1, combo2, combo3, combo4), fill = TRUE)
}

pt_3D_sharing <- function(df) {
  print("TABLE D")
  g1 <- df[Tissue == "PB" & HematoLineage == "CFCBULK" & TimePointType == "early"]
  g2 <- df[Tissue == "PB" & HematoLineage == "CFCBULK" & TimePointType == "late"]
  g3 <- df[Tissue == "PB" & HematoLineage == "CFCBULK" & TimePointType == "very late"]
  g4_1 <- df[Tissue == "BM" & HematoLineage %in% c("CFCBULK", "CD34") & TimePointType == "early"]
  g4_2 <- df[Tissue == "BM" & Cell.type == "HSPC" & TimePointType == "early"]
  g5_1 <- df[Tissue == "BM" & HematoLineage %in% c("CFCBULK", "CD34") & TimePointType == "late"]
  g5_2 <- df[Tissue == "BM" & Cell.type == "HSPC" & TimePointType == "late"]
  g6_1 <- df[Tissue == "BM" & HematoLineage %in% c("CFCBULK", "CD34") & TimePointType == "very late"]
  g6_2 <- df[Tissue == "BM" & Cell.type == "HSPC" & TimePointType == "very late"]
  
  tbl_1 <- NULL
  try({
    tbl_1 <- is_sharing(g1, g2, g3, g4_1, g5_1, g6_1,
               group_key = c("Tissue", "HematoLineage", "TimePointType"),
               minimal = TRUE)
  })
  
  tbl_2 <- NULL
  try({
    tbl_2 <- is_sharing(g1, g2, g3, g4_2, g5_2, g6_2,
               group_keys = list(
                 g1 = c("Tissue", "HematoLineage", "TimePointType"),
                 g2 = c("Tissue", "HematoLineage", "TimePointType"),
                 g3 = c("Tissue", "HematoLineage", "TimePointType"),
                 g4 = c("Tissue", "Cell.type", "TimePointType"),
                 g5 = c("Tissue", "Cell.type", "TimePointType"),
                 g6 = c("Tissue", "Cell.type", "TimePointType")
                ),
               minimal = TRUE)
  })
  
  g1[, group := "CFCBULK"]
  g2[, group := "CFCBULK"]
  g3[, group := "CFCBULK"]
  g4_1[HematoLineage == "CFCBULK", group := "CFCBULK"]
  g4_1[HematoLineage == "CD34", group := "CD34"]
  g5_1[HematoLineage == "CFCBULK", group := "CFCBULK"]
  g5_1[HematoLineage == "CD34", group := "CD34"]
  g6_1[HematoLineage == "CFCBULK", group := "CFCBULK"]
  g6_1[HematoLineage == "CD34", group := "CD34"]
  g4_2[, group := "HSPC"]
  g5_2[, group := "HSPC"]
  g6_2[, group := "HSPC"]
  
  print("TABLE D - combos 2")
  combo2 <- NULL
  try({
    combo2 <- is_sharing(
      rbindlist(list(g1, g2, g3, g4_1, g4_2, g5_1, g5_2, g6_1, g6_2)),
      n_comp = 2,
      group_key = c("Tissue", "group", "TimePointType")
    )
  })
  
  print("TABLE D - combos 3")
  combo3 <- NULL
  try({
    combo3 <- is_sharing(
      rbindlist(list(g1, g2, g3, g4_1, g4_2, g5_1, g5_2, g6_1, g6_2)),
      n_comp = 3,
      group_key = c("Tissue", "group", "TimePointType")
    )
  })
  
  print("TABLE D - combos 4")
  combo4 <- NULL
  try({
    combo4 <- is_sharing(
      rbindlist(list(g1, g2, g3, g4_1, g4_2, g5_1, g5_2, g6_1, g6_2)),
      n_comp = 4,
      group_key = c("Tissue", "group", "TimePointType")
    )
  })
  
  print("TABLE D - combos 5")
  combo5 <- NULL
  try({
    combo5 <- is_sharing(
      rbindlist(list(g1, g2, g3, g4_1, g4_2, g5_1, g5_2, g6_1, g6_2)),
      n_comp = 5,
      group_key = c("Tissue", "group", "TimePointType")
    )
  })
  
  rbindlist(list(tbl_1, tbl_2, combo2, combo3, combo4, combo5), fill = TRUE)
}

pt_3E_sharing <- function(df) {
  print("TABLE E")
  g1 <- df[Tissue == "PB" & HematoLineage == "CFCBULK" & TimePointType == "early"]
  g2 <- df[Tissue == "PB" & HematoLineage == "CFCBULK" & TimePointType == "late"]
  g3 <- df[Tissue == "PB" & HematoLineage == "CFCBULK" & TimePointType == "very late"]
  g4_1 <- df[Tissue == "BM" & AddedField2 == "DX" &
               HematoLineage %in% c("CFCBULK", "CD34") & TimePointType == "early"]
  g4_2 <- df[Tissue == "BM" & AddedField2 == "DX" &
               Cell.type == "HSPC" & TimePointType == "early"]
  g5_1 <- df[Tissue == "BM" & AddedField2 == "DX" &
               HematoLineage %in% c("CFCBULK", "CD34") & TimePointType == "late"]
  g5_2 <- df[Tissue == "BM" & AddedField2 == "DX" &
               Cell.type == "HSPC" & TimePointType == "late"]
  g6_1 <- df[Tissue == "BM" & AddedField2 == "DX" &
               HematoLineage %in% c("CFCBULK", "CD34") & TimePointType == "very late"]
  g6_2 <- df[Tissue == "BM" & AddedField2 == "DX" &
               Cell.type == "HSPC" & TimePointType == "very late"]
  g7_1 <- df[Tissue == "BM" & AddedField2 == "SX" &
               HematoLineage %in% c("CFCBULK", "CD34") & TimePointType == "early"]
  g7_2 <- df[Tissue == "BM" & AddedField2 == "SX" &
               Cell.type == "HSPC" & TimePointType == "early"]
  g8_1 <- df[Tissue == "BM" & AddedField2 == "SX" &
               HematoLineage %in% c("CFCBULK", "CD34") & TimePointType == "late"]
  g8_2 <- df[Tissue == "BM" & AddedField2 == "SX" &
               Cell.type == "HSPC" & TimePointType == "late"]
  g9_1 <- df[Tissue == "BM" & AddedField2 == "SX" &
               HematoLineage %in% c("CFCBULK", "CD34") & TimePointType == "very late"]
  g9_2 <- df[Tissue == "BM" & AddedField2 == "SX" &
               Cell.type == "HSPC" & TimePointType == "very late"]
  tbl_1 <- NULL
  try({
    tbl_1 <- is_sharing(g1, g2, g3, g4_1, g5_1, g6_1, g7_1, g8_1, g9_1,
                        group_keys = list(
                          g1 = c("Tissue", "HematoLineage", "TimePointType"),
                          g2 = c("Tissue", "HematoLineage", "TimePointType"),
                          g3 = c("Tissue", "HematoLineage", "TimePointType"),
                          g4 = c("Tissue", "AddedField2", "HematoLineage", "TimePointType"),
                          g5 = c("Tissue", "AddedField2", "HematoLineage", "TimePointType"),
                          g6 = c("Tissue", "AddedField2", "HematoLineage", "TimePointType"),
                          g7 = c("Tissue", "AddedField2", "HematoLineage", "TimePointType"),
                          g8 = c("Tissue", "AddedField2", "HematoLineage", "TimePointType"),
                          g9 = c("Tissue", "AddedField2", "HematoLineage", "TimePointType")
                        ),
                        minimal = TRUE)
  })
  
  tbl_2 <- NULL
  try({
    tbl_2 <- is_sharing(g1, g2, g3, g4_2, g5_2, g6_2, g7_2, g8_2, g9_2,
               group_keys = list(
                 g1 = c("Tissue", "HematoLineage", "TimePointType"),
                 g2 = c("Tissue", "HematoLineage", "TimePointType"),
                 g3 = c("Tissue", "HematoLineage", "TimePointType"),
                 g4 = c("Tissue", "AddedField2", "Cell.type", "TimePointType"),
                 g5 = c("Tissue", "AddedField2", "Cell.type", "TimePointType"),
                 g6 = c("Tissue", "AddedField2", "Cell.type", "TimePointType"),
                 g7 = c("Tissue", "AddedField2", "Cell.type", "TimePointType"),
                 g8 = c("Tissue", "AddedField2", "Cell.type", "TimePointType"),
                 g9 = c("Tissue", "AddedField2", "Cell.type", "TimePointType")
               ),
               minimal = TRUE)
  })
  
  g1[, group := "CFCBULK"]
  g2[, group := "CFCBULK"]
  g3[, group := "CFCBULK"]
  g1[, group2 := NA_character_]
  g2[, group2 := NA_character_]
  g3[, group2 := NA_character_]
  g4_1[HematoLineage == "CFCBULK", group := "CFCBULK"]
  g4_1[HematoLineage == "CD34", group := "CD34"]
  g5_1[HematoLineage == "CFCBULK", group := "CFCBULK"]
  g5_1[HematoLineage == "CD34", group := "CD34"]
  g6_1[HematoLineage == "CFCBULK", group := "CFCBULK"]
  g6_1[HematoLineage == "CD34", group := "CD34"]
  g7_1[HematoLineage == "CFCBULK", group := "CFCBULK"]
  g7_1[HematoLineage == "CD34", group := "CD34"]
  g8_1[HematoLineage == "CFCBULK", group := "CFCBULK"]
  g8_1[HematoLineage == "CD34", group := "CD34"]
  g9_1[HematoLineage == "CFCBULK", group := "CFCBULK"]
  g9_1[HematoLineage == "CD34", group := "CD34"]
  g4_2[, group := "HSPC"]
  g5_2[, group := "HSPC"]
  g6_2[, group := "HSPC"]
  g7_2[, group := "HSPC"]
  g8_2[, group := "HSPC"]
  g9_2[, group := "HSPC"]
  g4_1[, group2 := AddedField2]
  g4_2[, group2 := AddedField2]
  g5_1[, group2 := AddedField2]
  g5_2[, group2 := AddedField2]
  g6_1[, group2 := AddedField2]
  g6_2[, group2 := AddedField2]
  g7_1[, group2 := AddedField2]
  g7_2[, group2 := AddedField2]
  g8_1[, group2 := AddedField2]
  g8_2[, group2 := AddedField2]
  g9_1[, group2 := AddedField2]
  g9_2[, group2 := AddedField2]
  
  print("TABLE E - combo2")
  combo2 <- NULL
  try({
    combo2 <- is_sharing(
      rbindlist(list(g1, g2, g3, g4_1, g4_2, g5_1, g5_2, g6_1, g6_2,
                     g7_1, g7_2, g8_1, g8_2, g9_1, g9_2)),
      n_comp = 2,
      group_key = c("Tissue", "group2", "group", "TimePointType")
    )
  })
  
  print("TABLE E - combo3")
  combo3 <- NULL
  try({
    combo3 <- is_sharing(
      rbindlist(list(g1, g2, g3, g4_1, g4_2, g5_1, g5_2, g6_1, g6_2,
                     g7_1, g7_2, g8_1, g8_2, g9_1, g9_2)),
      n_comp = 3,
      group_key = c("Tissue", "group2", "group", "TimePointType")
    )
  })
  
  print("TABLE E - combo4")
  combo4 <- NULL
  try({
    combo4 <- is_sharing(
      rbindlist(list(g1, g2, g3, g4_1, g4_2, g5_1, g5_2, g6_1, g6_2,
                     g7_1, g7_2, g8_1, g8_2, g9_1, g9_2)),
      n_comp = 4,
      group_key = c("Tissue", "group2", "group", "TimePointType")
    )
  })
  
  print("TABLE E - combo5")
  combo5 <- NULL
  try({
    combo5 <- is_sharing(
      rbindlist(list(g1, g2, g3, g4_1, g4_2, g5_1, g5_2, g6_1, g6_2,
                     g7_1, g7_2, g8_1, g8_2, g9_1, g9_2)),
      n_comp = 5,
      group_key = c("Tissue", "group2", "group", "TimePointType")
    )
  })
  
  print("TABLE E - combo6")
  combo6 <- NULL
  try({
    combo6 <- is_sharing(
      rbindlist(list(g1, g2, g3, g4_1, g4_2, g5_1, g5_2, g6_1, g6_2,
                     g7_1, g7_2, g8_1, g8_2, g9_1, g9_2)),
      n_comp = 6,
      group_key = c("Tissue", "group2", "group", "TimePointType")
    )
  })
  
  print("TABLE E - combo7")
  combo7 <- NULL
  try({
    combo7 <- is_sharing(
      rbindlist(list(g1, g2, g3, g4_1, g4_2, g5_1, g5_2, g6_1, g6_2,
                     g7_1, g7_2, g8_1, g8_2, g9_1, g9_2)),
      n_comp = 7,
      group_key = c("Tissue", "group2", "group", "TimePointType")
    )
  })
  
  print("TABLE E - combo8")
  print("Calculating 8-way combos")
  combo8 <- NULL
  try({
    combo8 <- is_sharing(
      rbindlist(list(g1, g2, g3, g4_1, g4_2, g5_1, g5_2, g6_1, g6_2,
                     g7_1, g7_2, g8_1, g8_2, g9_1, g9_2)),
      n_comp = 8,
      group_key = c("Tissue", "group2", "group", "TimePointType")
    )
  })
  
  rbindlist(list(tbl_1, tbl_2, combo2, combo3, combo4, combo5, combo6, combo7, combo8),
            fill = TRUE)
}

pt_3_sharing <- function(df) {
  list(
    C = pt_3C_sharing(df),
    D = pt_3D_sharing(df),
    E = pt_3E_sharing(df)
  )
}

write_sheet_3 <- function(wb, sheet, tables) {
  curr_row <- 1
  curr_row <- write_section("C", tables, wb, sheet, curr_row)
  curr_row <- write_section("D", tables, wb, sheet, curr_row)
  curr_row <- write_section("E", tables, wb, sheet, curr_row)
  return(NULL)
}

pt_3_sharing_results <- pt_3_sharing(data_selection_1)

pt_3_sharing_results_per_patient <- purrr::map(pt1_patients, ~ {
  pt_3_sharing(data_selection_1[SubjectID == .x])
}) %>%
  purrr::set_names(pt1_patients)

#save.image()

wb3 <- createWorkbook()
addWorksheet(wb3, "Overall")
write_sheet_3(wb3, "Overall", pt_3_sharing_results)
for (pt in names(pt_3_sharing_results_per_patient)) {
  addWorksheet(wb3, pt)
  write_sheet_3(wb3, pt, pt_3_sharing_results_per_patient[[pt]])
}
saveWorkbook(wb3, file = fs::path(sharing_fold, paste(file_prefix, 
                                                      "sharing_analyses_pt3.xlsx", 
                                                      sep = "_")),
             overwrite = TRUE)

## SHARING PT4 ----
pt4_patients <- c("WASHE02", "WASCUP04", "WAS1001", "WAS1003", "WAS1004",
                  "WAS1007", "WAS1008", "WAS1009")
data_selection_2 <- final_matrix %>%
  left_join(final_af, by = "CompleteAmplificationID") %>%
  filter(SubjectID %in% pt4_patients, TimePointType == "very late")
setDT(data_selection_2)

### Define functions for the different tables
pt_4A_sharing <- function(df) {
  ## No need for selection on tp since data_selection_2 contains only very
  ## late tps
  g1_1 <- df[Tissue == "PB" & HematoLineage == "CFCBULK"]
  g2_1 <- df[Tissue == "BM" & HematoLineage %in% c("CFCBULK", "CD34", 
                                                   "PrimHSPC", "MyeloHSPC")]
  g2_2 <- df[Tissue == "BM" & Cell.type == "HSPC"]
  g2_3 <- df[Tissue == "BM" & HematoLineage == "LymphoHSPC" &
               Cell.type == "HSPC"]
  
  g2_4 <- df[Tissue == "BM" & AddedField2 %in% c("DX", "SX") & 
               HematoLineage %in% c("CFCBULK", "CD34", 
                                    "PrimHSPC", "MyeloHSPC")]
  g2_5 <- df[Tissue == "BM" & AddedField2 %in% c("DX", "SX") &
               Cell.type == "HSPC"]
  g2_6 <- df[Tissue == "BM" & AddedField2 %in% c("DX", "SX") & 
               HematoLineage == "LymphoHSPC" &
               Cell.type == "HSPC"]
  
  tbl_1 <- NULL
  try({
    tbl_1 <- is_sharing(
      g1_1, g2_1,
      group_key = c("Tissue", "HematoLineage", "TimePointType")
    )
  })
  tbl_2 <- NULL
  try({
    tbl_2 <- is_sharing(
      g1_1, g2_2,
      group_keys = list(
        g1 = c("Tissue", "HematoLineage", "TimePointType"),
        g2 = c("Tissue", "Cell.type", "TimePointType")
        )
    )
  })
  tbl_3 <- NULL
  try({
    tbl_3 <- is_sharing(
      g1_1, g2_3,
      group_keys = list(
        g1 = c("Tissue", "HematoLineage", "TimePointType"),
        g2 = c("Tissue", "HematoLineage", "Cell.type", "TimePointType")
      )
    )
  })
  tbl_4 <- NULL
  try({
    tbl_4 <- is_sharing(
      g1_1, g2_4,
      group_keys = list(
        g1 = c("Tissue", "HematoLineage", "TimePointType"),
        g2 = c("Tissue", "AddedField2", "HematoLineage", "TimePointType")
      ),
      minimal = TRUE
    )
  })
  tbl_5 <- NULL
  try({
    tbl_5 <- is_sharing(
      g1_1, g2_5,
      group_keys = list(
        g1 = c("Tissue", "HematoLineage", "TimePointType"),
        g2 = c("Tissue", "AddedField2", "Cell.type", "TimePointType")
      )
    )
  })
  tbl_6 <- NULL
  try({
    tbl_6 <- is_sharing(
      g1_1, g2_6,
      group_keys = list(
        g1 = c("Tissue", "HematoLineage", "TimePointType"),
        g2 = c("Tissue", "AddedField2", "HematoLineage", "Cell.type", "TimePointType")
      )
    )
  })
  tbl_7 <- NULL
  try({
    tbl_7 <- is_sharing(
      g2_1[HematoLineage == "CFCBULK"],
      g2_1[HematoLineage == "CD34"],
      group_key = c("Tissue", "HematoLineage", "TimePointType")
    )
  })
  tbl_8 <- NULL
  try({
    tbl_8 <- is_sharing(
      g2_4[HematoLineage == "CFCBULK"],
      g2_4[HematoLineage == "CD34"],
      group_key = c("Tissue", "AddedField2", "HematoLineage", "TimePointType")
    )
  })
  tbl_9 <- NULL
  try({
    tbl_9 <- is_sharing(
      g2_1[HematoLineage == "CFCBULK"],
      g2_2,
      group_keys = list(
        g1 = c("Tissue", "HematoLineage", "TimePointType"),
        g2 = c("Tissue", "Cell.type", "TimePointType")
      )
    )
  })
  tbl_10 <- NULL
  try({
    tbl_10 <- is_sharing(
      g2_4[HematoLineage == "CFCBULK"],
      g2_5,
      group_keys = list(
        g1 = c("Tissue", "HematoLineage", "TimePointType"),
        g2 = c("Tissue", "AddedField2", "Cell.type", "TimePointType")
      )
    )
  })
  tbl_11 <- NULL
  try({
    tbl_11 <- is_sharing(
      g2_2,
      g2_1[HematoLineage == "CD34"],
      group_keys = list(
        g1 = c("Tissue", "Cell.type", "TimePointType"),
        g2 = c("Tissue", "HematoLineage", "TimePointType")
      )
    )
  })
  tbl_12 <- NULL
  try({
    tbl_12 <- is_sharing(
      g2_5,
      g2_4[HematoLineage == "CD34"],
      group_keys = list(
        g1 = c("Tissue", "AddedField2", "Cell.type", "TimePointType"),
        g2 = c("Tissue", "AddedField2", "HematoLineage", "TimePointType")
      )
    )
  })
  
  rbindlist(list(tbl_1, tbl_2, tbl_3, tbl_4, tbl_5, tbl_6, tbl_7, tbl_8,
                 tbl_9, tbl_10, tbl_11, tbl_12))
}

pt_4B_sharing <- function(df) {
  ## No need for selection on tp since data_selection_2 contains only very
  ## late tps
  g1_1 <- df[Tissue %in% c("PB", "BM") & HematoLineage == "CFCBULK"]
  g1_2 <- df[Tissue == "BM" & AddedField2 %in% c("DX", "SX") & HematoLineage == "CFCBULK"]
  g1_3 <- df[Tissue == "BM" & Cell.type %in% c("CD34", "HSPC")]
  g2_1 <- df[Tissue == "PB" & HematoLineage %in% c("Myeloid", "Lymphoid")]
  
  tbl_1 <- NULL
  try({
    tbl_1 <- is_sharing(
      g1_1, g2_1,
      group_key = c("Tissue", "HematoLineage", "TimePointType")
    )
  })
  tbl_2 <- NULL
  try({
    tbl_2 <- is_sharing(
      g1_2, g2_1,
      group_keys = list(
        g1 = c("Tissue", "AddedField2", "HematoLineage", "TimePointType"),
        g2 = c("Tissue", "HematoLineage", "TimePointType")
      )
    )
  })
  tbl_3 <- NULL
  try({
    tbl_3 <- is_sharing(
      g1_1, g2_1,
      group_keys = list(
        g1 = c("Tissue", "HematoLineage", "TimePointType"),
        g2 = c("Tissue", "HematoLineage", "Cell.type", "TimePointType")
      )
    )
  })
  tbl_4 <- NULL
  try({
    tbl_4 <- is_sharing(
      g1_2, g2_1,
      group_keys = list(
        g1 = c("Tissue", "AddedField2", "HematoLineage", "TimePointType"),
        g2 = c("Tissue", "HematoLineage", "Cell.type", "TimePointType")
      )
    )
  })
  tbl_5 <- NULL
  try({
    tbl_5 <- is_sharing(
      g1_3, g2_1,
      group_keys = list(
        g1 = c("Tissue", "Cell.type", "TimePointType"),
        g2 = c("Tissue", "HematoLineage", "TimePointType")
      )
    )
  })
  tbl_6 <- NULL
  try({
    tbl_6 <- is_sharing(
      g1_3, g2_1,
      group_keys = list(
        g1 = c("Tissue", "Cell.type", "TimePointType"),
        g2 = c("Tissue", "HematoLineage", "Cell.type", "TimePointType")
      )
    )
  })
  
  rbindlist(list(tbl_1, tbl_2, tbl_3, tbl_4, tbl_5, tbl_6))
}

pt_4C_sharing <- function(df) {
  ## No need for selection on tp since data_selection_2 contains only very
  ## late tps
  g1 <- df[Tissue == "PB" & HematoLineage == "CFCBULK"]
  g2_1 <- df[Tissue == "BM" & AddedField2 == "DX" & 
               HematoLineage %in% c("CFCBULK", "CD34")] 
  g2_2 <- df[Tissue == "BM" & AddedField2 == "DX" & 
               Cell.type == "HSPC"] 
  g3_1 <- df[Tissue == "BM" & AddedField2 == "SX" & 
               HematoLineage %in% c("CFCBULK", "CD34")] 
  g3_2 <- df[Tissue == "BM" & AddedField2 == "SX" & 
               Cell.type == "HSPC"] 
  
  tbl_1 <- NULL
  try({
    tbl_1 <- is_sharing(
      g1, g2_1, g3_1,
      group_keys = list(
        g1 = c("Tissue", "HematoLineage", "TimePointType"),
        g2 = c("Tissue", "AddedField2", "HematoLineage", "TimePointType"),
        g3 = c("Tissue", "AddedField2", "HematoLineage", "TimePointType")
      )
    )
  })
  tbl_2 <- NULL
  try({
    tbl_2 <- is_sharing(
      g1, g2_2, g3_2,
      group_keys = list(
        g1 = c("Tissue", "HematoLineage", "TimePointType"),
        g2 = c("Tissue", "AddedField2", "Cell.type", "TimePointType"),
        g3 = c("Tissue", "AddedField2", "Cell.type", "TimePointType")
      )
    )
  })
  
  g1[, group := "CFCBULK"]
  g1[, group2 := NA_character_]
  g2_1[, group := HematoLineage]
  g2_1[, group2 := AddedField2]
  g2_2[, group := HematoLineage]
  g2_2[, group2 := AddedField2]
  g3_1[, group := "HSPC"]
  g3_1[, group2 := AddedField2]
  g3_2[, group := "HSPC"]
  g3_2[, group2 := AddedField2]
  
  combos <- NULL
  try({
    combos <- is_sharing(
      rbindlist(list(g1, g2_1, g2_2, g3_1, g3_2)),
      n_comp = 2,
      group_key = c("Tissue", "group2", "group", "TimePointType")
    )
  })
  rbindlist(list(tbl_1, tbl_2, combos), fill = TRUE)
}

pt_4_sharing <- function(df) {
  list(
    A = pt_4A_sharing(df),
    B = pt_4B_sharing(df),
    C = pt_4C_sharing(df)
  )
}

write_sheet_4 <- function(wb, sheet, tables) {
  curr_row <- 1
  curr_row <- write_section("A", tables, wb, sheet, curr_row)
  curr_row <- write_section("B", tables, wb, sheet, curr_row)
  curr_row <- write_section("C", tables, wb, sheet, curr_row)
  return(NULL)
}

pt_4_sharing_results <- pt_4_sharing(data_selection_2)

pt_4_sharing_results_per_patient <- purrr::map(pt4_patients, ~ {
  pt_4_sharing(data_selection_2[SubjectID == .x])
}) %>%
  purrr::set_names(pt4_patients)

wb4 <- createWorkbook()
addWorksheet(wb4, "Overall")
write_sheet_4(wb4, "Overall", pt_4_sharing_results)
for (pt in names(pt_4_sharing_results_per_patient)) {
  addWorksheet(wb4, pt)
  write_sheet_4(wb4, pt, pt_4_sharing_results_per_patient[[pt]])
}
saveWorkbook(wb4, file = fs::path(sharing_fold, paste(file_prefix, 
                                                      "sharing_analyses_pt4.xlsx", 
                                                      sep = "_")),
             overwrite = TRUE)
