## Min size plots


# ============================================================
# 1. Load packages and locate files
# ============================================================

library(tidyverse)

base_dir <- paste0(
  "/home/bio.aau.dk/hd95lp/Projects/RETHINK/",
  "p_nanopore_amplicon/6_snakemake_usearch_unoise_minsize_test/",
  "global_ASVs_min_size_complex/"
)

# Find map_raw.tsv for all primer regions
map_files <- list.files(
  base_dir,
  pattern = "^map_raw\\.tsv$",
  recursive = TRUE,
  full.names = TRUE
)


# ============================================================
# 1. Filter lenght of ASVs 
# ============================================================


f_OPR <- paste0(base_dir,"global_ASVs_OPR/global_ASVs.fa")
#f_V4 <-  paste0(base_dir, "global_ASVs_V4/global_ASVs.fa")
#f_V1V3 <- paste0(base_dir,"V1V3/global_ASVs.fa")
f_V1V8 <- paste0(base_dir, "global_ASVs_V1V8/global_ASVs.fa")


# Put in a named vector so we keep primer names
files <- c(OPR = f_OPR,
  #V4 = f_V4,
  #V1V3 = f_V1V3,
  V1V8 = f_V1V8
)


df_all <- map2_dfr(
  .x = files,
  .y = names(files),
  .f = ~ {
    seqs <- Biostrings::readDNAStringSet(.x)
    tibble(
      primerset = .y,
      seq_name  = names(seqs),
      length = Biostrings::width(seqs)
    )
  }
) 

# lenght_df_all %>% 
#   distinct(primerset, limit_max, limit_min)

length_limits <- tibble::tribble(
  ~primerset, ~limit_max, ~limit_min,
  "OPR",          5463,       3281,
  "V4",            316,        190,    # 320
  "V1V3",          624,        376,
  "V1V8",         1716,       1031
)

lenght_df <- df_all %>% 
  left_join(length_limits) %>% 
  filter((length > limit_min & length < limit_max)) %>% 
  mutate(global_asv = str_remove(
    seq_name,
    ";size=\\d+;$"
  ),)

df_all %>% distinct(seq_name)
lenght_df %>% distinct(global_asv)

# ============================================================
# 2. Read depths from earlier analysis
# ============================================================

# Exclude the hyper-accurate basecalling datasets
reads_per_barcode_filt <- reads_per_barcode %>%
  filter(basecalling != "hyper")


# ------------------------------------------------------------
# Calculate the actual depth of each all_reads dataset
# ------------------------------------------------------------

all_read_depths <- reads_per_barcode_filt %>%
  
  # Only retain rows representing all available filtered reads
  filter(
    subsampling_group == "sample_size_all_reads"
  ) %>%
  
  mutate(
    # Use the same platform names as in map_data
    platform = case_when(
      str_detect(
        basecalling,
        regex("^pacbio$", ignore_case = TRUE)
      ) ~ "PacBio",
      
      str_detect(
        basecalling,
        regex("^(super|sup|ont|nanopore)$", ignore_case = TRUE)
      ) ~ "ONT",
      
      TRUE ~ NA_character_
    ),
    
    # Rename metadata columns to match map_data
    primer = primerset,
    sample_type = sample_content
  ) %>%
  
  # Remove basecalling types that were not classified
  filter(
    !is.na(platform)
  ) %>%
  
  group_by(
    primer,
    sample_type,
    platform
  ) %>%
  
  summarise(
    # Actual total number of filtered reads across all barcodes
    total_reads = sum(
      reads,
      na.rm = TRUE
    ),
    
    # Retained for quality control
    n_barcodes = n_distinct(barcode),
    
    .groups = "drop"
  )


# Inspect the calculated all-read depths
# all_read_depths %>%
#   arrange(
#     primer,
#     sample_type,
#     platform
#   ) %>%
#   print(
#     n = Inf,
#     width = Inf
#   )


# ============================================================
# 3. Read and parse all map_raw.tsv files
# ============================================================


# Extract primer name from the immediate parent directory:
# global_ASVs_V1V3 -> V1V3
get_primer <- function(file) {
  
  basename(dirname(file)) %>%
    str_remove("^global_ASVs_")
}


# Read one map_raw.tsv file
read_map_file <- function(file) {
  
  read_tsv(
    file,
    col_names = c(
      "local_asv",
      "global_asv"
    ),
    col_types = "cc",
    show_col_types = FALSE
  ) %>%
    
    transmute(
      primer = get_primer(file),
      
      # Original dataset-specific ASV name
      local_asv,
      
      # Remove Zotu1_, Zotu2_, etc. to recover dataset name
      dataset = str_remove(
        local_asv,
        "^Zotu\\d+_"
      ),
      
      # Remove abundance annotation from the global ASV
      global_asv = str_remove(
        global_asv,
        ";size=\\d+;$"
      ),
      
      source_file = file
    )
}


# Read map_raw.tsv for all primer regions
map_data <- map_dfr(
  map_files,
  read_map_file
) %>% 
  filter(global_asv %in% unlist(lenght_df$global_asv))


# ============================================================
# 4. Parse dataset metadata and add actual read depth
# ============================================================

map_data <- map_data %>%
  
  # Remove duplicate all-read datasets
  #
  # Remove:
  #   ..._minsize_4_all_reads
  #
  # Retain:
  #   ..._minsize_4_sample_size_all_reads
  filter(
    !str_detect(
      dataset,
      "_minsize_[0-9]+_all_reads$"
    )
  ) %>%
  
  mutate(
    # --------------------------------------------------------
    # Sequencing platform
    # --------------------------------------------------------
    
    platform = case_when(
      str_detect(
        dataset,
        regex("^pacbio_", ignore_case = TRUE)
      ) ~ "PacBio",
      
      # The super datasets represent ONT
      str_detect(
        dataset,
        regex("^(super|ont|nanopore)_", ignore_case = TRUE)
      ) ~ "ONT",
      
      TRUE ~ NA_character_
    ),
    
    
    # --------------------------------------------------------
    # Sample type
    # --------------------------------------------------------
    #
    # Example:
    # super_all_OPR_Zmock_unoise_minsize_4_sample_size_100000
    #               ^^^^^
    
    sample_type = str_match(
      dataset,
      "^[^_]+_all_[^_]+_(.+?)_unoise_minsize_"
    )[, 2],
    
    
    # --------------------------------------------------------
    # UNOISE minsize
    # --------------------------------------------------------
    
    minsize = str_match(
      dataset,
      "_minsize_([0-9]+)_"
    )[, 2] %>%
      as.integer(),
    
    
    # --------------------------------------------------------
    # Read-depth label
    # --------------------------------------------------------
    #
    # Extract:
    #   sample_size_10000       -> 10000
    #   sample_size_100000      -> 100000
    #   sample_size_all_reads   -> all_reads
    
    read_token = str_match(
      dataset,
      "_sample_size_([0-9]+|all_reads)$"
    )[, 2],
    
    
    # Standardise all_reads to "all"
    read_label = case_when(
      read_token == "all_reads" ~ "all",
      TRUE ~ read_token
    ),
    
    
    # Numerical subset depth
    #
    # This is NA for sample_size_all_reads
    subset_depth = suppressWarnings(
      as.numeric(read_label)
    ),
    
    
    # Identify all_reads datasets
    is_all = read_label == "all"
  ) %>%
  
  select(
    -read_token
  ) %>%
  
  # ----------------------------------------------------------
# Join actual all_reads depths
# ----------------------------------------------------------

left_join(
  all_read_depths,
  by = c(
    "primer",
    "sample_type",
    "platform"
  )
) %>%
  
  mutate(
    # For numerical subsets, use the depth encoded in the name.
    #
    # For sample_size_all_reads, use the actual summed number
    # of filtered reads across barcodes.
    actual_depth = if_else(
      is_all,
      total_reads,
      subset_depth
    ),
    
    # Optional diagnostic label
    depth_type = if_else(
      is_all,
      "All available reads",
      "Numerical subset"
    )
  )


# ============================================================
# 5. Check metadata parsing
# ============================================================

parsing_problems <- map_data %>%
  filter(
    is.na(platform) |
      is.na(sample_type) |
      is.na(minsize) |
      is.na(read_label)
  ) %>%
  distinct(
    primer,
    dataset
  )

parsing_problems



# ============================================================
# 6. Check all_reads depth assignment
# ============================================================

missing_all_read_depths <- map_data %>%
  filter(
    is_all,
    is.na(actual_depth)
  ) %>%
  distinct(
    primer,
    sample_type,
    platform,
    dataset
  )

missing_all_read_depths


# map_data %>%
#   filter(is_all) %>%
#   distinct(
#     primer,
#     sample_type,
#     platform,
#     minsize,
#     dataset,
#     n_barcodes,
#     total_reads,
#     actual_depth
#   ) %>%
#   arrange(
#     primer,
#     sample_type,
#     platform,
#     minsize
#   ) %>%
#   print(
#     n = Inf,
#     width = Inf
#   )

# ============================================================
# Create one global-ASV set per dataset
# ============================================================

asv_sets <- map_data %>%
  filter(
    !is.na(global_asv),
    global_asv != ""
  ) %>%
  
  group_by(
    primer,
    sample_type,
    platform,
    minsize,
    read_label,
    subset_depth,
    actual_depth,
    is_all,
    depth_type
  ) %>%
  
  summarise(
    # Unique global ASVs recovered in this dataset
    asvs = list(
      sort(unique(global_asv))
    ),
    
    # Total number of unique global ASVs
    n_asvs = n_distinct(global_asv),
    
    # Retain the original dataset name for checking
    dataset = paste(
      sort(unique(dataset)),
      collapse = " | "
    ),
    
    .groups = "drop"
  )


duplicate_asv_sets <- asv_sets %>%
  count(
    primer,
    sample_type,
    platform,
    minsize,
    actual_depth
  ) %>%
  filter(n > 1)

duplicate_asv_sets


pacbio_reference_sets <- asv_sets %>%
  filter(
    platform == "PacBio",
    minsize == 8
  ) %>%
  arrange(
    primer,
    sample_type,
    actual_depth
  )

pacbio_reference_sets %>%
  select(
    primer,
    sample_type,
    minsize,
    read_label,
    actual_depth,
    depth_type,
    n_asvs
  ) %>%
  print(
    n = Inf,
    width = Inf
  )

# ============================================================
# Select an ASV set at a requested comparison depth
# ============================================================

pick_depth_set <- function(platform_data, requested_depth) {
  
  # ----------------------------------------------------------
  # First choice: dataset at the exact requested depth
  # ----------------------------------------------------------
  
  exact_set <- platform_data %>%
    filter(
      actual_depth == requested_depth
    ) %>%
    
    # In the unlikely event that both a numerical subset and
    # all_reads have the same depth, prefer the numerical subset
    arrange(is_all) %>%
    
    slice(1)
  
  if (nrow(exact_set) == 1) {
    
    return(
      list(
        asvs = exact_set$asvs[[1]],
        n_asvs = exact_set$n_asvs[[1]],
        
        actual_depth = exact_set$actual_depth[[1]],
        is_all = exact_set$is_all[[1]],
        
        capped = FALSE,
        
        dataset = exact_set$dataset[[1]]
      )
    )
  }
  
  
  # ----------------------------------------------------------
  # Second choice: all_reads when the requested depth is
  # greater than the total available depth
  # ----------------------------------------------------------
  
  all_set <- platform_data %>%
    filter(is_all) %>%
    slice(1)
  
  if (
    nrow(all_set) == 1 &&
    all_set$actual_depth[[1]] < requested_depth
  ) {
    
    return(
      list(
        asvs = all_set$asvs[[1]],
        n_asvs = all_set$n_asvs[[1]],
        
        actual_depth = all_set$actual_depth[[1]],
        is_all = TRUE,
        
        capped = TRUE,
        
        dataset = all_set$dataset[[1]]
      )
    )
  }
  
  
  # ----------------------------------------------------------
  # No valid depth-matched or capped dataset exists
  # ----------------------------------------------------------
  
  NULL
}



# ============================================================
# Compare one ONT minsize with PacBio minsize 8
# ============================================================




  ont_minsize_value = c(2,4,6,8)
  pacbio_minsize_value = 8
  
    # ----------------------------------------------------------
  # Find primer/sample combinations available for ONT
  # ----------------------------------------------------------
  
  ont_keys <- asv_sets %>%
    filter(
      platform == "ONT",
      minsize == ont_minsize_value
    ) %>%
    distinct(
      primer,
      sample_type
    )
  
  
  # ----------------------------------------------------------
  # Find primer/sample combinations available for the
  # PacBio minsize-8 reference
  # ----------------------------------------------------------
  
  pacbio_keys <- asv_sets %>%
    filter(
      platform == "PacBio",
      minsize == pacbio_minsize_value
    ) %>%
    distinct(
      primer,
      sample_type
    )
  
  
  # Only retain combinations present for both platforms
  comparison_keys <- inner_join(
    ont_keys,
    pacbio_keys,
    by = c(
      "primer",
      "sample_type"
    )
  )
  
  
  included_groups_vec = c(
    #"sample_size_100", "sample_size_1000", 
    #"sample_size_10000" , 
    #"sample_size_100000" , 
    "sample_size_1000000",     "sample_size_2000000",     "sample_size_5000000", 
    #"sample_size_500", "sample_size_5000", 
    "sample_size_50000" , "sample_size_500000" ,
    #"sample_size_200", "sample_size_2000", 
    "sample_size_20000" , "sample_size_200000" ,"sample_size_all_reads" 
    
  )
   
complete_asv_sets <- function(
    asv_sets,
    comparison_keys,
    included_groups_vec,
    ont_minsize_values = c(2, 4, 6, 8),
    pacbio_minsize_value = 8
) {
  
  # Convert group names to the read_label format used in asv_sets
  depth_labels <- tibble(
    included_group = included_groups_vec,
    read_label = case_when(
      included_group == "sample_size_all_reads" ~ "all",
      TRUE ~ stringr::str_remove(
        included_group,
        "^sample_size_"
      )
    )
  ) %>%
    mutate(
      expected_depth = as.numeric(
        stringr::str_extract(read_label, "\\d+")
      )
    )
  
  # Total available reads are independent of minsize
  all_read_depths <- asv_sets %>%
    filter(
      read_label == "all" |
        is_all %in% TRUE
    ) %>%
    group_by(
      primer,
      sample_type,
      platform
    ) %>%
    summarise(
      all_read_depth = if (
        all(is.na(actual_depth))
      ) {
        NA_real_
      } else {
        max(actual_depth, na.rm = TRUE)
      },
      .groups = "drop"
    )
  
  base_keys <- comparison_keys %>%
    distinct(primer, sample_type)
  
  # Define all required platform × minsize combinations
  analysis_keys <- bind_rows(
    tidyr::expand_grid(
      base_keys,
      platform = "ONT",
      minsize = ont_minsize_values
    ),
    tidyr::expand_grid(
      base_keys,
      platform = "PacBio",
      minsize = pacbio_minsize_value
    )
  )
  
  # Construct the complete expected dataset grid
  expected_grid <- tidyr::expand_grid(
    analysis_keys,
    depth_labels
  ) %>%
    left_join(
      all_read_depths,
      by = c(
        "primer",
        "sample_type",
        "platform"
      )
    ) %>%
    
    # Do not create an exact numerical subset above the number
    # of available reads. These comparisons should be capped.
    filter(
      read_label == "all" |
        expected_depth <= all_read_depth
    )
  
  expected_grid %>%
    left_join(
      asv_sets,
      by = c(
        "primer",
        "sample_type",
        "platform",
        "minsize",
        "read_label"
      )
    ) %>%
    mutate(
      subset_depth = coalesce(
        subset_depth,
        expected_depth
      ),
      
      actual_depth = coalesce(
        actual_depth,
        if_else(
          read_label == "all",
          all_read_depth,
          expected_depth
        )
      ),
      
      is_all = coalesce(
        is_all,
        read_label == "all"
      ),
      
      depth_type = coalesce(
        depth_type,
        if_else(
          read_label == "all",
          "All available reads",
          "Numerical subset"
        )
      ),
      
      # Missing datasets correspond to zero inferred ASVs
      asvs = purrr::map(
        asvs,
        \(x) {
          if (
            is.null(x) ||
            length(x) == 0 ||
            all(is.na(x))
          ) {
            character(0)
          } else {
            x
          }
        }
      ),
      
      n_asvs = lengths(asvs),
      
      dataset = coalesce(
        dataset,
        paste(
          primer,
          sample_type,
          platform,
          paste0("minsize_", minsize),
          included_group,
          sep = "__"
        )
      )
    ) %>%
    select(
      primer,
      sample_type,
      platform,
      minsize,
      read_label,
      subset_depth,
      actual_depth,
      is_all,
      depth_type,
      asvs,
      n_asvs,
      dataset
    ) %>%
    arrange(
      primer,
      sample_type,
      platform,
      minsize,
      actual_depth
    )
}


asv_sets_complete <- complete_asv_sets(
  asv_sets = asv_sets,
  comparison_keys = comparison_keys,
  included_groups_vec = included_groups_vec,
  ont_minsize_values = c(2, 4, 6, 8),
  pacbio_minsize_value = 8
)



build_minsize_comparison <- function(
    ont_minsize_value,
    asv_sets = asv_sets_complete,
    pacbio_minsize_value = 8
) {
  
  
  # ----------------------------------------------------------
  # Perform the comparisons
  # ----------------------------------------------------------
  
  purrr::pmap_dfr(
    comparison_keys,
    
    function(primer, sample_type) {
      
      primer_value <- primer
      sample_value <- sample_type
      
      
      # ONT datasets for the current minsize
      ont_data <- asv_sets %>%
        filter(
          .data$primer == .env$primer_value,
          .data$sample_type == .env$sample_value,
          platform == "ONT",
          minsize == ont_minsize_value
        )
      
      
      # PacBio reference datasets at minsize 8
      pacbio_data <- asv_sets %>%
        filter(
          .data$primer == .env$primer_value,
          .data$sample_type == .env$sample_value,
          platform == "PacBio",
          minsize == pacbio_minsize_value
        )
      
      
      # ------------------------------------------------------
      # Candidate comparison depths
      #
      # This is the union of all available ONT and PacBio
      # actual depths.
      # ------------------------------------------------------
      
      requested_depths <- sort(
        unique(
          c(
            ont_data$actual_depth,
            pacbio_data$actual_depth
          )
        )
      )
      
      
      # ------------------------------------------------------
      # Compare platforms at each candidate depth
      # ------------------------------------------------------
      
      purrr::map_dfr(
        requested_depths,
        
        function(requested_depth) {
          
          ont_selected <- pick_depth_set(
            platform_data = ont_data,
            requested_depth = requested_depth
          )
          
          pacbio_selected <- pick_depth_set(
            platform_data = pacbio_data,
            requested_depth = requested_depth
          )
          
          
          # Skip depths where one platform has neither:
          # 1. an exact-depth subset, nor
          # 2. a capped all_reads dataset
          if (
            is.null(ont_selected) ||
            is.null(pacbio_selected)
          ) {
            
            return(tibble())
          }
          
          
          # ASVs present in both ONT and PacBio
          shared_asvs <- intersect(
            ont_selected$asvs,
            pacbio_selected$asvs
          )
          
          
          # ASVs present in either platform
          union_asvs <- union(
            ont_selected$asvs,
            pacbio_selected$asvs
          )
          
          
          # Classify whether the comparison is exact or capped
          comparison_type <- case_when(
            ont_selected$capped &&
              pacbio_selected$capped ~ "Both capped",
            
            ont_selected$capped ~ "ONT capped",
            
            pacbio_selected$capped ~ "PacBio capped",
            
            TRUE ~ "Exact depth"
          )
          
          
          tibble(
            primer = primer_value,
            sample_type = sample_value,
            
            ont_minsize = ont_minsize_value,
            pacbio_minsize = pacbio_minsize_value,
            
            # Position used on the x-axis
            requested_depth = requested_depth,
            
            # Actual reads represented by the selected datasets
            ont_depth_used = ont_selected$actual_depth,
            pacbio_depth_used = pacbio_selected$actual_depth,
            
            ont_is_all = ont_selected$is_all,
            pacbio_is_all = pacbio_selected$is_all,
            
            ont_capped = ont_selected$capped,
            pacbio_capped = pacbio_selected$capped,
            
            comparison_type = comparison_type,
            
            # Total ASVs recovered by each platform
            n_ont_all = ont_selected$n_asvs,
            n_pacbio_reference = pacbio_selected$n_asvs,
            
            # ONT ASVs also found in PacBio
            n_shared = length(shared_asvs),
            
            # Percentage of PacBio ASVs recovered by ONT
            pacbio_recovery_pct = if_else(
              pacbio_selected$n_asvs > 0,
              100 * length(shared_asvs) /
                pacbio_selected$n_asvs,
              NA_real_
            ),
            
            # Percentage of ONT ASVs supported by PacBio
            ont_supported_pct = if_else(
              ont_selected$n_asvs > 0,
              100 * length(shared_asvs) /
                ont_selected$n_asvs,
              NA_real_
            ),
            
            
            # Retain original dataset names for diagnostics
            ont_dataset = ont_selected$dataset,
            pacbio_dataset = pacbio_selected$dataset
          )
        }
      )
    }
  )
}







# ============================================================
# Identify all available ONT minsize values
# ============================================================

ont_minsizes <- asv_sets_complete %>%
  filter(
    platform == "ONT"
  ) %>%
  distinct(minsize) %>%
  arrange(minsize) %>%
  pull(minsize)

ont_minsizes


# ============================================================
# Compare every ONT minsize with PacBio minsize 8
# ============================================================

comparison_results <- purrr::map_dfr(
  ont_minsizes,
  
  ~ build_minsize_comparison(
    ont_minsize_value = .x,
    asv_sets = asv_sets_complete,
    pacbio_minsize_value = 8
  )
)

comparison_results %>%
  select(
    primer,
    sample_type,
    ont_minsize,
    requested_depth,
    ont_depth_used,
    pacbio_depth_used,
    comparison_type,
    n_ont_all,
    n_pacbio_reference,
    n_shared,
    pacbio_recovery_pct
  ) %>%
  arrange(
    ont_minsize,
    primer,
    sample_type,
    requested_depth
  ) %>%
  print(
    n = 100,
    width = Inf
  )


# Number of exact and capped comparisons
comparison_results %>%
  count(
    ont_minsize,
    primer,
    comparison_type
  ) %>%
  arrange(
    ont_minsize,
    primer,
    comparison_type
  ) %>%
  print(n = Inf)





# 
# plot_data <- ASV_matches_complete %>% 
#   mutate(basecalling  = factor(basecalling, levels =c("pacbio", "super"), labels = c("Pacbio", "ONT")),
#          sample_content = factor(sample_content, levels =c("Zmock", "Zfecal", "AD", "WWTP", "Soil")), 
#          Basecalling_Matchtype = factor(paste0(basecalling, "_",ASV_type),
#                                         levels = c(
#                                           "ONT_shared", 
#                                           "Pacbio_total",
#                                           "ONT_total",
#                                           "Pacbio_shared"),
#                                         labels = c(
#                                           "ONT (exact match)", 
#                                           "Pacbio (all)",
#                                           "ONT (all)",
#                                           "Pacbio_shared")),
#          primerset = factor(primerset, levels = c("V4", "V1V3", "V1V8", "OPR"),
#                             labels = c("V4", "V1-V3", "V1-V8", "OPR")), 
#   ) %>% 
#   filter(Basecalling_Matchtype != "Pacbio_shared") 
# 
# 










# ============================================================
# Plot one ONT minsize against PacBio minsize 8
# ============================================================
custom_colors_fill = c(
  "white",
  #"#0077B6",
  "white","#0077B6",
  "#CC007A",
  "#0077B6"
  
)

custom_colors_color = c(
  "#0077B6",  
  "#CC007A",
  "#0077B6",
  "#CC007A",
  "#0077B6",
  "#0077B6"
)








library(tidyverse)
library(ggh4x)

# ============================================================================
# Transform comparison_results into the previous plot_data structure
# ============================================================================


custom_breaks_minimal=c(0, 10, 100, 1000,  10000,100000, 1000000) # Mortens suggestion


make_comparison_plot_data <- function(
    comparison_results,
    primer_select = "V1V8",
    ont_minsize_value = c(2,4,6,8),
    pacbio_minsize_value = 8
) {
  
  comparison_results %>%
    filter(
      ont_minsize %in% c(ont_minsize_value),
      pacbio_minsize == pacbio_minsize_value,
      sample_type %in% c("Zfecal", "AD", "WWTP", "Soil"),
      primer == primer_select
    ) %>%
    
    # Convert the three abundance measures to long format
    pivot_longer(
      cols = c(
        n_shared,
        n_pacbio_reference,
        n_ont_all
      ),
      names_to = "metric",
      values_to = "n_ASV"
    ) %>%
    
    mutate(
      # ----------------------------------------------------------------------
      # Platform and ASV category
      # ----------------------------------------------------------------------
      
      basecalling_raw = case_when(
        metric == "n_pacbio_reference" ~ "pacbio",
        metric %in% c("n_shared", "n_ont_all") ~ "super"
      ),
      
      ASV_type = case_when(
        metric == "n_shared" ~ "shared",
        metric %in% c(
          "n_pacbio_reference",
          "n_ont_all"
        ) ~ "total"
      ),
      
      # Use the actual number of reads used for each platform
      reads_used = case_when(
        basecalling_raw == "pacbio" ~ pacbio_depth_used,
        basecalling_raw == "super"  ~ ont_depth_used
      ),
      
      # ----------------------------------------------------------------------
      # Factors and labels
      # ----------------------------------------------------------------------
      
      basecalling = factor(
        basecalling_raw,
        levels = c("pacbio", "super"),
        labels = c("PacBio", "ONT")
      ),
      
      sample_content = factor(
        sample_type,
        levels = c("Zfecal", "AD", "WWTP", "Soil"),
        labels = c("Fecal", "AD", "WWTP", "Soil")
      ),
      
      Basecalling_Matchtype = factor(
        metric,
        levels = c(
          "n_shared",
          "n_pacbio_reference",
          "n_ont_all"
        ),
        labels = c(
          "ONT (exact match)",
          "PacBio (all)",
          "ONT (all)"
        )
      ),
      
      primerset = factor(
        primer,
        levels = c("V4", "V1V3", "V1V8", "OPR"),
        labels = c("V4", "V1-V3", "V1-V8", "OPR")
      )
    ) %>%
    
    # When one platform is capped, several requested depths may correspond
    # to the same actual depth and ASV count. Remove these duplicate points.
    distinct(
      primerset,
      sample_content,
      ont_minsize,
      pacbio_minsize,
      Basecalling_Matchtype,
      reads_used,
      n_ASV,
      .keep_all = TRUE
    ) %>%
    
    arrange(
      sample_content,
      primerset,
      Basecalling_Matchtype,
      reads_used
    )
}




# ============================================================================
# Plot function
# ============================================================================

plot_complex_communities <- function(
    plot_data
    ) {
  
  
  comparison_line_data <- comparison_line_data %>%
    group_by(sample_content, ont_minsize) %>%
    mutate(
      panel_max = max(plot_data$n_ASV[
        plot_data$sample_content == first(sample_content) &
          plot_data$ont_minsize == first(ont_minsize)
      ], na.rm = TRUE),
      
      label_x = sqrt(x_pacbio * x_ont),  # midpoint on log scale
      label_y = pmax(y_pacbio, y_ont) + 0.07 * panel_max,
      difference = y_pacbio - y_ont,
      label = paste0("\u0394ASVs = ", difference)
    ) %>%
    ungroup()
  
  
  line_colors <- c(
    "ONT (exact match)" = "#0077B6",
    "PacBio (all)"      = "#CC007A",
    "ONT (all)"         = "#56B4E9"
  )
  
  point_fills <- c(
    "ONT (exact match)" = "#0077B6",
    "PacBio (all)"      = "white",
    "ONT (all)"         = "white"
  )
  
  ggplot(
    plot_data,
    aes(
      x = reads_used,
      y = n_ASV,
      group = Basecalling_Matchtype,
      color = Basecalling_Matchtype
    )
  ) +
    
    # Lines
    geom_line(
      linewidth = 0.7,
      alpha = 0.9
    ) +
    
    # Filled point: ONT ASVs matching PacBio
    # Open points: total PacBio and total ONT ASVs
    geom_point(
      aes(fill = Basecalling_Matchtype),
      shape = 21,
      size = 2.8,
      stroke = 0.8
    ) +
    
    ggh4x::facet_grid2(
      rows = vars(sample_content),
      cols = vars(ont_minsize),
      scales = "free_y",
      labeller = labeller(
        ont_minsize = \(x) paste0("minsize: ", x)
      )
    ) +
    
    ############ Comparosic
    # geom_segment(
    #   data = comparison_line_data,
    #   aes(
    #     x = x_pacbio,
    #     y = y_pacbio,
    #     xend = x_ont,
    #     yend = y_ont
    #   ),
    #   inherit.aes = FALSE,
    #   colour = "#4D4D4D",
    #   linewidth = 0.8,
    #   linetype = "dashed"
    # ) +
    # geom_text_repel(
    #   data = comparison_line_data,
    #   aes(
    #     x = sqrt(x_pacbio * x_ont),   # midpoint on log scale
    #     y = (y_pacbio + y_ont) / 2,
    #     label = label
    #   ),
    #   inherit.aes = FALSE,
    #   color = "black",
    #   size = 3
    # ) +
    ###################
  scale_x_log10(
    breaks = c(custom_breaks_minimal),#custom_breaks_minimal,
    labels = c(format(custom_breaks_minimal, big.mark = ",", scientific = FALSE))
  ) +
    
    scale_color_manual(
      values = line_colors,
      drop = FALSE
    ) +
    
    scale_fill_manual(
      values = point_fills,
      drop = FALSE,
      guide = "none"
    ) +
    
    guides(
      color = guide_legend(
        title = NULL,
        override.aes = list(
          shape = 21,
          size = 3,
          fill = unname(point_fills)
        )
      )
    ) +
    
    labs(
      x = "Number of reads",
      y = "Number of ASVs"
    ) +
    
    theme_marie_clean(
      bg = "#ECECEC"
    ) +
    
    theme(
      axis.title = element_text(
        face = "bold"
      ),
      
      axis.text.x = element_text(
        angle = 45,
        hjust = 1,
        vjust = 1
      ),
      
      strip.text.x = element_text(hjust = 0,
        face = "bold"
      ),
      
      strip.text.y = element_text(
        face = "bold"
      ),
      
      # strip.background = element_rect(
      #   fill = "grey90",
      #   color = NA
      # ),
      
      legend.position = "bottom",
      
      panel.grid.minor = element_blank()
    )
}




# ============================================================================
# Make the plot
# ============================================================================
plot_data_V1V8 <- make_comparison_plot_data(
  comparison_results = comparison_results,
  primer_select = "V1V8",
  ont_minsize_value = c(2,4,6,8),
  pacbio_minsize_value = 8
)

comparison_line_data <- plot_data_V1V8 %>%
  filter(Basecalling_Matchtype %in% c("PacBio (all)", "ONT (exact match)")) %>%
  group_by(primerset, sample_content, ont_minsize, Basecalling_Matchtype) %>%
  slice_max(order_by = reads_used, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(primerset, sample_content, ont_minsize, Basecalling_Matchtype, reads_used, n_ASV) %>%
  pivot_wider(
    names_from = Basecalling_Matchtype,
    values_from = c(reads_used, n_ASV)
  ) %>%
  rename(
    x_pacbio = `reads_used_PacBio (all)`,
    y_pacbio = `n_ASV_PacBio (all)`,
    x_ont    = `reads_used_ONT (exact match)`,
    y_ont    = `n_ASV_ONT (exact match)`
  )

depth_plot_V1V8 <- plot_complex_communities(
  plot_data = plot_data_V1V8
)

depth_plot_V1V8
plot_name="Fig_S10_minsize_complex_V1V8"
plot_folder=folder_dir_supp
ggsave(paste0(plot_dir,"/" ,plot_folder, "/", plot_name, ".png"), width = 12, height = 9)
ggsave(paste0(plot_dir,"/" ,plot_folder, "/", plot_name, ".svg"), width = 12, height = 9)


plot_data_OPR <- make_comparison_plot_data(
  comparison_results = comparison_results,
  primer_select = "OPR",
  ont_minsize_value = c(2,4,6,8),
  pacbio_minsize_value = 8
)


comparison_line_data <- plot_data_OPR %>%
  filter(Basecalling_Matchtype %in% c("PacBio (all)", "ONT (exact match)")) %>%
  group_by(primerset, sample_content, ont_minsize, Basecalling_Matchtype) %>%
  slice_max(order_by = reads_used, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(primerset, sample_content, ont_minsize, Basecalling_Matchtype, reads_used, n_ASV) %>%
  pivot_wider(
    names_from = Basecalling_Matchtype,
    values_from = c(reads_used, n_ASV)
  ) %>%
  rename(
    x_pacbio = `reads_used_PacBio (all)`,
    y_pacbio = `n_ASV_PacBio (all)`,
    x_ont    = `reads_used_ONT (exact match)`,
    y_ont    = `n_ASV_ONT (exact match)`
  )

depth_plot <- plot_complex_communities(
  plot_data = plot_data_OPR
)

depth_plot
plot_name="Fig_S11_minsize_complex_OPR"
plot_folder=folder_dir_supp
ggsave(paste0(plot_dir,"/" ,plot_folder, "/", plot_name, ".png"), width = 12, height = 9)
ggsave(paste0(plot_dir,"/" ,plot_folder, "/", plot_name, ".svg"), width = 12, height = 9)



