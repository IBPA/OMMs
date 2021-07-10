#' Commonly used constants and functions
  
# Constants
MONO_COLUMNS <- c('Glucose', 'Galactose', 'Fructose', 'Xylose', 'Arabinose', 'Fucose', 'Rhamnose', 'GlcA', 'GalA', 'GlcNAc', 'GalNAc', 'Mannose', 'Allose', 'Ribose')

# Load the food data. The 'data_file' contains the food descriptions along with 
#  their glycan content. The 'metadata_file' contains the %water content.
load_food_data <- function(data_file, metadata_file){
  # A) Load data
  df <- read.table(data_file, sep = ",", quote = '"', header = TRUE)
  df_data <- df[,MONO_COLUMNS]
  df_data[is.na(df_data)] <- 0
  df[,MONO_COLUMNS] <- df_data
  colnames(df)[4] <- "uid"
  
  # B) Load water content metadata
  df_meta <- read_excel(metadata_file)
  df_meta <- df_meta[,c("Unique Identifier Code","% water content")]
  colnames(df_meta) <- c("uid", "water")
  df_meta$uid <- as.numeric(str_remove(df_meta$uid, "^0+"))
  
  # C) Validate water content data
  uids_invalid_water_content <- unlist(df_meta[is.na(df_meta$water) | df_meta$water<0 | df_meta$water>100, c("uid") ], use.names = FALSE)
  if (length(uids_invalid_water_content)>0){
    uids_str <- paste(uids_invalid_water_content, collapse = ",")
    warning(sprintf("Some foods have invalid water content (not between 0-100). Corresponding records with missing values will be removed and the water content of the rest will be clamped to be in range. The food uids are: \n%s", uids_str))
  }
  
  # D) Preprocess the water content data
  df_meta$dry <- (100-df_meta$water)/100
  df_meta$dry <- pmin(pmax(df_meta$dry, 0.0), 1.0) # Clamp the values
  df_meta <- df_meta[!is.na(df_meta$dry),] # Remove missing value records
  df_meta$water <- NULL
  
  # E) Incorporate dry content % of each food into the corresponding glycan content vector
  df_comb <- merge(df, df_meta, by = "uid")
  df_comb[,MONO_COLUMNS] <- df_comb[,MONO_COLUMNS] * df_comb$dry
  
  return(df_comb[,c("uid", MONO_COLUMNS)])
}

# Perform minimax normalization on each column
minmax_normalize <- function(df_data){
  scale_minmax <- function(x){
    (x- min(x)) /(max(x)-min(x))
  }
  df_data_scaled <- as.data.frame(lapply(df_data, scale_minmax))
  df_data_scaled[is.na(df_data_scaled)] <- 0
  
  return(df_data_scaled)
}

# Save mixed-meals into the filename.
save_MMs <- function(df_food_proportions, filename) {
  # A) Find the nonzero food proportions in each mixed-meal
  nz_proportions <- apply(df_food_proportions, 1, FUN = function(x){x[x!=0]})
  
  # B) Create One line for each mixed meal
  lines <- c()
  for (meal_idx in 1:length(nz_proportions)) {
    cur_meal_props <- nz_proportions[[meal_idx]]
    line_str <- paste(sprintf("%s:%.2f", names(cur_meal_props), unname(cur_meal_props)), collapse = ",")
    lines <- c(lines, line_str)
  }
  
  # C) Save the lines
  writeLines(lines, filename)
}

# Return the friendly names of the provided food 'uids' from the 'data_file'
get_friendly_food_names <- function(uids, data_file){
  df <- read.table(data_file, sep = ",", quote = '"', header = TRUE)
  df_names <- df[df$Unique.Identifier.Code %in% uids, c("Unique.Identifier.Code", "Simple.name")]
  rownames(df_names) <- df_names$Unique.Identifier.Code
  return(unname(df_names[uids, c("Simple.name")]))
}

# Load mixed-meals from 'filename'.
load_MMs <- function(filename) {
  # A) Read all mixed-meals (one in each line) and extract all the food uids used.
  all_lines <- str_split(readLines(filename), ",")
  uids <- unique(unlist(unname(data.frame(str_split(unlist(all_lines), ":"))[1,])))
  num_meals <- length(all_lines)
  
  # B) Construct the matrix, parse each line, and fill the matrix
  MMs <- matrix(0, nrow = num_meals, ncol = length(uids), dimnames = list(1:num_meals, uids))
  for (i in 1:num_meals){
    line_cur <- all_lines[[i]]
    df_cur <- data.frame(str_split(line_cur, ":"), stringsAsFactors = FALSE)
    fids_cur <- unlist(unname(df_cur[1,]))
    proportions_cur <- unname(df_cur[2,])
    MMs[i,fids_cur] <- as.numeric(unlist(proportions_cur))
  }
  
  return(MMs)
}