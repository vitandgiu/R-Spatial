# ============================================================
# Functions in GeometryOperations
# Author: Andrea Giuseppe Vitali
# Date: 13 Nov 2025
# check_and_fix_validity(geom_sf)
# Description: Check and Fix Topological Geometries irregularities
#
# AdList2Admat(Adlist, type_distance)
# Description: Convert adjacet list into adhjacent matrix
#
# calcola_spatWeighNetMat(pnts, ID = NULL, nb_neighbors, distMetric = "duration",
# osrmServer = "https://routing.openstreetmap.de/",
# osrmProfile = "car")
# Description: From points list and connection matrix to travel time OD weighted matrix
# 
# lw_to_graph(lw, sf_obj = NULL, save_path = NULL, directed = FALSE)
# Description: Converts a spatial weights object (listw)
#              into an igraph network, with coordinates
#              automatically extracted from an sf object.
# ============================================================

library(sf)
library(osrm)
library(spdep)

### RENAME ACTIVE GEOMETRY
rename_geometry <- function(x, name) {
  current <- attr(x, "sf_column")
  names(x)[names(x) == current] <- name
  st_geometry(x) <- name
  x
}

### Validity CHECK
# Check validity of geometries
check_and_fix_validity <- function(geom_sf) {
  # Check validity of each geometry with reasons
  valid_checks     <- st_is_valid(geom_sf)
  validity_reasons <- st_is_valid(geom_sf, reason = TRUE)
  
  cat("Checking topological validity of geometries...\n")
  
  if (all(valid_checks)) {
    cat("✅ All geometries are topologically valid!\n")
    return(geom_sf)  # Return original if all valid
  } else {
    # Identify invalid geometries and their reasons
    invalid_indices <- which(!valid_checks)
    cat("❌ The following geometries are NOT valid (by index):\n")
    print(invalid_indices)
    
    cat("Reasons for invalidity:\n")
    print(validity_reasons[invalid_indices])
    
    cat("Summary of invalid geometries:\n")
    print(geom_sf[invalid_indices, ])
    
    # Attempt to fix invalid geometries
    cat("Attempting to fix invalid geometries with st_make_valid()...\n")
    geom_fixed <- st_make_valid(geom_sf)
    
    # Check validity again after fixing
    fixed_checks <- st_is_valid(geom_fixed)
    if (all(fixed_checks)) {
      cat("✔ All geometries are valid after applying st_make_valid()\n")
    } else {
      still_invalid <- which(!fixed_checks)
      cat("⚠️ Warning! The following geometries are still invalid after fixing:\n")
      print(still_invalid)
    }
    
    return(geom_fixed)  # Return fixed geometries (partially or fully fixed)
  }
}

compare_crs <- function(obj1, obj2, verbose = TRUE) {
  # Extract CRS objects
  crs1 <- st_crs(obj1)
  crs2 <- st_crs(obj2)
  
  # Initialize report
  report <- list(
    objects = c(deparse(substitute(obj1)), deparse(substitute(obj2))),
    identical_crs = identical(crs1, crs2),
    sf_equal = (crs1 == crs2),
    epsg_match = !is.null(crs1$epsg) && !is.null(crs2$epsg) && crs1$epsg == crs2$epsg,
    wkt_match = identical(crs1$wkt, crs2$wkt),
    can_transform = st_can_transform(crs1, crs2),
    transform_identity = FALSE  # Will compute if needed
  )
  
  # If not identical, check if transformation is identity (expensive but definitive)
  if (!report$identical_crs && report$sf_equal) {
    report$transform_identity <- st_is_longlat(st_crs(crs1)) == st_is_longlat(st_crs(crs2))
  }
  
  # Pretty report
  if (verbose) {
    cat(sprintf("\nDetails:\n  Obj1 input: %s\n  Obj2 input: %s\n", 
                crs1$input, crs2$input), "\n")
    
    cat("CRS Comparison Report\n")
    cat("=====================\n")
    cat(sprintf("Objects: %s vs %s\n\n", report$objects[1], report$objects[2]))
    
    cat("Strict Checks:\n")
    cat(sprintf("  identical(): %s\n", ifelse(report$identical_crs, "✓ PASS", "✗ FAIL")))
    cat(sprintf("  sf::==     : %s\n", ifelse(report$sf_equal, "✓ PASS", "✗ FAIL")))
    
    cat("\nSemantic Checks:\n")
    cat(sprintf("  EPSG match : %s (%s vs %s)\n", 
                ifelse(report$epsg_match, "✓", "✗"), 
                crs1$epsg %||% "NULL", crs2$epsg %||% "NULL"))
    cat(sprintf("  WKT match  : %s\n", ifelse(report$wkt_match, "✓ PASS", "✗ FAIL")))
    cat(sprintf("  Transform? : %s\n", ifelse(report$can_transform, "✓ (no reproj needed)", "✗ DIFFERENT")))
    
    cat("\nRecommendation:\n")
    if (report$sf_equal && report$can_transform) {
      cat("✅ CRS are **FUNCTIONALLY IDENTICAL**. Safe for overlay/join.\n")
    } else if (report$can_transform) {
      cat("⚠️  CRS differ but transformable. Use st_transform() before analysis.\n")
    } else {
      cat("❌ CRS fundamentally incompatible.\n")
    }
  }
  
  return(invisible(report))
}

# AdList2Admat <- function(Adlist, type_distance) {
#   # Convert to data frame if not already
#   df <- as.data.frame(Adlist)
# 
#   # Check if required columns exist
#   if (!("Idx" %in% colnames(df))) {
#     stop("Input Adlist must contain 'Idx' column")
#   }
#   if (!(type_distance %in% colnames(df))) {
#     stop(paste("Input Adlist must contain", type_distance, "column"))
#   }
# 
#   # Extract origin and destination indices from 'Idx' column
#   # Assuming 'Idx' has format "origin_destination" such as "12_34"
#   split_idx <- strsplit(as.character(df$Idx), "_")
#   origins <- as.numeric(sapply(split_idx, `[`, 1))
#   destinations <- as.numeric(sapply(split_idx, `[`, 2))
# 
#   # Determine matrix size as max index found
#   max_idx <- max(c(origins, destinations))
# 
#   # Initialize empty square matrix
#   adjacency_matrix <- matrix(0, nrow = max_idx, ncol = max_idx)
# 
#   # Assign distances to matrix
#   for (i in seq_along(origins)) {
#     adjacency_matrix[origins[i], destinations[i]] <- df[[type_distance]][i]
#   }
# 
#   return(adjacency_matrix)
# }
# 
# AdList2Admat <- function(Adlist,
#                          type_distance,
#                          idx_col = "Idx",
#                          fill_value = 0,
#                          allow_duplicate = c("last", "min", "max", "mean", "first"),
#                          symmetric = FALSE,
#                          ids = NULL) {
#   allow_duplicate <- match.arg(allow_duplicate)
#   
#   df <- st_drop_geometry(Adlist) |> as.data.frame()
#   
#   if (!idx_col %in% names(df)) {
#     stop(sprintf("Input must contain '%s' column.", idx_col))
#   }
#   if (!type_distance %in% names(df)) {
#     stop(sprintf("Input must contain '%s' column.", type_distance))
#   }
#   
#   idx_chr <- as.character(df[[idx_col]])
#   dist_val <- suppressWarnings(as.numeric(df[[type_distance]]))
#   
#   bad_dist <- is.na(dist_val)
#   if (any(bad_dist)) {
#     warning(sprintf("Dropping %d rows with NA/non-numeric '%s'.",
#                     sum(bad_dist), type_distance))
#   }
#   
#   mat_idx <- stringr::str_match(idx_chr, "^(\\d+)_(\\d+)$")
#   bad_idx <- is.na(mat_idx[, 1])
#   
#   if (any(bad_idx)) {
#     warning(sprintf("Dropping %d malformed '%s' values; expected 'origin_destination'.",
#                     sum(bad_idx), idx_col))
#   }
#   
#   keep <- !(bad_idx | bad_dist)
#   if (!any(keep)) {
#     stop("No valid origin-destination rows found after validation.")
#   }
#   
#   origins <- as.integer(mat_idx[keep, 2])
#   destinations <- as.integer(mat_idx[keep, 3])
#   values <- dist_val[keep]
#   
#   od_df <- data.frame(
#     origin = origins,
#     destination = destinations,
#     value = values
#   )
#   
#   od_df <- od_df |>
#     group_by(origin, destination) |>
#     summarise(
#       value = dplyr::case_when(
#         allow_duplicate == "last"  ~ dplyr::last(value),
#         allow_duplicate == "first" ~ dplyr::first(value),
#         allow_duplicate == "min"   ~ min(value, na.rm = TRUE),
#         allow_duplicate == "max"   ~ max(value, na.rm = TRUE),
#         allow_duplicate == "mean"  ~ mean(value, na.rm = TRUE)
#       ),
#       .groups = "drop"
#     )
#   
#   if (is.null(ids)) {
#     ids <- sort(unique(c(od_df$origin, od_df$destination)))
#   } else {
#     ids <- sort(unique(as.integer(ids)))
#   }
#   
#   id_map <- setNames(seq_along(ids), ids)
#   i <- unname(id_map[as.character(od_df$origin)])
#   j <- unname(id_map[as.character(od_df$destination)])
#   
#   valid_map <- !(is.na(i) | is.na(j))
#   if (!all(valid_map)) {
#     warning(sprintf("Dropping %d OD rows whose IDs are not in the supplied 'ids'.",
#                     sum(!valid_map)))
#     i <- i[valid_map]
#     j <- j[valid_map]
#     od_df <- od_df[valid_map, , drop = FALSE]
#   }
#   
#   n <- length(ids)
#   adjacency_matrix <- matrix(fill_value, nrow = n, ncol = n,
#                              dimnames = list(as.character(ids), as.character(ids)))
#   adjacency_matrix[cbind(i, j)] <- od_df$value
#   
#   if (symmetric) {
#     adjacency_matrix[cbind(j, i)] <- od_df$value
#   }
#   
#   adjacency_matrix
# }

# AdList2Admat <- function(Adlist,
#                               type_distance,
#                               idx_col = "Idx",
#                               fill_value = 0,
#                               allow_duplicate = c("last", "min", "max", "mean", "first"),
#                               symmetric = FALSE,
#                               ids = NULL) {
#   
#   allow_duplicate <- match.arg(allow_duplicate)
#   
#   # ---- drop geometry if present (fast base check) ----
#   if (inherits(Adlist, "sf")) {
#     df <- sf::st_drop_geometry(Adlist)
#   } else {
#     df <- as.data.frame(Adlist)
#   }
#   
#   # ---- column checks ----
#   if (!idx_col %in% names(df)) stop("Missing idx_col")
#   if (!type_distance %in% names(df)) stop("Missing type_distance")
#   
#   idx_chr <- df[[idx_col]]
#   val <- as.numeric(df[[type_distance]])
#   
#   # ---- fast validity mask ----
#   ok <- !is.na(val)
#   
#   # fast parse "origin_destination" WITHOUT regex
#   split <- strsplit(idx_chr[ok], "_", fixed = TRUE)
#   
#   origins <- as.integer(vapply(split, `[`, 1, FUN.VALUE = character(1)))
#   dests   <- as.integer(vapply(split, `[`, 2, FUN.VALUE = character(1)))
#   
#   # ---- drop malformed ----
#   good <- !is.na(origins) & !is.na(dests)
#   origins <- origins[good]
#   dests   <- dests[good]
#   values  <- val[ok][good]
#   
#   if (length(values) == 0) stop("No valid OD pairs")
#   
#   # ---- FAST duplicate handling (no dplyr) ----
#   key <- paste(origins, dests, sep = "_")
#   
#   if (allow_duplicate != "first") {
#     if (allow_duplicate == "last") {
#       ord <- order(key)  # keep last by overwriting
#       origins <- origins[ord]
#       dests   <- dests[ord]
#       values  <- values[ord]
#       key     <- key[ord]
#     }
#     
#     # aggregate via split-index trick (VERY fast base R)
#     agg <- switch(allow_duplicate,
#                   
#                   min  = tapply(values, key, min),
#                   max  = tapply(values, key, max),
#                   mean = tapply(values, key, mean),
#                   first = values  # already ordered case handled
#     )
#     
#     if (allow_duplicate != "first") {
#       key <- names(agg)
#       values <- as.numeric(agg)
#       
#       split_k <- strsplit(key, "_", fixed = TRUE)
#       origins <- as.integer(vapply(split_k, `[`, 1, FUN.VALUE = character(1)))
#       dests   <- as.integer(vapply(split_k, `[`, 2, FUN.VALUE = character(1)))
#     }
#   }
#   
#   # ---- ID handling ----
#   if (is.null(ids)) {
#     ids <- sort(unique(c(origins, dests)))
#   } else {
#     ids <- sort(as.integer(ids))
#   }
#   
#   id_map <- integer(max(ids))
#   id_map[ids] <- seq_along(ids)
#   
#   i <- id_map[origins]
#   j <- id_map[dests]
#   
#   ok2 <- i > 0 & j > 0
#   i <- i[ok2]
#   j <- j[ok2]
#   values <- values[ok2]
#   
#   # ---- matrix build (fast indexing) ----
#   n <- length(ids)
#   mat <- matrix(fill_value, n, n,
#                 dimnames = list(as.character(ids), as.character(ids)))
#   
#   mat[cbind(i, j)] <- values
#   
#   if (symmetric) {
#     mat[cbind(j, i)] <- values
#   }
#   
#   mat
# }

AdList2Admat <- function(Adlist,
                         type_distance = "duration",
                         fill_value = 0,
                         allow_duplicate = c("last", "min", "max", "mean", "first"),
                         symmetric = FALSE,
                         ids = NULL) {
  
  allow_duplicate <- match.arg(allow_duplicate)
  
  # ---- drop geometry if sf ----
  if (inherits(Adlist, "sf")) {
    df <- sf::st_drop_geometry(Adlist)
  } else {
    df <- as.data.frame(Adlist)
  }
  
  # ---- checks ----
  if (!all(c("i", "j") %in% names(df))) {
    stop("Columns 'i' and 'j' must be present in Adlist")
  }
  if (!type_distance %in% names(df)) {
    stop("Missing type_distance column")
  }
  
  origins <- as.integer(df$i)
  dests   <- as.integer(df$j)
  values  <- as.numeric(df[[type_distance]])
  
  # ---- remove NA ----
  ok <- !is.na(origins) & !is.na(dests) & !is.na(values)
  origins <- origins[ok]
  dests   <- dests[ok]
  values  <- values[ok]
  
  if (length(values) == 0) stop("No valid OD pairs")
  
  # ---- duplicates handling (FAST, no strings) ----
  key <- origins * 1e9 + dests  # fast unique key (no paste)
  
  if (allow_duplicate == "last") {
    keep <- !duplicated(key, fromLast = TRUE)
    origins <- origins[keep]
    dests   <- dests[keep]
    values  <- values[keep]
    
  } else if (allow_duplicate == "first") {
    keep <- !duplicated(key)
    origins <- origins[keep]
    dests   <- dests[keep]
    values  <- values[keep]
    
  } else {
    # aggregation
    agg <- switch(allow_duplicate,
                  min  = tapply(values, key, min),
                  max  = tapply(values, key, max),
                  mean = tapply(values, key, mean))
    
    key_u <- as.numeric(names(agg))
    values <- as.numeric(agg)
    
    origins <- key_u %/% 1e9
    dests   <- key_u %% 1e9
  }
  
  # ---- ID handling ----
  if (is.null(ids)) {
    ids <- sort(unique(c(origins, dests)))
  } else {
    ids <- sort(as.integer(ids))
  }
  
  id_map <- integer(max(ids))
  id_map[ids] <- seq_along(ids)
  
  i <- id_map[origins]
  j <- id_map[dests]
  
  ok2 <- i > 0 & j > 0
  i <- i[ok2]
  j <- j[ok2]
  values <- values[ok2]
  
  # ---- build matrix ----
  n <- length(ids)
  mat <- matrix(fill_value, n, n,
                dimnames = list(as.character(ids), as.character(ids)))
  
  mat[cbind(i, j)] <- values
  
  if (symmetric) {
    mat[cbind(j, i)] <- values
  }
  
  return(mat)
}

calcola_spatWeighNetMat <- function(pnts, ID = NULL, nb_neighbors, distMetric = "duration",
                                    osrmServer = "https://routing.openstreetmap.de/",
                                    osrmProfile = "car") {
  if (is.null(ID)) {
    pnts$ID <- as.character(1:nrow(pnts))
    ID <- "ID"
  }
  
  # Check CRS of points and transform if necessary
  if (is.na(st_crs(pnts))) {
    stop("Input points 'pnts' have no CRS defined.")
  }
  if (st_crs(pnts)$epsg != 4326) {
    warning("Transforming 'pnts' CRS to EPSG:4326 for routing compatibility.")
    pnts <- st_transform(pnts, 4326)
  }
  
  # Convert nb list to adjacency matrix if needed
  if (inherits(nb_neighbors, "nb")) {
    orig_nb_neighbors <- nb_neighbors
    nnb <- length(orig_nb_neighbors)
    nb_matrix <- matrix(0, nrow = nnb, ncol = nnb)
    for (i in seq_len(nnb)) {
      nb_list <- orig_nb_neighbors[[i]]
      if (length(nb_list) > 0) {
        nb_matrix[i, nb_list] <- 1
      }
    }
    nb_neighbors <- nb_matrix
  }
  
  stopifnot(is.matrix(nb_neighbors),
            nrow(nb_neighbors) == ncol(nb_neighbors),
            all(nb_neighbors %in% c(0, 1)))
  
  n <- nrow(pnts)
  spatWeighNetMat <- matrix(0, n, n)
  rownames(spatWeighNetMat) <- colnames(spatWeighNetMat) <- pnts[[ID]]
  
  # Initialize travelDist SF object with CRS to avoid CRS mismatch on rbind
  travelDist <- st_sf(id_from = character(),
                      id_to = character(),
                      geometry = st_sfc(crs = 4326))
  
  estrai_valore <- function(route) {
    if (distMetric == "duration") return(route$duration)
    else return(route$distance)
  }
  
  for (i in seq_len(n)) {
    for (j in seq_len(n)) {
      if (i != j && nb_neighbors[i, j] == 1) {
        route <- tryCatch({
          osrmRoute(as.numeric(st_coordinates(pnts)[i, ]), as.numeric(st_coordinates(pnts)[j, ]),
                    osrm.server = osrmServer,
                    osrm.profile = osrmProfile)
        }, error = function(e) {
          warning(sprintf("Routing failed from %s to %s: %s", pnts[[ID]][i], pnts[[ID]][j], e$message))
          return(NULL)
        })
        
        if (!is.null(route)) {
          spatWeighNetMat[i, j] <- estrai_valore(route)
          # Create sf object for this route with proper CRS
          new_route_sf <- st_sf(id_from = pnts[[ID]][i],
                                id_to = pnts[[ID]][j],
                                geometry = route$geometry,
                                crs = 4326)
          travelDist <- rbind(travelDist, new_route_sf)
        }
      }
    }
  }
  
  list(spatWeighNetMat = spatWeighNetMat, travelDist = travelDist)
}


lw_to_graph <- function(lw, sf_obj = NULL, save_path = NULL, directed = FALSE) {
  # Load required packages safely
  if (!requireNamespace("sf", quietly = TRUE) ||
      !requireNamespace("igraph", quietly = TRUE) ||
      !requireNamespace("dplyr", quietly = TRUE) ||
      !requireNamespace("purrr", quietly = TRUE)) {
    stop("Please install 'sf', 'igraph', 'dplyr', and 'purrr' packages.")
  }
  
  # --- Validation ---
  if (!inherits(lw, "listw")) stop("`lw` must be a spatial weights (listw) object.")
  
  # Extract neighbors and weights
  nb <- lw$neighbours
  weights <- lw$weights
  
  # --- Build edge list ---
  edges <- purrr::map2_dfr(seq_along(nb), nb, ~ {
    if (length(.y) == 0) return(NULL)
    tibble::tibble(
      from = .x,
      to = .y,
      weight = weights[[.x]][seq_along(.y)]
    )
  })
  
  # --- Node data ---
  nodes <- tibble::tibble(id = seq_along(nb))
  
  # --- Extract coordinates from sf object ---
  if (!is.null(sf_obj)) {
    if (!inherits(sf_obj, "sf")) stop("`sf_obj` must be an sf object with geometry column.")
    
    geom <- sf::st_geometry(sf_obj)
    geom_type <- unique(sf::st_geometry_type(geom))
    
    # Use centroids if geometry is polygonal
    if (any(grepl("POLYGON", geom_type, ignore.case = TRUE))) {
      coords <- sf::st_coordinates(sf::st_centroid(geom))
    } else {
      coords <- sf::st_coordinates(geom)
    }
    
    # Add coordinates to node table
    nodes <- nodes |>
      dplyr::mutate(x = coords[,1], y = coords[,2])
    
    # Also attach sf attributes if desired
    attrs <- sf_obj |> sf::st_drop_geometry()
    nodes <- dplyr::bind_cols(nodes, attrs)
  }
  
  # --- Build igraph object ---
  g <- igraph::graph_from_data_frame(d = edges, vertices = nodes, directed = directed)
  
  # --- Optional save to file ---
  if (!is.null(save_path)) {
    igraph::write_graph(g, save_path, format = "graphml")
    message("Graph saved to: ", normalizePath(save_path))
  }
  
  return(g)
}

#### It returns a tidy table of geometry type frequencies from a sf collection  
geom_freq_table <- function(sf_obj) {
  # Validate input
  if (!inherits(sf_obj, "sf")) stop("Input must be an sf object")
  
  # Get types per feature
  geom_types <- st_geometry_type(sf_obj, by_geometry = TRUE)
  
  # Tidy table
  tbl <- geom_types %>%
    droplevels() %>%
    table() %>%
    as.data.frame() %>%
    setNames(c("Geometry_Type", "Count")) %>%
    mutate(
      Percentage = round(Count / sum(Count) * 100, 1),
      Total_Features = sum(Count)
    ) %>%
    arrange(desc(Count))
  
  return(tbl)
}

##### It splits a Spatial Feature Collection into a named list of sf objects, one per unique geometry type
split_by_geom <- function(sf_obj) {
  # Validate input
  if (!inherits(sf_obj, "sf")) stop("Input must be an sf object")
  
  # Get geometry types per feature
  geom_types <- st_geometry_type(sf_obj, by_geometry = TRUE)
  unique_geoms <- unique(droplevels(geom_types))
  
  # Split into list by type
  split_list <- lapply(as.character(unique_geoms), function(g_type) {
    idx <- which(geom_types == g_type)
    sf_subset <- sf_obj[idx, , drop = FALSE]  # Preserve sf class
    attr(sf_subset, "geom_type") <- g_type
    attr(sf_subset, "n_features") <- length(idx)
    return(sf_subset)
  })
  
  # Name list and return
  names(split_list) <- paste0("geom_", as.character(unique_geoms))
  return(split_list)
}

######## CUMULATIVE ISOCHRONES AVOID MANY CALLS TO osrmIsochrone function
cumulative_isochrones <- function(sf_obj) {
  
  sf_obj <- st_make_valid(sf_obj) %>%
    arrange(ID, rangeIso)
  
  result <- sf_obj
  
  split_ids <- split(sf_obj, sf_obj$ID)
  
  out_all <- list()
  
  for (id in names(split_ids)) {
    
    df <- split_ids[[id]]
    
    geoms <- st_geometry(df)
    cum_geoms <- vector("list", length(geoms))
    
    cum_geoms[[1]] <- geoms[[1]]
    
    for (k in 2:length(geoms)) {
      cum_geoms[[k]] <- st_union(cum_geoms[[k-1]], geoms[[k]])
    }
    
    df$geometry <- st_sfc(cum_geoms, crs = st_crs(df))
    
    out_all[[id]] <- df
  }
  
  result <- do.call(rbind, out_all)
  return(result)
}

#### CLIP ISOCHRONES
clip_isochrones <- function(iso_sf, region_sf) {
  library(sf)
  
  # 1. CRS alignment
  if (st_crs(iso_sf) != st_crs(region_sf)) {
    region_sf <- st_transform(region_sf, st_crs(iso_sf))
  }
  
  # 2. Fix geometries (important)
  iso_sf    <- st_make_valid(iso_sf)
  region_sf <- st_make_valid(region_sf)
  
  # 3. Intersection (NO union)
  iso_sf <- iso_sf[st_intersects(iso_sf, region_sf, sparse = FALSE), ]
  clipped <- st_intersection(iso_sf, region_sf)
  
  # 4. Remove empty geometries
  clipped <- clipped[!st_is_empty(clipped), ]
  
  # 5. Keep polygonal geometries only
  clipped <- st_collection_extract(clipped, "POLYGON")
  
  return(clipped)
}

plot_iso_catchment_tmap <- function(icat_row, supermarkets, malta, outdir) {
  
  # Ensure geometry columns
  st_geometry(icat_row) <- "geom"
  st_geometry(supermarkets) <- "geom"
  st_geometry(malta) <- "geom"
  
  # Create output directory
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  
  # Safe filename helper
  safe_name <- function(x) str_replace_all(as.character(x), "[^A-Za-z0-9_-]", "_")
  
  # Extract info
  target_id  <- icat_row$ID[1]
  map_engine <- icat_row$Map_engine[1]
  range_iso  <- icat_row$rangeIso[1]
  
  map_engine_safe <- safe_name(map_engine)
  range_iso_safe  <- safe_name(range_iso)
  sup_id_safe     <- safe_name(target_id)
  
  # Split points
  target_superm <- supermarkets %>% filter(ID_superm == target_id)
  competitors   <- supermarkets %>% filter(ID_superm != target_id)
  
  if (nrow(target_superm) == 0) stop("Target supermarket not found")
  
  # ---- Build map ----
  p <- 
    tm_shape(malta) +
    tm_polygons(fill = "white", col = "grey70", lwd = 0.6) +
    
    tm_shape(icat_row) +
    tm_polygons(fill = "red", fill_alpha = 0.35, col = "black", lwd = 1) +
    
    tm_shape(competitors) +
    tm_dots(fill = "grey75", size = 0.4, alpha = 0.8) +
    
    tm_shape(target_superm) +
    tm_symbols(shape = 8, col = "red", size = 1.2)
  
  # ---- Legends ----
  p <- p +
    tm_add_legend(
      type = "symbols",
      col = c("red", "grey75"),
      labels = c("Target Supermarket", "Competitors"),
      shape = c(8, 16),
      size = c(1.2, 1),
      title = "Legend"
    ) +
    tm_add_legend(
      type = "polygons",
      fill = "red",
      col = "black",
      labels = paste0("Isochrone ≤ ", range_iso, " min"),
      title = ""
    )
  
  # ---- Title & subtitle ----
  subtitle_text <- paste(
    "Map engine:", map_engine,
    "| rangeIso:", range_iso,
    "| supermarket ID:", target_id
  )
  
  p <- p +
    tm_title("Catchment Area - Isochrone Method") +
    tm_layout(
      title = subtitle_text,
      title.size = 0.8,
      legend.outside = TRUE,
      legend.outside.position = "left",
      legend.outside.x = -0.03,
      legend.outside.y = 0.01,
      legend.outside.size = 0.22,
      outer.margins = c(0.05, 0.01, 0.05, 0.17),
      frame = FALSE
    )
  
  # ---- Save ----
  outfile <- file.path(
    outdir,
    paste0("PlotCatIso_", map_engine_safe, "_", range_iso_safe, "_sup_", sup_id_safe, ".png")
  )
  
  tmap_save(p, outfile, width = 9, height = 7, dpi = 300)
  
  return(p)
}

###INVERSE DISTANCE DECAY
inverseDistance <- function(W, alpha=1){
  # W=ED; alpha=1
  ID <- as.matrix(1/(W^alpha))
  ID[!is.finite(ID)] <- 0
  return(ID)
}

### 
build_pairwise_intersections <- function(df) {
  
  # spatial index (VERY fast, C++)
  nb <- st_intersects(df)
  
  results <- list()
  k <- 1
  
  for (i in seq_along(nb)) {
    
    # only likely neighbors
    js <- nb[[i]]
    js <- js[js > i]
    
    if (length(js) == 0) next
    
    for (j in js) {
      
      # intersection ONLY here (minimal calls)
      inter <- st_intersection(df[i, ], df[j, ])
      
      if (nrow(inter) == 0) next
      
      inter <- st_make_valid(inter)
      
      # handle geometry collections safely
      if (any(grepl("GEOMETRYCOLLECTION", st_geometry_type(inter)))) {
        inter <- st_collection_extract(inter, "POLYGON")
      }
      
      if (nrow(inter) == 0) next
      
      inter <- st_union(inter)
      
      if (st_is_empty(inter)) next
      
      results[[k]] <- st_sf(
        ID_i = df$ID[i],
        ID_j = df$ID[j],
        Map_engine = df$Map_engine[1],
        rangeIso = df$rangeIso[1],
        geometry = st_sfc(inter, crs = st_crs(df))
      )
      
      k <- k + 1
    }
  }
  
  if (length(results) == 0) return(st_sf())
  
  do.call(rbind, results)
}

clean_polygonal_sf <- function(x) {
  
  # 1. Fix invalid geometries
  x <- sf::st_make_valid(x)
  
  # 2. Drop empty geometries
  x <- x[!sf::st_is_empty(x), ]
  
  # 3. Keep only polygonal + geometry collections
  x <- x[sf::st_geometry_type(x) %in% c("POLYGON", "MULTIPOLYGON", "GEOMETRYCOLLECTION"), ]
  
  # 4. Extract polygons ONLY where needed (vectorized over whole object)
  x <- sf::st_collection_extract(x, "POLYGON")
  
  # 5. Ensure consistent output type
  x <- sf::st_cast(x, "MULTIPOLYGON", warn = FALSE)
  
  return(x)
}

########################################################################################################################################################################
########################################################################################################################################################################
add_int_weights <- function(intcat, icat) {
  stopifnot(inherits(intcat, "sf"), inherits(icat, "sf"))
  
  # i side (origin)
  icat_i <- icat %>%
    st_drop_geometry() %>%
    select(ID, rangeIso, Map_engine, IsoPopulation) %>%
    rename(ID_i = ID, IsoPopulation_i = IsoPopulation)
  
  # j side (destination)
  icat_j <- icat %>%
    st_drop_geometry() %>%
    select(ID, rangeIso, Map_engine, IsoPopulation) %>%
    rename(ID_j = ID, IsoPopulation_j = IsoPopulation)
  
  # Attach BOTH populations to EACH intersection (edge)
  out <- intcat %>%
    left_join(icat_i, by = c("ID_i", "rangeIso", "Map_engine")) %>%
    left_join(icat_j, by = c("ID_j", "rangeIso", "Map_engine")) %>%
    mutate(
      w_ij = ifelse(IsoPopulation_i > 0, IntPopulation / IsoPopulation_i, 0),
      w_ji = ifelse(IsoPopulation_j > 0, IntPopulation / IsoPopulation_j, 0)
    )
  
  return(out)
}

build_binary_weights_matrix <- function(df, all_ids) {
  
  all_ids <- as.character(all_ids)
  
  W <- matrix(0,
              nrow = length(all_ids),
              ncol = length(all_ids),
              dimnames = list(all_ids, all_ids))
  
  i_idx <- match(as.character(df$ID_i), all_ids)
  j_idx <- match(as.character(df$ID_j), all_ids)
  
  W[cbind(i_idx, j_idx)] <- 1
  W[cbind(j_idx, i_idx)] <- 1   # symmetry
  
  diag(W) <- 0
  
  return(W)
}

build_weights_matrix <- function(df, all_ids) {
  
  all_ids <- as.character(all_ids)
  
  # aggregate duplicates (important!)
  df_agg <- df %>%
    dplyr::group_by(ID_i, ID_j) %>%
    dplyr::summarise(
      w_ij = sum(w_ij, na.rm = TRUE),
      w_ji = sum(w_ji, na.rm = TRUE),
      .groups = "drop"
    )
  
  W <- matrix(0,
              nrow = length(all_ids),
              ncol = length(all_ids),
              dimnames = list(all_ids, all_ids))
  
  i_idx <- match(as.character(df_agg$ID_i), all_ids)
  j_idx <- match(as.character(df_agg$ID_j), all_ids)
  
  W[cbind(i_idx, j_idx)] <- df_agg$w_ij
  W[cbind(j_idx, i_idx)] <- df_agg$w_ji
  
  diag(W) <- 0
  
  return(W)
}

########################################################################################################################################################################
########################################################################################################################################################################

# library(stringr)
# library(ggplot2)
# plot_iso_catchment <- function(icat_row, supermarkets, malta, outdir) {
#   st_geometry(icat_row) <- "geom"
#   st_geometry(supermarkets) <- "geom"
#   st_geometry(malta) <- "geom"
#   
#   dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
#   
#   safe_name <- function(x) str_replace_all(as.character(x), "[^A-Za-z0-9_-]", "_")
#   
#   target_id <- icat_row$ID[1]
#   map_engine <- safe_name(icat_row$Map_engine[1])
#   range_iso  <- safe_name(icat_row$rangeIso[1])
#   sup_id     <- safe_name(target_id)
#   
#   target_superm <- supermarkets %>% filter(ID_superm == target_id)
#   competitors   <- supermarkets %>% filter(ID_superm != target_id)
#   
#   if (nrow(target_superm) == 0) stop("Target supermarket not found in supermarkets")
#   
#   p <- ggplot() +
#     geom_sf(data = malta, fill = "white", color = "grey70", linewidth = 0.4) +
#     geom_sf(data = competitors, color = "grey75", size = 2, alpha = 0.8) +
#     geom_sf(data = icat_row, fill = "red", color = "black", alpha = 0.35, linewidth = 0.6) +
#     geom_sf(data = target_superm, color = "red", size = 3.2, shape = 8, stroke = 1.1) +
#     coord_sf(expand = FALSE) +
#     labs(
#       title = "Catchment Area - Isochrone Method",
#       subtitle = paste("Map engine:", icat_row$Map_engine[1],
#                        "| rangeIso:", icat_row$rangeIso[1],
#                        "| supermarket ID:", target_id),
#       x = NULL, y = NULL
#     ) +
#     theme_minimal(base_size = 12) +
#     theme(
#       panel.grid = element_blank(),
#       legend.position = "none",
#       plot.title = element_text(size = 18),
#       plot.subtitle = element_text(size = 11)
#     )
#   
#   outfile <- file.path(
#     outdir,
#     paste0("PlotCatIso_", map_engine, "_", range_iso, "_sup_", sup_id, ".png")
#   )
#   
#   ggsave(outfile, p, width = 9, height = 7, dpi = 300)
#   return(p)
# }

# compute_intersections <- function(df_group) {
#   # function(df_group) {
#     
#     df_group <- st_make_valid(df_group)
#     
#     # spatial index once (IMPORTANT)
#     idx <- st_intersects(df_group)
#     
#     results <- list()
#     k <- 1
#     
#     for (i in seq_along(idx)) {
#       
#       # only candidates that actually intersect
#       j_vals <- idx[[i]]
#       
#       for (j in j_vals[j_vals > i]) {
#         
#         inter <- st_intersection(df_group[i, ], df_group[j, ])
#         
#         if (nrow(inter) == 0) next
#         
#         inter <- st_collection_extract(inter, "POLYGON")
#         if (nrow(inter) == 0) next
#         
#         inter <- st_union(inter)
#         
#         results[[k]] <- st_sf(
#           i = df_group$ID[i],
#           j = df_group$ID[j],
#           rangeIso = df_group$rangeIso[1],
#           geometry = st_sfc(inter, crs = st_crs(df_group))
#         )
#         
#         k <- k + 1
#       }
#     }
#     
#     if (length(results) == 0) return(st_sf())
#     
#     do.call(rbind, results)
# }
# 
# build_intersection_sf <- function(df) {
# 
#   df <- st_make_valid(df)
#   
#   idx <- st_intersects(df)
#   
#   out <- list()
#   k <- 1
#   
#   for (i in seq_along(idx)) {
#     for (j in idx[[i]]) {
#       
#       if (j <= i) next
#       
#       inter <- st_intersection(df[i, ], df[j, ])
#       
#       if (nrow(inter) == 0) next
#       
#       inter <- st_collection_extract(inter, "POLYGON")
#       if (nrow(inter) == 0) next
#       
#       out[[k]] <- st_sf(
#         ID_i = df$ID[i],
#         ID_j = df$ID[j],
#         Map_engine = df$Map_engine[1],
#         rangeIso = df$rangeIso[1],
#         geometry = st_union(inter)
#       )
#       
#       k <- k + 1
#     }
#   }
#   
#   if (length(out) == 0) return(st_sf())
#   
#   do.call(rbind, out)
# }

#   library(sf)
#   
#   df_group <- st_make_valid(df_group)
#   
#   n <- nrow(df_group)
#   results <- list()
#   count <- 1
#   
#   for (i in 1:(n - 1)) {
#     for (j in (i + 1):n) {
#       
#       # Pre-filter (BIG speedup)
#       if (!st_intersects(df_group[i, ], df_group[j, ], sparse = FALSE)) next
#       
#       inter <- tryCatch(
#         st_intersection(df_group[i, ], df_group[j, ]),
#         error = function(e) NULL
#       )
#       
#       if (is.null(inter) || nrow(inter) == 0 || all(st_is_empty(inter))) next
#       
#       inter <- st_make_valid(inter)
#       
#       # Extract polygons only
#       inter <- st_collection_extract(inter, "POLYGON")
#       if (nrow(inter) == 0 || all(st_is_empty(inter))) next
#       
#       # Dissolve
#       inter <- st_union(inter)
#       
#       if (st_is_empty(inter)) next
#       
#       results[[count]] <- st_sf(
#         i = df_group$ID[i],
#         j = df_group$ID[j],
#         rangeIso = unique(df_group$rangeIso),
#         geometry = st_sfc(inter, crs = st_crs(df_group))
#       )
#       
#       count <- count + 1
#     }
#   }
#   
#   if (length(results) == 0) {
#     return(st_sf(
#       i = integer(0),
#       j = integer(0),
#       rangeIso = integer(0),
#       geometry = st_sfc(crs = st_crs(df_group))
#     ))
#   }
#   
#   do.call(rbind, results)
# }

# extractPop <- function(catchment, pop_raster) {
#   
#   # CRS check (IMPORTANT)
#   if (st_crs(catchment) != st_crs(pop_raster)) {
#     catchment <- st_transform(catchment, st_crs(pop_raster))
#   }
#   
#   overlaps <- catchment %>%
#     group_by(rangeIso) %>%
#     group_map(~ compute_intersections(.x), .keep = TRUE) %>%
#     bind_rows()
#   
#   if (nrow(overlaps) > 0) {
#     overlaps <- overlaps[!st_is_empty(overlaps), ]
#     
#     overlaps$IntPopulation <- exact_extract(pop_raster, overlaps, "sum")
#   }
#   
#   return(overlaps)
# }

