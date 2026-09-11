#!/usr/bin/env Rscript
# maligner_dp multi-rank strip figure.
# One image per molecule: the query, then ranks 1-10, each drawn as a smoothed and an
# unsmoothed reference strip. Every reference pair carries the number of fragments that
# smoothing removed.
#
# Usage: Rscript maligner_multirank_strip.R <parent_dir> <smoothed_ref.maps> <unsmoothed_ref.maps>

library(ggplot2)

# ============== Command-line arguments ==============
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  cat("Usage: Rscript maligner_multirank_strip.R <dir> <smoothed_ref.maps> <unsmoothed_ref.maps>\n")
  quit(status = 1)
}
PARENT_DIR <- args[1]
SMOOTHED_REF_FILE <- args[2]
UNSMOOTHED_REF_FILE <- args[3]

cat("Parent folder:", PARENT_DIR, "\n")
cat("Smoothed Reference:", SMOOTHED_REF_FILE, "\n")
cat("Unsmoothed Reference:", UNSMOOTHED_REF_FILE, "\n")
cat("========================================\n")

# ============== Parameters ==============
MAX_RANKS       <- 10
SMOOTHING_THRESHOLD <- 100
LABEL_LW        <- 0.6
LABEL_COLOR     <- "green"

# ============== Reading .maps files (same as maligner_heatstrip.R) ==============

read_query_labels <- function(map_file) {
  line <- readLines(map_file, n = 1)
  parts <- strsplit(line, "\t")[[1]]
  frag_lengths <- as.numeric(parts[4:length(parts)])
  return(c(0, cumsum(frag_lengths)))
}

read_query_fragments <- function(map_file) {
  line <- readLines(map_file, n = 1)
  parts <- strsplit(line, "\t")[[1]]
  return(as.numeric(parts[4:length(parts)]))
}

read_ref_labels <- function(map_file, map_id) {
  lines <- readLines(map_file)
  for (line in lines) {
    parts <- strsplit(line, "\t")[[1]]
    if (parts[1] == map_id) {
      frag_lengths <- as.numeric(parts[4:length(parts)])
      return(c(0, cumsum(frag_lengths)))
    }
  }
  warning(paste("Map ID not found:", map_id))
  return(NULL)
}

read_ref_fragments <- function(map_file, map_id) {
  lines <- readLines(map_file)
  for (line in lines) {
    parts <- strsplit(line, "\t")[[1]]
    if (parts[1] == map_id) {
      return(as.numeric(parts[4:length(parts)]))
    }
  }
  return(NULL)
}


# ============== Fragment matching (same as maligner_heatstrip.R) ==============

match_fragments <- function(smoothed_frags_all, unsmoothed_frags_all,
                             ref_start_idx, ref_end_idx) {
  # Deterministic smoothed -> unsmoothed mapping.
  # merge: exact lookup of the cumsum boundary (R match())
  # delete: keep the fragment indices above the threshold
  n_smooth <- ref_end_idx - ref_start_idx
  if (n_smooth == 0) return(NULL)

  is_merge <- abs(sum(smoothed_frags_all) - sum(unsmoothed_frags_all)) /
              sum(unsmoothed_frags_all) < 0.01

  if (is_merge) {
    sm_labels  <- c(0, cumsum(smoothed_frags_all))
    unsm_labels <- c(0, cumsum(unsmoothed_frags_all))

    mapping <- integer(n_smooth)
    for (s in 1:n_smooth) {
      target_bp <- sm_labels[ref_start_idx + s + 1]
      idx <- match(target_bp, unsm_labels)
      if (!is.na(idx)) {
        mapping[s] <- idx - 1
      }
    }

    start_bp <- sm_labels[ref_start_idx + 1]
    start_idx <- match(start_bp, unsm_labels)
    if (is.na(start_idx)) start_idx <- which.min(abs(unsm_labels - start_bp))
    u_start_idx <- start_idx - 1

  } else {
    above_thr <- which(unsmoothed_frags_all > SMOOTHING_THRESHOLD)

    mapping <- integer(n_smooth)
    for (s in 1:n_smooth) {
      sm_idx <- ref_start_idx + s
      if (sm_idx <= length(above_thr)) {
        mapping[s] <- above_thr[sm_idx]
      }
    }

    if (mapping[1] > 0) {
      u_start_idx <- mapping[1] - 1
    } else {
      return(NULL)
    }
  }

  matched <- sum(mapping > 0)
  if (matched == 0) return(NULL)

  list(
    u_start_idx = u_start_idx,
    mapping = mapping,
    u_end_idx = max(mapping[mapping > 0])
  )
}


# ============== Unsmoothed track alignment (same as maligner_heatstrip.R) ==============

align_unsmoothed_track <- function(smoothed_frags, smoothed_labels_scaled,
                                    unsmoothed_labels, match_result) {
  if (is.null(match_result)) return(NULL)

  u_start_idx <- match_result$u_start_idx
  mapping     <- match_result$mapping
  u_end_idx   <- match_result$u_end_idx

  unsmoothed_range <- unsmoothed_labels[(u_start_idx + 1):(u_end_idx + 1)]
  unsmoothed_range <- unsmoothed_range - unsmoothed_range[1]
  n_unsmooth_range <- length(unsmoothed_range)

  anchors_smooth   <- c(smoothed_labels_scaled[1])
  anchors_unsmooth <- c(0)

  for (s_idx in seq_along(mapping)) {
    if (mapping[s_idx] > 0) {
      s_label_pos <- smoothed_labels_scaled[s_idx + 1]
      u_frag_idx  <- mapping[s_idx]
      u_label_idx_in_range <- u_frag_idx - u_start_idx + 1
      if (u_label_idx_in_range >= 1 && u_label_idx_in_range <= n_unsmooth_range) {
        anchors_smooth   <- c(anchors_smooth,   s_label_pos)
        anchors_unsmooth <- c(anchors_unsmooth, unsmoothed_range[u_label_idx_in_range])
      }
    }
  }

  result <- numeric(n_unsmooth_range)
  for (i in 1:n_unsmooth_range) {
    u_pos <- unsmoothed_range[i]
    seg_idx <- max(which(anchors_unsmooth <= u_pos))
    if (seg_idx < length(anchors_unsmooth)) {
      u_start_seg <- anchors_unsmooth[seg_idx]
      u_end_seg   <- anchors_unsmooth[seg_idx + 1]
      s_start_seg <- anchors_smooth[seg_idx]
      s_end_seg   <- anchors_smooth[seg_idx + 1]
      if (u_end_seg > u_start_seg) {
        ratio <- (u_pos - u_start_seg) / (u_end_seg - u_start_seg)
        result[i] <- s_start_seg + ratio * (s_end_seg - s_start_seg)
      } else {
        result[i] <- s_start_seg
      }
    } else {
      if (length(anchors_smooth) >= 2 && length(anchors_unsmooth) >= 2) {
        last_s_len <- anchors_smooth[length(anchors_smooth)] - anchors_smooth[length(anchors_smooth) - 1]
        last_u_len <- anchors_unsmooth[length(anchors_unsmooth)] - anchors_unsmooth[length(anchors_unsmooth) - 1]
        if (last_u_len > 0) {
          scale <- last_s_len / last_u_len
          result[i] <- anchors_smooth[length(anchors_smooth)] +
                       (u_pos - anchors_unsmooth[length(anchors_unsmooth)]) * scale
        } else {
          result[i] <- anchors_smooth[length(anchors_smooth)]
        }
      } else {
        result[i] <- u_pos
      }
    }
  }
  return(result)
}


# ============== Per-rank data extraction ==============

prepare_rank_data <- function(df, row_index, query_folder,
                               smoothed_ref_file, unsmoothed_ref_file,
                               query_labels_full) {
  ref_id   <- basename(as.character(df[row_index, "ref_map"]))
  is_reverse <- (df[row_index, "is_forward"] == "R")
  direction  <- if (is_reverse) "-" else "+"
  score    <- round(as.numeric(df[row_index, "total_rescaled_score"]), 2)
  m_sc     <- round(as.numeric(df[row_index, "m_score"]), 2)
  fp_pct   <- round(as.numeric(df[row_index, "query_miss_rate"]) * 100, 1)
  fn_pct   <- round(as.numeric(df[row_index, "ref_miss_rate"]) * 100, 1)

  query_start_idx <- as.numeric(df[row_index, "query_start"])
  query_end_idx   <- as.numeric(df[row_index, "query_end"])
  ref_start_idx   <- as.numeric(df[row_index, "ref_start"])
  ref_end_idx     <- as.numeric(df[row_index, "ref_end"])

  q_offset    <- query_labels_full[query_start_idx + 1]
  n_sm_frags  <- ref_end_idx - ref_start_idx

  # Smoothed ref labels
  sm_ref_labels <- read_ref_labels(smoothed_ref_file, ref_id)
  if (is.null(sm_ref_labels)) return(NULL)

  smoothed_aligned <- sm_ref_labels[(ref_start_idx + 1):(ref_end_idx + 1)]
  smoothed_aligned <- smoothed_aligned - smoothed_aligned[1]
  if (max(smoothed_aligned) == 0) return(NULL)

  # Chunks, needed for matched_query_positions and for the inner anchor
  chunk_string <- as.character(df[row_index, "chunk_string"])
  chunks <- strsplit(chunk_string, ";")[[1]]
  chunks <- chunks[chunks != ""]
  if (length(chunks) == 0) return(NULL)

  query_bp <- 0
  matched_query_positions <- c(0)
  first_chunk_r_end  <- NULL
  last_chunk_r_start <- NULL

  for (j in seq_along(chunks)) {
    parts <- as.numeric(strsplit(chunks[j], ",")[[1]])
    q_len <- parts[3]; r_start <- parts[4]; r_end <- parts[5]
    if (j == 1) first_chunk_r_end <- r_end
    last_chunk_r_start <- r_start
    query_bp <- query_bp + q_len
    matched_query_positions <- c(matched_query_positions, query_bp)
  }

  # Inner anchor scaling (same logic as maligner_heatstrip.R)
  n_mq <- length(matched_query_positions)
  r_li <- as.integer(first_chunk_r_end  - ref_start_idx + 1)
  r_ri <- as.integer(last_chunk_r_start - ref_start_idx + 1)
  use_inner <- !is.na(r_li) && !is.na(r_ri) && (n_mq >= 4) &&
               (r_li >= 1) && (r_li <= length(smoothed_aligned)) &&
               (r_ri >= 1) && (r_ri <= length(smoothed_aligned)) &&
               (r_ri > r_li) && (smoothed_aligned[r_ri] > smoothed_aligned[r_li])

  if (use_inner) {
    q_left  <- matched_query_positions[2]
    q_right <- matched_query_positions[n_mq - 1]
    r_left  <- smoothed_aligned[r_li]
    r_right <- smoothed_aligned[r_ri]
    inner_scale  <- (q_right - q_left) / (r_right - r_left)
    inner_offset <- q_left - r_left * inner_scale
    smoothed_scaled <- smoothed_aligned * inner_scale + inner_offset
  } else {
    smoothed_scaled <- smoothed_aligned * (query_bp / max(smoothed_aligned))
  }

  # Unsmoothed ref
  unsm_ref_labels <- read_ref_labels(unsmoothed_ref_file, ref_id)
  unsm_frags_all  <- read_ref_fragments(unsmoothed_ref_file, ref_id)
  sm_frags_all    <- read_ref_fragments(smoothed_ref_file, ref_id)

  unsmoothed_scaled <- c()
  n_unsm_frags <- NA

  if (!is.null(unsm_frags_all) && !is.null(sm_frags_all)) {
    sm_aligned_frags <- sm_frags_all[(ref_start_idx + 1):ref_end_idx]
    match_result <- match_fragments(sm_frags_all, unsm_frags_all,
                                     ref_start_idx, ref_end_idx)

    if (!is.null(match_result)) {
      n_unsm_frags <- match_result$u_end_idx - match_result$u_start_idx

      unsmoothed_scaled <- align_unsmoothed_track(
        sm_aligned_frags, smoothed_scaled, unsm_ref_labels, match_result)

      if (is.null(unsmoothed_scaled)) {
        u_range <- unsm_ref_labels[(match_result$u_start_idx + 1):(match_result$u_end_idx + 1)]
        u_range <- u_range - u_range[1]
        if (max(u_range) > 0) {
          sm_lo <- min(smoothed_scaled); sm_hi <- max(smoothed_scaled)
          unsmoothed_scaled <- u_range / max(u_range) * (sm_hi - sm_lo) + sm_lo
        } else {
          unsmoothed_scaled <- u_range
        }
      }
    }
  }

  # Offset into full-molecule coordinates
  smoothed_scaled <- smoothed_scaled + q_offset
  if (length(unsmoothed_scaled) > 0)
    unsmoothed_scaled <- unsmoothed_scaled + q_offset

  # A reverse alignment mirrors the reference strip left to right
  if (is_reverse) {
    flip_sum <- 2 * q_offset + query_bp
    smoothed_scaled <- flip_sum - smoothed_scaled
    if (length(unsmoothed_scaled) > 0)
      unsmoothed_scaled <- flip_sum - unsmoothed_scaled
  }

  delta <- if (!is.na(n_unsm_frags)) n_unsm_frags - n_sm_frags else 0L

  list(
    smoothed_labels  = smoothed_scaled,
    unsmoothed_labels = unsmoothed_scaled,
    n_sm_frags  = n_sm_frags,
    n_unsm_frags = if (is.na(n_unsm_frags)) n_sm_frags else n_unsm_frags,
    delta     = delta,
    ref_id    = ref_id,
    direction = direction,
    score     = score,
    m_score   = m_sc,
    fp_pct    = fp_pct,
    fn_pct    = fn_pct
  )
}


# ============== Multi-rank plot ==============

create_multirank_plot <- function(query_name, query_labels_full, ranks, output_path) {
  n_ranks  <- length(ranks)
  total_bp <- max(query_labels_full)
  n_frags  <- length(query_labels_full) - 1

  # Track layout
  TRACK_H     <- 0.035
  GAP_WITHIN  <- 0.005   # smoothed to raw, inside one rank
  GAP_BETWEEN <- 0.018   # between rank groups
  GAP_QUERY   <- 0.025   # data ↔ rank1

  total_y <- TRACK_H + GAP_QUERY +
             n_ranks * (2 * TRACK_H + GAP_WITHIN) +
             max(0, n_ranks - 1) * GAP_BETWEEN

  y_cur <- total_y
  query_y <- c(y_cur - TRACK_H, y_cur)
  y_cur <- y_cur - TRACK_H - GAP_QUERY

  track_ys <- list()
  for (i in 1:n_ranks) {
    sm_y  <- c(y_cur - TRACK_H, y_cur)
    y_cur <- y_cur - TRACK_H - GAP_WITHIN
    raw_y <- c(y_cur - TRACK_H, y_cur)
    if (i < n_ranks) y_cur <- y_cur - TRACK_H - GAP_BETWEEN
    track_ys[[i]] <- list(sm = sm_y, raw = raw_y)
  }

  # Pad, so a reference label at a negative offset still fits
  pad_bp <- 0
  for (rd in ranks) {
    if (length(rd$smoothed_labels) > 0)
      pad_bp <- max(pad_bp, abs(min(0, min(rd$smoothed_labels))))
    if (length(rd$unsmoothed_labels) > 0)
      pad_bp <- max(pad_bp, abs(min(0, min(rd$unsmoothed_labels))))
  }

  label_x <- -pad_bp - total_bp * 0.01
  info_x  <- total_bp + pad_bp + total_bp * 0.01

  # --- Data track ---
  q_df <- data.frame(x = query_labels_full, ymin = query_y[1], ymax = query_y[2])
  p <- ggplot() +
    annotate("rect", xmin = -pad_bp, xmax = total_bp + pad_bp,
             ymin = query_y[1], ymax = query_y[2],
             fill = "black", color = "black", linewidth = 0.3) +
    geom_segment(data = q_df, aes(x = x, xend = x, y = ymin, yend = ymax),
                 color = LABEL_COLOR, linewidth = LABEL_LW, lineend = "butt") +
    annotate("text", x = label_x, y = mean(query_y),
             label = "Data", hjust = 1, size = 3.5, fontface = "bold",
             color = "black") +
    annotate("text", x = info_x, y = mean(query_y),
             label = paste0(n_frags, " frags, ", round(total_bp / 1000, 1), " kb"),
             hjust = 0, size = 2.5, color = "black")

  # --- Rank tracks ---
  for (i in 1:n_ranks) {
    rd    <- ranks[[i]]
    sm_y  <- track_ys[[i]]$sm
    raw_y <- track_ys[[i]]$raw

    # Smoothed ref
    if (length(rd$smoothed_labels) > 0) {
      sm_df <- data.frame(x = rd$smoothed_labels,
                          ymin = sm_y[1], ymax = sm_y[2])
      sm_info <- paste0(rd$ref_id, " (", rd$direction, ") | ",
                        rd$score, " | m:", rd$m_score,
                        " | FN:", rd$fn_pct, "% FP:", rd$fp_pct, "%")
      p <- p +
        annotate("rect", xmin = -pad_bp, xmax = total_bp + pad_bp,
                 ymin = sm_y[1], ymax = sm_y[2],
                 fill = "black", color = "black", linewidth = 0.3) +
        geom_segment(data = sm_df,
                     aes(x = x, xend = x, y = ymin, yend = ymax),
                     color = LABEL_COLOR, linewidth = LABEL_LW,
                     lineend = "butt") +
        annotate("text", x = label_x, y = mean(sm_y),
                 label = paste0("R", i, " Ref (sm)"),
                 hjust = 1, size = 3, color = "blue") +
        annotate("text", x = info_x, y = mean(sm_y),
                 label = sm_info, hjust = 0, size = 2.5, color = "blue")
    }

    # Unsmoothed ref
    has_unsm <- length(rd$unsmoothed_labels) > 2
    if (has_unsm) {
      unsm_df <- data.frame(x = rd$unsmoothed_labels,
                            ymin = raw_y[1], ymax = raw_y[2])
      delta_str <- if (rd$delta > 0) paste0("+", rd$delta, " by smoothing")
                   else if (rd$delta == 0) "same"
                   else as.character(rd$delta)
      raw_info <- paste0(rd$n_unsm_frags, " frags [", delta_str, "]")
      p <- p +
        annotate("rect", xmin = -pad_bp, xmax = total_bp + pad_bp,
                 ymin = raw_y[1], ymax = raw_y[2],
                 fill = "black", color = "black", linewidth = 0.3) +
        geom_segment(data = unsm_df,
                     aes(x = x, xend = x, y = ymin, yend = ymax),
                     color = LABEL_COLOR, linewidth = LABEL_LW,
                     lineend = "butt") +
        annotate("text", x = label_x, y = mean(raw_y),
                 label = paste0("R", i, " Ref (raw)"),
                 hjust = 1, size = 3, color = "purple") +
        annotate("text", x = info_x, y = mean(raw_y),
                 label = raw_info, hjust = 0, size = 2.5, color = "purple")
    } else {
      p <- p +
        annotate("rect", xmin = -pad_bp, xmax = total_bp + pad_bp,
                 ymin = raw_y[1], ymax = raw_y[2],
                 fill = "#222222", color = "#444444", linewidth = 0.3) +
        annotate("text", x = label_x, y = mean(raw_y),
                 label = paste0("R", i, " Ref (raw)"),
                 hjust = 1, size = 3, color = "purple")
    }
  }

  # --- Axes / theme (same as maligner_heatstrip.R) ---
  p <- p +
    scale_x_continuous(
      breaks = seq(0, total_bp, by = round(total_bp / 5, -3)),
      labels = function(x) paste0(round(x / 1000), "k")
    ) +
    coord_cartesian(
      xlim = c(min(-total_bp * 0.05, -pad_bp - total_bp * 0.02),
               max(total_bp * 1.01, total_bp + pad_bp + total_bp * 0.01)),
      ylim = c(-0.01, total_y + 0.01),
      clip = "off"
    ) +
    labs(
      title = paste0(query_name, " — Multi-rank comparison (", n_ranks, " ranks)"),
      x = "Position (bp)", y = NULL
    ) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title    = element_text(face = "bold", size = 14),
      panel.grid.major.y = element_blank(),
      panel.grid.minor   = element_blank(),
      axis.text.y   = element_blank(),
      panel.grid.major.x = element_blank(),
      axis.ticks.x  = element_line(color = "grey50", linewidth = 0.3),
      axis.ticks.length.x = unit(1.5, "mm"),
      plot.background  = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      plot.margin = margin(10, 120, 2, 40)
    )

  plot_h <- 2 + n_ranks * 0.8
  ggsave(output_path, plot = p, width = 16, height = plot_h, dpi = 1000)
}


# ============== Main ==============

direct_mal <- list.files(PARENT_DIR, pattern = "_mal\\.txt$",
                         full.names = TRUE, recursive = FALSE)
if (length(direct_mal) > 0) {
  mal_files <- direct_mal
  cat("Single-folder mode\n")
} else {
  subdirs <- list.dirs(PARENT_DIR, recursive = FALSE)
  mal_files <- unlist(lapply(subdirs, function(d) {
    list.files(d, pattern = "_mal\\.txt$", full.names = TRUE, recursive = FALSE)
  }))
}

if (length(mal_files) == 0) {
  cat("Error: no _mal.txt file found.\n")
  quit(status = 1)
}

cat("Files found:", length(mal_files), "\n\n")

for (mal_file in mal_files) {
  cat("Processing:", mal_file, "\n")
  data_folder <- dirname(mal_file)
  folder_name <- basename(data_folder)

  tryCatch({
    df <- read.delim(mal_file, header = TRUE, stringsAsFactors = FALSE)
    if (nrow(df) > 0 && df[1, 1] == "query_map") df <- df[-1, ]

    df$total_rescaled_score <- as.numeric(df$total_rescaled_score)
    df <- df[order(df$total_rescaled_score, decreasing = FALSE), ]
    top_n <- min(MAX_RANKS, nrow(df))
    df_top <- df[1:top_n, ]
    rownames(df_top) <- 1:nrow(df_top)

    # Locate the query .maps file
    query_map_name <- basename(as.character(df_top[1, "query_map"]))
    query_map_path <- file.path(data_folder, paste0(query_map_name, ".maps"))
    if (!file.exists(query_map_path)) {
      maps_files <- list.files(data_folder, pattern = "\\.maps$", full.names = TRUE)
      if (length(maps_files) > 0) {
        query_map_path <- maps_files[1]
      } else {
        cat("  no .maps file, skipping\n\n"); next
      }
    }
    query_labels_full <- read_query_labels(query_map_path)

    cat("  rows:", nrow(df), "/ taking top", top_n, "\n")

    # One entry per rank
    ranks <- list()
    for (i in 1:top_n) {
      tryCatch({
        rd <- prepare_rank_data(df_top, i, data_folder,
                                SMOOTHED_REF_FILE, UNSMOOTHED_REF_FILE,
                                query_labels_full)
        if (!is.null(rd)) {
          ranks[[length(ranks) + 1]] <- rd
          cat(sprintf("  R%d: %s (%s) score=%.1f sm=%d unsm=%d delta=%+d\n",
                      i, rd$ref_id, rd$direction, rd$score,
                      rd$n_sm_frags, rd$n_unsm_frags, rd$delta))
        } else {
          cat("  R", i, ": bad data, skipping\n", sep = "")
        }
      }, error = function(e) {
        cat("  R", i, " error: ", e$message, "\n", sep = "")
      })
    }

    if (length(ranks) > 0) {
      output_path <- file.path(data_folder,
                               paste0(folder_name, "_multirank_strip.png"))
      create_multirank_plot(folder_name, query_labels_full, ranks, output_path)
      cat("  saved:", basename(output_path), "\n")
    }

  }, error = function(e) {
    cat("  error:", e$message, "\n")
  })
  cat("\n")
}

cat("========================================\n")
cat("Batch complete\n")
