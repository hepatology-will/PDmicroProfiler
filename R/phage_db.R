#' Query Phage Database
#'
#' Functions for querying and visualizing phage database information.

# Internal function: get votu_db data
.get_votu_db <- function() {
  if (exists("votu_db")) {
    return(votu_db)
  }
  stop("votu_db data not found. Please ensure the package is properly installed.")
}

#' Query phages by host taxonomy and phage features
#'
#' Retrieve phage records that infect a specific host at a given taxonomic level,
#' with optional filtering by phage characteristics.
#'
#' @param host_level Host taxonomic level. One of: "Host.Phylum", "Host.Class", "Host.Order",
#'   "Host.Family", "Host.Genus", "Host.Species". If NULL (default), host filtering is skipped.
#' @param host_value The name of the host taxon (e.g., "Firmicutes", "p__Firmicutes"). Required if host_level is provided.
#' @param ignore_case Logical. Whether to ignore case in host matching (default: TRUE).
#' @param remove_prefix Logical. Whether to remove common taxonomic prefixes (e.g., "p__", "c__")
#'   from the host value before matching (default: TRUE).
#' @param family Character vector of phage families (e.g., "Siphoviridae"). Default NULL.
#' @param lifestyle Character vector of lifestyles (e.g., "temperate", "virulent"). Default NULL.
#' @param specialization Character vector of host specialization (e.g., "specialist", "generalist"). Default NULL.
#' @param novelty Character vector of novelty status (e.g., "known", "novel"). Default NULL.
#' @param ... Additional arguments (currently unused).
#'
#' @return A data frame with matching phage records, or NULL if none found.
#' @export
#'
#' @examples
#' \dontrun{
#' result1 <- query_by_host(host_level = "Host.Phylum", host_value = "Firmicutes")
#' result2 <- query_by_host(
#'   host_level = "Host.Phylum",
#'   host_value = "Actinobacteria",
#'   lifestyle = "temperate"
#' )
#' }
query_by_host <- function(host_level = NULL,
                          host_value = NULL,
                          ignore_case = TRUE,
                          remove_prefix = TRUE,
                          family = NULL,
                          lifestyle = NULL,
                          specialization = NULL,
                          novelty = NULL,
                          ...) {

  data <- .get_votu_db()

  if (!is.null(host_level)) {
    valid_levels <- c("Host.Phylum", "Host.Class", "Host.Order",
                      "Host.Family", "Host.Genus", "Host.Species")
    if (!host_level %in% valid_levels) {
      stop("Invalid host classification level. Valid options: ",
           paste(valid_levels, collapse = ", "))
    }
    if (is.null(host_value)) {
      stop("host_value must be provided when host_level is specified.")
    }

    if (remove_prefix) {
      host_value_clean <- gsub("^[a-z]__", "", host_value)
    } else {
      host_value_clean <- host_value
    }

    col_vals <- unique(data[[host_level]])

    if (ignore_case) {
      pattern <- tolower(host_value_clean)
      matches_idx <- grep(pattern, tolower(col_vals))
    } else {
      pattern <- host_value_clean
      matches_idx <- grep(pattern, col_vals)
    }

    matched_vals <- col_vals[matches_idx]

    if (length(matched_vals) == 0) {
      message("No matching host taxa found.")
      return(NULL)
    }

    result <- data[data[[host_level]] %in% matched_vals, ]

    if (nrow(result) == 0) {
      message("No matching records found.")
      return(NULL)
    }

    message(sprintf("Found %d matching host taxa: %s",
                    length(matched_vals),
                    paste(matched_vals, collapse = ", ")))
    message(sprintf("Returning %d phage records after host filter.", nrow(result)))

  } else {
    result <- data
  }

  if (!is.null(family)) result <- result[result$Family %in% family, ]
  if (!is.null(lifestyle)) result <- result[result$lifestyle %in% lifestyle, ]
  if (!is.null(specialization)) result <- result[result$Specialization %in% specialization, ]
  if (!is.null(novelty)) result <- result[result$novelty_status %in% novelty, ]

  if (nrow(result) == 0) {
    message("No records match the combined criteria.")
    return(NULL)
  }

  return(result)
}

#' Query phages by their own features
#'
#' Retrieve phage records based on phage characteristics such as family,
#' lifestyle, specialization, and novelty status.
#'
#' @param family Character vector of phage families (e.g., "Siphoviridae").
#' @param lifestyle Character vector of lifestyles (e.g., "temperate", "virulent").
#' @param specialization Character vector of host specialization (e.g., "specialist", "generalist").
#' @param novelty Character vector of novelty status (e.g., "known", "novel").
#'
#' @return A data frame with matching phage records, or NULL if none found.
#' @export
#'
#' @examples
#' \dontrun{
#' result <- query_by_phage_features(family = "Siphoviridae", lifestyle = "temperate")
#' }
query_by_phage_features <- function(family = NULL,
                                    lifestyle = NULL,
                                    specialization = NULL,
                                    novelty = NULL) {
  data <- .get_votu_db()
  result <- data
  if (!is.null(family)) result <- result[result$Family %in% family, ]
  if (!is.null(lifestyle)) result <- result[result$lifestyle %in% lifestyle, ]
  if (!is.null(specialization)) result <- result[result$Specialization %in% specialization, ]
  if (!is.null(novelty)) result <- result[result$novelty_status %in% novelty, ]
  if (nrow(result) == 0) { message("No matching records found"); return(NULL) }
  return(result)
}

#' Query a specific phage by its ID
#'
#' Retrieve a phage record using its unique identifier (qseqid).
#'
#' @param id The phage ID (qseqid) to search for.
#'
#' @return A data frame with the matching phage record, or NULL if not found.
#' @export
#'
#' @examples
#' \dontrun{
#' result <- query_by_id("vOTU_12345")
#' }
query_by_id <- function(id) {
  data <- .get_votu_db()
  result <- data[data$qseqid == id, ]
  if (nrow(result) == 0) { message("No matching phage ID found"); return(NULL) }
  return(result)
}

#' Dashboard for a specific phage family
#'
#' Generate a multi-plot dashboard summarizing a phage family, including
#' genome length distribution, host distribution, lifestyle, specialization, and novelty.
#'
#' @param family_name Name of the phage family (e.g., "Siphoviridae").
#' @param host_level Host taxonomic level to display in the host distribution plot.
#'   Default is "Host.Phylum".
#' @param top_n Number of top host categories to show (others are grouped as "Others").
#'   Default is 10.
#'
#' @return A grid of ggplot2 plots arranged by gridExtra.
#' @export
#'
#' @examples
#' \dontrun{
#' phage_family_dashboard("Siphoviridae", host_level = "Host.Genus", top_n = 8)
#' }
phage_family_dashboard <- function(family_name, host_level = "Host.Phylum", top_n = 10) {
  if (!requireNamespace("gridExtra", quietly = TRUE)) stop("Please install 'gridExtra'")
  
  data <- .get_votu_db()
  
  unique_data <- data %>%
    dplyr::distinct(qseqid, .keep_all = TRUE) %>%
    dplyr::filter(Family == family_name)
  
  family_data <- data %>%
    dplyr::filter(Family == family_name)
  
  if (nrow(family_data) == 0 || nrow(unique_data) == 0) {
    message("No records found for this phage family")
    return(NULL)
  }
  
  morandi_colors <- c(
    "#e7a6a6", "#F0A388", "#FDEBE3", "#E39077",
    "#D04E70", "#C8859F", "#D89090", "#70D0C3", "#d9d9d9",
    "#847151", "#5B7F75", "#FC6E53", "#73755e", "#F5DEB1",
    "#c5dfbf", "#aec6d7", "#f1c8db", "#C091B5",
    "#C785D8", "#D1B0C5", "#865348", "#E5B965",
    "#AB8A63", "#ecab62", "#eedf72", "#B893B4"
  )
  
  # 指定颜色：只固定你要求修改的几个子图
  specialization_colors <- c(
    "Specialist" = "#3e427d",
    "Generalist" = "#c2bdd2"
  )
  
  lifestyle_colors <- c(
    "Virulent" = "#e6d7ac",
    "Uncertain Virulent" = "#d5be88",
    "Uncertain Temperate" = "#c4a461",
    "Temperate" = "#d09f3b"
  )
  
  novelty_status_colors <- c(
    "novel" = "#f7e094",
    "partial" = "#a8c1cf",
    "identical" = "#6e93a7"
  )
  
  genome_length_fill <- "#e7c5af"
  
  full_host_data <- family_data %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(host_level))) %>%
    dplyr::summarise(Count = dplyr::n(), .groups = "drop") %>%
    dplyr::arrange(dplyr::desc(Count))
  
  stats_text <- paste(
    paste("Family:", family_name),
    paste("Total vOTUs:", nrow(unique_data)),
    paste("Average Length:", format(round(mean(unique_data$Length)), big.mark = ","), "bp"),
    paste(
      "Top Host:",
      ifelse(
        nrow(full_host_data) > 0,
        paste(full_host_data[[1]][1], "(", full_host_data$Count[1], ")"),
        "None"
      )
    ),
    paste("Total Hosts:", nrow(full_host_data)),
    paste(
      "Dominant Lifestyle:",
      ifelse(
        nrow(unique_data) > 0,
        names(sort(table(unique_data$lifestyle), decreasing = TRUE)[1]),
        "None"
      )
    ),
    paste(
      "Common Specialization:",
      ifelse(
        nrow(unique_data) > 0,
        names(sort(table(unique_data$Specialization), decreasing = TRUE)[1]),
        "None"
      )
    ),
    sep = " | "
  )
  
  stats_plot <- ggplot2::ggplot() +
    ggplot2::theme_void() +
    ggplot2::geom_text(
      ggplot2::aes(x = 0.5, y = 0.5, label = stats_text),
      size = 4,
      hjust = 0.5,
      vjust = 0.5
    ) +
    ggplot2::theme(
      panel.border = ggplot2::element_rect(color = "black", fill = NA, linewidth = 1),
      plot.margin = ggplot2::margin(10, 10, 10, 10)
    )
  
  p_length <- ggplot2::ggplot(unique_data, ggplot2::aes(x = Length)) +
    ggplot2::geom_histogram(
      fill = genome_length_fill,
      binwidth = 2000,
      alpha = 0.85
    ) +
    ggplot2::geom_density(
      ggplot2::aes(y = ggplot2::after_stat(count) * 2000),
      color = morandi_colors[9],
      linewidth = 1
    ) +
    ggplot2::labs(
      title = "Genome Length Distribution",
      x = "Genome Length (bp)",
      y = "Number of vOTUs"
    ) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
      panel.grid.minor = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      panel.border = ggplot2::element_rect(color = "black", fill = NA, linewidth = 1)
    ) +
    ggplot2::scale_x_continuous(labels = scales::comma)
  
  host_data <- family_data %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(host_level))) %>%
    dplyr::summarise(Count = dplyr::n(), .groups = "drop") %>%
    dplyr::arrange(dplyr::desc(Count))
  
  if (nrow(host_data) > top_n) {
    others <- data.frame(
      Category = "Others",
      Count = sum(host_data$Count[(top_n + 1):nrow(host_data)])
    )
    colnames(others)[1] <- host_level
    
    host_data_top <- host_data[1:top_n, ]
    host_data <- dplyr::bind_rows(host_data_top, others)
  }
  
  host_data[[1]] <- factor(host_data[[1]], levels = host_data[[1]])
  
  n_hosts <- nrow(host_data)
  host_colors <- grDevices::colorRampPalette(morandi_colors)(n_hosts)
  
  if ("Others" %in% host_data[[1]]) {
    host_colors[host_data[[1]] == "Others"] <- "#888888"
  }
  
  p_host <- ggplot2::ggplot(
    host_data,
    ggplot2::aes(
      x = stats::reorder(.data[[host_level]], Count),
      y = Count,
      fill = .data[[host_level]]
    )
  ) +
    ggplot2::geom_bar(stat = "identity", alpha = 0.85) +
    ggplot2::coord_flip() +
    ggplot2::labs(
      title = paste(
        "Top",
        min(top_n, nrow(host_data)),
        "Host Distribution:",
        gsub("Host\\.", "", host_level)
      ),
      x = "",
      y = "Number of Associations"
    ) +
    ggplot2::scale_fill_manual(values = host_colors) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
      legend.position = "none",
      panel.grid.major.y = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(color = "black", fill = NA, linewidth = 1)
    )
  
  spec_data <- unique_data %>%
    dplyr::group_by(Specialization) %>%
    dplyr::summarise(Count = dplyr::n(), .groups = "drop") %>%
    dplyr::arrange(dplyr::desc(Count))
  
  p_specialization <- ggplot2::ggplot(
    spec_data,
    ggplot2::aes(
      x = stats::reorder(Specialization, Count),
      y = Count,
      fill = Specialization
    )
  ) +
    ggplot2::geom_bar(stat = "identity", alpha = 0.85) +
    ggplot2::geom_text(
      ggplot2::aes(label = Count),
      hjust = -0.3,
      size = 4
    ) +
    ggplot2::coord_flip() +
    ggplot2::labs(
      title = "Host Specialization",
      x = "",
      y = ""
    ) +
    ggplot2::scale_fill_manual(values = specialization_colors) +
    ggplot2::expand_limits(y = max(spec_data$Count) * 1.2) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
      legend.position = "none",
      panel.grid.major.y = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(color = "black", fill = NA, linewidth = 1)
    )
  
  life_data <- unique_data %>%
    dplyr::group_by(lifestyle) %>%
    dplyr::summarise(Count = dplyr::n(), .groups = "drop") %>%
    dplyr::arrange(dplyr::desc(Count))
  
  p_lifestyle <- ggplot2::ggplot(
    life_data,
    ggplot2::aes(
      x = stats::reorder(lifestyle, Count),
      y = Count,
      fill = lifestyle
    )
  ) +
    ggplot2::geom_bar(stat = "identity", alpha = 0.85) +
    ggplot2::geom_text(
      ggplot2::aes(label = Count),
      vjust = -0.3,
      size = 4
    ) +
    ggplot2::labs(
      title = "Lifestyle Distribution",
      x = "",
      y = "Number of vOTUs"
    ) +
    ggplot2::scale_fill_manual(values = lifestyle_colors) +
    ggplot2::expand_limits(y = max(life_data$Count) * 1.2) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
      legend.position = "none",
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      panel.border = ggplot2::element_rect(color = "black", fill = NA, linewidth = 1)
    )
  
  novelty_data <- unique_data %>%
    dplyr::group_by(novelty_status) %>%
    dplyr::summarise(Count = dplyr::n(), .groups = "drop") %>%
    dplyr::arrange(dplyr::desc(Count))
  
  p_novelty <- ggplot2::ggplot(
    novelty_data,
    ggplot2::aes(
      x = stats::reorder(novelty_status, Count),
      y = Count,
      fill = novelty_status
    )
  ) +
    ggplot2::geom_bar(stat = "identity", alpha = 0.85) +
    ggplot2::geom_text(
      ggplot2::aes(label = Count),
      vjust = -0.3,
      size = 4
    ) +
    ggplot2::labs(
      title = "Novelty Status",
      x = "",
      y = ""
    ) +
    ggplot2::scale_fill_manual(values = novelty_status_colors) +
    ggplot2::expand_limits(y = max(novelty_data$Count) * 1.2) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
      legend.position = "none",
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      panel.border = ggplot2::element_rect(color = "black", fill = NA, linewidth = 1)
    )
  
  gridExtra::grid.arrange(
    stats_plot,
    gridExtra::arrangeGrob(p_length, p_host, ncol = 2, widths = c(1.5, 1)),
    gridExtra::arrangeGrob(p_specialization, p_lifestyle, p_novelty, ncol = 3),
    nrow = 3,
    heights = c(0.3, 1.2, 0.8),
    top = grid::textGrob(
      paste("Phage Family Dashboard:", family_name),
      gp = grid::gpar(fontsize = 18, fontface = "bold", fontfamily = "sans")
    )
  )
}

#' Visualize summary of phage query results
#'
#' Create bar plots or pie charts summarizing the distribution of
#' phage families, lifestyles, specializations, and novelty status from a query result.
#'
#' @param data A data frame returned by query functions (e.g., query_by_host).
#' @param plot_type Type of plot: "bar" (default) or "pie".
#' @param color_palette Optional vector of colors. If NULL, a default Morandi color palette is used.
#' @param n_top Number of top categories to show for family distribution (others grouped as "Other").
#' @return A combined ggplot object (patchwork) with four subplots.
#' @export
#'
#' @examples
#' \dontrun{
#' result <- query_by_host(host_level = "Host.Phylum", host_value = "Firmicutes")
#' plot_phage_summary(result, plot_type = "pie")
#' }
plot_phage_summary <- function(data, plot_type = "bar", color_palette = NULL, n_top = 10) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Please install ggplot2.")
  if (!requireNamespace("patchwork", quietly = TRUE)) stop("Please install patchwork.")
  if (!requireNamespace("dplyr", quietly = TRUE)) stop("Please install dplyr.")
  if (!requireNamespace("rlang", quietly = TRUE)) stop("Please install rlang.")
  if (!requireNamespace("scales", quietly = TRUE)) stop("Please install scales.")
  
  morandi_colors <- c(
    "#a5c2d3", "#a6435d", "#824421", "#c63179", "#E39077",
    "#6b67a0", "#246d9f", "#a55629", "#339cd1", "#399139",
    "#9d7225", "#5B7F75", "#c82621", "#7E9161", "#edb66d",
    "#3ABD9C", "#82bed7", "#E8D0DC", "#a7cb84", "#3e427d",
    "#5f3887", "#e19eaa", "#92545A", "#209170", "#E5B965",
    "#AB8A63", "#588979", "#c55d1a", "#d4a220", "#B893B4"
  )
  
  if (is.null(color_palette)) {
    color_palette <- morandi_colors
  }
  
  lifestyle_colors <- c(
    "Virulent" = "#e6d7ac",
    "Uncertain Virulent" = "#d5be88",
    "Uncertain Temperate" = "#c4a461",
    "Temperate" = "#d09f3b"
  )
  
  specialization_colors <- c(
    "Specialist" = "#3e427d",
    "Generalist" = "#c2bdd2"
  )
  
  novelty_colors <- c(
    "novel" = "#f7e094",
    "partial" = "#a8c1cf",
    "identical" = "#6e93a7"
  )
  
  required_cols <- c("Family", "lifestyle", "Specialization", "novelty_status")
  missing <- setdiff(required_cols, colnames(data))
  if (length(missing) > 0) {
    stop("Input data missing required columns: ", paste(missing, collapse = ", "))
  }
  
  if (!"qseqid" %in% colnames(data)) {
    stop("Input data missing required column: qseqid")
  }
  
  plot_type <- match.arg(plot_type, choices = c("bar", "pie"))
  
  data_unique <- data %>%
    dplyr::distinct(qseqid, .keep_all = TRUE)
  
  .get_colors_for_col <- function(counts, col, default_palette) {
    categories <- as.character(counts[[col]])
    
    if (col == "lifestyle") {
      color_map <- lifestyle_colors
    } else if (col == "Specialization") {
      color_map <- specialization_colors
    } else if (col == "novelty_status") {
      color_map <- novelty_colors
    } else {
      return(stats::setNames(rep(default_palette, length.out = length(categories)), categories))
    }
    
    colors <- color_map[categories]
    
    if (any(is.na(colors))) {
      missing_categories <- categories[is.na(colors)]
      fallback_colors <- rep(default_palette, length.out = length(missing_categories))
      colors[is.na(colors)] <- fallback_colors
    }
    
    stats::setNames(colors, categories)
  }
  
  .plot_dist <- function(df, col, title, plot_type, default_palette) {
    counts <- df %>%
      dplyr::count(!!rlang::sym(col)) %>%
      dplyr::arrange(dplyr::desc(n))
    
    if (col == "Family" && nrow(counts) > n_top) {
      top <- counts[1:n_top, ]
      other_count <- sum(counts$n[(n_top + 1):nrow(counts)])
      other <- data.frame("Other", other_count, stringsAsFactors = FALSE)
      colnames(other) <- c(col, "n")
      counts <- dplyr::bind_rows(top, other)
    }
    
    counts[[col]] <- as.character(counts[[col]])
    counts[[col]] <- factor(counts[[col]], levels = counts[[col]])
    
    plot_colors <- .get_colors_for_col(counts, col, default_palette)
    
    if (plot_type == "bar") {
      p <- ggplot2::ggplot(
        counts,
        ggplot2::aes(
          x = stats::reorder(.data[[col]], n),
          y = n,
          fill = .data[[col]]
        )
      ) +
        ggplot2::geom_bar(stat = "identity") +
        ggplot2::coord_flip() +
        ggplot2::labs(x = "", y = "Count", title = title) +
        ggplot2::theme_bw() +
        ggplot2::theme(
          legend.position = "none",
          plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
          panel.grid.major.y = ggplot2::element_blank()
        ) +
        ggplot2::scale_fill_manual(values = plot_colors)
    } else {
      counts <- counts %>%
        dplyr::mutate(
          prop = n / sum(n),
          label = paste0(
            .data[[col]], "\n",
            n, " (", scales::percent(prop, accuracy = 0.1), ")"
          )
        )
      
      p <- ggplot2::ggplot(
        counts,
        ggplot2::aes(x = "", y = n, fill = .data[[col]])
      ) +
        ggplot2::geom_bar(stat = "identity", width = 1) +
        ggplot2::coord_polar("y", start = 0) +
        ggplot2::labs(title = title, fill = col) +
        ggplot2::theme_void() +
        ggplot2::theme(
          legend.position = "right",
          plot.title = ggplot2::element_text(hjust = 0.5, face = "bold")
        ) +
        ggplot2::scale_fill_manual(values = plot_colors)
    }
    
    return(p)
  }
  
  get_colors <- function(n) {
    rep(color_palette, length.out = n)
  }
  
  p1 <- .plot_dist(
    data_unique,
    "Family",
    "Phage Family",
    plot_type,
    get_colors(length(unique(data_unique$Family)))
  )
  
  p2 <- .plot_dist(
    data_unique,
    "lifestyle",
    "Lifestyle",
    plot_type,
    get_colors(length(unique(data_unique$lifestyle)))
  )
  
  p3 <- .plot_dist(
    data_unique,
    "Specialization",
    "Specialization",
    plot_type,
    get_colors(length(unique(data_unique$Specialization)))
  )
  
  p4 <- .plot_dist(
    data_unique,
    "novelty_status",
    "Novelty",
    plot_type,
    get_colors(length(unique(data_unique$novelty_status)))
  )
  
  p <- p1 + p2 + p3 + p4 + patchwork::plot_layout(ncol = 2)
  
  return(p)
}
