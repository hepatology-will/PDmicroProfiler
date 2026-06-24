# ============================================================================
# vOTU Functional Analysis Module
# ============================================================================

# Internal function: get significant_ko_data
.get_ko_data <- function(data = NULL) {
  if (!is.null(data)) return(data)
  if (exists("significant_ko_data")) return(significant_ko_data)
  stop("significant_ko_data not found. Please ensure the package is properly installed.")
}

# Internal function: get differential vOTUs list (not exported)
.get_diff_votus <- function(data, p_threshold = 0.05, adj = TRUE,
                            fc_threshold = 1, direction = "both",
                            min_genes = NULL) {
  p_col <- if (adj) "adj.p.value" else "p.value"

  votu_summary <- data %>%
    dplyr::group_by(vOTU) %>%
    dplyr::summarise(
      Family = dplyr::first(Family),
      direction = dplyr::first(direction),
      min_p = min(.data[[p_col]], na.rm = TRUE),
      median_fc = stats::median(fold_change, na.rm = TRUE),
      median_log2fc = stats::median(log2_fold_change, na.rm = TRUE),
      n_KO = dplyr::n(),
      .groups = "drop"
    )

  votu_summary <- votu_summary %>% dplyr::filter(min_p < p_threshold)

  if (direction == "PD") {
    votu_summary <- votu_summary %>%
      dplyr::filter(direction == "Enriched_in_PD" & median_log2fc > log2(fc_threshold))
  } else if (direction == "HC") {
    votu_summary <- votu_summary %>%
      dplyr::filter(direction == "Enriched_in_HC" & median_log2fc < -log2(fc_threshold))
  } else {
    votu_summary <- votu_summary %>%
      dplyr::filter(abs(median_log2fc) > log2(fc_threshold))
  }

  if (!is.null(min_genes)) {
    votu_summary <- votu_summary %>% dplyr::filter(n_KO >= min_genes)
  }

  if (nrow(votu_summary) == 0) return(NULL)

  return(votu_summary$vOTU)
}

# ----------------------------------------------------------------------------
# Query functions
# ----------------------------------------------------------------------------

#' Query all records for a specific vOTU
#'
#' @param votu_id vOTU ID
#' @param data Optional data frame
#' @return Data frame
#' @export
#'
#' @examples
#' \dontrun{
#' get_votu_info("vOTU_123")
#' }
get_votu_info <- function(votu_id, data = NULL) {
  data <- .get_ko_data(data)
  result <- data[data$vOTU == votu_id, ]
  if (nrow(result) == 0) warning("No records found for this vOTU.")
  result
}

#' Query records containing a specific KO number
#'
#' @param ko_id KO number
#' @param data Optional data frame
#' @param exact Logical; if TRUE, exact match; if FALSE, fuzzy match (default)
#' @return Data frame
#' @export
#'
#' @examples
#' \dontrun{
#' get_ko_info("K00001")
#' }
get_ko_info <- function(ko_id, data = NULL, exact = FALSE) {
  data <- .get_ko_data(data)
  if (exact) result <- data[data$KO == ko_id, ]
  else result <- data[grep(ko_id, data$KO, fixed = TRUE), ]
  if (nrow(result) == 0) warning("No records containing this KO found.")
  result
}

#' Query records containing a specific pathway description
#'
#' @param pathway_name Pathway description
#' @param data Optional data frame
#' @param exact Logical; if TRUE, exact match; if FALSE, fuzzy match (default)
#' @return Data frame
#' @export
#'
#' @examples
#' \dontrun{
#' get_pathway_info("Glycolysis")
#' }
get_pathway_info <- function(pathway_name, data = NULL, exact = FALSE) {
  data <- .get_ko_data(data)
  if (exact) result <- data[data$pathway == pathway_name, ]
  else result <- data[grep(pathway_name, data$pathway, ignore.case = TRUE), ]
  if (nrow(result) == 0) warning("No records containing this pathway found.")
  result
}

#' Query all records for a specific phage family
#'
#' @param family_name Family name
#' @param data Optional data frame
#' @return Data frame
#' @export
#'
#' @examples
#' \dontrun{
#' get_family_info("Siphoviridae")
#' }
get_family_info <- function(family_name, data = NULL) {
  data <- .get_ko_data(data)
  result <- data[data$Family == family_name, ]
  if (nrow(result) == 0) warning("No records found for this family.")
  result
}

#' Query records for a specific enrichment direction
#'
#' @param direction "Enriched_in_PD" or "Enriched_in_HC"
#' @param data Optional data frame
#' @return Data frame
#' @export
#'
#' @examples
#' \dontrun{
#' get_direction_info("Enriched_in_PD")
#' }
get_direction_info <- function(direction = c("Enriched_in_PD", "Enriched_in_HC"), data = NULL) {
  direction <- match.arg(direction)
  data <- .get_ko_data(data)
  result <- data[data$direction == direction, ]
  if (nrow(result) == 0) warning("No records found for this direction.")
  result
}

#' Get list of significantly differential vOTUs
#'
#' @param p_threshold p-value threshold
#' @param data Optional data frame
#' @param adj Logical; use adjusted p-value (default TRUE)
#' @param return_records Logical; if TRUE, return full data frame; if FALSE, return vector of vOTU IDs
#' @return Character vector or data frame
#' @export
#'
#' @examples
#' \dontrun{
#' get_significant_votus(p_threshold = 0.01)
#' }
get_significant_votus <- function(p_threshold = 0.05, data = NULL, adj = TRUE, return_records = FALSE) {
  data <- .get_ko_data(data)
  p_col <- if (adj) "adj.p.value" else "p.value"
  sig <- data[data[[p_col]] < p_threshold, ]
  if (nrow(sig) == 0) { warning("No significant vOTUs found."); return(NULL) }
  if (return_records) sig else unique(sig$vOTU)
}

#' Summarize functional information for a single vOTU
#'
#' @param votu_id vOTU ID
#' @param data Optional data frame
#' @return List with summary statistics
#' @export
#'
#' @examples
#' \dontrun{
#' summarize_votu_function("vOTU_123")
#' }
summarize_votu_function <- function(votu_id, data = NULL) {
  data <- .get_ko_data(data)
  sub <- data[data$vOTU == votu_id, ]
  if (nrow(sub) == 0) stop("vOTU not found.")
  list(
    vOTU = votu_id,
    Family = unique(sub$Family),
    direction = unique(sub$direction),
    min_pvalue = min(sub$p.value),
    min_adj_pvalue = min(sub$adj.p.value),
    KO_count = length(unique(sub$KO)),
    KO_list = unique(sub$KO),
    pathway_count = length(unique(sub$pathway[!is.na(sub$pathway) & sub$pathway != ""])),
    pathway_list = unique(sub$pathway[!is.na(sub$pathway) & sub$pathway != ""]),
    proteins = sub$Protein
  )
}

# ----------------------------------------------------------------------------
# Visualization functions
# ----------------------------------------------------------------------------

# Morandi-style palette for bar plots
.bar_morandi_palette <- c(
  "#cae1ca", "#7fa9d1", "#ed8c58", "#e6dad2", "#bfbbdb",
  "#89c7c1", "#aacedd", "#347aae", "#f7ba79", "#ee8474",
  "#d8c6a3", "#b7c9a8", "#9fb7cf", "#d6b4a7", "#c8b8d8",
  "#8fb6a7", "#d9c9b8", "#b6a6a0", "#a9c3c7", "#c9b27c"
)

#' Plot KO or pathway distribution for a single vOTU
#'
#' @param votu_id vOTU ID
#' @param data Optional data frame
#' @param by "KO" or "pathway"
#' @param top_n Number of top categories to show
#' @param plot_type "bar" or "bubble"
#' @param title Custom title
#' @return ggplot2 object
#' @export
#'
#' @examples
#' \dontrun{
#' plot_votu_function("vOTU_123", by = "KO", top_n = 5)
#' }
plot_votu_function <- function(votu_id, data = NULL, by = c("KO", "pathway"),
                               top_n = 10, plot_type = c("bar", "bubble"),
                               title = NULL) {
  by <- match.arg(by)
  plot_type <- match.arg(plot_type)
  data <- .get_ko_data(data)
  
  sub <- data[data$vOTU == votu_id, ]
  if (nrow(sub) == 0) stop("No records found for this vOTU.")
  
  if (by == "KO") {
    col <- "KO"
    col_name <- "KO"
  } else {
    col <- "pathway"
    col_name <- "Pathway"
    sub <- sub[!is.na(sub$pathway) & sub$pathway != "", ]
    if (nrow(sub) == 0) stop("This vOTU has no pathway annotations.")
  }
  
  counts <- sub %>%
    dplyr::count(!!rlang::sym(col)) %>%
    dplyr::arrange(dplyr::desc(n))
  
  if (nrow(counts) > top_n) {
    counts <- counts[1:top_n, ]
  }
  
  counts[[col]] <- factor(counts[[col]], levels = counts[[col]])
  
  if (plot_type == "bar") {
    bar_colors <- stats::setNames(
      rep(.bar_morandi_palette, length.out = length(levels(counts[[col]]))),
      levels(counts[[col]])
    )
    
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
      ggplot2::scale_fill_manual(values = bar_colors) +
      ggplot2::labs(
        x = col_name,
        y = "Count",
        title = if (is.null(title)) {
          paste("vOTU:", votu_id, "-", by, "distribution")
        } else {
          title
        }
      ) +
      ggplot2::theme_bw() +
      ggplot2::theme(
        legend.position = "none",
        plot.title = ggplot2::element_text(hjust = 0.5, face = "bold")
      )
  } else {
    sub_agg <- sub %>%
      dplyr::group_by(!!rlang::sym(col)) %>%
      dplyr::summarise(
        count = dplyr::n(),
        avg_logp = mean(-log10(p), na.rm = TRUE),
        .groups = "drop"
      ) %>%
      dplyr::arrange(dplyr::desc(count)) %>%
      dplyr::slice_head(n = top_n)
    
    p <- ggplot2::ggplot(
      sub_agg,
      ggplot2::aes(
        x = stats::reorder(.data[[col]], count),
        y = avg_logp,
        size = count,
        color = .data[[col]]
      )
    ) +
      ggplot2::geom_point(alpha = 0.7) +
      ggplot2::scale_size_continuous(range = c(3, 10), name = "Count") +
      ggplot2::labs(
        x = col_name,
        y = "-log10(avg p-value)",
        title = if (is.null(title)) {
          paste("vOTU:", votu_id, "-", by, "bubble")
        } else {
          title
        }
      ) +
      ggplot2::coord_flip() +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        legend.position = "right",
        axis.text.y = ggplot2::element_text(size = 10),
        plot.title = ggplot2::element_text(hjust = 0.5, face = "bold")
      )
  }
  
  # Add border and axis ticks
  p <- p + ggplot2::theme(
    panel.border = ggplot2::element_rect(color = "black", fill = NA),
    axis.ticks = ggplot2::element_line()
  )
  
  p
}


#' Plot KO or pathway distribution for a phage family
#'
#' @param family_name Family name
#' @param data Optional data frame
#' @param by "KO" or "pathway"
#' @param top_n Number of top categories to show
#' @param plot_type "bar" or "bubble"
#' @param title Custom title
#' @return ggplot2 object
#' @export
#'
#' @examples
#' \dontrun{
#' plot_family_function("Siphoviridae", by = "pathway", top_n = 5)
#' }
plot_family_function <- function(family_name, data = NULL, by = c("KO", "pathway"),
                                 top_n = 10, plot_type = c("bar", "bubble"),
                                 title = NULL) {
  by <- match.arg(by)
  plot_type <- match.arg(plot_type)
  data <- .get_ko_data(data)
  
  sub <- data[data$Family == family_name, ]
  if (nrow(sub) == 0) stop("No records found for this family.")
  
  if (by == "KO") {
    col <- "KO"
    col_name <- "KO"
  } else {
    col <- "pathway"
    col_name <- "Pathway"
    sub <- sub[!is.na(sub$pathway) & sub$pathway != "", ]
    if (nrow(sub) == 0) stop("This family has no pathway annotations.")
  }
  
  counts <- sub %>%
    dplyr::count(!!rlang::sym(col)) %>%
    dplyr::arrange(dplyr::desc(n))
  
  if (nrow(counts) > top_n) {
    counts <- counts[1:top_n, ]
  }
  
  counts[[col]] <- factor(counts[[col]], levels = counts[[col]])
  
  if (plot_type == "bar") {
    bar_colors <- stats::setNames(
      rep(.bar_morandi_palette, length.out = length(levels(counts[[col]]))),
      levels(counts[[col]])
    )
    
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
      ggplot2::scale_fill_manual(values = bar_colors) +
      ggplot2::labs(
        x = col_name,
        y = "Count",
        title = if (is.null(title)) {
          paste("Family:", family_name, "-", by, "distribution")
        } else {
          title
        }
      ) +
      ggplot2::theme_bw() +
      ggplot2::theme(
        legend.position = "none",
        plot.title = ggplot2::element_text(hjust = 0.5, face = "bold")
      )
  } else {
    sub_agg <- sub %>%
      dplyr::group_by(!!rlang::sym(col)) %>%
      dplyr::summarise(
        count = dplyr::n(),
        avg_logp = mean(-log10(p), na.rm = TRUE),
        .groups = "drop"
      ) %>%
      dplyr::arrange(dplyr::desc(count)) %>%
      dplyr::slice_head(n = top_n)
    
    p <- ggplot2::ggplot(
      sub_agg,
      ggplot2::aes(
        x = stats::reorder(.data[[col]], count),
        y = avg_logp,
        size = count,
        color = .data[[col]]
      )
    ) +
      ggplot2::geom_point(alpha = 0.7) +
      ggplot2::scale_size_continuous(range = c(3, 10), name = "Count") +
      ggplot2::labs(
        x = col_name,
        y = "-log10(avg p-value)",
        title = if (is.null(title)) {
          paste("Family:", family_name, "-", by, "bubble")
        } else {
          title
        }
      ) +
      ggplot2::coord_flip() +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        legend.position = "right",
        axis.text.y = ggplot2::element_text(size = 10),
        plot.title = ggplot2::element_text(hjust = 0.5, face = "bold")
      )
  }
  
  # Add border and axis ticks
  p <- p + ggplot2::theme(
    panel.border = ggplot2::element_rect(color = "black", fill = NA),
    axis.ticks = ggplot2::element_line()
  )
  
  p
}

#' Plot pathway enrichment comparison between PD and HC groups for differential vOTUs
#'
#' @param data Optional data frame (significant_ko_data)
#' @param p_threshold p-value threshold for pathway enrichment (default 0.05)
#' @param adj Logical; use adjusted p-value for vOTU filtering (default TRUE)
#' @param fc_threshold Fold change threshold for vOTU filtering (default 1)
#' @param direction Differential direction: "both", "PD", "HC" (default "both")
#' @param top_n Number of top significant pathways to show (default 20)
#' @param min_genes Minimum number of KO genes per vOTU (optional)
#' @return ggplot2 bar plot object
#' @export
#'
#' @examples
#' \dontrun{
#' plot_diff_pathway(p_threshold = 0.01, top_n = 15)
#' }
plot_diff_pathway <- function(data = NULL,
                              p_threshold = 0.05,
                              adj = TRUE,
                              fc_threshold = 1,
                              direction = c("both", "PD", "HC"),
                              top_n = 20,
                              min_genes = NULL) {
  direction <- match.arg(direction)
  data <- .get_ko_data(data)

  diff_votus <- .get_diff_votus(data,
                                p_threshold = p_threshold,
                                adj = adj,
                                fc_threshold = fc_threshold,
                                direction = direction,
                                min_genes = min_genes)

  if (is.null(diff_votus) || length(diff_votus) == 0) {
    stop("No differential vOTUs found with the given criteria.")
  }

  diff_data <- data %>% dplyr::filter(vOTU %in% diff_votus)
  diff_data <- diff_data[!is.na(diff_data$pathway) & diff_data$pathway != "", ]
  if (nrow(diff_data) == 0) stop("No pathway data for differential vOTUs.")

  pathway_counts <- diff_data %>%
    dplyr::group_by(pathway, direction) %>%
    dplyr::summarise(count = dplyr::n(), .groups = "drop") %>%
    tidyr::pivot_wider(names_from = direction, values_from = count, values_fill = 0)

  if (!"Enriched_in_PD" %in% colnames(pathway_counts)) pathway_counts$Enriched_in_PD <- 0
  if (!"Enriched_in_HC" %in% colnames(pathway_counts)) pathway_counts$Enriched_in_HC <- 0

  total_PD <- sum(pathway_counts$Enriched_in_PD)
  total_HC <- sum(pathway_counts$Enriched_in_HC)

  pathway_stats <- pathway_counts %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      a = Enriched_in_PD, b = total_PD - Enriched_in_PD,
      c = Enriched_in_HC, d = total_HC - Enriched_in_HC,
      cont_mat = list(matrix(c(a, b, c, d), nrow = 2, byrow = TRUE)),
      fisher = list(stats::fisher.test(cont_mat)),
      p.value = fisher$p.value,
      odds_ratio = fisher$estimate,
      log2_odds_ratio = log2(odds_ratio)
    ) %>%
    dplyr::ungroup() %>%
    dplyr::select(pathway, Enriched_in_PD, Enriched_in_HC, p.value, odds_ratio, log2_odds_ratio)

  pathway_stats$adj.p.value <- stats::p.adjust(pathway_stats$p.value, method = "BH")

  sig_pathways <- pathway_stats %>%
    dplyr::filter(adj.p.value < p_threshold) %>%
    dplyr::arrange(dplyr::desc(abs(log2_odds_ratio))) %>%
    dplyr::slice_head(n = top_n)

  if (nrow(sig_pathways) == 0) {
    warning("No significant pathways found at adjusted p < ", p_threshold)
    return(NULL)
  }

  plot_data <- sig_pathways %>%
    tidyr::pivot_longer(cols = c(Enriched_in_PD, Enriched_in_HC),
                        names_to = "group", values_to = "gene_count") %>%
    dplyr::mutate(group = ifelse(group == "Enriched_in_PD", "PD", "HC"))

  plot_data <- plot_data %>%
    dplyr::group_by(pathway) %>%
    dplyr::mutate(total = sum(gene_count)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(pathway = stats::reorder(pathway, total))

  p <- ggplot2::ggplot(plot_data,
                       ggplot2::aes(x = pathway,
                                    y = gene_count,
                                    fill = group)) +
    ggplot2::geom_bar(stat = "identity", position = "dodge") +
    ggplot2::scale_fill_manual(values = c("PD" = "#b35761", "HC" = "#a9c7da")) +
    ggplot2::labs(x = "", y = "Gene Count",
                  title = "Pathway enrichment in differential vOTUs (PD vs HC)") +
    ggplot2::coord_flip() +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "top",
                   axis.text.y = ggplot2::element_text(size = 10),
                   plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"))

  # Add border and axis ticks
  p <- p + ggplot2::theme(
    panel.border = ggplot2::element_rect(color = "black", fill = NA),
    axis.ticks = ggplot2::element_line()
  )
  return(p)
}
