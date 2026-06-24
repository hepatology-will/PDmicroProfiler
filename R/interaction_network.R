#' Find bacteria-phage interactions
#'
#' Search for phage-host interactions based on bacterial taxonomy.
#'
#' @param bacteria_name Bacteria name (can be species, genus, family, order, class, phylum level).
#' @param groups Groups to search, options: "HC", "PD", "Shared". Default is all three.
#' @param tax_level Taxonomic level, options: "species", "genus", "family", "order", "class", "phylum".
#' @param correlation_threshold Correlation threshold, default 0.
#' @return A data frame with interaction information.
#' @export
#'
#' @examples
#' \dontrun{
#' species_phages <- find_phage_host_interactions(
#'   bacteria_name = "Achromobacter_xylosoxidans",
#'   groups = c("HC", "PD", "Shared"),
#'   tax_level = "species",
#'   correlation_threshold = 0.5
#' )
#' }
find_phage_host_interactions <- function(bacteria_name,
                                         groups = c("HC", "PD", "Shared"),
                                         tax_level = "species",
                                         correlation_threshold = 0) {

  all_results <- data.frame()

  for (group in groups) {
    # Select edges and nodes based on group
    edges <- switch(group,
                    "HC" = interaction_edges_HC,
                    "PD" = interaction_edges_PD,
                    "Shared" = interaction_edges_shared)

    nodes <- switch(group,
                    "HC" = interaction_nodes_HC,
                    "PD" = interaction_nodes_PD,
                    "Shared" = interaction_nodes_shared)

    # Match bacterial name based on taxonomic level
    if (tax_level == "species") {
      pattern <- paste0("s__.*", bacteria_name, ".*")
      matched_bacteria <- nodes$id[grepl(pattern, nodes$id, ignore.case = TRUE) & nodes$group == "bacteria"]
    } else if (tax_level == "genus") {
      pattern <- paste0(".*", bacteria_name, ".*")
      matched_genus <- bac_tax$Genus[grepl(pattern, bac_tax$Genus, ignore.case = TRUE)]
      matched_bacteria <- bac_tax$id[bac_tax$Genus %in% matched_genus]
    } else if (tax_level == "family") {
      pattern <- paste0(".*", bacteria_name, ".*")
      matched_families <- bac_tax$Family[grepl(pattern, bac_tax$Family, ignore.case = TRUE)]
      matched_bacteria <- bac_tax$id[bac_tax$Family %in% matched_families]
    } else if (tax_level == "phylum") {
      pattern <- paste0(".*", bacteria_name, ".*")
      matched_phyla <- bac_tax$Phylum[grepl(pattern, bac_tax$Phylum, ignore.case = TRUE)]
      matched_bacteria <- bac_tax$id[bac_tax$Phylum %in% matched_phyla]
    } else if (tax_level == "class") {
      pattern <- paste0(".*", bacteria_name, ".*")
      matched_classes <- bac_tax$Class[grepl(pattern, bac_tax$Class, ignore.case = TRUE)]
      matched_bacteria <- bac_tax$id[bac_tax$Class %in% matched_classes]
    } else if (tax_level == "order") {
      pattern <- paste0(".*", bacteria_name, ".*")
      matched_orders <- bac_tax$Order[grepl(pattern, bac_tax$Order, ignore.case = TRUE)]
      matched_bacteria <- bac_tax$id[bac_tax$Order %in% matched_orders]
    }

    if (length(matched_bacteria) == 0) {
      message(paste("No matching bacteria found for", bacteria_name, "at", tax_level, "level in", group, "group"))
      next
    }

    interactions <- edges[edges$Source %in% matched_bacteria, ]
    interactions <- interactions[abs(interactions$Correlation) >= correlation_threshold, ]

    if (nrow(interactions) == 0) {
      message(paste("No interactions found for", bacteria_name, "in", group, "group"))
      next
    }

    interactions$Group <- group

    interactions <- merge(interactions, bac_tax, by.x = "Source", by.y = "id", all.x = TRUE)
    interactions <- merge(interactions, votu_tax, by.x = "Target", by.y = "id", all.x = TRUE,
                          suffixes = c("_bacteria", "_virus"))

    all_results <- rbind(all_results, interactions)
  }

  if (nrow(all_results) == 0) {
    return(NULL)
  }

  return(all_results)
}

#' Visualize bacteria-phage interaction network
#'
#' Creates an interactive network plot using visNetwork.
#'
#' @param interaction_data Output from find_phage_host_interactions.
#' @param correlation_threshold Correlation threshold, default 0.3.
#' @param title Network title.
#' @return An interactive visNetwork object.
#' @export
#'
#' @examples
#' \dontrun{
#' plot_interaction_network(species_phages, correlation_threshold = 0.5)
#' }
plot_interaction_network <- function(interaction_data,
                                     correlation_threshold = 0.3,
                                     title = "Bacteria-phage interaction network") {

  if (!requireNamespace("visNetwork", quietly = TRUE)) {
    stop("Please install visNetwork package: install.packages('visNetwork')")
  }

  filtered_data <- interaction_data[abs(interaction_data$Correlation) >= correlation_threshold, ]

  if (nrow(filtered_data) == 0) {
    message(paste("No interactions with correlation above", correlation_threshold))
    return(NULL)
  }

  bacteria_nodes <- data.frame(
    id = unique(filtered_data$Source),
    label = unique(filtered_data$Source),
    group = "Bacteria",
    title = paste0(
      "Phylum: ", filtered_data$Phylum_bacteria[match(unique(filtered_data$Source), filtered_data$Source)], "<br>",
      "Class: ", filtered_data$Class_bacteria[match(unique(filtered_data$Source), filtered_data$Source)], "<br>",
      "Order: ", filtered_data$Order_bacteria[match(unique(filtered_data$Source), filtered_data$Source)], "<br>",
      "Family: ", filtered_data$Family_bacteria[match(unique(filtered_data$Source), filtered_data$Source)], "<br>",
      "Genus: ", filtered_data$Genus_bacteria[match(unique(filtered_data$Source), filtered_data$Source)], "<br>",
      "Species: ", filtered_data$Species_bacteria[match(unique(filtered_data$Source), filtered_data$Source)]
    ),
    stringsAsFactors = FALSE
  )

  phage_nodes <- data.frame(
    id = unique(filtered_data$Target),
    label = unique(filtered_data$Target),
    group = "Phage",
    title = paste0(
      "Virus ID: ", unique(filtered_data$Target), "<br>",
      "Phylum: ", filtered_data$Phylum_virus[match(unique(filtered_data$Target), filtered_data$Target)], "<br>",
      "Class: ", filtered_data$Class_virus[match(unique(filtered_data$Target), filtered_data$Target)], "<br>",
      "Order: ", filtered_data$Order_virus[match(unique(filtered_data$Target), filtered_data$Target)], "<br>",
      "Family: ", filtered_data$Family_virus[match(unique(filtered_data$Target), filtered_data$Target)]
    ),
    stringsAsFactors = FALSE
  )

  nodes <- rbind(bacteria_nodes, phage_nodes)

  edges <- data.frame(
    from = filtered_data$Source,
    to = filtered_data$Target,
    value = abs(filtered_data$Correlation),
    title = paste0("Correlation: ", round(filtered_data$Correlation, 3)),
    color = ifelse(filtered_data$Correlation > 0, "#EECBC9", "#D1E9EF"),
    stringsAsFactors = FALSE
  )

  visNetwork::visNetwork(nodes, edges, main = title) %>%
    visNetwork::visGroups(groupname = "Bacteria", color = "#D88090", shape = "dot") %>%
    visNetwork::visGroups(groupname = "Phage", color = "#87DBCB", shape = "triangle") %>%
    visNetwork::visLegend() %>%
    visNetwork::visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE) %>%
    visNetwork::visPhysics(enabled = FALSE) %>%
    visNetwork::visLayout(randomSeed = 123) %>%
    visNetwork::visInteraction(navigationButtons = TRUE, dragNodes = TRUE, dragView = TRUE, zoomView = TRUE)
}

#' Visualize bacteria-phage interaction network (static)
#'
#' Creates a static network plot using ggplot2 and igraph, suitable for saving
#' as PNG, PDF, or other raster/vector formats via \code{ggsave()}. Bacteria
#' nodes are placed at the centre; phage nodes with positive correlations are
#' arranged on the right arc and those with negative correlations on the left
#' arc. Edge width is proportional to the absolute correlation value.
#'
#' @param interaction_data Output from \code{find_phage_host_interactions()}.
#' @param correlation_threshold Numeric. Only edges with absolute correlation
#'   greater than or equal to this value are shown. Default: \code{0.3}.
#' @param top_n_phage Integer. Maximum number of phage nodes to display (ranked
#'   by absolute correlation, descending). Default: \code{30}.
#' @param title Character. Custom plot title. If \code{NULL} (default), a title
#'   is generated automatically from the data.
#' @param bacteria_color Character. Fill colour for bacteria nodes (hex or named
#'   colour). Default: \code{"#D88090"}.
#' @param phage_color Character. Colour for phage nodes (hex or named colour).
#'   Default: \code{"#87DBCB"}.
#'
#' @return A \code{ggplot2} object that can be further modified or saved with
#'   \code{ggsave()}.
#' @export
#'
#' @examples
#' \dontrun{
#' interactions <- find_phage_host_interactions(
#'   bacteria_name = "Lachnospiraceae",
#'   groups = "PD",
#'   tax_level = "family",
#'   correlation_threshold = 0.3
#' )
#' p <- plot_interaction_network_static(interactions, correlation_threshold = 0.3)
#' ggsave("network_PD.png", p, width = 10, height = 10, dpi = 300)
#' }
plot_interaction_network_static <- function(interaction_data,
                                            correlation_threshold = 0.3,
                                            top_n_phage = 30,
                                            title = NULL,
                                            bacteria_color = "#D88090",
                                            phage_color = "#87DBCB") {

  if (!requireNamespace("igraph", quietly = TRUE)) {
    stop("Package 'igraph' is required. Please install it: install.packages('igraph')")
  }

  # Filter by correlation threshold
  filtered_data <- interaction_data[abs(interaction_data$Correlation) >= correlation_threshold, ]

  if (is.null(filtered_data) || nrow(filtered_data) == 0) {
    message(paste("No interactions with absolute correlation >=", correlation_threshold))
    return(NULL)
  }

  # Keep only top_n_phage phages ranked by absolute correlation
  filtered_data <- filtered_data[order(-abs(filtered_data$Correlation)), ]
  if (nrow(filtered_data) > top_n_phage) {
    filtered_data <- filtered_data[1:top_n_phage, ]
  }

  bac_nodes   <- unique(filtered_data$Source)
  phage_nodes <- unique(filtered_data$Target)
  n_phage     <- length(phage_nodes)

  # Build igraph object (used only for the node name list)
  g <- igraph::graph_from_data_frame(
    d        = filtered_data[, c("Source", "Target", "Correlation")],
    directed = FALSE
  )
  node_names <- igraph::V(g)$name

  # ── Layout ─────────────────────────────────────────────────────────────────
  # Bacteria: centre (or small inner circle when multiple bacteria exist)
  # Phages with positive correlation: right semicircle
  # Phages with negative correlation: left semicircle
  coords <- matrix(0, nrow = length(node_names), ncol = 2,
                   dimnames = list(node_names, c("x", "y")))

  pos_phages <- unique(filtered_data$Target[filtered_data$Correlation >  0])
  neg_phages <- unique(filtered_data$Target[filtered_data$Correlation <= 0])
  r <- 2.8  # circle radius

  # Place bacteria nodes: single bacteria stays at (0,0);
  # multiple bacteria are spread on a small inner circle so labels don't overlap
  n_bac <- length(bac_nodes)
  if (n_bac == 1) {
    coords[bac_nodes[1], ] <- c(0, 0)
  } else {
    r_bac <- 0.7  # inner circle radius for bacteria
    bac_angles <- seq(0, 2 * pi, length.out = n_bac + 1)[-( n_bac + 1)]
    for (i in seq_along(bac_nodes)) {
      nm <- bac_nodes[i]
      if (nm %in% rownames(coords))
        coords[nm, ] <- c(r_bac * cos(bac_angles[i]), r_bac * sin(bac_angles[i]))
    }
  }

  if (length(pos_phages) > 0) {
    angles <- seq(pi / 2, -pi / 2, length.out = length(pos_phages) + 2)
    angles <- angles[-c(1, length(angles))]
    for (i in seq_along(pos_phages)) {
      nm <- pos_phages[i]
      if (nm %in% rownames(coords))
        coords[nm, ] <- c(r * cos(angles[i]), r * sin(angles[i]))
    }
  }

  if (length(neg_phages) > 0) {
    angles <- seq(pi / 2, 3 * pi / 2, length.out = length(neg_phages) + 2)
    angles <- angles[-c(1, length(angles))]
    for (i in seq_along(neg_phages)) {
      nm <- neg_phages[i]
      if (nm %in% rownames(coords))
        coords[nm, ] <- c(r * cos(angles[i]), r * sin(angles[i]))
    }
  }

  # ── Node data frame ─────────────────────────────────────────────────────────
  node_df <- data.frame(
    name  = node_names,
    x     = coords[node_names, 1],
    y     = coords[node_names, 2],
    type  = ifelse(node_names %in% bac_nodes, "Bacteria", "Phage"),
    label = ifelse(node_names %in% bac_nodes,
                   gsub("s__", "", node_names),
                   gsub(".*__", "", node_names)),
    stringsAsFactors = FALSE
  )

  # ── Edge data frame ─────────────────────────────────────────────────────────
  edge_df <- filtered_data[, c("Source", "Target", "Correlation")]
  colnames(edge_df) <- c("from", "to", "Correlation")
  edge_df <- merge(edge_df, node_df[, c("name", "x", "y")],
                   by.x = "from", by.y = "name")
  colnames(edge_df)[c(ncol(edge_df) - 1, ncol(edge_df))] <- c("x_from", "y_from")
  edge_df <- merge(edge_df, node_df[, c("name", "x", "y")],
                   by.x = "to", by.y = "name")
  colnames(edge_df)[c(ncol(edge_df) - 1, ncol(edge_df))] <- c("x_to", "y_to")
  edge_df$Correlation <- as.numeric(edge_df$Correlation)
  edge_df$direction   <- ifelse(edge_df$Correlation > 0, "Positive", "Negative")

  # ── Auto title ──────────────────────────────────────────────────────────────
  if (is.null(title)) {
    grp_label <- if ("Group" %in% colnames(filtered_data))
      paste(unique(filtered_data$Group), collapse = "/") else ""
    bac_label <- gsub("s__", "", bac_nodes[1])
    title <- paste0("Bacteria-Phage Interaction Network",
                    if (nchar(grp_label) > 0) paste0(" (", grp_label, ")") else "")
  }

  # ── Plot ────────────────────────────────────────────────────────────────────
  phage_df    <- subset(node_df, type == "Phage")
  bacteria_df <- subset(node_df, type == "Bacteria")

  p <- ggplot2::ggplot() +
    # Edges
    ggplot2::geom_segment(
      data = edge_df,
      ggplot2::aes(x = x_from, y = y_from, xend = x_to, yend = y_to,
                   color = direction, linewidth = abs(Correlation)),
      alpha = 0.55
    ) +
    ggplot2::scale_linewidth_continuous(range = c(0.3, 2.2),
                                        name  = "|Correlation|") +
    ggplot2::scale_color_manual(
      values = c("Positive" = "#C94F6D", "Negative" = "#4A90A4"),
      name   = "Correlation"
    ) +
    # Phage nodes
    ggplot2::geom_point(
      data  = phage_df,
      ggplot2::aes(x = x, y = y),
      shape = 17, color = phage_color, size = 3.5, alpha = 0.9
    ) +
    # Bacteria nodes
    ggplot2::geom_point(
      data  = bacteria_df,
      ggplot2::aes(x = x, y = y),
      shape = 21, fill = bacteria_color, color = "white",
      size = 9, stroke = 1.5
    ) +
    # Phage labels (small, positioned away from centre)
    ggplot2::geom_text(
      data  = phage_df,
      ggplot2::aes(x = x, y = y, label = label,
                   hjust = ifelse(x >= 0, -0.12, 1.12),
                   vjust = ifelse(y >= 0, -0.7,  1.5)),
      size = 2.4, color = "grey30"
    ) +
    # Bacteria labels (bold, black, positioned outward from centre)
    ggplot2::geom_text(
      data  = bacteria_df,
      ggplot2::aes(x = x, y = y, label = label,
                   hjust = ifelse(x > 0.1, -0.15, ifelse(x < -0.1, 1.15, 0.5)),
                   vjust = ifelse(y >= 0, -2.0, 2.8)),
      size = 3.0, fontface = "bold", color = "black"
    ) +
    ggplot2::labs(
      title    = title,
      subtitle = paste0("Correlation threshold: ", correlation_threshold,
                        "  |  Phage nodes displayed: ", n_phage),
      x = NULL, y = NULL
    ) +
    ggplot2::theme_void(base_size = 12) +
    ggplot2::theme(
      plot.title      = ggplot2::element_text(size = 14, face = "bold",
                                              hjust = 0.5, margin = ggplot2::margin(b = 4)),
      plot.subtitle   = ggplot2::element_text(size = 9, color = "grey55",
                                              hjust = 0.5, margin = ggplot2::margin(b = 8)),
      legend.position = "right",
      legend.title    = ggplot2::element_text(size = 9, face = "bold"),
      legend.text     = ggplot2::element_text(size = 8),
      plot.background = ggplot2::element_rect(fill = "white", color = NA),
      plot.margin     = ggplot2::margin(15, 20, 15, 20)
    ) +
    ggplot2::coord_fixed(xlim = c(-(r + 1.2), r + 1.2),
                         ylim = c(-(r + 1.0), r + 1.0))

  return(p)
}
