#' Get Microbe Abundance Data
#'
#' Retrieves abundance data for a single microbe with built-in fuzzy matching.
#' Supports both bacteria and viruses. Handles prefixes like 's__' in bacterial names.
#'
#' @param name Single microbe name for query (required). Supports fuzzy matching.
#' @param type Microbe type: "bacteria", "virus", or "auto" (automatic detection, default).
#' @param group Group to subset: "PD" or "HC". Default is NULL (return all groups).
#' @param area Area(s) to subset: e.g., "USA", "Finland", "China1", "Japan", "China2".
#'             Can be a single value or vector. Default is NULL (return all areas).
#' @param exact Logical. If TRUE, perform exact matching; if FALSE (default), fuzzy matching.
#' @param ignore_case Logical. Whether to ignore case in matching (default: TRUE).
#' @param remove_prefix Logical. Whether to ignore taxonomic prefixes (e.g., 's__') when matching (default: TRUE).
#' @return A data frame with columns: sample, abundance, group, area, microbe, type.
#' @export
#'
#' @examples
#' \dontrun{
#' # Get abundance for a bacterium
#' get_microbe_abundance("Achromobacter_xylosoxidans")
#'
#' # Get abundance for a virus in PD group in Finland
#' get_microbe_abundance("vOTU_123", type = "virus", group = "PD", area = "Finland")
#' }
get_microbe_abundance <- function(name, type = "auto", group = NULL, area = NULL,
                                  exact = FALSE, ignore_case = TRUE,
                                  remove_prefix = TRUE) {

  if (missing(name) || length(name) == 0) {
    stop("A single microbe name must be provided")
  }

  if (length(name) > 1) {
    name <- name[1]
    warning("Only single microbe name is supported. Using the first name.")
  }

  # 直接使用内部数据（已在 globalVariables 中声明）
  if (!exists("group_area") || !exists("bacteria_data") || !exists("virus_data")) {
    stop("Required data not loaded. Please ensure the package is properly installed.")
  }

  # 匹配函数
  find_matches <- function(query, dataset, is_bacteria = FALSE) {
    candidates <- rownames(dataset)

    if (exact) {
      if (ignore_case) {
        matches <- candidates[tolower(candidates) == tolower(query)]
      } else {
        matches <- candidates[candidates == query]
      }
    } else {
      if (ignore_case) {
        pattern <- tolower(query)
        search_pool <- tolower(candidates)
      } else {
        pattern <- query
        search_pool <- candidates
      }

      matches <- candidates[grep(pattern, search_pool)]

      # 如果需要，去除前缀后再次匹配
      if (remove_prefix && is_bacteria) {
        no_prefix <- gsub("^[a-z]__", "", candidates)
        if (ignore_case) {
          additional_matches <- grep(tolower(query), tolower(no_prefix))
        } else {
          additional_matches <- grep(query, no_prefix)
        }
        matches <- unique(c(matches, candidates[additional_matches]))
      }
    }

    return(matches)
  }

  # 确定数据集和匹配
  if (type == "auto") {
    bacteria_matches <- find_matches(name, bacteria_data, TRUE)
    virus_matches <- find_matches(name, virus_data, FALSE)

    if (length(bacteria_matches) > 0) {
      dataset <- bacteria_data
      matches <- bacteria_matches
      dataset_name <- "bacteria"
    } else if (length(virus_matches) > 0) {
      dataset <- virus_data
      matches <- virus_matches
      dataset_name <- "virus"
    } else {
      stop(sprintf("No matches found for '%s' in any dataset.", name))
    }
  } else {
    type <- match.arg(type, c("bacteria", "virus"))
    dataset <- if (type == "bacteria") bacteria_data else virus_data
    matches <- find_matches(name, dataset, type == "bacteria")
    dataset_name <- type

    if (length(matches) == 0) {
      stop(sprintf("No matches found for '%s' in %s dataset.", name, dataset_name))
    }
  }

  # 报告匹配
  if (length(matches) > 1) {
    message(sprintf("Found %d matches for '%s': %s",
                    length(matches), name, paste(matches, collapse = ", ")))
  }

  # 提取数据
  extract_data <- function(microbe_name, dataset) {
    if (!microbe_name %in% rownames(dataset)) {
      return(NULL)
    }

    abundance <- as.numeric(dataset[microbe_name, ])
    result <- data.frame(
      sample = colnames(dataset),
      abundance = abundance,
      stringsAsFactors = FALSE
    )

    result <- merge(result, group_area, by = "sample", all.x = TRUE)

    if (!is.null(group)) {
      group <- match.arg(group, c("PD", "HC"))
      result <- result[result$group == group, ]
    }

    if (!is.null(area)) {
      valid_areas <- unique(group_area$area)
      if (all(area %in% valid_areas)) {
        result <- result[result$area %in% area, ]
      } else {
        invalid <- setdiff(area, valid_areas)
        warning(sprintf("Invalid area(s): %s. Using all areas.", paste(invalid, collapse = ", ")))
      }
    }

    result$microbe <- microbe_name
    result$type <- dataset_name

    return(result)
  }

  # 处理多个匹配
  if (length(matches) == 1) {
    result <- extract_data(matches, dataset)
  } else {
    result_list <- lapply(matches, function(m) extract_data(m, dataset))
    result_list <- result_list[!sapply(result_list, is.null)]

    if (length(result_list) == 0) {
      stop("No valid matches found.")
    }

    result <- do.call(rbind, result_list)
  }

  # 检查NA值
  if (any(is.na(result$abundance))) {
    warning(sprintf("Found %d NA values in abundance data.", sum(is.na(result$abundance))))
  }

  rownames(result) <- NULL
  return(result)
}

#' Visualize Microbe Abundance
#'
#' Creates boxplots comparing microbe abundance between PD and HC groups.
#' Supports single microbe visualization with built-in fuzzy matching.
#' Can compare across multiple areas (countries).
#'
#' @param name Microbe name (required). Supports fuzzy matching.
#' @param type Microbe type: "bacteria", "virus", or "auto" (default).
#' @param group Group to subset: "PD", "HC", or NULL (all groups, default).
#' @param area Area(s) to compare: e.g., "USA", "Finland", etc.
#'             Can be a single value or vector. Default is NULL (use all areas).
#' @param log_transform Apply log10 transformation? Default: TRUE.
#' @param significance Control significance annotation:
#'                     - "none": no significance annotation (default)
#'                     - "all": show all significance levels (including "ns")
#'                     - "sig_only": only show significant results (p < 0.05)
#' @param group_colors Custom colors for groups.
#' @param show_points Show individual data points? Default: FALSE.
#' @param point_size Size of points if shown. Default: 1.5.
#' @param point_alpha Transparency of points if shown. Default: 0.5.
#' @param ... Additional arguments passed to get_microbe_abundance().
#' @return A ggplot2 object.
#' @export
#'
#' @examples
#' \dontrun{
#' plot_microbe_abundance("Achromobacter_xylosoxidans", log_transform = TRUE)
#' }
plot_microbe_abundance <- function(name, type = "auto", group = NULL, area = NULL,
                                   log_transform = TRUE, significance = "none",
                                   group_colors = c("PD" = "#a7434f", "HC" = "#ABCBDF"),
                                   show_points = FALSE, point_size = 1.5,
                                   point_alpha = 0.5, ...) {

  # 获取数据
  data <- get_microbe_abundance(name = name, type = type, group = group, area = area, ...)

  if (is.null(data) || nrow(data) == 0) {
    stop("No data to plot.")
  }

  # 验证significance参数
  significance <- match.arg(significance, c("none", "all", "sig_only"))

  # 清理数据
  data <- data[!is.na(data$abundance), ]

  # 准备绘图数据
  if (log_transform) {
    data$plot_value <- log10(data$abundance + 1e-10)
    y_label <- "log10(Abundance)"
  } else {
    data$plot_value <- data$abundance
    y_label <- "Abundance"
  }

  # 确定分面
  n_microbes <- length(unique(data$microbe))
  n_areas <- length(unique(data$area))

  # 自动生成标题
  prefix <- if (unique(data$type)[1] == "bacteria") "Bacteria: " else "Virus: "
  microbe_name <- unique(data$microbe)[1]

  if (n_areas == 1) {
    title <- sprintf("%s%s in %s", prefix, microbe_name, unique(data$area)[1])
  } else if (n_areas > 1) {
    if (n_areas <= 3) {
      title <- sprintf("%s%s in %s", prefix, microbe_name, paste(unique(data$area), collapse = ", "))
    } else {
      title <- sprintf("%s%s in %d areas", prefix, microbe_name, n_areas)
    }
  } else {
    title <- sprintf("%s%s", prefix, microbe_name)
  }

  # 构建基础图形
    p <- ggplot2::ggplot(
    data,
    ggplot2::aes(x = group, y = plot_value, fill = group, color = group)
  ) +
    ggplot2::geom_violin(
      width = 0.65,
      alpha = 0.16,
      trim = TRUE,
      linewidth = 0.35
    ) +
    ggplot2::geom_jitter(
      width = 0.08,
      size = 0.9,
      alpha = 0.35,
      shape = 16,
      show.legend = FALSE
    ) +
    ggplot2::geom_boxplot(
      width = 0.12,
      alpha = 0.55,
      outlier.shape = NA,
      color = "black",
      linewidth = 0.35
    ) +
    ggplot2::scale_fill_manual(values = group_colors) +
    ggplot2::scale_color_manual(values = group_colors) +
    ggplot2::labs(title = title, x = "Group", y = y_label) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      panel.grid.major = ggplot2::element_line(color = "grey92"),
      panel.grid.minor = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(color = "black", linewidth = 0.8),
      axis.text = ggplot2::element_text(size = 12, color = "black"),
      axis.title = ggplot2::element_text(size = 14, face = "bold"),
      legend.position = "top",
      legend.title = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(size = 16, face = "bold", hjust = 0.5)
    )


  # 添加分面
  if (n_areas > 1) {
    p <- p + ggplot2::facet_wrap(~ area, scales = "free_y") +
      ggplot2::theme(strip.text = ggplot2::element_text(size = 10, face = "bold"))
  }

  # 显著性检验
  if (significance != "none" && n_microbes == 1 && requireNamespace("ggsignif", quietly = TRUE)) {

    # 辅助函数：计算p值符号和标签位置
    calculate_significance <- function(data_subset) {
      if (is.null(data_subset) || nrow(data_subset) < 3) {
        return(NULL)
      }

      groups <- unique(data_subset$group)
      if (length(groups) < 2) {
        return(NULL)
      }

      group_counts <- table(data_subset$group)
      if (any(group_counts < 2)) {
        return(NULL)
      }

      if (stats::var(data_subset$plot_value, na.rm = TRUE) == 0) {
        return(NULL)
      }

      tryCatch({
        test_result <- stats::wilcox.test(plot_value ~ group, data = data_subset,
                                          exact = FALSE, na.action = na.omit)
        p_value <- test_result$p.value

        if (is.na(p_value) || is.nan(p_value)) {
          return(NULL)
        }

        if (p_value < 0.001) {
          symbol <- "***"
        } else if (p_value < 0.01) {
          symbol <- "**"
        } else if (p_value < 0.05) {
          symbol <- "*"
        } else {
          symbol <- "ns"
        }

        if (significance == "sig_only" && symbol == "ns") {
          return(NULL)
        }

        y_values <- data_subset$plot_value
        y_min <- min(y_values, na.rm = TRUE)
        y_max <- max(y_values, na.rm = TRUE)
        y_range <- y_max - y_min

        if (y_range < 0.1) {
          label_position <- y_max + 0.1
        } else {
          label_position <- y_max + y_range * 0.1
        }

        return(list(
          symbol = symbol,
          p_value = p_value,
          position = label_position,
          data = data_subset
        ))
      }, error = function(e) {
        return(NULL)
      })
    }

    # 处理单个或多个区域的情况
    if (n_areas > 1) {
      for (area_name in unique(data$area)) {
        area_data <- data[data$area == area_name, ]
        sig_info <- calculate_significance(area_data)

        if (!is.null(sig_info)) {
          p <- p + ggsignif::geom_signif(
            data = area_data,
            comparisons = list(c("HC", "PD")),
            annotations = sig_info$symbol,
            y_position = sig_info$position,
            tip_length = 0.02,
            textsize = 4,
            vjust = 0.5,
            map_signif_level = FALSE
          )
        }
      }
    } else {
      sig_info <- calculate_significance(data)

      if (!is.null(sig_info)) {
        p <- p + ggsignif::geom_signif(
          comparisons = list(c("HC", "PD")),
          annotations = sig_info$symbol,
          y_position = sig_info$position,
          tip_length = 0.02,
          textsize = 4,
          vjust = 0.5,
          map_signif_level = FALSE
        )
      }
    }
  } else if (significance != "none" && !requireNamespace("ggsignif", quietly = TRUE)) {
    warning("Package 'ggsignif' is required for significance annotation. Please install it.")
  }

  return(p)
}
