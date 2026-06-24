#' Significant KO pathway data for vOTUs
#'
#' A dataset containing KO and pathway annotations for significant vOTUs.
#'
#' @docType data
#' @usage data(significant_ko_data)
#' @format A data frame with columns:
#' \describe{
#'   \item{vOTU}{vOTU ID}
#'   \item{Family}{Phage family}
#'   \item{KO}{KO number}
#'   \item{pathway}{Pathway description}
#'   \item{direction}{Enrichment direction: "Enriched_in_PD" or "Enriched_in_HC"}
#'   \item{p.value}{P-value}
#'   \item{adj.p.value}{Adjusted p-value}
#'   \item{fold_change}{Fold change}
#'   \item{log2_fold_change}{Log2 fold change}
#'   \item{Protein}{Protein annotation}
#'   \item{median_CL}{Median in control group (if applicable)}
#'   \item{median_PD}{Median in PD group (if applicable)}
#'   \item{p}{Raw p-value (if different from p.value)}
#' }
#' @source \url{https://example.com}
"significant_ko_data"
