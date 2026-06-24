#' Bacteria-phage interaction edges for HC group
#'
#' A data frame of edges (interactions) between bacteria and phages in the HC group.
#'
#' @docType data
#' @usage data(interaction_edges_HC)
#' @format A data frame with columns:
#' \describe{
#'   \item{Source}{Bacterial OTU ID}
#'   \item{Target}{vOTU ID}
#'   \item{Correlation}{Correlation value}
#'   \item{interaction}{Interaction type (if applicable)}
#' }
#' @source \url{https://example.com}
"interaction_edges_HC"

#' Bacteria-phage interaction edges for PD group
#'
#' A data frame of edges (interactions) between bacteria and phages in the PD group.
#'
#' @docType data
#' @usage data(interaction_edges_PD)
#' @format A data frame with columns:
#' \describe{
#'   \item{Source}{Bacterial OTU ID}
#'   \item{Target}{vOTU ID}
#'   \item{Correlation}{Correlation value}
#'   \item{interaction}{Interaction type (if applicable)}
#' }
#' @source \url{https://example.com}
"interaction_edges_PD"

#' Bacteria-phage interaction edges shared between groups
#'
#' A data frame of edges (interactions) between bacteria and phages that are shared.
#'
#' @docType data
#' @usage data(interaction_edges_shared)
#' @format A data frame with columns:
#' \describe{
#'   \item{Source}{Bacterial OTU ID}
#'   \item{Target}{vOTU ID}
#'   \item{Correlation}{Correlation value}
#'   \item{interaction}{Interaction type (if applicable)}
#' }
#' @source \url{https://example.com}
"interaction_edges_shared"
