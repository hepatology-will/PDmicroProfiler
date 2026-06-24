#' Bacteria-phage interaction nodes for HC group
#'
#' A data frame of nodes (bacteria and phages) in the HC interaction network.
#'
#' @docType data
#' @usage data(interaction_nodes_HC)
#' @format A data frame with columns:
#' \describe{
#'   \item{id}{Node ID}
#'   \item{group}{Node type: "bacteria" or "phage"}
#'   \item{Family}{Phage family or bacterial taxonomy (if available)}
#'   \item{enrichment}{Enrichment information (optional)}
#' }
#' @source \url{https://example.com}
"interaction_nodes_HC"

#' Bacteria-phage interaction nodes for PD group
#'
#' A data frame of nodes (bacteria and phages) in the PD interaction network.
#'
#' @docType data
#' @usage data(interaction_nodes_PD)
#' @format A data frame with columns:
#' \describe{
#'   \item{id}{Node ID}
#'   \item{group}{Node type: "bacteria" or "phage"}
#'   \item{Family}{Phage family or bacterial taxonomy (if available)}
#'   \item{enrichment}{Enrichment information (optional)}
#' }
#' @source \url{https://example.com}
"interaction_nodes_PD"

#' Bacteria-phage interaction nodes shared between groups
#'
#' A data frame of nodes (bacteria and phages) in the shared interaction network.
#'
#' @docType data
#' @usage data(interaction_nodes_shared)
#' @format A data frame with columns:
#' \describe{
#'   \item{id}{Node ID}
#'   \item{group}{Node type: "bacteria" or "phage"}
#'   \item{Family}{Phage family or bacterial taxonomy (if available)}
#'   \item{enrichment}{Enrichment information (optional)}
#' }
#' @source \url{https://example.com}
"interaction_nodes_shared"
