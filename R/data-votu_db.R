#' Phage database
#'
#' A dataset containing information about phages, including host taxonomy,
#' lifestyle, specialization, novelty status, and genome length.
#'
#' @docType data
#' @usage data(votu_db)
#' @format A data frame with rows for each phage and columns:
#' \describe{
#'   \item{qseqid}{Phage ID}
#'   \item{Family}{Phage family}
#'   \item{lifestyle}{Lifestyle (temperate/virulent)}
#'   \item{Specialization}{Host specialization}
#'   \item{novelty_status}{Novelty status}
#'   \item{Length}{Genome length}
#'   \item{Host.Domain}{Host domain}
#'   \item{Host.Phylum}{Host phylum}
#'   \item{Host.Class}{Host class}
#'   \item{Host.Order}{Host order}
#'   \item{Host.Family}{Host family}
#'   \item{Host.Genus}{Host genus}
#'   \item{Host.Species}{Host species}
#'   \item{Lineage}{Full lineage (if available)}
#'   \item{Species_rep}{Species representative}
#' }
#' @source \url{https://example.com}
"votu_db"
