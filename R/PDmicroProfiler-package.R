#' @importFrom magrittr %>%
#' @importFrom rlang .data :=
#' @importFrom stats na.omit
#' @keywords internal
"_PACKAGE"

utils::globalVariables(c(
  ".", "sample", "group", "area", "abundance", "microbe", "type",
  "qseqid", "Family", "lifestyle", "Specialization", "novelty_status", "Length",
  "Host.Phylum", "Host.Class", "Host.Order", "Host.Family", "Host.Genus", "Host.Species",
  "Count", "Pathway_Name", "Gene_Count", "Avg_pvalue", "KO", "pathway", "direction",
  "Enriched_in_PD", "Enriched_in_HC", "rate_PD", "rate_HC", "significance",
  "pathway_list", "occurrence_rate", "KO_label", "Source", "Target", "Correlation",
  "Phylum_bacteria", "Class_bacteria", "Order_bacteria", "Family_bacteria",
  "Genus_bacteria", "Species", "Phylum_virus", "Class_virus", "Order_virus",
  "Family_virus", "id", "label", "title", "from", "to", "value", "color",
  "contingency_table", "fisher_test", "p.value", "adj.p.value", "rate_diff",
  "vOTU", "p", "group_area", "bacteria_data", "virus_data", "bacteria_tax", "virus_tax",
  "votu_db", "interaction_edges_HC", "interaction_edges_PD", "interaction_edges_shared",
  "interaction_nodes_HC", "interaction_nodes_PD", "interaction_nodes_shared",
  "bac_tax", "votu_tax", "significant_ko_data",
  # 新增变量
  "fold_change", "log2_fold_change", "min_p", "median_log2fc", "n_KO",
  "count", "a", "b", "d", "cont_mat", "fisher", "odds_ratio", "log2_odds_ratio",
  "gene_count", "total", "avg_logp", "plot_value", "prop", "na.omit",
  "n"  # 新增
))
