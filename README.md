# PDmicroProfiler

<p align="center">
  <img src="PDmicroProfiler_README_workflow.png" width="1000">
</p>

<p align="center">
  <strong>An R package for exploring Parkinson's disease-associated gut bacteriome, virome, and phage–bacteria interactions</strong>
</p>

---

## Overview

**PDmicroProfiler** is an R package designed for standardized exploration and visualization of Parkinson's disease (PD)-associated gut bacterial and viral signatures.

The package integrates processed multi-cohort metagenomic resources and provides a unified framework for:

* querying bacterial and viral abundance profiles;
* comparing microbial features between PD and healthy control groups;
* exploring viral operational taxonomic units (vOTUs);
* summarizing phage family characteristics;
* visualizing viral functional annotations;
* investigating phage-host and bacteria–phage interaction networks.

PDmicroProfiler was developed to make PD-associated gut microbiome and virome resources easier to access, visualize, and interpret.

---

## Background

Parkinson's disease has been increasingly linked to alterations in the gut microbial ecosystem. However, compared with bacterial communities, the viral component of the gut microbiome, especially bacteriophages, remains relatively underexplored.

PDmicroProfiler was developed from an integrated multi-cohort metagenomic analysis of PD and healthy control fecal samples. It supports reproducible investigation of gut bacterial taxa, viral vOTUs, viral functional features, and inferred bacteria–phage interaction patterns.

---

## Installation

You can install the development version of **PDmicroProfiler** from GitHub using `remotes`:

```r
# install.packages("remotes")
remotes::install_github("hepatology-will/PDmicroProfiler")
```

Alternatively, you can use `devtools`:

```r
# install.packages("devtools")
devtools::install_github("hepatology-will/PDmicroProfiler")
```

Load the package:

```r
library(PDmicroProfiler)
```

---

## Main features

PDmicroProfiler provides four major functional modules for exploring PD-associated gut bacteriome, virome, viral functions, and bacteria–phage interactions.

| Module                          | Description                                                                                                     | Representative functions                                                                                                                    |
| ------------------------------- | --------------------------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------- |
| Microbial abundance profiling   | Query and visualize bacterial taxa, viral families, or viral vOTUs across PD and healthy control samples.       | `get_microbe_abundance()`<br>`plot_microbe_abundance()`                                                                                     |
| Phage database exploration      | Explore phages by host, family, lifestyle, specialization, and related phage characteristics.                   | `query_by_host()`<br>`query_by_id()`<br>`phage_family_dashboard()`<br>`plot_phage_summary()`                                                |
| Viral functional annotation     | Retrieve and visualize vOTU-level and family-level KEGG Ortholog annotations and differential pathway profiles. | `get_votu_info()`<br>`get_ko_info()`<br>`get_pathway_info()`<br>`plot_votu_function()`<br>`plot_family_function()`<br>`plot_diff_pathway()` |
| Bacteria–phage network analysis | Identify and visualize phage-host or bacteria–phage interaction networks in PD, HC, or shared networks.         | `find_phage_host_interactions()`<br>`plot_interaction_network()`<br>`plot_interaction_network_static()`                                     |

---

## Data resources

PDmicroProfiler includes processed data resources generated from multi-cohort PD gut metagenomic datasets.

| Data type                   | Description                                                             |
| --------------------------- | ----------------------------------------------------------------------- |
| Bacterial abundance         | Bacterial abundance profiles for PD and healthy control samples         |
| Viral vOTU abundance        | Viral operational taxonomic unit abundance profiles                     |
| Sample metadata             | Group and cohort information                                            |
| Viral taxonomy              | Family-level and vOTU-level viral taxonomic annotations                 |
| Viral functional annotation | KEGG Ortholog and pathway annotations for vOTUs                         |
| Phage-host links            | Predicted associations between phages and bacterial hosts               |
| Interaction networks        | Bacteria–phage association networks for PD, HC, and shared interactions |

---

## Quick start

The following examples demonstrate the four major modules of **PDmicroProfiler**.

### 1. Abundance profiling

Retrieve sample-level abundance data for a selected bacterial species and visualize its distribution between PD and healthy control samples.

```r
library(PDmicroProfiler)

abund_df <- get_microbe_abundance(
  name = "Bifidobacterium_catenulatum",
  type = "bacteria",
  area = "USA"
)

head(abund_df)

plot_microbe_abundance(
  name = "Bifidobacterium_catenulatum",
  type = "bacteria",
  area = "USA"
)
```

---

### 2. Phage family dashboard

Summarize and visualize the characteristics of a selected phage family.

```r
family_info <- get_family_info(
  family_name = "Siphoviridae"
)

head(family_info)

phage_family_dashboard(
  family_name = "Siphoviridae",
  host_level  = "Host.Phylum",
  top_n       = 10
)
```

---

### 3. Viral functional annotation

Retrieve functional annotation information for a selected phage family and visualize the most frequent KEGG Ortholog annotations.

```r
family_info <- get_family_info(
  family_name = "Myoviridae"
)

head(family_info)

plot_family_function(
  family_name = "Myoviridae",
  top_n = 10,
  plot_type = "bar"
)
```

---

### 4. Bacteria–phage network analysis

Load the built-in bacteria–phage interaction network, extract Streptococcus-associated edges, and visualize the corresponding interaction network.

```r
data("interaction_edges_HC")

strepto_edges <- subset(
  interaction_edges_HC,
  grepl("Streptococcus", Source)
)

head(strepto_edges)

plot_interaction_network(
  strepto_edges,
  title = "Streptococcus phage interaction network"
)
```

---

## Function reference

### Abundance profiling

| Function                   | Description                                                                 |
| -------------------------- | --------------------------------------------------------------------------- |
| `get_microbe_abundance()`  | Retrieve sample-level abundance of a bacterial taxon, viral family, or vOTU |
| `plot_microbe_abundance()` | Visualize microbial abundance differences between groups                    |

### Phage database exploration

| Function                    | Description                                                      |
| --------------------------- | ---------------------------------------------------------------- |
| `query_by_host()`           | Query phages associated with a selected host                     |
| `query_by_id()`             | Query phage or vOTU information by identifier                    |
| `query_by_phage_features()` | Query phages by family, lifestyle, host range, or novelty status |
| `phage_family_dashboard()`  | Generate a summary dashboard for a selected phage family         |
| `plot_phage_summary()`      | Visualize family-level composition of queried phages             |

### Viral functional annotation

| Function                 | Description                                                              |
| ------------------------ | ------------------------------------------------------------------------ |
| `get_votu_info()`        | Retrieve taxonomic, host, and functional information for a selected vOTU |
| `get_ko_info()`          | Retrieve information for a selected KEGG Ortholog                        |
| `get_pathway_info()`     | Retrieve pathway-level annotation information                            |
| `plot_votu_function()`   | Visualize KO annotations for a selected vOTU                             |
| `plot_family_function()` | Visualize frequent KO annotations in a viral family                      |
| `plot_diff_pathway()`    | Compare pathway enrichment patterns between PD- and HC-enriched vOTUs    |

### Network analysis

| Function                            | Description                                                |
| ----------------------------------- | ---------------------------------------------------------- |
| `find_phage_host_interactions()`    | Identify bacteria–phage interactions for a selected target |
| `plot_interaction_network()`        | Generate an interactive bacteria–phage network             |
| `plot_interaction_network_static()` | Generate a static bacteria–phage network                   |

---

## Example outputs

PDmicroProfiler can generate several types of visualization outputs:

* abundance comparison plots;
* multi-cohort abundance plots;
* phage family composition pie charts;
* phage family dashboards;
* vOTU-level functional bubble plots;
* family-level KO bar plots;
* differential pathway plots;
* interactive bacteria–phage networks;
* static network plots.

These outputs are designed to support exploratory analysis and publication-ready visualization of PD-associated gut microbiome and virome features.

---

## Repository structure

```text
PDmicroProfiler/
├── R/                         # R functions
├── data/                      # Built-in package datasets
├── man/                       # Function documentation
├── DESCRIPTION                # Package metadata
├── NAMESPACE                  # Exported functions
├── LICENSE                    # License information
├── README.md                  # GitHub README
└── PDmicroProfiler.Rproj      # RStudio project file
```

---

## Citation

If you use **PDmicroProfiler** in your research, please cite the R package:

```bibtex
@software{PDmicroProfiler,
  author = {Jiale Liu and Wei Hou},
  title = {PDmicroProfiler: an integrated platform for multi-cohort profiling of gut bacterial and viral signatures in Parkinson's disease},
  year = {2026},
  url = {https://github.com/hepatology-will/PDmicroProfiler},
  version = {0.1.0}
}
```

## License

This project is licensed under the MIT License.
