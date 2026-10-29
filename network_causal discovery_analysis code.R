#-------network analysis and causal discovery-------

# load packages

library(haven)
library(qgraph)
library(igraph)
library(pcalg)
library(dplyr)
library(IsingFit)
library(bootnet)
library(reshape2)
library(ggplot2)
library(ggrepel)


# data preparation ------------------------------------------------------
# SHARE data are not included in this repository.
# Authorized users should download the SHARE Wave 2 data
# and place the file in the local "data" directory.

data_path <- "data"
results_path <- "results"

# Create output directory
dir.create(results_path, recursive = TRUE, showWarnings = FALSE)

# Load SHARE Wave 2 data
data_file <- file.path(data_path, "sharew2_rel9_health.sav")

if (!file.exists(data_file)) {
  stop(
    "SHARE data file not found. ",
    "Please place sharew2_rel9_health.sav in the data/ directory."
  )
}

df_share <- haven::read_sav(data_file)
df_share$eurod <- as.numeric(zap_labels(df_share$eurod))


# Keep participants with non-zero EURO-D score
df_depr <- df_share %>%
  filter(eurod != 0) %>%
  select(1, starts_with("euro"))


# Undirected Network ------------------------------------------------------

# Select the 12 EURO-D symptom items
df_analysis <- df_depr[, 2:13]

# Remove SPSS value labels
df_network <- as.data.frame(
  lapply(df_analysis, haven::zap_labels)
)

# Estimate network
fit_network <- IsingFit(
  df_network,
  gamma = 0.25,
  plot = FALSE
)

# set labels
symptom_labels <- c(
  "depr",      # euro1
  "pess",           # euro2
  "suic",        # euro3
  "guil",               # euro4
  "slep",      # euro5
  "intr",    # euro6
  "irri",             # euro7
  "appt",    # euro8
  "fatg",        # euro9
  "conc", # euro10
  "enjy",    # euro11
  "tear"         # euro12
)

# plot network
pdf(file.path(results_path, "network.pdf"), width = 8, height = 8)
q_network = qgraph(fit_network$weiadj,
                   layout = "circular",
                   repulsion = 0.8,
                   edge.color = fit_network$color,
                   color = rep("#AEDFF780", 12),
                   labels = symptom_labels,
                   vsize = 8,
                   label.cex =1.2,
                   label.scale = FALSE)
dev.off()

# save layout
layout_fixed <- q_network$layout

# possible edges
sum(fit_network$weiadj != 0) / 2
n_nodes <- ncol(fit_network$weiadj)
n_possible_edges <- n_nodes * (n_nodes - 1) / 2


# stability
colnames(df_network) <- symptom_labels
boot1 <- bootnet(df_network, nBoots = 2000, 
                 default = "IsingFit", 
                 type = "nonparametric")

boot2 <- bootnet(df_network, nBoots = 2000, 
                 default = "IsingFit", 
                 type = "case",
                 statistics = c("strength", "closeness", "betweenness"))

plot(boot1, labels = TRUE, order = "sample")
plot(boot1, plot = "difference", onlyNonZero = TRUE, order = "sample", labels = TRUE)
plot(boot1, statistics = "strength", plot = "difference", order = "sample", labels = TRUE)
plot(boot2, labels = TRUE, statistics = c("strength", "closeness", "betweenness"), plot = "area")

corStability(boot2)


# Causal network ------------------------------------------------------
# estimate PAG
fit_fci <- fci(suffStat = list(C = cor(df_analysis), n = nrow(df_analysis)), 
               indepTest = gaussCItest, 
               selectionBias=FALSE,
               alpha = 0.1, 
               verbose=TRUE,
               labels = symptom_labels)

plot.fci.igraph <- function(fci){
  require(igraph) # using igraph
  amat <- fci@amat
  amat.dir <- amat
  amat.dir[amat.dir==1 | amat.dir==3] <- 0
  amat.dir[amat.dir==2 & t(amat.dir)==3] <- 1
  amat.dbl <- amat
  amat.dbl[amat.dbl==2 | amat.dbl==3] <- 0
  amat.tot <- amat.dir+amat.dbl
  amat.tot[amat.tot>1] <- 1
  amat.col <- amat.tot
  amat.col[amat.col==1 & t(amat.col)==0] <- 3
  graph <- graph_from_adjacency_matrix(amat.tot,weighted=TRUE)
  graph.col <- graph_from_adjacency_matrix(amat.col,weighted=TRUE)
  E(graph)$color <- ifelse(E(graph.col)$weight==3, "red", "blue")
  l <- layout_fixed
  plot(graph,
       edge.arrow.size=0.5,
       vertex.size=20,
       vertex.label.cex = 1,
       vertex.frame.color="black",
       vertex.color= "#AEDFF780",
       vertex.label.color="black", 
       edge.curved=FALSE, 
       rescale=TRUE,
       layout=l*1.1)
}

# save plot
pdf(file.path(results_path, "fci_unstable.pdf"), width = 12, height = 8)
plot.fci.igraph(fit_fci)
dev.off()

# FCI subsampling
# create 500 subsamples at 80% of total sample size
# estimate causal graph
# strength of evidence > 0.75 proportion of subsamples edge present

set.seed(587)
nsamples <- 500
n <- nrow(df_analysis)
size <- 0.8
subsize <- round(size*n,0)
nvar <- ncol(df_analysis)
fci.subs <- array(0,dim=c(nvar,nvar,nsamples))

for(i in 1:nsamples){
  subsample <- sample(1:n,subsize,replace=FALSE)
  fci.sub <- df_analysis[subsample,]
  
  fci.output.subs <- fci(
    suffStat=list(C = cov2cor(cov(fci.sub)),
                  n = nrow(fci.sub)),
    indepTest=gaussCItest,
    selectionBias=FALSE,
    labels=symptom_labels,
    alpha=0.1,
    verbose=FALSE)
  
  amat <- fci.output.subs@amat
  amat.dir <- amat
  amat.dir[amat==1] <- 0 # remove uncertain connections
  amat.dir[amat==2 & t(amat==3)] <- 1 # if i->j 
  amat.dir[amat.dir!=1] <- 0
  
  amat.bid <- amat
  amat.bid[amat==2 & t(amat==3)] <- 0 # remove if i->j 
  amat.bid[amat==3 & t(amat==2)] <- 0 # remove if j->i 	
  amat.bid[amat==2 & t(amat==2)] <- 1 # if i<->j 
  amat.bid[amat==2 & t(amat==1)] <- 1 # if io->j 
  amat.bid[amat==3 & t(amat==1)] <- 1 # if io--j 
  
  amat <- amat.dir + amat.bid
  fci.subs[,,i] <- amat
}

res <- apply(fci.subs,c(1,2),sum)/nsamples

colnames(res) <- rownames(res) <- symptom_labels

alpha = 0.75 # threshold for subsamples

fci.stable <- res
fci.stable[fci.stable >= alpha] <- 1
fci.stable[fci.stable < alpha] <- 0


# possible edges
n_possible_directed <- 12 * 11
directed_mat <- (fci.stable == 1) & (t(fci.stable) == 0)
diag(directed_mat) <- FALSE

n_detected_directed <- sum(directed_mat)

n_possible_directed
n_detected_directed


#FCI plot
amat_stable <- fci.stable
amat_stable[amat_stable==1 & t(amat_stable)==0] <- 3
graph_stable <- graph_from_adjacency_matrix(amat_stable,weighted=TRUE)
graph_stable.col <- graph_from_adjacency_matrix(amat_stable,weighted=TRUE)
E(graph_stable)$color <- ifelse(E(graph_stable.col)$weight==3, "red", "blue")

pdf(file.path(results_path, "fci_stable.pdf"),width = 12,height = 8)

plot(graph_stable,
     edge.arrow.size=0.5,
     vertex.size=20,
     vertex.label.cex = 1,
     vertex.frame.color="black",
     vertex.color= "#AEDFF780",
     vertex.label.color="black", 
     edge.curved=FALSE, 
     rescale=TRUE,
     layout=layout_fixed*1.1)

dev.off()


# centrality comparison ---------------------------------------------

# network
cent_net <- centrality_auto(q_network)
Strength <- scale(cent_net$node.centrality$Strength)

#PAG
M <- fci.stable
dir_mat <- (M == 1) & (t(M) == 0)
diag(dir_mat) <- FALSE
storage.mode(dir_mat) <- "integer"

out_deg <- as.numeric(scale(rowSums(dir_mat)))
in_deg  <- as.numeric(scale(colSums(dir_mat)))

# Combine
df_cent_c <- data.frame(
  Node = symptom_labels,
  InDegree = in_deg,
  OutDegree = out_deg,
  Strength = Strength
)

# comparison plot
df_melt <- melt(df_cent_c, id.vars = "Node")
df_melt$Node <- factor(df_melt$Node, levels = df_cent_c$Node)

plot_com <- ggplot(df_melt, aes(x = Node, y = value, color = variable, group = variable)) +
  geom_line() +
  geom_point() +
  labs(y = "Z-score", x = "Node", color = NULL) +
  theme_classic() +
  theme(panel.grid.major.y = element_line(linetype = "dashed", color = "grey80"),
        plot.title = element_text(size = 10, hjust = 0)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

plot_com

# correlation plot
cor_fci <- cor.test(df_cent_c$OutDegree,df_cent_c$Strength,
                    method = "spearman",exact = FALSE)

cor_fci

plot_cor_fci <- ggplot(df_cent_c, aes(x = OutDegree, y = Strength)) +
  geom_point(size = 2.5) +
  ggrepel::geom_text_repel(
    aes(label = Node),
    size = 3.5
  ) +
  annotate(
    "text",
    x = -Inf, y = Inf,
    hjust = -0.05, vjust = 1.2,
    label = "Spearman \u03c1 = .41, p = .182"
  ) +
  labs(
    x = "FCI outdegree",
    y = "Strength centrality"
  ) +
  theme_classic()

#save plot
ggsave(
  file.path(results_path, "centrality_outdegree.pdf"),
  plot = plot_cor_fci,
  width = 7,
  height = 5,
  units = "in"
)

ggsave(
  file.path(results_path, "centrality_degree_comparison.pdf"),
  plot = plot_com,
  width = 7,
  height = 5,
  units = "in"
)



# Generalized Adjustment Criterion (GAC) ------------------------------
# Exploratory analysis of the full-sample FCI PAGFCI causal effect

amat <- fit_fci@amat

# Extract directed edges: X -> Y
# pcalg: 2 = arrowhead, 3 = tail
dir_mat <- (amat == 2) & (t(amat) == 3)
diag(dir_mat) <- FALSE

edge_idx <- which(dir_mat, arr.ind = TRUE)

full_edges <- data.frame(
  From = symptom_labels[edge_idx[, 1]],
  To = symptom_labels[edge_idx[, 2]],
  stringsAsFactors = FALSE
)

cat("Number of directed edges:", nrow(full_edges), "\n")
print(full_edges)


# Find valid adjustment sets
get_adjustment_sets <- function(amat, x, y) {
  
  # Explicitly test the empty adjustment set
  empty_test <- pcalg::gac(
    amat = amat,
    x = x,
    y = y,
    z = integer(0),
    type = "pag"
  )
  
  # Empty set {} is valid
  if (isTRUE(empty_test$gac)) {
    return(list(integer(0)))
  }
  
  # Otherwise search for minimal non-empty adjustment sets
  adj <- pcalg::adjustment(
    amat = amat,
    amat.type = "pag",
    x = x,
    y = y,
    set.type = "minimal"
  )
  
  adj
}


# Apply GAC to every directed edge
adj_results <- lapply(seq_len(nrow(full_edges)), function(i) {
  
  x <- edge_idx[i, 1]
  y <- edge_idx[i, 2]
  
  sets <- get_adjustment_sets(amat, x, y)
  
  identifiable <- length(sets) > 0
  
  if (!identifiable) {
    
    n_sets <- 0L
    min_size <- NA_integer_
    set_text <- NA_character_
    
  } else {
    
    n_sets <- length(sets)
    min_size <- min(lengths(sets))
    
    set_text <- paste(
      vapply(sets, function(z) {
        
        if (length(z) == 0) {
          "{}"
        } else {
          paste0(
            "{",
            paste(symptom_labels[z], collapse = ", "),
            "}"
          )
        }
        
      }, character(1)),
      collapse = " | "
    )
  }
  
  data.frame(
    From = full_edges$From[i],
    To = full_edges$To[i],
    Valid_adjustment = identifiable,
    N_minimal_sets = n_sets,
    Min_set_size = min_size,
    Adjustment_sets = set_text,
    stringsAsFactors = FALSE
  )
})


# results
adjustment_table <- do.call(rbind, adj_results)

print(adjustment_table, row.names = FALSE)