# Visualization function for cell interactions
visualization <- function(input_file, output_dir, tr, ir, hr) {
  
  # Load cell data
  cell_data <- read.csv(input_file)
  cell_type_list <<- sort(unique(cell_data[[3]]))
  cell_type_num <- length(cell_type_list)
  input_file_name <- basename(tools::file_path_sans_ext(input_file))
  
  filename <- file.path(output_dir, paste0(input_file_name, "_model_TR_", tr, "_IR_", ir, "_HR_", hr, ".Rda"))
  load(file = filename)
  #the model file has variables "coef"    "confid"  "family"  "formula"
  #confid contains the columns "Coefficient" "SE"  "CI95_Lower"  "CI95_Upper"   

  # Process coefficient intensities
  coef_intensity <- coef[1:cell_type_num]
  coef_intensity <- ifelse(1:cell_type_num == 1, 
                           coef_intensity[1:cell_type_num], 
                           coef_intensity[1:cell_type_num] + coef_intensity[1])
  
  coef_interaction <- coef[(cell_type_num + 1):length(coef)]

  # Process cell type links and nodes
  cell_type_links = process_links(confid, cell_type_num, tr, ir)
  cell_type_nodes_list = process_nodes(coef_interaction, cell_type_num)
  cell_type_nodes = cell_type_nodes_list[[1]]
  cell_type_self_interaction_links = cell_type_nodes_list[[2]]
  
  # Further process nodes and save visualizations
  process_and_visualize("pairwise", cell_type_links, cell_type_nodes, input_file_name, output_dir, tr, ir, hr)
  process_and_visualize("self", cell_type_self_interaction_links, cell_type_nodes, input_file_name, output_dir, tr, ir, hr)
}

# Initializes an empty data frame with given column names and number of columns
initialize_dataframe <- function(ncols, col_names, nrow=NULL) {
  df <- data.frame(matrix(ncol = ncols))
  colnames(df) <- col_names
  return(df)
}

# Process interactions between different cell types
process_links <- function(confid, cell_type_num, tr, ir) {
  
  cell_type_links <- initialize_dataframe(4, c('from', 'to', 'weight', 'sign'))
  idx <- 1
  
  for (i in 1:(cell_type_num-1)){
    for (j in (i+1):(cell_type_num)){
      for (r in seq(ir, tr, ir)){
        
        current_interaction_name <- paste("InteractionmarkX", cell_type_list[i], "xX", cell_type_list[j], "x", r, sep = '')
        current_coef <- confid[current_interaction_name,"Coefficient"]

        if (confid[current_interaction_name,"CI95_Lower"] > 0 | confid[current_interaction_name,"CI95_Upper"] < 0){
           current_sign <- ifelse(sign(current_coef) == 1, 'darkblue', 'darkred')
        
           cell_type_links[idx,] <- c(cell_type_list[i], cell_type_list[j], abs(current_coef), current_sign)
           idx <- idx + 1
        }
      }
    }
  }
  return(cell_type_links)
}

# Process individual cell types' interactions and calculates influence
process_nodes <- function(coef_interaction, cell_type_num) {
  
  cell_type_nodes <- initialize_dataframe(2, c('cell_type', 'influence'), cell_type_num)
  cell_type_self_interaction_links <- initialize_dataframe(4, c('from', 'to', 'weight', 'sign'))
  

  interaction_names <- names(coef_interaction)
  idx <- 1
  
  for (i in 1:cell_type_num){
    current_interaction_name <- paste('X', cell_type_list[i], sep = '')
    current_influence <- 0
    
    for (j in 1:length(coef_interaction)){
      if (grepl(current_interaction_name, interaction_names[j], fixed=TRUE) 
          & grepl(paste('InteractionmarkX', cell_type_list[i], 'xX', cell_type_list[i], sep = ''), interaction_names[j], fixed=TRUE)){
        
        current_sign <- ifelse(sign(coef_interaction[j]) == 1, 'darkblue', 'darkred')
        
        cell_type_self_interaction_links[idx,] <- c(cell_type_list[i], paste(cell_type_list[i], '2', sep='_'), abs(coef_interaction[j]), current_sign)
        current_influence <- current_influence + coef_interaction[j]
        idx <- idx + 1
      }
    }
    cell_type_nodes[i,] <- c(cell_type_list[i], current_influence)
  }
  return(list(cell_type_nodes, cell_type_self_interaction_links))
}


process_and_visualize <- function(interaction_type, cell_type_links, cell_type_nodes, input_file_name, output_dir, tr, ir, hr) {
  if (interaction_type == 'pairwise'){
    x = c(0, -0.951, 0.951, -0.588, 0.588)
    y = c(1, 0.309, 0.309, -0.809, -0.809)
    cell_type_nodes = cbind(cell_type_nodes, x, y)
  } else if (interaction_type == 'self'){
    x = c(0, 0, 0, 0, 0)
    y = c(4, 3, 2, 1, 0)
    cell_type_nodes = cbind(cell_type_nodes, x, y)
    cell_type_nodes = cell_type_nodes
    cell_type_nodes_copy = cell_type_nodes
    cell_type_nodes_copy$cell_type = paste(cell_type_nodes_copy$cell_type,2,sep='_')
    cell_type_nodes_copy$x = cell_type_nodes_copy$x + 1
    cell_type_nodes = rbind(cell_type_nodes, cell_type_nodes_copy)
  }
  
  
  # Scaling functions
  scaling_node_weight <- function(weight) {
    positive_weight <- weight - min(weight)
    scaled_weight <- scale(positive_weight, center = FALSE, scale = max(positive_weight)/40) + 10
    return(scaled_weight)
  }
  
  scaling_edge_weight <- function(weight) {
    positive_weight <- weight - min(weight)
    scaled_weight <- scale(positive_weight, center = FALSE, scale = max(positive_weight)/19.9) + 0.1
    return(scaled_weight)
  }
  
  # Scale node weights
  cell_type_nodes$influence <- scaling_node_weight(as.numeric(cell_type_nodes$influence))

  # Scale edge weights
  cell_type_links$weight <- scaling_edge_weight(as.numeric(cell_type_links$weight))
  
  # Generate and save pairwise interaction plot
  filename <- file.path(output_dir, paste0(input_file_name, "_", interaction_type,"_interaction_TR_", tr, "_IR_", ir, "_HR_", hr, ".png"))
  png(filename, width = 3000, height = 2500, res = 300)
  
  net <- graph_from_data_frame(d=cell_type_links, vertices=cell_type_nodes, directed=F)
  V(net)$size <- as.numeric(V(net)$influence)

  if (interaction_type == "self"){
    num_of_vertices <- length(V(net)) / 2
    distinct_colors <- rainbow(num_of_vertices)
    V(net)$color <- c(distinct_colors, distinct_colors)
  } else{
    num_of_vertices <- length(V(net))
    distinct_colors <- rainbow(num_of_vertices)
    V(net)$color <- distinct_colors
  }
  
  E(net)$width <- as.numeric(E(net)$weight)
  E(net)$edge.color <- E(net)$sign
  
  plot(net, edge.color = E(net)$sign, vertex.label=NA)
  legend("topright", legend = V(net)$name[1:num_of_vertices], fill = distinct_colors, cex = 1.5, box.lwd = 2)
  
  garbage <- dev.off()
  cat(paste0(interaction_type, " interactions visualized!\n"))
}



