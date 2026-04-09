simulation <- function(input_file, output_dir, tr, ir, hr) {
  tr <<- tr
  cat("Simulating synthetic image...\n")
  
  cell_data <- read.csv(input_file)
  input_file_name <- basename(tools::file_path_sans_ext(input_file))
  filename <- file.path(output_dir, paste0(input_file_name, "_model_TR_", tr, "_IR_", ir, "_HR_", hr, ".Rda"))
  load(file = filename)

  W <<- owin(c(min(cell_data[[1]]), max(cell_data[[1]])),
             c(min(cell_data[[2]]), max(cell_data[[2]])))

  cell_type_list <<- sort(unique(cell_data[[3]]))
  cell_type_num <<- length(cell_type_list)

  G <- list()
  coef_names <- names(coef)

  for (i in seq_len(cell_type_num)) {
    coef_current <- c()
    for (r in seq(ir, tr, ir)) {
      for (j in seq_len(cell_type_num)) {
        mark1 <- cell_type_list[i]
        mark2 <- cell_type_list[j]
        interaction_name <- if (i <= j) {
          paste('X', mark1, 'x', 'X', mark2, 'x', r, sep = '')
        } else {
          paste('X', mark2, 'x', 'X', mark1, 'x', r, sep = '')
        }
        coef_current <- c(coef_current, coef[grep(interaction_name, coef_names, value = FALSE)])
      }
    }
    G[[i]] <- coef_current
  }
  G <<- G

  B <- c()
  for (i in seq_len(cell_type_num)) {
    B_current <- if (i == 1) {
      exp(coef[i])
    } else {
      exp(coef[i] + coef[1])
    }
    B <- c(B, B_current)
  }
  B <<- B

  resample_percent <- 100
  model_mark_prob <- (dplyr::count(cell_data, marks)$n) / nrow(cell_data)

  intensity <- nrow(cell_data) / area(W)
  cell_pattern_random <- rpoispp(intensity, win = W)

  pattern_exist <- data.frame(cell_pattern_random)

  pattern_exist[, 3] <- as.factor(sample(cell_type_list, nrow(pattern_exist), replace = TRUE, model_mark_prob))
  cell_pattern_num <- nrow(pattern_exist)
  resample_num <- ceiling(cell_pattern_num * resample_percent / 100)

  if (resample_num != 0) {
    for (i in seq_len(resample_num)) {
      point_current_idx <- sample(1:cell_pattern_num, 1)
      point_current <- pattern_exist[point_current_idx, ]
      counts_current <- t(apply(point_current[,1:2], 1, calc_dist, pattern_exist))
      mark_prob_current <- apply(counts_current, 1, calc_mark_prob)
      point_current[, 3] <- mark_prob_current
      pattern_exist[point_current_idx, ] <- point_current
    }
  }

  marked_simulated_pattern <- ppp(c(pattern_exist[,1]), c(pattern_exist[,2]), window = W, marks = factor(pattern_exist[,3], levels = cell_type_list))
  filename <- file.path(output_dir, paste0(input_file_name, "_synthetic_pattern_TR_", tr, "_IR_", ir, "_HR_", hr, ".Rda"))
  save(marked_simulated_pattern, file = filename)
  # load(filename)
  synthetic_image_vis(data.frame(marked_simulated_pattern), input_file_name, output_dir, tr, ir, hr)
  cat("Simulation completed and synthetic image saved!\n")
}

truncate_polygon <- function(polygon, center, cutoff_distance) {
  set.seed(3)
  
  if (is.list(polygon)) {
    polygon <- cbind(polygon$x, polygon$y)
  }
  
  if (!all(polygon[1, ] == polygon[nrow(polygon), ])) {
    polygon <- rbind(polygon, polygon[1, ])  
  }
  
  center <- as.numeric(unname(center))
  pol <- st_polygon(list(polygon)) %>% st_sfc() %>% st_sf()
  center_sf <- st_point(center) %>% st_sfc() %>% st_sf()
  buffered_point <- st_buffer(center_sf, dist = cutoff_distance)
  intersection <- st_intersection(pol, buffered_point)
  
  if (!is.null(intersection) && nrow(intersection) > 0) {
    return(intersection)  # Returns an sf object
  } else {
    return(NULL)
  }
}




synthetic_image_vis <- function(pattern_df, input_file_name, output_dir, tr, ir, hr){
  
  colors <- rainbow(cell_type_num)
  
  plot_base <- ggplot() + 
    theme_void() +
    theme(legend.position = "none", plot.background = element_rect(fill = "black"))+ scale_fill_identity()
  
  buffer_distance <- 20
  points <- as.matrix(pattern_df[, 1:2])
  marks <- as.integer(pattern_df$mark)
  
  vor <- deldir(points)
  tile_list <- tile.list(vor)
  plot_data_list <- list() 
  suppressWarnings({
    for (i in seq_along(tile_list)) {
      tile <- tile_list[[i]]
      polygon <- tile[c("x", "y")]
      
      polygon_matrix <- cbind(polygon$x, polygon$y)
      
      if (!all(polygon_matrix[1, ] == polygon_matrix[nrow(polygon_matrix), ])) {
        polygon_matrix <- rbind(polygon_matrix, polygon_matrix[1, ])
      }
      center_coords <- as.numeric(unname(points[i,]))  
      truncated <- truncate_polygon(polygon_matrix, center_coords, buffer_distance)
      if (!is.null(truncated)) {
        temp_data <- st_coordinates(truncated)  
        temp_data <- data.frame(long = temp_data[,1], lat = temp_data[,2])
        temp_data$color <- as.character(colors[marks[i] + 1])
        temp_data$group <- i  
        plot_data_list[[i]] <- temp_data
      }
    }
  })
  
  plot_data <- do.call(rbind, plot_data_list)
  plot_base <- plot_base + 
    geom_polygon(data = plot_data, aes(x = long, y = lat, fill = color, group = group), color = "white", linewidth = 0.3) 
  filename <- file.path(output_dir, paste0(input_file_name, "_synthetic_image_TR_", tr, "_IR_", ir, "_HR_", hr, ".png"))
  ggsave(filename = filename, dpi = 500, bg = "black", width=10, height=10)
}


get_random_sampled_marks <- function(types, num, prob) {
  sampled_marks <- sample(types, num, TRUE, prob)
  return(sampled_marks)
}

calc_dist <- function(new_points, existing_points) {
  new_points <- data.frame(t(new_points))
  x <- new_points[1, 1]
  y <- new_points[1, 2]
  xmin <- max(0, x-tr)
  ymin <- max(0, y-tr)
  xmax <- min(W$xrange[2], x+tr)
  ymax <- min(W$yrange[2], y+tr)
  
  existing_points_within_range <- existing_points[which(existing_points[,1] > xmin & existing_points[,1] < xmax & existing_points[,2] > ymin & existing_points[,2] < ymax),]
  eu_dist <- proxy::dist(new_points, existing_points_within_range[,1:2])
  
  existing_points_within_range <- cbind(existing_points_within_range, c(eu_dist))
  existing_points_within_range <- existing_points_within_range[which(existing_points_within_range[,4] < 500 & existing_points_within_range[,4] > 0), 3:4]
  
  mark_count <- matrix(0, nrow = 5, ncol = 5)
  for (i in seq_len(nrow(existing_points_within_range))) {
    mark_count_row <- as.integer(existing_points_within_range[i, 1])
    mark_count_col <- existing_points_within_range[i, 2] %/% 100 + 1
    mark_count[mark_count_row, mark_count_col] <- mark_count[mark_count_row, mark_count_col] + 1
  }
  return(mark_count)
}

calc_mark_prob <- function(counts) {
  m_prob <- c()
  for (m in cell_type_list) {
    m_int <- as.integer(m) + 1
    G_type <- G[[m_int]]
    # print(counts)
    
    counts_type <- counts
    interact_term <- exp(sum(G_type * counts_type))
    intensity_term <- B[m_int]
    m_prob <- c(m_prob, unname(intensity_term * interact_term))
  }
  
  m_prob <- m_prob / sum(m_prob)
  return(factor(sample(cell_type_list, 1, prob = m_prob), levels = cell_type_list))
}