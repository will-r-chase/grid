library(tidyverse)
library(EnvStats)
library(ambient)
library(geos)
library(wkutils)
library(mgcv)

regon_pts <- function(cx = 0, cy = 0, r = 10, edges = 3, start_angle = 0, points = 10000) {
  a = start_angle * (pi / 180)
  n = floor(points / edges)
  
  corners <- tibble(angle = seq(0 + a, 2*pi + a, length.out = edges + 1), x = r*cos(angle) + cx, y = r*sin(angle) + cy)
  
  edges_list <- list()
  
  for(i in 1:edges) {
    edges_list[[i]] <- tibble(id = 1:n) %>%
      mutate(d = id / n, 
             x = (1 - d) * corners$x[i] + d * corners$x[i + 1],
             y = (1 - d) * corners$y[i] + d * corners$y[i + 1])
  }
  
  bind_rows(edges_list)
  
}

gen_regon <- function(cx = 0, cy = 0, r = 10, edges = 3, start_angle = 0) {
  a = start_angle * (pi / 180)
  tibble(angle = seq(0 + a, 2*pi + a, length.out = edges + 1), x = r*cos(angle) + cx, y = r*sin(angle) + cy)
}

gen_rect <- function(x = 0, y = 0, l = 10, w = 5) {
  tibble(x = c(x, x + w, x + w, x), y = c(y, y, y - l, y - l))
}

# a function that you define a number of lines in X direction with x-gap
# same for Y, with y-gap
# define how many points are in each line
# give each line a unique group
# to do later: give the points some properties, like inertia
remove_overlaps <- function(polygons) {
  poly_geos <- geos_make_polygon(
    polygons$x, polygons$y, 
    feature_id = polygons$id
  )
  
  # set up loop variables
  poly_geos_iter <- poly_geos
  
  # st_intersects will give you the same object
  poly_geos_intersections <- geos_intersects_matrix(
    poly_geos_iter, 
    poly_geos_iter
  )
  max_iter <- 1000
  
  for (i in seq_len(max_iter)) {
    n_intersections <- vapply(poly_geos_intersections, length, integer(1))
    if (all(n_intersections == 1)) {
      break;
    }
    
    poly_geos_iter <- poly_geos_iter[-which.max(n_intersections)]
    poly_geos_intersections <- geos_intersects_matrix(
      poly_geos_iter, 
      poly_geos_iter
    )
  }
  
  geos_write_wkb(poly_geos_iter) %>%
    wkb_coords() %>%
    select(id = feature_id, x, y)
}


lerp <- function(start_x, start_y, end_x, end_y, points) {
  tibble(id = 1:points) %>%
    mutate(d = id / points,
           x = (1 - d) * start_x + d * end_x,
           y = (1 - d) * start_y + d * end_y)
}

# grid_gen(h_lines = 30, h_xstart = -0.9, h_ystart = 0, 
#          h_xend = 20, h_gap = 20 / 30, h_points = 100,
#          v_lines = 30, v_xstart = 0, v_ystart = -0.9, 
#          v_xend = 20, v_yend = 20, v_gap = 20 / 30, v_points = 100) %>%

grid_gen <- function(h_lines = NULL, v_lines = NULL, h_xstart, h_xend, h_ystart, 
                      h_yend, v_xstart, v_ystart, v_xend, v_yend, 
                      h_gap, v_gap, h_points, v_points) {
  
  if(is.null(h_gap)) {
    h_gap <- h_xend / h_lines
  }
  if(is.null(v_gap)) {
    v_gap <- v_yend / v_lines
  }
  
  if(!is.null(h_lines)) {
    #make a df for the start and end points for the horizontal lines
    h_df <- tibble(lines = 1:h_lines,
                   x_start = h_xstart,
                   y_start = h_ystart + h_gap * (lines - 1),
                   x_end = h_xend,
                   y_end = y_start)
  
    #interpolate points between the horizontal endpoints
    h_interpolated <-
      h_df %>%
      group_by(lines) %>%
      group_map( ~lerp(.x$x_start, .x$y_start, .x$x_end, .x$y_end, h_points)) %>%
      bind_rows(.id = 'line') %>%
      mutate(line_direction = 'horizontal')
  }
  
  if(!is.null(v_lines)) {
    #make a df for the start and end points for the vertical lines
    v_df <- tibble(lines = 1:v_lines, 
                   x_start = v_xstart + v_gap * (lines - 1),
                   y_start = v_ystart,
                   x_end = x_start,
                   y_end = v_yend)
    
    #interpolate points between the vertical endpoints
    v_interpolated <-
      v_df %>%
      group_by(lines) %>%
      group_map( ~lerp(.x$x_start, .x$y_start, .x$x_end, .x$y_end, v_points)) %>%
      bind_rows(.id = 'line') %>%
      mutate(line_direction = 'vertical')
  }
  
  obj <- list()
  class(obj) <- "grid"
  
  #if we have both vertical and horizontal, adjust line id and bind them into 1 df
  #else just return the vertical or horizontal df
  if(!is.null(h_lines) & !is.null(v_lines)) {
    v_interpolated <-
      v_interpolated %>%
      mutate(line = as.character(as.numeric(line) + h_lines))
    
    obj$grid <- rbind(h_interpolated, v_interpolated)
  } else if (!is.null(h_lines) & is.null(v_lines)) {
    obj$grid <- h_interpolated
  } else {
    obj$grid <- v_interpolated
  }
  
  obj$bounds <- c(min(obj$grid$x), max(obj$grid$y), max(obj$grid$x), min(obj$grid$y))
  return(obj)
  
}

grid_regons <- function(grid, n = NULL, edges = NULL, cx = NULL, cy = NULL, r = NULL, color = NULL, overlap = FALSE) {
  #randomize and set up parameters 
  if(is.null(n)) {n <- sample(3:15, 1)}
  if(is.null(edges)) {
    edges <- sample(3:8, n, replace = TRUE)
  } else if(length(edges) == 1) {
    edges <- rep(edges, n)
  }
  if(is.null(cx)) {
    cx <- sample(seq(grid$bounds[1], grid$bounds[3], by = 0.1), n, replace = TRUE)
  }
  if(is.null(cy)) {
    cy <- sample(seq(grid$bounds[4], grid$bounds[2], by = 0.1), n, replace = TRUE)
  }
  if(is.null(r)) {
    r <- rpareto(n, location = sample(seq(max(grid$bounds) / 20, max(grid$bounds) / 10, by = 0.1), 1), shape = sample(3:4, 1))
  }
  if(is.null(color)) {
    colors <- c("#666666")
  } else {
    colors <- color
  }
  
  regons <- 
    tibble(edges = edges, cx = cx, cy = cy, r = r) %>%
    rowwise() %>%
    mutate(data = list(gen_regon(cx = cx, cy = cy, edges = edges, r = r, start_angle = ifelse(edges %% 2 == 0, 180 / edges, 90))))
  
  if(!overlap) {
    regons_df <- bind_rows(regons$data, .id = "id") %>%
      remove_overlaps() %>%
      group_by(id) %>%
      mutate(color = sample(colors, 1, replace = TRUE))
  } else {
    regons_df <- bind_rows(regons$data, .id = "id") %>%
      group_by(id) %>%
      mutate(color = sample(colors, 1, replace = TRUE))
  }
  
  grid$regons <- regons_df 
  return(grid)
}

grid_rects <- function(grid, name = NULL, n = NULL, x = NULL, y = NULL, l = NULL, w = NULL, color = NULL, overlap = FALSE) {
  #randomize and set up parameters 
  if(is.null(n)) {n <- sample(5:50, 1)} 
  
  if(is.null(x)) {
    x <- sample(seq(grid$bounds[1], grid$bounds[3], by = 0.1), n, replace = TRUE)
  }
  if(is.null(y)) {
    y <- sample(seq(grid$bounds[4], grid$bounds[2], by = 0.1), n, replace = TRUE)
  }
  if(is.null(l)) {
    #l <- rpareto(n, location = sample(seq(max(grid$bounds) / 30, max(grid$bounds) / 4, by = 0.1), 1), shape = sample(3:4, 1))
    l <- sample(seq(max(grid$bounds) / 30, max(grid$bounds) / 4, by = 0.1), n, replace = TRUE)
  }
  if(is.null(w)) {
    #w <- rpareto(n, location = sample(seq(max(grid$bounds) / 30, max(grid$bounds) / 4, by = 0.1), 1), shape = sample(3:4, 1))
    w <- sample(seq(max(grid$bounds) / 30, max(grid$bounds) / 4, by = 0.1), n, replace = TRUE)
  }
  if(is.null(color)) {
    colors <- c("#666666")
  } else {
    colors <- color
  }
  
  rects <- 
    tibble(x = x, y = y, l = l, w = w) %>%
    rowwise() %>%
    mutate(data = list(gen_rect(x = x, y = y, l = l, w = w)))
  
  if(!overlap) {
    rects_df <- bind_rows(rects$data, .id = "id") %>%
      remove_overlaps() %>%
      group_by(id) %>%
      mutate(color = sample(colors, 1, replace = TRUE))
  } else {
    rects_df <- bind_rows(rects$data, .id = "id") %>%
      group_by(id) %>%
      mutate(color = sample(colors, 1, replace = TRUE))
  }
  
  grid$rects <- rects_df 
  return(grid)
}

point_in_polygon <- function(polys, points) {
  pts <- points %>% select(x, y)
  
  inout_pts <- 
    polys %>%
    select(x, y) %>%
    split(., polys$id) %>%
    map(., ~in.out(as.matrix(.x), as.matrix(pts))) %>%
    map(., ~cbind(.x, pts)) %>%
    map(., ~rename(.x, "inout" = ".x")) %>%
    #map(., ~filter(.x, inout == TRUE)) %>%
    bind_rows(.id = "polygon_id") %>%
    mutate(id = rep(1:nrow(pts), max(polys$id)))
  
  in_pts <- 
    inout_pts %>%
    filter(inout == TRUE)
  
  out_pts <-
    inout_pts %>%
    filter(inout == FALSE) %>%
    distinct(id, .keep_all = TRUE) %>%
    anti_join(in_pts, by = "id")
  
  test_inout <- rbind(in_pts, out_pts) %>%
    arrange(id)
}

assign_inertia <- function(grid, in_intertia = 1, out_inertia = 0) {
  points <- grid$grid
  rects <- grid$rects
  regons <- grid$regons
  if(is.null(rects)) {
    polys <- regons %>% ungroup()
  } else {
    polys <- rects %>% ungroup()
  }
  
  pip <- point_in_polygon(polys, points)
  
  points <-
    points %>%
    mutate(inout = pip$inout,
           intertia = ifelse(inout, in_intertia, out_inertia))
  
  grid$grid <- points
  return(grid)
}



test <- 
  grid_gen(h_lines = 30, h_xstart = -0.9, h_ystart = 0, 
                 h_xend = 20, h_gap = 20 / 30, h_points = 100,
                 v_lines = 30, v_xstart = 0, v_ystart = -0.9, 
                 v_xend = 20, v_yend = 20, v_gap = 20 / 30, v_points = 100) %>%
  #grid_regons(n = 15, color = neon2) %>%
  grid_rects(n = 20, color = neon1) %>%
  assign_inertia()



ggplot() +
  #geom_path(data = test$grid, aes(x = x, y = y, group = line), size = 0.4) +
  #geom_polygon(data = test$regons, aes(x = x, y = y, group = id, fill = color), alpha = 0.8) +
  geom_polygon(data = test$rects, aes(x = x, y = y, group = id, fill = color), alpha = 0.7) +
  geom_point(data = test$grid, aes(x = x, y = y, color = inout), size = 0.4) +
  scale_color_manual(values = c("grey", "red")) +
  scale_fill_identity() +
  theme_void()
  
ggsave("inout_test.jpg", width = 8, height = 8)

##function to add polygons to a grid scene
##maybe later add gradient or other random thing
##should vector fields be added here too?


grid_vectors <- function(grid, ) {
  
}

##function to add inertia to points based on some rule
##random, point in polygon, point based on gradient, point based on ID or index
grid_inertia <- function() {}

##function to project points based on inertia and given vector field
##add option to allow different vector fields for different points
grid_project <- function() {}


