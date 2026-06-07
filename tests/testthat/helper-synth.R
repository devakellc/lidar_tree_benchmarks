# Synthetic fixtures shared across tests. No real NEON data is ever required.
suppressMessages({ library(lidR); library(data.table) })

# A labelled point table: two well-separated "trees", each a small vertical
# cluster. Columns X,Y,Z (metres) + a known instance id. Tree A apex (10,10,18),
# tree B apex (40,12,12).
synth_labelled_points <- function() {
  a <- data.table(X = c(10, 10.2, 9.8, 10.1), Y = c(10, 10.1, 9.9, 10.0),
                  Z = c(18, 14, 9, 4),  id = 1L)
  b <- data.table(X = c(40, 40.1, 39.9),       Y = c(12, 12.1, 11.9),
                  Z = c(12, 8, 3),      id = 2L)
  noise <- data.table(X = 25, Y = 25, Z = 1.0,  id = NA_integer_)  # unassigned
  rbind(a, b, noise)
}

# A minimal NORMALIZED LAS (ground at 0) with two clusters, for arm smoke tests.
synth_las_normalized <- function() {
  set.seed(1)
  mk <- function(cx, cy, top, n = 60) {
    data.frame(X = rnorm(n, cx, 0.7), Y = rnorm(n, cy, 0.7),
               Z = runif(n, 0, top))
  }
  df <- rbind(mk(10, 10, 18), mk(40, 12, 12))
  df$Z[1]  <- 18; df$Y[1]  <- 10; df$X[1]  <- 10   # guarantee a tree-A apex
  df$Z[61] <- 12; df$Y[61] <- 12; df$X[61] <- 40   # guarantee a tree-B apex
  las <- LAS(df)
  st_crs(las) <- 32611L
  las
}

