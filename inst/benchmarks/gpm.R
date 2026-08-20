# Tiny matrix-free GPM harness. Not part of the package namespace.
# Run from the package root:
#   Rscript inst/benchmarks/gpm.R

if (!requireNamespace("gprocrustes", quietly = TRUE)) {
  pkgload::load_all()
} else {
  library(gprocrustes)
}

rot2 <- function(theta) {
  matrix(c(cos(theta), sin(theta), -sin(theta), cos(theta)), 2, 2)
}

make_problem <- function(n, k, d = 2) {
  X <- scale(matrix(rnorm(n * d), n, d), scale = FALSE)
  views <- vector("list", k)
  names(views) <- paste0("V", seq_len(k))
  views[[1L]] <- X
  for (i in seq.int(2L, k)) {
    Q <- qr.Q(qr(matrix(rnorm(d * d), d, d)))
    views[[i]] <- X %*% Q
  }
  views
}

bench_one <- function(n, k) {
  views <- make_problem(n, k)
  t_gpm <- system.time({
    gpm <- gpa(views, transform = proc_orthogonal("O"),
               control = gpa_control(backend = "matrix_free", init = "spectral"))
  })[["elapsed"]]
  t_gower <- system.time({
    gw <- gpa(views, transform = proc_orthogonal("O"), solver = "gower_bcd",
              control = gpa_control(accelerate = FALSE))
  })[["elapsed"]]
  data.frame(
    n = n, k = k,
    gpm_sec = t_gpm, gower_sec = t_gower,
    gpm_obj = gpm$objective, gower_obj = gw$objective,
    gpm_status = gpm$optimality_status,
    backend = gpm$backend
  )
}

grid <- expand.grid(n = c(40, 80), k = c(4, 8), KEEP.OUT.ATTRS = FALSE)
set.seed(1)
out <- do.call(rbind, lapply(seq_len(nrow(grid)), function(i) {
  bench_one(grid$n[[i]], grid$k[[i]])
}))
print(out)
