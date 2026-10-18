#' Observed-data optimization on the product of rotation/scale manifolds.
#'
#' The public `mm` engine uses local Cayley charts and limited-memory BFGS
#' with an analytic gradient, without a dense Hessian over all landmarks.
#' Every accepted step is checked against the true landmark-vector loss.
#' Numerical stagnation is never promoted to first-order stationarity.
#' @noRd
.mm_gpa <- function(data, spec, gauge, metric, alphas, control, anchor, loss) {
  n <- length(data$global_ids)
  d <- ncol(data$views[[1L]])
  if (any(vapply(data$views, ncol, integer(1)) != d)) {
    .gproc_stop("dimension_mismatch", "Observed-data fitting requires equal coordinate dimensions.")
  }
  tol <- control$tolerance
  if (length(tol) != 1L || !is.finite(tol) || tol <= 0 ||
      length(control$max_iterations) != 1L || is.na(control$max_iterations) || control$max_iterations < 1) {
    .gproc_stop("invalid_problem", "Positive tolerance and iteration count are required.")
  }
  mode <- .gower_scale_mode(spec, gauge)
  prepared <- .gower_prepare_views(data, spec, mode$mode)
  work <- data
  work$views <- lapply(prepared, `[[`, "X")
  chan <- .gproc_cell_channel(work, metric)
  covariance <- .gproc_fitting_covariance(metric)
  A <- covariance$coordinate
  if (!is.null(A) && (!identical(dim(A), c(d, d)) || any(!is.finite(A)) ||
      max(abs(A - t(A))) > 1e-10 * max(1, max(abs(A))) ||
      min(.gproc_eigen_sym(A)$values) < -1e-12 * max(1, max(abs(A))))) {
    .gproc_stop("invalid_problem", "Coordinate covariance must be finite, symmetric and positive semidefinite.")
  }
  Q <- .gproc_coordinate_precision(covariance)
  if (is.null(Q)) Q <- diag(d)
  if (!identical(dim(Q), c(d, d)) || any(!is.finite(Q)) ||
      max(abs(Q - t(Q))) > 1e-10 * max(1, max(abs(Q)))) {
    .gproc_stop("invalid_problem", "Coordinate precision must be a finite symmetric d by d matrix.")
  }
  if (min(.gproc_eigen_sym(Q)$values) < -1e-12 * max(1, max(abs(Q)))) {
    .gproc_stop("invalid_problem", "Coordinate precision must be positive semidefinite.")
  }
  model <- list(data = work, Q = Q, weights = chan$weight, cells = chan$mask,
                alpha = alphas, loss = loss, spec = spec, mode = mode$mode,
                anchor = anchor, landmark = list(), D = numeric())
  model$translate <- .gproc_eliminates_translation(spec) && mode$mode != "preshape"
  model$scale <- identical(spec$scaling, "isotropic") && mode$mode != "preshape"
  for (nm in names(work$views)) {
    X <- work$views[[nm]]
    w <- .gproc_local_landmark_weights(metric$landmark, work$row_map[[nm]], n, nm)
    model$landmark[[nm]] <- .gproc_row_weights(w, nrow(X))
    keep <- work$observed[[nm]]
    Xc <- X[keep, , drop = FALSE]
    if (model$translate && nrow(Xc)) Xc <- sweep(Xc, 2, colMeans(Xc), "-")
    model$D[nm] <- alphas[[nm]] * sum(Xc^2)
  }
  model$energy <- sum(model$D)
  if (model$scale && mode$mode == "gower" &&
      (model$energy <= 0 || any(model$D[alphas > 0] <= 0))) {
    .gproc_stop("invalid_problem", "Collective scale needs positive centered size in each weighted view.")
  }
  inner <- control
  inner$max_iterations <- min(control$max_iterations, 40L)
  inner$nstart <- 1L
  raw <- .gower_gpa(data, spec, gauge, metric, alphas, inner, anchor)
  state <- list(transforms = raw$work_transforms, M = raw$consensus)
  if (!model$translate) state$transforms <- lapply(state$transforms, function(tr) {tr$t[] <- 0; tr})
  if (!model$scale) state$transforms <- lapply(state$transforms, function(tr) {tr$s <- 1; tr})
  if (!is.null(anchor)) state$transforms[[anchor]] <- proc_identity(spec, d)
  # Decode once to establish the declared scale constraints before any objective.
  chart <- .mm_chart(state, model)
  state <- .mm_decode(chart$par, chart, model)$state
  current <- .mm_evaluate(state, model)
  history <- list()
  status <- "maximum_iterations"
  t0 <- proc.time()[["elapsed"]]
  for (it in seq_len(control$max_iterations)) {
    old <- current$value
    chart <- .mm_chart(state, model)
    cache <- NULL
    evaluate <- function(par) {
      if (!is.null(cache) && identical(par, cache$par)) return(cache)
      decoded <- .mm_decode(par, chart, model)
      val <- .mm_evaluate(decoded$state, model)
      cache <<- list(par = par, value = val$value, state = decoded$state,
                     gradient = .mm_chart_gradient(val, decoded, chart, model))
      cache
    }
    opt <- stats::optim(chart$par, function(p) evaluate(p)$value,
                        function(p) evaluate(p)$gradient, method = "L-BFGS-B",
                        control = list(maxit = 50L, factr = 1, pgtol = min(1e-8, sqrt(tol) / 10)))
    proposal <- evaluate(opt$par)
    accepted <- is.finite(proposal$value) && proposal$value <= old
    if (accepted) state <- proposal$state
    current <- .mm_evaluate(state, model)
    stat <- .mm_stationarity(current, state, model)
    rel <- abs(old - current$value) / max(1, abs(old))
    orth <- .gower_orthogonality(state$transforms)
    history[[it]] <- data.frame(iteration = it, objective = current$value,
      relative_change = rel, accepted_mm = accepted, stationarity = stat,
      orthogonality_residual = orth, elapsed_time = proc.time()[["elapsed"]] - t0)
    if (rel <= tol && stat <= sqrt(tol) && orth <= 1e-8) {
      status <- "converged"
      break
    }
    if (!accepted) break
  }
  raw$work_transforms <- state$transforms
  raw$transforms <- .gower_original_transforms(state$transforms, prepared, spec, d)
  raw$consensus <- state$M
  raw$aligned <- .gproc_apply_all(work, state$transforms, n, d)
  raw$objective <- current$value
  raw$history <- do.call(rbind, history)
  raw$stationarity <- .mm_stationarity(current, state, model)
  raw$numerical_status <- status
  raw$optimality_status <- if (status == "converged") "first_order_stationary" else "not_converged"
  raw$cell_masks <- chan$mask
  raw$cell_weights <- chan$weight
  raw$landmark_weights <- lapply(names(work$views), function(nm) {
    current$q[[nm]] * model$landmark[[nm]]
  })
  names(raw$landmark_weights) <- names(work$views)
  raw
}

#' A local chart, retaining the determinant component of each fitted rotation.
#' @noRd
.mm_chart <- function(state, model) {
  d <- ncol(state$M)
  pairs <- which(upper.tri(matrix(0, d, d)), arr.ind = TRUE)
  par <- numeric()
  fields <- list()
  append <- function(value) {
    idx <- length(par) + seq_along(value)
    par <<- c(par, value)
    idx
  }
  for (nm in names(state$transforms)) {
    if (identical(nm, model$anchor) || model$alpha[[nm]] == 0) next
    tr <- state$transforms[[nm]]
    fields[[nm]] <- list(rotation = append(rep(0, nrow(pairs))),
      scale = if (model$scale) append(log(max(tr$s, 1e-8))) else integer(),
      translation = if (model$translate) append(tr$t) else integer())
  }
  m_idx <- append(as.vector(state$M))
  list(par = par, fields = fields, M = m_idx, base = state, pairs = pairs)
}

#' @noRd
.mm_decode <- function(par, chart, model) {
  state <- chart$base
  d <- ncol(state$M)
  inverse <- list()
  for (nm in names(chart$fields)) {
    fld <- chart$fields[[nm]]
    A <- matrix(0, d, d)
    A[chart$pairs] <- par[fld$rotation]
    A <- A - t(A)
    Li <- solve(diag(d) - A / 2)
    state$transforms[[nm]]$R <- chart$base$transforms[[nm]]$R %*% Li %*% (diag(d) + A / 2)
    inverse[[nm]] <- Li
    if (model$translate) state$transforms[[nm]]$t <- par[fld$translation]
    if (model$scale) state$transforms[[nm]]$s <- exp(par[fld$scale])
  }
  if (model$scale && model$mode == "gower" && length(chart$fields)) {
    free <- names(chart$fields)
    fixed <- setdiff(names(state$transforms), free)
    remaining <- model$energy - sum(vapply(fixed, function(nm) {
      model$D[[nm]] * state$transforms[[nm]]$s^2
    }, numeric(1)))
    eta <- vapply(free, function(nm) par[chart$fields[[nm]]$scale], numeric(1))
    s <- exp(eta - max(eta))
    s <- s * sqrt(remaining / sum(model$D[free] * s^2))
    for (nm in free) state$transforms[[nm]]$s <- s[[nm]]
  }
  M <- matrix(par[chart$M], nrow(state$M), d)
  mnorm <- 1
  if (model$mode == "fixed_consensus") {
    if (model$translate) M <- sweep(M, 2, colMeans(M), "-")
    mnorm <- sqrt(sum(M^2))
    if (mnorm <= .Machine$double.eps) .gproc_stop("invalid_problem", "A fixed-scale consensus needs positive size.")
    M <- M / mnorm
  }
  state$M <- M
  list(state = state, inverse = inverse, mnorm = mnorm)
}

#' True vector loss and its Euclidean scores. Cell weights enter as
#' diag(sqrt(w)) Q diag(sqrt(w)); landmark reliability multiplies rho.
#' @noRd
.mm_evaluate <- function(state, model) {
  value <- 0
  out <- list(R = list(), s = numeric(), t = list(), M = state$M * 0, q = list())
  for (nm in names(state$transforms)) {
    tr <- state$transforms[[nm]]
    X <- as.matrix(model$data$views[[nm]])
    map <- model$data$row_map[[nm]]
    XR <- X %*% tr$R
    E <- sweep(tr$s * XR, 2, tr$t, "+") - state$M[map, , drop = FALSE]
    sw <- sqrt(model$weights[[nm]][map, , drop = FALSE])
    Z <- E * sw
    ZQ <- Z %*% model$Q
    r <- sqrt(pmax(rowSums(Z * ZQ), 0))
    q <- .gproc_irls_weight(r, model$loss)
    l2 <- identical(model$loss$family, "squared_l2")
    a <- model$alpha[[nm]] * model$landmark[[nm]]
    value <- value + sum(a * if (l2) r^2 else .gproc_rho(r, model$loss))
    G <- ((if (l2) 2 else 1) * a * q) * (ZQ * sw)
    out$R[[nm]] <- tr$s * crossprod(X, G)
    out$s[nm] <- sum(XR * G)
    out$t[[nm]] <- colSums(G)
    out$M[map, ] <- out$M[map, , drop = FALSE] - G
    out$q[[nm]] <- q
  }
  out$value <- value
  out
}

#' @noRd
.mm_chart_gradient <- function(grad, decoded, chart, model) {
  state <- decoded$state
  ans <- numeric(length(chart$par))
  free <- names(chart$fields)
  sg <- grad$s[free] * vapply(state$transforms[free], `[[`, numeric(1), "s")
  if (model$scale && model$mode == "gower" && length(free)) {
    energy <- model$D[free] * vapply(state$transforms[free], function(tr) tr$s^2, numeric(1))
    sg <- sg - energy / sum(energy) * sum(sg)
  }
  for (nm in free) {
    fld <- chart$fields[[nm]]
    Li <- decoded$inverse[[nm]]
    H <- t(Li) %*% crossprod(chart$base$transforms[[nm]]$R, grad$R[[nm]]) %*% t(Li)
    ans[fld$rotation] <- (H - t(H))[chart$pairs]
    if (model$scale) ans[fld$scale] <- sg[[nm]]
    if (model$translate) ans[fld$translation] <- grad$t[[nm]]
  }
  GM <- grad$M
  if (model$mode == "fixed_consensus") {
    if (model$translate) GM <- sweep(GM, 2, colMeans(GM), "-")
    GM <- (GM - sum(GM * state$M) * state$M) / decoded$mnorm
  }
  ans[chart$M] <- as.vector(GM)
  ans
}

#' First-order residual on all free rotation, translation, scale and consensus
#' blocks, using the true loss rather than a surrogate's stationarity.
#' @noRd
.mm_stationarity <- function(grad, state, model) {
  residual <- 0
  free <- names(state$transforms)[names(state$transforms) != model$anchor]
  if (is.null(model$anchor)) free <- names(state$transforms)
  free <- free[model$alpha[free] > 0]
  for (nm in free) {
    H <- crossprod(state$transforms[[nm]]$R, grad$R[[nm]])
    residual <- max(residual, sqrt(sum(((H - t(H)) / 2)^2)))
    if (model$translate) residual <- max(residual, sqrt(sum(grad$t[[nm]]^2)))
  }
  if (model$scale && length(free)) {
    gs <- grad$s[free]
    if (model$mode == "gower") {
      q <- sqrt(model$D[free]) * vapply(state$transforms[free], `[[`, numeric(1), "s")
      gs <- gs / sqrt(model$D[free])
      gs <- gs - sum(gs * q) / sum(q^2) * q
    }
    residual <- max(residual, sqrt(sum(gs^2)))
  }
  GM <- grad$M
  if (model$mode == "fixed_consensus") {
    if (model$translate) GM <- sweep(GM, 2, colMeans(GM), "-")
    GM <- GM - sum(GM * state$M) * state$M
  }
  max(residual, sqrt(sum(GM^2)))
}
