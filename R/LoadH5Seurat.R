#' @include zzz.R
#' @include h5Seurat.R
#' @include GetObject.R
#' @include AssembleObject.R
#'
NULL

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Generics
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Load a saved \code{Seurat} object from an h5Seurat file
#'
#' @param file,x Name of h5Seurat or connected h5Seurat file to load
#' @param assays One of:
#' \itemize{
#'  \item A character vector with names of assays
#'  \item A character vector with one or more of \code{counts}, \code{data},
#'  \code{scale.data} describing which slots of \strong{all assays} to load
#'  \item A named list where each entry is either the name of an assay or a vector
#'  describing which slots (described above) to take from which assay
#'  \item \code{NULL} for all assays
#' }
#' @param reductions One of:
#' \itemize{
#'  \item A character vector with names of reductions
#'  \item \code{NULL} for all reductions
#'  \item \code{NA} for \link[Seurat:IsGlobal]{global} reductions
#'  \item \code{FALSE} for no reductions
#' }
#' \strong{Note}: Only reductions associated with an assay loaded in
#' \code{assays} or marked as \link[Seurat:IsGlobal]{global} will be loaded
#' @param graphs One of:
#' \itemize{
#'  \item A character vector with names of graphs
#'  \item \code{NULL} for all graphs
#'  \item \code{FALSE} for no graphs
#' }
#' \strong{Note}: Only graphs associated with an assay loaded in \code{assays}
#' will be loaded
#' @param neighbors One of:
#' \itemize{
#'  \item A character vector with the names of neighbors
#'  \item \code{NULL} for all neighbors
#'  \item \code{FALSE} for no neighbors
#' }
#' @param images One of:
#' \itemize{
#'  \item A character vector with names of images
#'  \item \code{NULL} for all images
#'  \item \code{NA} for \link[Seurat:IsGlobal]{global} images
#'  \item \code{FALSE} for no images
#' }
#' @param meta.data Load object metadata
#' @param commands Load command information \cr
#' \strong{Note}: only commands associated with an assay loaded in
#' \code{assays} will be loaded
#' @param misc Load miscellaneous data
#' @param tools Load tool-specific information
#' @param verbose Show progress updates
#' @param ... Arguments passed to other methods
#'
#' @return A \code{Seurat} object with the data requested
#'
#' @export
#'
LoadH5Seurat <- function(file, ...) {
  UseMethod(generic = 'LoadH5Seurat', object = file)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Methods
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @rdname LoadH5Seurat
#' @method LoadH5Seurat character
#' @export
#'
LoadH5Seurat.character <- function(
  file,
  assays = NULL,
  reductions = NULL,
  graphs = NULL,
  neighbors = NULL,
  images = NULL,
  meta.data = TRUE,
  commands = TRUE,
  misc = is.null(x = assays),
  tools = is.null(x = assays),
  verbose = TRUE,
  ...
) {
  hfile <- h5Seurat$new(filename = file, mode = 'r')
  on.exit(expr = hfile$close_all())
  return(LoadH5Seurat(
    file = hfile,
    assays = assays,
    reductions = reductions,
    graphs = graphs,
    neighbors = neighbors,
    images = images,
    meta.data = meta.data,
    commands = commands,
    misc = misc,
    tools = tools,
    verbose = verbose,
    ...
  ))
}

#' @rdname LoadH5Seurat
#' @method LoadH5Seurat H5File
#' @export
#'
LoadH5Seurat.H5File <- function(
  file,
  assays = NULL,
  reductions = NULL,
  graphs = NULL,
  neighbors = NULL,
  images = NULL,
  meta.data = TRUE,
  commands = TRUE,
  misc = is.null(x = assays),
  tools = is.null(x = assays),
  verbose = TRUE,
  ...
) {
  return(LoadH5Seurat(
    file = as.h5Seurat(x = file),
    assays = assays,
    reductions = reductions,
    graphs = graphs,
    neighbors = neighbors,
    images = images,
    meta.data = meta.data,
    commands = commands,
    misc = misc,
    tools = tools,
    verbose = verbose,
    ...
  ))
}

#' @importFrom Seurat as.Seurat
#'
#' @rdname LoadH5Seurat
#' @method LoadH5Seurat h5Seurat
#' @export
#'
LoadH5Seurat.h5Seurat <- function(
  file,
  assays = NULL,
  reductions = NULL,
  graphs = NULL,
  neighbors = NULL,
  images = NULL,
  meta.data = TRUE,
  commands = TRUE,
  misc = is.null(x = assays),
  tools = is.null(x = assays),
  verbose = TRUE,
  ...
) {
  return(as.Seurat(
    x = file,
    assays = assays,
    reductions = reductions,
    graphs = graphs,
    neighbors = neighbors,
    images = images,
    meta.data = meta.data,
    commands = commands,
    misc = misc,
    tools = tools,
    verbose = verbose,
    ...
  ))
}

#' @importClassesFrom Seurat Seurat
#' @importFrom methods slot<-
#' @importFrom Seurat as.Seurat DefaultAssay Cells
#' Idents<- Idents Project<- Project
#' AddMetaData
#'
#' @aliases as.Seurat
#'
#' @rdname LoadH5Seurat
#' @method as.Seurat h5Seurat
#' @export
#'
as.Seurat.h5Seurat <- function(
  x,
  assays = NULL,
  reductions = NULL,
  graphs = NULL,
  neighbors = NULL,
  images = NULL,
  meta.data = TRUE,
  commands = TRUE,
  misc = TRUE,
  tools = TRUE,
  verbose = TRUE,
  ...
) {
  index <- x$index()
  obj.all <- all(vapply(
    X = c(assays, reductions, graphs),
    FUN = is.null,
    FUN.VALUE = logical(length = 1L)
  ))
  # Load Assays
  assays <- GetAssays(assays = assays, index = index)
  if (!DefaultAssay(object = index) %in% names(x = assays)) {
    active.assay <- names(x = assays)[1]
    warning(
      "Default assay not requested, using ",
      active.assay,
      " instead",
      call. = FALSE,
      immediate. = TRUE
    )
  } else {
    active.assay <- DefaultAssay(object = index)
  }
  assay.objects <- vector(mode = 'list', length = length(x = assays))
  names(x = assay.objects) <- names(x = assays)
  for (assay in names(x = assays)) {
    assay.objects[[assay]] <- AssembleAssay(
      assay = assay,
      file = x,
      slots = assays[[assay]],
      verbose = verbose
    )
  }
  default.assay <- list(assay.objects[[active.assay]])
  names(x = default.assay) <- active.assay
  object <- new(
    Class = 'Seurat',
    assays = default.assay,
    active.assay = active.assay,
    meta.data = data.frame(row.names = Cells(x = x)),
    version = package_version(x = x$version())
  )
  for (assay in names(x = assay.objects)) {
    if (assay != active.assay) {
      object[[assay]] <- assay.objects[[assay]]
    }
  }
  # Load DimReducs
  reductions <- GetDimReducs(
    reductions = reductions,
    index = index,
    assays = assays
  )
  for (reduc in reductions) {
    if (verbose) {
      message("Adding reduction ", reduc)
    }
    reduction <- AssembleDimReduc(
      reduction = reduc,
      file = x,
      verbose = verbose
    )
    if (isTRUE(x = getOption(x = 'SeuratDisk.dimreducs.allglobal', default = FALSE))) {
      slot(object = reduction, name = 'global') <- TRUE
    }
    object[[reduc]] <- reduction
  }
  # Load Graphs
  graphs <- GetGraphs(graphs = graphs, index = index, assays = assays)
  for (graph in graphs) {
    if (verbose) {
      message("Adding graph ", graph)
    }
    object[[graph]] <- AssembleGraph(graph = graph, file = x, verbose = verbose)
  }
  # Load Neighbors
  neighbors <- GetNeighbors(neighbors = neighbors, index = index)
  for (neighbor in neighbors) {
    if (verbose) {
      message("Adding neighbors ", neighbor)
    }
    object[[neighbor]] <- AssembleNeighbor(
      neighbor = neighbor,
      file = x,
      verbose = verbose
      )
  }
  # Load SpatialImages
  if (packageVersion(pkg = 'Seurat') >= numeric_version(x = spatial.version)) {
    images <- GetImages(images = images, index = index, assays = assays)
    for (image in images) {
      if (verbose) {
        message("Adding image ", image)
      }
      object[[image]] <- AssembleImage(
        image = image,
        file = x,
        verbose = verbose
      )
    }
  }
  # Load SeuratCommands
  if (commands) {
    if (verbose) {
      message("Adding command information")
    }
    cmds <- GetCommands(index = index, assays = assays)
    cmdlogs <- vector(mode = 'list', length = length(x = cmds))
    names(x = cmdlogs) <- cmds
    for (cmd in cmds) {
      cmdlogs[[cmd]] <- AssembleSeuratCommand(
        cmd = cmd,
        file = x,
        verbose = verbose
      )
    }
    slot(object = object, name = 'commands') <- cmdlogs
  }
  # Load meta.data
  if (meta.data) {
    if (verbose) {
      message("Adding cell-level metadata")
    }
    md <- as.data.frame(x = x[['meta.data']], row.names = Cells(x = x))
    if (ncol(x = md)) {
      object <- AddMetaData(object = object, metadata = md)
    }
  }
  # Set cell identities and object project
  Idents(object = object) <- Idents(object = x)
  Project(object = object) <- Project(object = x)
  # Load misc
  if (misc) {
    if (verbose) {
      message("Adding miscellaneous information")
    }
    slot(object = object, name = 'misc') <- as.list(x = x[['misc']])
  }
  # Load tools
  if (tools) {
    if (verbose) {
      message("Adding tool-specific results")
    }
    slot(object = object, name = 'tools') <- as.list(x = x[['tools']])
  }
  # Load no.assay information
  if (obj.all && !is.null(x = index$no.assay)) {
    if (verbose) {
      message("Adding data that was not associated with an assay")
    }
    for (graph in index$no.assay$graphs) {
      object[[graph]] <- AssembleGraph(
        graph = graph,
        file = x,
        verbose = verbose
      )
    }
    for (cmd in index$no.assay$commands) {
      object[[cmd]] <- AssembleSeuratCommand(
        cmd = cmd,
        file = x,
        verbose = verbose
      )
    }
  }
  return(object)
}

# Too slow
.dgCMatrixGroupSubsetV1 <- function(x, i, j) {
  if (!x$exists(name = 'data') || !x$exists(name = 'indices')
        || !x$exists(name = 'indptr')) {
    stop("x is not a HDF5 group of dgCMatrix", call. = FALSE)
  }

  dims <- h5attr(x, "dims")
  data <- x[["data"]]
  indices <- x[["indices"]]
  indptr <- x[["indptr"]]

  if (missing(i))
    i <- seq_len(dims[1])
  if (missing(j))
    j <- seq_len(dims[2])

  coldata <- lapply(
    X = j,
    FUN = function(col) {
      # Get non-zero element number of current column
      col.n <- indptr[col + 1] - indptr[col]
      # It's the pos of start and end element of current column.
      # Conver coords to 1-based.
      start <- indptr[col] + 1
      end <- start + col.n - 1
      # Get pos of non-zero element, that match supplied rows.
      col.pos <- seq(start, end)
      data.pos <- match(i - 1, indices[col.pos]) + start - 1
      data.return <- data.pos
      data.return[!is.na(data.pos)] <- data[data.pos[!is.na(data.pos)]]
      data.return[is.na(data.pos)] <- 0
      data.return
    }
  )

  mtx <- do.call(cbind, coldata)
  mtx[is.na(mtx)] <- 0
  mtx
}

.dgCMatrixGroupSubset <- function(x, i, j) {
  if (!x$exists(name = 'data') || !x$exists(name = 'indices')
        || !x$exists(name = 'indptr')) {
    stop("x is not a HDF5 group of dgCMatrix", call. = FALSE)
  }

  dims <- h5attr(x, "dims")
  data <- x[["data"]]
  indices <- x[["indices"]]
  indptr <- x[["indptr"]]

  if (missing(i))
    i <- seq_len(dims[1])
  if (missing(j))
    j <- seq_len(dims[2])

  # Get column information
  newcol.pos <- indptr[j]
  col.n <- indptr[j + 1] - newcol.pos
  start <- newcol.pos + 1
  end <- start + col.n - 1

  # Positions in non-zero data vector
  pos <- integer(sum(col.n))
  n.sum <- 0
  for (idx in seq_along(col.n)) {
    pos[n.sum + seq_len(col.n[idx])] <- seq(start[idx], end[idx])
    n.sum <- n.sum + col.n[idx]
  }

  # Get row index of non-zero data in original matrix
  row.index <- indices[pos]

  # Get position of queried data
  matched.pos <- integer(length(i) * length(j))
  n.sum.1 <- n.sum.2 <- 0
  for (idx in seq_along(col.n)) {
    rng.1 <- n.sum.1 + seq_len(col.n[idx])
    rng.2 <- n.sum.2 + seq_len(length(i))
    matched.pos[rng.2] <- match(i - 1, row.index[rng.1]) + start[idx] - 1
    n.sum.1 <- n.sum.1 + col.n[idx]
    n.sum.2 <- n.sum.2 + length(i)
  }
  matched.pos.na <- is.na(matched.pos)

  # Extract actual data
  matched.data <- matched.pos
  matched.data[!matched.pos.na] <- data[matched.pos[!matched.pos.na]]
  matched.data[matched.pos.na] <- 0

  mtx <- matrix(
    matched.data,
    nrow = length(i),
    ncol = length(j),
    byrow = FALSE
  )

  mtx
}

#' Access cellular data
#'
#' Retrieves data (feature expression, PCA scores, metrics, etc.) for a set
#' of cells in a Seurat object
#'
#' @inheritParams Seurat::FetchData
#' @param object h5Seurat object
#'
#' @return A data frame with cells as rows and cellular data as columns
#'
#' @export
FetchCellData <- function(object, vars, cells = NULL, slot = 'data') {
  # Convert cell names to positions
  cells <- cells %||% seq(1, object[["cell.names"]]$dims)
  if (is.character(x = cells)) {
    cells <- match(cells, Cells(object))
  }
  # Get a list of all objects to search through and their keys
  keyed.field <- c('assays', 'reductions')
  object.keys <- lapply(
    X = keyed.field,
    FUN = function(field) {
      k <- vapply(
        X = names(object[[field]]),
        FUN = function(x) Key(object[[field]][[x]]),
        FUN.VALUE = character(length = 1L),
        USE.NAMES = FALSE
      )
      names(k) <- names(object[[field]])
      return(k)
    }
  )
  names(object.keys) <- keyed.field

  # Find all vars that are keyed
  keyed.types <- lapply(
    X = object.keys,
    FUN = function(keys) {
      lapply(
        X = keys,
        FUN = function(key) {
          if (length(x = key) == 0 || nchar(x = key) == 0) {
            return(integer(length = 0L))
          }
          return(grep(pattern = paste0('^', key), x = vars))
        }
      )

    }
  )

  keyed.types <- Filter(
    f = length,
    x = lapply(
      X = keyed.types,
      FUN = function(x) Filter(f = length, x = x)
    )
  )

  data.fetched <- lapply(
    X = names(keyed.types),
    FUN = function(type) {
      keyed.vars <- keyed.types[[type]]
      lapply(
        X = names(keyed.vars),
        FUN = function(x) {
          vars.use <- vars[keyed.vars[[x]]]
          key.use <- object.keys[[type]][x]
          data.return <- if (type == 'reductions') {
            vars.use <- grep(
              pattern = paste0('^', key.use, '[[:digit:]]+$'),
              x = vars.use,
              value = TRUE
            )
            if (length(x = vars.use) > 0) {
              tryCatch(
                expr = {
                  dim.pos <- as.integer(sub(key.use, '', vars.use))
                  embed.data <- object[[type]][[x]][["cell.embeddings"]]
                  embeddings <- embed.data[cells, dim.pos, drop = FALSE]
                  colnames(embeddings) <- vars.use
                  embeddings
                },
                error = function(...) {
                  return(NULL)
                }
              )
            } else {
              NULL
            }
          } else if (type == 'assays') {
            vars.use <- gsub(
              pattern = paste0('^', key.use),
              replacement = '',
              x = vars.use
            )
            features <- object[[type]][[x]][["features"]][]
            #TODO: try not to load whole matrix
            data.assay <- object[[type]][[x]][[slot]]
            feature.pos <- match(vars.use, features)
            if (sum(is.na(feature.pos)) > 0) {
              warning(
                "features not found: ",
                paste(vars.use[is.na(feature.pos)], collapse = ", "),
                immediate. = TRUE, call. = FALSE)
              vars.use <- vars.use[!is.na(feature.pos)]
              feature.pos <- feature.pos[!is.na(feature.pos)]
            }
            data.vars <- t(x = .dgCMatrixGroupSubset(
              x = data.assay, i = feature.pos, j = cells))
            if (ncol(data.vars) > 0) {
              colnames(x = data.vars) <- paste0(key.use, vars.use)
            }
            data.vars
          }
          data.return <- as.list(x = as.data.frame(x = data.return))
          return(data.return)
        }
      )
    }
  )

  # Unlist twice
  data.fetched <- unlist(x = data.fetched, recursive = FALSE)
  data.fetched <- unlist(x = data.fetched, recursive = FALSE)

  # Pull vars from object metadata
  meta.vars <- vars[vars %in% names(object[["meta.data"]]) & !(vars %in% names(x = data.fetched))]
  meta.list <- lapply(
    X = meta.vars,
    FUN = function(var) object[["meta.data"]][[var]][cells]
  )
  names(meta.list) <- meta.vars
  data.fetched <- c(data.fetched, meta.list)

  # Pull vars from the default assay
  default.assay.features <- object[["assays"]][[DefaultAssay(object)]][["features"]][]
  default.vars <- vars[vars %in% default.assay.features & !(vars %in% names(x = data.fetched))]
  # TODO: warning not found features.
  data.fetched <- c(
    data.fetched,
    tryCatch(
      expr = {
        #TODO: try not to load whole matrix
        data.assay <- object[["assays"]][[DefaultAssay(object)]][[slot]]
        feature.pos <- match(default.vars, default.assay.features)
        data.vars <- t(x = .dgCMatrixGroupSubset(
          x = data.assay, i = feature.pos, j = cells))
        if (ncol(data.vars) > 0) {
          colnames(x = data.vars) <- default.vars
        }
        as.list(x = as.data.frame(x = data.vars))
      },
      error = function(...) {
        return(NULL)
      }
    )
  )

  # Pull identities
  if ('ident' %in% vars && !'ident' %in% names(object[["meta.data"]])) {
    data.fetched[['ident']] <- structure(
      Idents(object = object)[cells],
      names = Cells(x = object)[cells]
    )
  }

  fetched <- names(x = data.fetched)
  vars.missing <- setdiff(x = vars, y = fetched)

  # Assembled fetched vars in a data frame
  data.fetched <- as.data.frame(
    x = data.fetched,
    row.names = Cells(x = object)[cells],
    stringsAsFactors = FALSE
  )

  data.order <- na.omit(object = pmatch(
    x = vars,
    table = fetched
  ))
  if (length(x = data.order) > 1) {
    data.fetched <- data.fetched[, data.order]
  }
  colnames(x = data.fetched) <- vars[vars %in% fetched]

  data.fetched
}
