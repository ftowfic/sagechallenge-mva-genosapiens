# Read only an explicitly supplied, existing local Level 5 extract. No hub resolver.
suppressPackageStartupMessages({library(signatureSearch); library(rhdf5)})
stopifnot(packageVersion('signatureSearch') == '1.24.0', getRversion() == '4.5.2')
local_reference <- function(refdb) {
  if (length(refdb) != 1L || !grepl('\\.h5$', refdb) || !file.exists(refdb))
    stop('An existing local HDF5 extract is required; remote references are disabled')
  normalizePath(refdb, mustWork=TRUE)
}
assignInNamespace('determine_refdb', local_reference, ns='signatureSearch')
args <- commandArgs(trailingOnly=TRUE)
load_indexed_matrix <- function(request) {
  stopifnot(file.exists(request$matrix))
  genes <- as.character(h5read(request$matrix, '0/META/ROW/id'))
  signatures <- as.character(h5read(request$matrix, '0/META/COL/id'))
  stopifnot(!anyDuplicated(genes), !anyDuplicated(signatures),
            length(genes) == request$expected_dimensions$rows,
            length(signatures) == request$expected_dimensions$columns)
  info <- h5ls(request$matrix)
  dimensions <- as.integer(strsplit(info$dim[info$group == '/0/DATA/0' & info$name == 'matrix'], ' x ', fixed=TRUE)[[1]])
  normal <- identical(dimensions, as.integer(c(length(genes), length(signatures))))
  transposed <- identical(dimensions, as.integer(c(length(signatures), length(genes))))
  stopifnot(normal || transposed)
  mapping <- function(ids) list2env(setNames(as.list(seq_along(ids)), ids), hash=TRUE, parent=emptyenv())
  list(matrix=normalizePath(request$matrix, mustWork=TRUE), normal=normal,
       rows=mapping(genes), columns=mapping(signatures), dimensions=request$expected_dimensions)
}

extract_indexed <- function(request, session, notify=function(block) NULL) {
  stopifnot(request$schema == 'genosapiens.track2-lincs-discovery-source/v1',
            identical(normalizePath(request$matrix, mustWork=TRUE), session$matrix),
            request$expected_dimensions$rows == session$dimensions$rows,
            request$expected_dimensions$columns == session$dimensions$columns,
            !file.exists(request$output), request$block_size >= 1L, request$block_size <= 256L,
            !anyDuplicated(request$row_ids), !anyDuplicated(request$signature_ids),
            length(request$row_ids) > 0L, length(request$signature_ids) > 0L)
  positions <- function(ids, mapping) as.integer(unlist(mget(as.character(ids), mapping, inherits=FALSE), use.names=FALSE))
  ri <- positions(request$row_ids, session$rows)
  ci <- positions(request$signature_ids, session$columns)
  stopifnot(!anyNA(ri), !anyNA(ci))
  stopifnot(dir.create(request$output, recursive=TRUE, showWarnings=FALSE, mode='0700'))
  blocks <- list()
  for (start in seq.int(1L, length(ci), by=request$block_size)) {
    positions <- seq.int(start, min(start + request$block_size - 1L, length(ci)))
    local <- ci[positions]
    rr <- sort(ri)
    cc <- sort(local)
    if (session$normal) {
      x <- h5read(session$matrix, '0/DATA/0/matrix', index=list(rr, cc), drop=FALSE)
    } else {
      x <- t(h5read(session$matrix, '0/DATA/0/matrix', index=list(cc, rr), drop=FALSE))
    }
    x <- x[match(ri, rr), match(local, cc), drop=FALSE]
    stopifnot(nrow(x) == length(ri), ncol(x) == length(local))
    colnames(x) <- request$signature_ids[positions]
    name <- sprintf('block-%06d.csv', length(blocks)+1L)
    write.csv(data.frame(gene_id=request$row_ids, x, check.names=FALSE),
              file.path(request$output, name), row.names=FALSE, na='')
    block <- list(path=name, row_ids=I(request$row_ids),
      signature_ids=I(request$signature_ids[positions]), rows=length(ri), columns=length(local))
    blocks[[length(blocks)+1L]] <- block
    notify(block)
  }
  completion <- list(schema='genosapiens.track2-lincs-extraction/v1',status='complete',
    orientation=if (session$normal) 'rows_by_columns' else 'columns_by_rows', blocks=blocks)
  temporary <- file.path(request$output, '.completion.json')
  jsonlite::write_json(completion, temporary, auto_unbox=TRUE, pretty=TRUE)
  stopifnot(file.rename(temporary, file.path(request$output, 'completion.json')))
  completion
}

if (length(args) == 2L && args[1] == '--session') {
  protocol <- 'genosapiens.track2-reader-session/v1'
  send <- function(type, request_id, ...) {
    cat(as.character(jsonlite::toJSON(list(protocol=protocol, type=type, request_id=request_id, ...),
                                     auto_unbox=TRUE, null='null')), '\n', sep='')
    flush(stdout())
  }
  request_id <- 'initialize'
  status <- tryCatch({
    init <- jsonlite::fromJSON(args[2])
    stopifnot(init$protocol == protocol, !is.null(init$matrix_identity))
    session <- load_indexed_matrix(init)
    send('ready', request_id, matrix_identity=init$matrix_identity,
         identifier_loads=1L, orientation=if (session$normal) 'rows_by_columns' else 'columns_by_rows')
    input <- file('stdin', open='r')
    repeat {
      line <- readLines(input, n=1L, warn=FALSE)
      stopifnot(length(line) == 1L, nchar(line, type='bytes') <= 2097152L)
      message <- jsonlite::fromJSON(line)
      stopifnot(message$protocol == protocol, is.character(message$request_id), length(message$request_id) == 1L)
      request_id <- message$request_id
      if (message$type == 'shutdown') {
        h5closeAll()
        send('closed', request_id)
        break
      }
      stopifnot(message$type == 'extract', is.numeric(message$timeout_seconds),
                message$timeout_seconds > 0, message$timeout_seconds <= 900)
      setTimeLimit(elapsed=message$timeout_seconds, transient=TRUE)
      completion <- extract_indexed(message$request, session, function(block) send('block', request_id, block=block))
      setTimeLimit(cpu=Inf, elapsed=Inf, transient=FALSE)
      send('complete', request_id, blocks=length(completion$blocks))
    }
    close(input)
    0L
  }, error=function(error) {
    setTimeLimit(cpu=Inf, elapsed=Inf, transient=FALSE)
    send('error', request_id, message=conditionMessage(error))
    1L
  }, finally=h5closeAll())
  quit(save='no', status=status)
}

if (length(args) == 2L && args[1] == '--indexed') {
  request <- jsonlite::fromJSON(args[2])
  session <- load_indexed_matrix(request)
  extract_indexed(request, session)
  h5closeAll()
  quit(save='no',status=0L)
}
stopifnot(length(args) == 2L, !file.exists(args[2]))
db <- local_reference(args[1])
genes <- as.character(h5read(db, 'rownames'))
conditions <- as.character(h5read(db, 'colnames'))
x <- h5read(db, 'assay')
stopifnot(!anyDuplicated(genes), !anyDuplicated(conditions),
          nrow(x) == length(genes), ncol(x) == length(conditions))
colnames(x) <- conditions
write.csv(data.frame(gene_id=genes, x, check.names=FALSE), args[2], row.names=FALSE, na='')
