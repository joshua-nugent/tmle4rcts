`%nin%` <- negate(`%in%`)

mean.harmonic <- function(a){
  1/mean(1/a) #compute the harmonic mean
}

aggregate_IC <- function(recordIC, id) {
  if (is.null(id)) return(recordIC)
  aggregatedIC <- as.matrix(aggregate(recordIC, list(id=id), sum)[, -1, drop = FALSE])
  num.records <- nrow(recordIC)
  num.clusters <- nrow(aggregatedIC)
  aggregatedIC <- aggregatedIC * num.clusters / num.records
  return(aggregatedIC)
}