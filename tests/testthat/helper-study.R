fixture <- function() {
  root <- tempfile("chronomicrobiome test ")
  path <- make_demo(root, seed=42L)
  list(root=root,path=path)
}
edit_config <- function(f, ...) {
  cfg <- yaml::read_yaml(f$path)
  edits <- list(...)
  for (nm in names(edits)) cfg[nm] <- edits[nm]
  yaml::write_yaml(cfg,f$path)
}
write_input <- function(x,path) utils::write.table(x,path,sep="\t",quote=TRUE,row.names=FALSE,na="NA")
