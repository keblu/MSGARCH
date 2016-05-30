f.error = function(message){
  cat(paste0("try/catch ", message, "\n"))
  return(FALSE)
}