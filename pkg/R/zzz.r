.First.lib <- function(lib, pkg){
    library.dynam("EBS", pkg, lib)
    cat("EBS Loaded \n")
   
    
}
