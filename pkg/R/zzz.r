.onAttach <- function(lib, pkg) {
    library.dynam("EBS", pkg, lib)
    packageStartupMessage("EBS Loaded \n")
   
    
}


