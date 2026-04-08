.onLoad <- function(libname, pkgname) {
    s3_register <- function(generic, class, method = NULL) {
        stopifnot(is.character(generic), length(generic) == 1)
        stopifnot(is.character(class), length(class) == 1)
        pieces <- strsplit(generic, "::")[[1]]
        stopifnot(length(pieces) == 2)
        package <- pieces[[1]]
        generic_name <- pieces[[2]]
        if (is.null(method))
            method <- get(paste0(generic_name, ".", class),
                          envir = parent.frame())
        if (package %in% loadedNamespaces())
            registerS3method(generic_name, class, method,
                             envir = asNamespace(package))
        setHook(packageEvent(package, "onLoad"), function(...) {
            registerS3method(generic_name, class, method,
                             envir = asNamespace(package))
        })
    }

    s3_register("ggplot2::autoplot", "bcaboot")
}
