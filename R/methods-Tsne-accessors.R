## Setters for Tsne class

setReplaceMethod(
    f = "setCoordinates",
    signature = "Tsne",
    definition = function(theObject, value){
        if (!is.data.frame(value)) stop("Coordinates should be data frame")
        validObject(theObject)
        theObject@coordinates[[1]] <- value
        return(theObject)
    })


setMethod(
    f = "getCoordinates",
    signature = "Tsne",
    definition = function(theObject){
        return(theObject@coordinates[[1]])
    })