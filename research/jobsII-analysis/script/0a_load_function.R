..loadFunction <- function(file,function.name) {

    eval(parse(text=paste(function.name," <- function(x) {0}",sep="")),envir = .GlobalEnv)
    suppressMessages(insertSource(file, functions=function.name))
    eval(parse(text=paste(function.name," <- ",function.name,"@.Data",sep="")),envir = .GlobalEnv)

}

..unloadFunction <- function(function.name) {

    eval(parse(text=paste("rm(",function.name,",envir = .GlobalEnv)",sep="")))

}
