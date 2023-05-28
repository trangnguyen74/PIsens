####%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### declare_data ####

#' Declare data to the package
#'
#' This function processes the data frame with user-provided variable names into a data frame that is used by the functions that implement the main analysis and sensitivity analysis methods.
#'
#' @param data The input dataset.
#' @param z.var Name of the assigned treatment variable.
#' @param c.var Name of the compliance type variable.
#' @param y.var Name of the outcome variable.
#' @param x.vars A character vector containing names of the covariates.
#' @param samp.wt Name of sampling weight variable. Default is \code{NULL}.
#' @param id Name of id variable. Default is \code{NULL}.
#' @param other.vars A vector of names of other variables to be kept in the dataset. Default is \code{NULL}.
#'
#' @return A data frame with just analysis variables, where assigned treatment is renamed \code{z}, compliance type is renamed \code{c}, outcome variable is renamed \code{y}, sampling weight (if used) is renamed \code{s.wt}, and id variable (if used) is renames \code{id}. Covariates (named in \code{x.vars}) and other variables (named in \code{other.vars}) are retained and not renamed.
#'
#' @export

declare_data <- function(data,
                         z.var,
                         c.var,
                         y.var,
                         x.vars,
                         samp.wt = NULL,
                         id = NULL,
                         other.vars = NULL) {

    not.names <- setdiff(c(z.var, c.var, y.var, x.vars, samp.wt, id, other.vars),
                         names(data))
    if (length(not.names) > 0) {
        stop(paste("variable(s)", paste(not.names, collapse = ", "),
                   "not available in data."))
    }


    data <- data[c(id, samp.wt, z.var, c.var, y.var, x.vars, other.vars)]


    names(data)[names(data) %in% c(z.var, c.var, y.var)] <- c("z", "c", "y")

    if (!.is.01(data$z)) stop(paste("z.var must be binary 0/1."))
    if (!.is.01(data$c)) stop(paste("c.var must be binary 0/1."))

    if (!all(is.na(data$c[data$z==0])))
        stop(paste("c.var should be NA for all z.var==0 units."))


    if (is.null(id)) { data$id <- 1:nrow(data)
    } else           { names(data)[names(data)==id] <- "id"
    }


    if (is.null(samp.wt)) { data$s.wt <- 1
    } else                {
        names(data)[names(data)==samp.wt] <- "s.wt"

        if (!is.numeric(data$s.wt))
            stop("samp.wt variable must be numeric.")

        if (anyNA(data$s.wt) || any(data$s.wt < 0))
            stop("samp.wt variable must not contain NA or negative values.")
    }


    data[c("id", "s.wt", "z", "c", "y", x.vars, other.vars)]
}




####%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### .is.declared ####

#' Check if data frame has been processed by function \code{declare_data()}
#'
#' @param data A data frame
#' @keywords internal

.is.declared <- function(data) {
    all(c("id", "s.wt", "z", "c", "y") %in% names(data))
}
