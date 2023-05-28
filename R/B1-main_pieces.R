##############################
#### estimate_weights

#' Estimate weights
#'
#' Estimate inverse propensity score weights (\code{z.wt}) or principal score weights (\code{c.wt} and \code{n.wt}), or both
#'
#' @param targeted If TRUE, principal score model is fit to weighted data that target the principal stratum covariate space. If FALSE (default), this weighting is not used. NOTE: if TRUE, \code{Z.form} must be specified.
#' @inheritParams .estimate_nuisance
#' @inheritParams .do_all_main
#' @export

estimate_weights <- function(data,
                             Z.form = NULL,
                             C.form = NULL,
                             targeted = FALSE) {

    if (!.is.declared(data))
        stop("Data needs to be declared to the package through function declare_data() first.")

    if (!is.null(C.form) && targeted && is.null(Z.form))
        stop("To estimate principal score weights with 'targeted' option, must specify Z.form.")

    nuis <- .estimate_nuisance(data = data,
                               Z.form = Z.form,
                               C.form = C.form,
                               targeted = targeted)

    nuis[!(names(nuis) %in% c("s.wt", "z", "c", "n", "y"))]

}




# TODO: make dataset for Y.FUN so user can build their Y.FUN properly

# TODO: make synth dataset


