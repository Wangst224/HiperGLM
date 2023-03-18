expit = function(design, coeff) {
    1 / (1 + exp(-design %*% coeff))
}

take_one_newton_step = function(x0, func_gradient, func_hessian, solver = "qr") {
    supported_solvers = c("qr", "lu")

    if (!(solver %in% supported_solvers)) {
        stop("Solver not supported. Available solvers: 'qr', 'lu'")
    }

    if (solver == "qr") {
        x1 = x0 - qr_solve(func_hessian(x0), func_gradient(x0))
    }

    else if (solver == "lu") {
        x1 = x0 - solve(func_hessian(x0), func_gradient(x0))
    }
}
