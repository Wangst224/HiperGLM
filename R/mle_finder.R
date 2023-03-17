find_mle_bfgs = function(design, outcome, func_likelihood, func_gradient) {

    n_pred = ncol(design)
    coeff_init = rep(1, n_pred)

    op_result = optim(coeff_init,
                      func_likelihood,
                      func_gradient,
                      design = design,
                      outcome = outcome,
                      method = "BFGS",
                      control = list(fnscale = -1))

    return(op_result$par)
}

find_mle_linear_pseudo_inv = function(design, outcome) {

    A = crossprod(design)
    b = crossprod(design, outcome)

    as.vector(qr.solve(A, b))
}

find_mle_logit_newton = function(design, outcome, max_iter = 100, init = rep(0, ncol(design))) {

    coeff = init
    iter = 0
    abs_tol = 1e-6
    rel_tol = 1e-6

    while (iter < max_iter) {
        pre_loglikelihood = log_likelihood_logit(design, outcome, coeff)

        coeff = take_one_newton_step(coeff,
                                     function(x) log_likelihood_logit_gradient(design, outcome, x),
                                     function(x) log_likelihood_logit_hessian(design, outcome, x))

        post_loglikelihood = log_likelihood_logit(design, outcome, coeff)

        if (are_all_close(pre_loglikelihood, post_loglikelihood, abs_tol = abs_tol, rel_tol =  rel_tol)) {
            break
        }

        iter = iter + 1

    }

    if (iter == max_iter) {
        warning("Max iteration reached without convergence.")
    }

    return(as.vector(coeff))
}
