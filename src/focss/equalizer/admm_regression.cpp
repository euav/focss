#include "admm_regression.h"
#include <armadillo>
#include "focss/functions.h"

namespace focss {
arma::cx_vec cx_soft_threshold(const arma::cx_vec& values,
                               const double& threshold) {
    arma::cx_vec new_values(values.n_rows);
    for (unsigned long i = 0; i < values.n_rows; ++i) {
        double magnitude = std::sqrt(norm(values(i)));
        double phase = arg(values(i));

        if (magnitude > threshold)
            new_values(i) = i_exp(phase) * (magnitude - threshold);
        else
            new_values(i) = 0;
    }

    return new_values;
}

arma::cx_vec cx_lasso_admm(const arma::cx_mat& A,
                           const arma::cx_vec& b,
                           const double& lambda) {
    using namespace arma;

    unsigned long MAX_ITERATIONS = 2000;
    double ABSTOL = 1e-10;
    double RELTOL = 1e-10;
    double rho = 25;
    double eps_primal;
    double eps_dual;

    unsigned long dim = A.n_cols;

    cx_vec Atb = A.t() * b;
    cx_mat lagrangian = A.t() * A + rho * eye<cx_mat>(dim, dim);
    cx_mat U = trimatu(chol(lagrangian));
    cx_mat Uh = trimatl(U.t());

    cx_vec u(dim, fill::zeros);
    cx_vec x(dim, fill::zeros);
    cx_vec z(dim, fill::zeros);
    cx_vec z_prev(dim, fill::zeros);

    unsigned long iteration = 0;
    while (iteration++ < MAX_ITERATIONS) {
        u = u + x - z;
        x = solve(Uh, Atb + rho * (z - u), solve_opts::fast);
        x = solve(U, x, solve_opts::fast);
        z_prev = z;
        z = cx_soft_threshold(x + u, lambda / rho);

        eps_primal =
            std::min(ABSTOL * sqrt(dim), RELTOL * std::max(norm(x), norm(z)));
        eps_dual = std::min(ABSTOL * sqrt(dim), RELTOL * rho * norm(u));
        if (norm(z - x) < eps_primal && norm(z - z_prev) < eps_dual) break;
    }

    return z;
}

arma::cx_vec cx_lad_admm(const arma::cx_mat& A, const arma::cx_vec& b) {
    using namespace arma;

    unsigned long MAX_ITERATIONS = 1000;
    double ABSTOL = 1e-6;
    double RELTOL = 1e-6;
    double rho = 100;
    double eps_primal;
    double eps_dual;

    unsigned long rows = A.n_rows;
    unsigned long cols = A.n_cols;

    cx_mat Ah = A.t();
    cx_mat AhA = A.t() * A + 1e-3;

    cx_vec u(rows, fill::zeros);
    cx_vec x(cols, fill::zeros);
    cx_vec z(rows, fill::zeros);
    cx_vec z_prev(rows, fill::zeros);

    unsigned long iteration = 0;
    while (iteration++ < MAX_ITERATIONS) {
        u = u + A * x - z - b;
        x = solve(AhA, Ah * (b + z - u), solve_opts::fast);
        z_prev = z;
        z = cx_soft_threshold(A * x - b + u, 1 / rho);

        eps_primal = ABSTOL * sqrt(rows) +
                     RELTOL * std::max(norm(A * x), std::max(norm(z), norm(b)));
        eps_dual = ABSTOL * sqrt(cols) + RELTOL * rho * norm(Ah * u);
        if (norm(A * x - z - b) < eps_primal &&
            rho * norm(Ah * (z - z_prev)) < eps_dual)
            break;
    }

    return z;
}
}  // namespace focss
