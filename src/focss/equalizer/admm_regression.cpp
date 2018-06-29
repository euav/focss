#include "admm_regression.h"
#include <armadillo>
#include "focss/functions.h"

namespace focss {
arma::cx_vec soft_threshold(const arma::cx_vec& values,
                            const double& threshold) {
    arma::cx_vec new_values(values.n_elem);
    for (unsigned long i = 0; i < values.n_elem; ++i) {
        double magnitude = std::sqrt(std::norm(values(i)));

        if (magnitude > threshold)
            new_values(i) *= 1 - threshold / magnitude;
        else
            new_values(i) = 0;
    }

    return new_values;
}

arma::cx_vec proximal_l2(const double threshold, const arma::cx_vec& values) {
    double length = arma::norm(values, 2);
    if (length > threshold)
        return (1 - threshold / length) * values;
    else
        return arma::cx_vec(values.n_elem, arma::fill::zeros);
}

void soft_inplace(const double& threshold, arma::cx_vec& values) {
    for (arma::uword i = 0; i < values.n_elem; ++i) {
        double magnitude = std::abs(values(i));

        if (magnitude > threshold)
            values(i) *= 1 - threshold / magnitude;
        else
            values(i) = 0;
    }
}

arma::cx_vec cx_lasso_admm(const arma::cx_mat& A,
                           const arma::cx_vec& b,
                           const double& lambda) {
    using namespace arma;

    double rho = 1;
    int MAX_ITERATIONS = 1000;
    uword n_features = A.n_cols;

    cx_vec Atb = A.t() * b;
    cx_mat lagrangian = A.t() * A + rho * eye<cx_mat>(n_features, n_features);
    cx_mat U = trimatu(chol(lagrangian));
    cx_mat Uh = trimatl(U.t());

    cx_vec u(n_features, fill::zeros);
    cx_vec x(n_features, fill::zeros);
    cx_vec z(n_features, fill::zeros);

    int iteration = 0;
    std::cout << "> iterating..." << std::endl;
    while (iteration++ < MAX_ITERATIONS) {
        x = solve(Uh, Atb + rho * (z - u));
        x = solve(U, x);
        z = soft_threshold(x + u, lambda / rho);
        u = u + x - z;

        if (iteration % 10 == 0) {
            std::cout << iteration << ' ' << norm(A * x - b, 2) << std::endl;
        }
    }

    return x;
}

arma::cx_vec cx_lad_admm(const arma::cx_mat& A, const arma::cx_vec& b) {
    using namespace arma;

    double rho = 1;
    uword MAX_ITERATIONS = 50;
    uword rows = A.n_rows;
    uword cols = A.n_cols;

    cx_mat At = A.t();
    cx_mat AtA = A.t() * A + 1e-10 * eye<cx_mat>(cols, cols);
    cx_mat U = trimatu(chol(AtA));
    cx_mat Ut = trimatl(U.t());

    cx_vec u(rows, fill::zeros);
    cx_vec x(cols, fill::zeros);
    cx_vec z(rows, fill::zeros);

    uword iteration = 0;
    std::cout << "> iterating..." << std::endl;
    while (iteration++ < MAX_ITERATIONS) {
        x = solve(Ut, At * (b + z - u));
        x = solve(U, x);
        z = soft_threshold(A * x - b + u, 1 / rho);
        u = u + A * x - z - b;

        if (iteration % 10 == 0) {
            std::cout << iteration << ' ' << norm(A * x - b, 1);
            std::cout << ' ' << norm(A * x - b, 2) << std::endl;
        }
    }

    return x;
}

arma::cx_vec cx_ridge_cholesky(const arma::cx_mat& A,
                               const arma::cx_vec& b,
                               const double& alpha) {
    using namespace arma;

    uword cols = A.n_cols;
    uword rows = A.n_rows;
    std::cout << "> decomposing matrix..." << std::endl;
    cx_vec Atb = A.t() * b;
    cx_mat tikhonov = A.t() * A;  // + alpha * eye<cx_mat>(cols, cols);
    tikhonov(cols - 1, cols - 1) += alpha * alpha * rows;
    cx_mat U = chol(tikhonov);
    cx_mat Ut = U.t();

    std::cout << "> solving triangular..." << std::endl;
    cx_vec y = solve(trimatl(Ut), Atb);
    cx_vec x = solve(trimatu(U), y);

    std::cout << "> residual: " << norm(A * x - b, 2) / rows << std::endl;
    return x;
}

arma::cx_vec cx_modulus_cholesky(const arma::cx_mat& A,
                                 const arma::cx_vec& b,
                                 const double& needed) {
    arma::uword index = A.n_cols - 1;
    arma::cx_vec x;
    double alpha_prev = 0;
    double alpha = 1e-2;
    double alpha_next;

    x = cx_ridge_cholesky(A, b, alpha_prev);
    double f_prev = std::abs(x(index)) - needed;

    x = cx_ridge_cholesky(A, b, alpha);
    double f = std::abs(x(index)) - needed;

    std::cout << "> needed: " << needed << std::endl;
    while (std::abs(f) > 1e-3) {
        std::cout << "> targets: " << f_prev << " " << f << std::endl;
        alpha_next = alpha - f * (alpha - alpha_prev) / (f - f_prev);
        alpha_prev = alpha;
        alpha = alpha_next;

        std::cout << "> alpha: " << alpha << std::endl;
        x = cx_ridge_cholesky(A, b, alpha);

        f_prev = f;
        f = std::abs(x(index)) - needed;
        std::cout << "> unity: " << f + needed << std::endl;
    }

    return x;
}

arma::cx_vec cx_ridge_qr(const arma::cx_mat& A,
                         const arma::cx_vec& b,
                         const double& alpha) {
    using namespace arma;

    uword rows = A.n_rows;
    uword cols = A.n_cols;
    std::cout << "> concatenating matrix..." << std::endl;
    cx_mat regularizer(1, cols, fill::zeros);
    regularizer(0, cols - 1) = alpha * std::sqrt(double(rows));
    cx_mat A_reg = join_vert(A, regularizer);
    cx_vec b_reg = join_vert(b, cx_vec{0});

    std::cout << "> solving with qr..." << std::endl;
    cx_vec x = solve(A_reg, b_reg);
    std::cout << "> residual: " << norm(A * x - b, 2) / rows << std::endl;
    return x;
}

arma::cx_vec cx_sag(const arma::cx_mat& A,
                    const arma::cx_vec& b,
                    const double& lambda) {
    using namespace arma;

    uword iteration = 0;
    uword MAX_ITERATIONS = 3000000;
    uword n_samples = A.n_rows;
    uword n_features = A.n_cols;
    uword n_seen = 0;
    uword sample;

    uvec seen(n_samples, fill::zeros);
    cx_vec x(n_features, fill::zeros);
    cx_vec direction(n_features, fill::zeros);
    cx_vec gradients(n_samples, fill::zeros);
    cx_double update;

    vec squared_sum(n_samples, fill::zeros);
    for (uword j = 0; j < n_features; ++j)
        for (uword i = 0; i < n_samples; ++i)
            squared_sum(i) += norm(A(i, j));

    double max_squared_sum = squared_sum.max();
    double step = 1.0 / (max_squared_sum + lambda * lambda);

    std::uniform_int_distribution<uword> dist(0, n_samples - 1);
    std::cout << "> iterating..." << std::endl;
    while (iteration++ < MAX_ITERATIONS) {
        sample = dist(global_urng());
        if (seen(sample) == 0) {
            seen(sample) = 1;
            n_seen++;
        }

        update = as_scalar(A.row(sample) * x - b(sample));
        direction += A.row(sample).t() * (update - gradients(sample));
        gradients(sample) = update;

        x -= (step / n_seen) * direction;
        // soft_inplace(step * lambda * n_samples, x);

        if (iteration % 100000 == 0) {
            std::cout << "  * " << iteration << ", res: ";
            std::cout << norm(A * x - b, 2) / n_samples << " (" << n_seen;
            std::cout << "); step: " << step << std::endl;
        }
    }

    return x;
}
}  // namespace focss
