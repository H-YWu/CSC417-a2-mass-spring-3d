#include <mass_matrix_particles.h>

void mass_matrix_particles(Eigen::SparseMatrixd &M, Eigen::Ref<const Eigen::VectorXd> q, double mass) {
    int n = q.size() / 3;
    M.resize(n, n);
    for (int i = 0; i < n; i ++) {
        M.insert(i, i) = mass;
    }
}
