#include <mass_matrix_particles.h>

void mass_matrix_particles(Eigen::SparseMatrixd &M, Eigen::Ref<const Eigen::VectorXd> q, double mass) {
    int nx3 = q.size();
    M.resize(nx3, nx3);
    for (int i = 0; i < nx3; i ++) {
        M.insert(i, i) = mass;
    }
}
