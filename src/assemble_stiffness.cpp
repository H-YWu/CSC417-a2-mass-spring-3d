#include <assemble_stiffness.h>

void assemble_stiffness(Eigen::SparseMatrixd &K, Eigen::Ref<const Eigen::VectorXd> q, Eigen::Ref<const Eigen::VectorXd> qdot, 
                     Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::MatrixXi> E, Eigen::Ref<const Eigen::VectorXd> l0, 
                     double k) { 
    // K = - (\sum_{i=0}^{m-1} S_i^T H_i S_i )
    //  [q0; q1] = S_i q
    //  S_i(x, p0 + x) = 1 and S_i(x + 3, p1 + x) = 1 for x = 0, 1, 2
    K.resize(q.size(), q.size());
    std::vector<Eigen::Triplet<double>> triples;
    for (int i = 0; i < E.rows(); i ++) {
        int p0 = E(i, 0) * 3;
        int p1 = E(i, 1) * 3;
        Eigen::Vector3d q0 = q.segment(p0, 3);
        Eigen::Vector3d q1 = q.segment(p1, 3);
        Eigen::Matrix66d Hi;
        d2V_spring_particle_particle_dq2(Hi, q0, q1, l0(i), k);
        for (int a = 0; a < 3; a ++ ) {
            for (int b = 0; b < 3; b ++ ) {
                triples.push_back(Eigen::Triplet<double>(p0 + a, p0 + b, -Hi(a, b)));
                triples.push_back(Eigen::Triplet<double>(p0 + a, p1 + b, -Hi(a, b + 3)));
                triples.push_back(Eigen::Triplet<double>(p1 + a, p0 + b, -Hi(a + 3, b)));
                triples.push_back(Eigen::Triplet<double>(p1 + a, p1 + b, -Hi(a + 3, b + 3)));
            }
        }
    }
    K.setFromTriplets(triples.begin(), triples.end());
}