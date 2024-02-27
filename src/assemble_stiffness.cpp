#include <assemble_stiffness.h>

void assemble_stiffness(Eigen::SparseMatrixd &K, Eigen::Ref<const Eigen::VectorXd> q, Eigen::Ref<const Eigen::VectorXd> qdot, 
                     Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::MatrixXi> E, Eigen::Ref<const Eigen::VectorXd> l0, 
                     double k) { 
    // K = - (\sum_{i=0}^{m-1} S_i^T H_i S_i )
    //  [q0; q1] = S_i q
    //  S_i(3a+b, p(a)+b) = 1 for a=0..1, b=0..2
    //  S_i^T(p(a)+b, 3a+b) = 1 for a=0..1, b=0..2
    //  H(x, y) for x=0..5, y=0..5
    //  HS(x, p(a)+b) -= H(x, 3a+b)
    //  S^THS(p(a0)+b0, p(a1)+b1) += HS(3a0+b0, p(a1)+b1) -= H(3a0+b0, 3a1+b1)
    K.resize(q.size(), q.size());
    std::vector<Eigen::Triplet<double>> triples;
    for (int i = 0; i < E.rows(); i ++) {
        int p0 = E(i, 0) * 3;
        int p1 = E(i, 1) * 3;
        Eigen::Vector3d q0 = q.segment(p0, 3);
        Eigen::Vector3d q1 = q.segment(p1, 3);
        Eigen::Matrix66d Hi;
        d2V_spring_particle_particle_dq2(Hi, q0, q1, l0(i), k);
        for (int b0 = 0; b0 < 3; b0 ++) {
            for (int b1 = 0; b1 < 3; b1 ++) {
                triples.push_back(Eigen::Triplet<double>(p0+b0, p0+b1, -Hi(b0,b1)));
                triples.push_back(Eigen::Triplet<double>(p1+b0, p1+b1, -Hi(b0+3,b1+3)));
                triples.push_back(Eigen::Triplet<double>(p0+b0, p1+b1, -Hi(b0,b1+3)));
                triples.push_back(Eigen::Triplet<double>(p1+b0, p0+b1, -Hi(b0+3,b1)));
            }
        }
        
    }
    K.setFromTriplets(triples.begin(), triples.end());
}