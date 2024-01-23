#include <assemble_stiffness.h>

void assemble_stiffness(Eigen::SparseMatrixd &K, Eigen::Ref<const Eigen::VectorXd> q, Eigen::Ref<const Eigen::VectorXd> qdot, 
                     Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::MatrixXi> E, Eigen::Ref<const Eigen::VectorXd> l0, 
                     double k) { 
    K.resize(q.size() * 3, q.size() * 3);
    std::vector<Eigen::Triplet<double>> triple;
    for (int i = 0; i < E.rows(); i ++) {
        int p0 = E(i, 0) * 3;
        int p1 = E(i, 1) * 3;
        Eigen::Vector3d q0 = q.segment(p0, 3);
        Eigen::Vector3d q1 = q.segment(p1, 3);
        Eigen::Vector3d rq0 = V.row(p0);
        Eigen::Vector3d rq1 = V.row(p1);
        Eigen::Vector3d dr = rq1 - rq0;
        double l0 = dr.dot(dr);
        Eigen::Matrix66d H;
        d2V_spring_particle_particle_dq2(H, q0, q1, l0, k);
        for (int a = 0; a < 3; a ++ ) {
            for (int b = 0; b < 3; b ++ ) {
                triple.push_back(Eigen::Triplet<double>(p0 + a, p0 +b, H(a, b)));
                triple.push_back(Eigen::Triplet<double>(p1 + a, p1 +b, H(a + 3, b + 3)));
            }
        }
    }
    K.setFromTriplets(triple.begin(), triple.end());
}