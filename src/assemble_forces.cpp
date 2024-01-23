#include <assemble_forces.h>
#include <iostream>

void assemble_forces(Eigen::VectorXd &f, Eigen::Ref<const Eigen::VectorXd> q, Eigen::Ref<const Eigen::VectorXd> qdot, 
                     Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::MatrixXi> E, Eigen::Ref<const Eigen::VectorXd> l0, 
                     double mass, double k) { 
    f.resize(q.size() * 3);
    f.setZero();
    for (int i = 0; i < E.rows(); i ++) {
        int p0 = E(i, 0) * 3;
        int p1 = E(i, 1) * 3;
        Eigen::Vector3d q0 = q.segment(p0, 3);
        Eigen::Vector3d q1 = q.segment(p1, 3);
        Eigen::Vector3d dx = q1 - q0;
        Eigen::Vector3d rq0 = V.row(p0);
        Eigen::Vector3d rq1 = V.row(p1);
        Eigen::Vector3d dr = rq1 - rq0;
        double lx = dx.dot(dx);
        double l0 = dr.dot(dr);
        Eigen::Vector6d force;
        dV_spring_particle_particle_dq(force, q0, q1, l0, k);   // -dx; dx 
        if (lx < l0) {  // contract
            f.segment(3 * p0, 3) += force.head(3); 
            f.segment(3 * p1, 3) += force.tail(3); 
        } else {    // strech
            f.segment(3 * p0, 3) += -force.head(3); 
            f.segment(3 * p1, 3) += -force.tail(3); 
        }
    }
}