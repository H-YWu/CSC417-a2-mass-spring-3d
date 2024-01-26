#include <assemble_forces.h>
#include <iostream>

void assemble_forces(Eigen::VectorXd &f, Eigen::Ref<const Eigen::VectorXd> q, Eigen::Ref<const Eigen::VectorXd> qdot, 
                     Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::MatrixXi> E, Eigen::Ref<const Eigen::VectorXd> l0, 
                     double mass, double k) { 
    f.resize(q.size());
    f.setZero();
    for (int i = 0; i < q.size(); i += 3) {
        Eigen::Vector3d gravityForce;
        dV_gravity_particle_dq(gravityForce, mass, Eigen::Vector3d(0., -9.8, 0.));
        f.segment(i, 3) += gravityForce;
    }
    for (int i = 0; i < E.rows(); i ++) {
        int p0 = E(i, 0) * 3;
        int p1 = E(i, 1) * 3;
        Eigen::Vector3d q0 = q.segment(p0, 3);
        Eigen::Vector3d q1 = q.segment(p1, 3);
        Eigen::Vector6d springForce;
        dV_spring_particle_particle_dq(springForce, q0, q1, l0(i), k);
        f.segment(p0, 3) += springForce.segment(0, 3); 
        f.segment(p1, 3) += springForce.segment(3, 3); 
    }
}