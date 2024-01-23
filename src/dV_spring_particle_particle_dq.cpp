#include <dV_spring_particle_particle_dq.h>

void dV_spring_particle_particle_dq(Eigen::Ref<Eigen::Vector6d> f, Eigen::Ref<const Eigen::Vector3d> q0,  Eigen::Ref<const Eigen::Vector3d>     q1, double l0, double stiffness) {
    Eigen::Vector3d dx = q1 - q0;
    Eigen::Vector6d BTdx;
    BTdx << -dx, dx;
    double dxdot = dx.dot(dx);
    f = 0.5 * stiffness * (sqrt(dxdot) - l0) * pow(dxdot, -0.5) * BTdx;

}