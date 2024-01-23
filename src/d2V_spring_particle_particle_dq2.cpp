#include <d2V_spring_particle_particle_dq2.h>

void d2V_spring_particle_particle_dq2(Eigen::Ref<Eigen::Matrix66d> H, Eigen::Ref<const Eigen::Vector3d> q0,  Eigen::Ref<const Eigen::Vector3d> q1, double l0, double stiffness) {
    Eigen::Matrix3d I = Eigen::Matrix3d::Identity();
    Eigen::Matrix66d BTB;
    BTB << I, -I,
            -I, I;
    Eigen::Vector3d dx = q1 - q0;
    Eigen::Vector6d BTdx;
    BTdx << -dx, dx;
    double dxdot = dx.dot(dx);
    H = 0.5 * stiffness * 
        (0.5 * l0 * pow(dxdot, -1.5) * BTdx.transpose() * BTdx +
        (pow(dxdot, 0.5) - l0) * pow(dxdot, -0.5) * BTB);
}