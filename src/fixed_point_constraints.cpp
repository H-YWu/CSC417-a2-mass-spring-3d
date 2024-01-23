#include <fixed_point_constraints.h>
#include <algorithm>
void fixed_point_constraints(Eigen::SparseMatrixd &P, unsigned int q_size, const std::vector<unsigned int> indices) {
    P.resize(3 * (q_size - indices.size()), 3 * q_size);
    int id = 0, cnt = 0;
    std::vector<Eigen::Triplet<int>> triple;
    for (int i = 0; i < q_size; i ++) {
        if (i == indices[id]) {
            id ++;
            continue;
        } else {
            for (int j = 0; j < 3; j ++) {
                triple.push_back(Eigen::Triplet<int>(i + j, cnt, 1));
            }
        }
    }
    P.setFromTriplets(triple.begin(), triple.end());
}