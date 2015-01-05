#ifndef _Deform_H
#define _Deform_H

#define ARAP_DLL_EXPORTS

#include <iostream>
#include <vector>
#include <cmath>

#include <Eigen\Eigen>
#include <Eigen\Sparse>

#include "WunderSVD3x3.h"

#include "arap_wrapper.h"

class ARAP_DLL_API Deform
{
public:

    enum ArapType
    {
        ORIGIN_HARD,
        HARD_SOFT
    };

    typedef std::vector< std::vector<int> > AdjList;
    typedef std::vector< std::vector<int> > TriangleList;
    typedef std::vector< float > VectorF;
    typedef std::vector< int > VectorI;

    typedef std::vector< Eigen::Vector3i > FaceList;

    struct Constraint
    {
        Constraint(Eigen::Vector3f _ctr, int _cid){ ctr = _ctr; cid = _cid; }

        Eigen::Vector3f ctr;
        int cid;
    };

    struct ConstraintCompare
    {
        bool operator()(const Constraint& ctr1, const Constraint& ctr2)
        {
            return ctr1.cid < ctr2.cid;
        }
    };

    typedef std::vector< Constraint > ConstraintSet;

public:
    Deform(){};
    ~Deform(){};
    Deform(const float *P_data, int P_Num, const AdjList &adj_list, const TriangleList &triangle_list);
    
    // complete deform call
    float *do_Deform();

    // deform call by iteration number
    float *do_Deform(int iter_num);
    
    // get current P' point i stores from P_Prime[3*i + 0]
    float *get_P_Prime();
    
    // do deform one iterate, returned current P'
    float *do_Deform_Iter(float &delta);

    // energy value
    float energy();
    
    // set energy tolerance
    void set_tolerance(float tol);

    // set max iteration
    void set_max_iteration(int iter_num);

    // set lambda. Default is 5
    void set_lambda(float lambda);

    // initial hard constraints
    void set_hard_ctrs(const VectorF &T, const VectorI &idx_T);

    // initial soft constraints
    void set_soft_ctrs(const VectorF &T, const VectorI &idx_T);

    // arap type
    void set_arap_type(ArapType at);
    
private:
    // pre-build weight matrix
    void build_weight_matrix();

    // pre-build laplacian matrix
    void build_laplacian_matrix();

    // build right-hand d
    void build_rh_d();

    // set linear equation system
    void set_linear_sys();

    // update R
    void update_Ri();

    // update P
    float update_P_Prime();

    // pre-compute weights
    float compute_wij(const float *p1, const float *p2, const float *p3, const float *p4 = nullptr);

    // find share vertex of an edge
    void find_share_vertex(int i, int j, VectorI &share_vertex);

    // hard constraint or not
    int is_hard_ctrs(int id);

    // soft constraint or not
    int is_soft_ctrs(int id);

private:
    ArapType at;

    const float* P_data;
    Eigen::Matrix3Xf P_Prime;
    Eigen::Matrix3Xf P;
    AdjList adj_list;
    FaceList face_list;

    ConstraintSet hard_ctrs;
    ConstraintSet soft_ctrs;

    std::vector<Eigen::Matrix3f> R;

    Eigen::SparseMatrix<float> Weight;
    Eigen::SparseMatrix<float> L;
    Eigen::SparseLU<Eigen::SparseMatrix<float>> chol;
//    Eigen::SimplicialCholesky<Eigen::SparseMatrix<float>> chol;
    Eigen::MatrixX3f d;

    int P_Num;
    int max_iter;
    double min_tol;
    float lambda;
};

#endif