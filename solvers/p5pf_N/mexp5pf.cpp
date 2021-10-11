#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <complex>

#include "mex.h"

#include <Eigen/Dense>
#include "p5pf.h"

namespace colmap {
inline double sign(const double z) {
    return z < 0 ? -1.0 : 1.0;
}
// Stolen from PoseLib implementation
/* Solves the quadratic equation a*x^2 + b*x + c = 0 */
int solve_quadratic_real(double a, double b, double c, double roots[2]) {

    double b2m4ac = b * b - 4 * a * c;
    if (b2m4ac < 0)
        return 0;

    double sq = std::sqrt(b2m4ac);

    // Choose sign to avoid cancellations
    roots[0] = (b > 0) ? (2 * c) / (-b - sq) : (2 * c) / (-b + sq);
    roots[1] = c / (a * roots[0]);

    return 2;
}
void solve_cubic_single_real(double c2, double c1, double c0, double &root) {
    double a = c1 - c2 * c2 / 3.0;
    double b = (2.0 * c2 * c2 * c2 - 9.0 * c2 * c1) / 27.0 + c0;
    double c = b * b / 4.0 + a * a * a / 27.0;
    if (c > 0) {
        c = std::sqrt(c);
        b *= -0.5;
        root = std::cbrt(b + c) + std::cbrt(b - c) - c2 / 3.0;
    } else {
        c = 3.0 * b / (2.0 * a) * std::sqrt(-3.0 / a);
        root = 2.0 * std::sqrt(-a / 3.0) * std::cos(std::acos(c) / 3.0) - c2 / 3.0;
    }
}
/* Solves the quartic equation x^4 + b*x^3 + c*x^2 + d*x + e = 0 */
int solve_quartic_real(double b, double c, double d, double e, double roots[4]) {

    // Find depressed quartic
    double p = c - 3.0 * b * b / 8.0;
    double q = b * b * b / 8.0 - 0.5 * b * c + d;
    double r = (-3.0 * b * b * b * b + 256.0 * e - 64.0 * b * d + 16.0 * b * b * c) / 256.0;

    // Resolvent cubic is now
    // U^3 + 2*p U^2 + (p^2 - 4*r) * U - q^2
    double bb = 2.0 * p;
    double cc = p * p - 4.0 * r;
    double dd = -q * q;

    // Solve resolvent cubic
    double u2;
    solve_cubic_single_real(bb, cc, dd, u2);

    if (u2 < 0)
        return 0;

    double u = sqrt(u2);

    double s = -u;
    double t = (p + u * u + q / u) / 2.0;
    double v = (p + u * u - q / u) / 2.0;

    int sols = 0;
    double disc = u * u - 4.0 * v;
    if (disc > 0) {
        roots[0] = (-u - sign(u) * std::sqrt(disc)) / 2.0;
        roots[1] = v / roots[0];
        sols += 2;
    }
    disc = s * s - 4.0 * t;
    if (disc > 0) {
        roots[sols] = (-s - sign(s) * std::sqrt(disc)) / 2.0;
        roots[sols + 1] = t / roots[sols];
        sols += 2;
    }

    for (int i = 0; i < sols; i++) {
        roots[i] = roots[i] - b / 4.0;

        // do one step of newton refinement
        double x = roots[i];
        double x2 = x * x;
        double x3 = x * x2;
        double dx = -(x2 * x2 + b * x3 + c * x2 + d * x + e) / (4.0 * x3 + 3.0 * b * x2 + 2.0 * c * x + d);
        roots[i] = x + dx;
    }
    return sols;
}

int p5pf(const std::vector<Eigen::Vector2d>& points2D,
    const std::vector<Eigen::Vector3d>& points3D,
    std::vector<Eigen::Matrix<double,3,4>>* poses,
    std::vector<double>* focal_lengths,
    bool normalize_input) {


    if (normalize_input) {
        double focal0 = 0.0;
        for (int i = 0; i < 5; ++i) {
            focal0 += points2D[i].norm();
        }
        focal0 /= 5;

        std::vector<Eigen::Vector2d> scaled_points2D(5);
        for (int i = 0; i < 5; ++i) {
            scaled_points2D[i] = points2D[i] / focal0;
        }

        int n_sols = p5pf(scaled_points2D, points3D, poses, focal_lengths, false);

        for (int i = 0; i < n_sols; ++i) {
            (*focal_lengths)[i] *= focal0;
        }
        return n_sols;
    }

    // Setup nullspace
    Eigen::Matrix<double, 8, 5> cc;
    for (int i = 0; i < 5; i++) {
        cc(0, i) = -points2D[i](1) * points3D[i](0);
        cc(1, i) = -points2D[i](1) * points3D[i](1);
        cc(2, i) = -points2D[i](1) * points3D[i](2);
        cc(3, i) = -points2D[i](1);
        cc(4, i) = points2D[i](0) * points3D[i](0);
        cc(5, i) = points2D[i](0) * points3D[i](1);
        cc(6, i) = points2D[i](0) * points3D[i](2);
        cc(7, i) = points2D[i](0);
    }

    Eigen::Matrix<double, 8, 8> Q = cc.householderQr().householderQ();
    Eigen::Matrix<double, 8, 3> N = Q.rightCols(3);

    // Compute coefficients for sylvester resultant
    double c11_1 = N(0, 1) * N(4, 1) + N(1, 1) * N(5, 1) + N(2, 1) * N(6, 1);
    double c12_1 = N(0, 1) * N(4, 2) + N(0, 2) * N(4, 1) + N(1, 1) * N(5, 2) + N(1, 2) * N(5, 1) + N(2, 1) * N(6, 2) + N(2, 2) * N(6, 1);
    double c12_2 = N(0, 0) * N(4, 1) + N(0, 1) * N(4, 0) + N(1, 0) * N(5, 1) + N(1, 1) * N(5, 0) + N(2, 0) * N(6, 1) + N(2, 1) * N(6, 0);
    double c13_1 = N(0, 2) * N(4, 2) + N(1, 2) * N(5, 2) + N(2, 2) * N(6, 2);
    double c13_2 = N(0, 0) * N(4, 2) + N(0, 2) * N(4, 0) + N(1, 0) * N(5, 2) + N(1, 2) * N(5, 0) + N(2, 0) * N(6, 2) + N(2, 2) * N(6, 0);
    double c13_3 = N(0, 0) * N(4, 0) + N(1, 0) * N(5, 0) + N(2, 0) * N(6, 0);
    double c21_1 = N(0, 1) * N(0, 1) + N(1, 1) * N(1, 1) + N(2, 1) * N(2, 1) - N(4, 1) * N(4, 1) - N(5, 1) * N(5, 1) - N(6, 1) * N(6, 1);
    double c22_1 = 2 * N(0, 1) * N(0, 2) + 2 * N(1, 1) * N(1, 2) + 2 * N(2, 1) * N(2, 2) - 2 * N(4, 1) * N(4, 2) - 2 * N(5, 1) * N(5, 2) - 2 * N(6, 1) * N(6, 2);
    double c22_2 = 2 * N(0, 0) * N(0, 1) + 2 * N(1, 0) * N(1, 1) + 2 * N(2, 0) * N(2, 1) - 2 * N(4, 0) * N(4, 1) - 2 * N(5, 0) * N(5, 1) - 2 * N(6, 0) * N(6, 1);
    double c23_1 = N(0, 2) * N(0, 2) + N(1, 2) * N(1, 2) + N(2, 2) * N(2, 2) - N(4, 2) * N(4, 2) - N(5, 2) * N(5, 2) - N(6, 2) * N(6, 2);
    double c23_2 = 2 * N(0, 0) * N(0, 2) + 2 * N(1, 0) * N(1, 2) + 2 * N(2, 0) * N(2, 2) - 2 * N(4, 0) * N(4, 2) - 2 * N(5, 0) * N(5, 2) - 2 * N(6, 0) * N(6, 2);
    double c23_3 = N(0, 0) * N(0, 0) + N(1, 0) * N(1, 0) + N(2, 0) * N(2, 0) - N(4, 0) * N(4, 0) - N(5, 0) * N(5, 0) - N(6, 0) * N(6, 0);

    double a4 = c11_1 * c11_1 * c23_3 * c23_3 - c11_1 * c12_2 * c22_2 * c23_3 - 2 * c11_1 * c13_3 * c21_1 * c23_3 + c11_1 * c13_3 * c22_2 * c22_2 + c12_2 * c12_2 * c21_1 * c23_3 - c12_2 * c13_3 * c21_1 * c22_2 + c13_3 * c13_3 * c21_1 * c21_1;
    double a3 = c11_1 * c13_2 * c22_2 * c22_2 + 2 * c13_2 * c13_3 * c21_1 * c21_1 + c12_2 * c12_2 * c21_1 * c23_2 + 2 * c11_1 * c11_1 * c23_2 * c23_3 - c11_1 * c12_1 * c22_2 * c23_3 - c11_1 * c12_2 * c22_1 * c23_3 - c11_1 * c12_2 * c22_2 * c23_2 - 2 * c11_1 * c13_2 * c21_1 * c23_3 - 2 * c11_1 * c13_3 * c21_1 * c23_2 + 2 * c11_1 * c13_3 * c22_1 * c22_2 + 2 * c12_1 * c12_2 * c21_1 * c23_3 - c12_1 * c13_3 * c21_1 * c22_2 - c12_2 * c13_2 * c21_1 * c22_2 - c12_2 * c13_3 * c21_1 * c22_1;
    double a2 = c11_1 * c11_1 * c23_2 * c23_2 + c13_2 * c13_2 * c21_1 * c21_1 + c11_1 * c13_1 * c22_2 * c22_2 + c11_1 * c13_3 * c22_1 * c22_1 + 2 * c13_1 * c13_3 * c21_1 * c21_1 + c12_2 * c12_2 * c21_1 * c23_1 + c12_1 * c12_1 * c21_1 * c23_3 + 2 * c11_1 * c11_1 * c23_1 * c23_3 - c11_1 * c12_1 * c22_1 * c23_3 - c11_1 * c12_1 * c22_2 * c23_2 - c11_1 * c12_2 * c22_1 * c23_2 - c11_1 * c12_2 * c22_2 * c23_1 - 2 * c11_1 * c13_1 * c21_1 * c23_3 - 2 * c11_1 * c13_2 * c21_1 * c23_2 + 2 * c11_1 * c13_2 * c22_1 * c22_2 - 2 * c11_1 * c13_3 * c21_1 * c23_1 + 2 * c12_1 * c12_2 * c21_1 * c23_2 - c12_1 * c13_2 * c21_1 * c22_2 - c12_1 * c13_3 * c21_1 * c22_1 - c12_2 * c13_1 * c21_1 * c22_2 - c12_2 * c13_2 * c21_1 * c22_1;
    double a1 = c11_1 * c13_2 * c22_1 * c22_1 + 2 * c13_1 * c13_2 * c21_1 * c21_1 + c12_1 * c12_1 * c21_1 * c23_2 + 2 * c11_1 * c11_1 * c23_1 * c23_2 - c11_1 * c12_1 * c22_1 * c23_2 - c11_1 * c12_1 * c22_2 * c23_1 - c11_1 * c12_2 * c22_1 * c23_1 - 2 * c11_1 * c13_1 * c21_1 * c23_2 + 2 * c11_1 * c13_1 * c22_1 * c22_2 - 2 * c11_1 * c13_2 * c21_1 * c23_1 + 2 * c12_1 * c12_2 * c21_1 * c23_1 - c12_1 * c13_1 * c21_1 * c22_2 - c12_1 * c13_2 * c21_1 * c22_1 - c12_2 * c13_1 * c21_1 * c22_1;
    double a0 = c11_1 * c11_1 * c23_1 * c23_1 - c11_1 * c12_1 * c22_1 * c23_1 - 2 * c11_1 * c13_1 * c21_1 * c23_1 + c11_1 * c13_1 * c22_1 * c22_1 + c12_1 * c12_1 * c21_1 * c23_1 - c12_1 * c13_1 * c21_1 * c22_1 + c13_1 * c13_1 * c21_1 * c21_1;

    a4 = 1.0 / a4;
    a3 *= a4;
    a2 *= a4;
    a1 *= a4;
    a0 *= a4;

    double roots[4];

    int n_roots = solve_quartic_real(a3, a2, a1, a0, roots);
    
    poses->resize(n_roots);
    focal_lengths->resize(n_roots);
    for (int i = 0; i < n_roots; i++) {
        // We have two quadratic polynomials in y after substituting x
        double a = roots[i];
        double c1a = c11_1;
        double c1b = c12_1 + c12_2 * a;
        double c1c = c13_1 + c13_2 * a + c13_3 * a * a;

        double c2a = c21_1;
        double c2b = c22_1 + c22_2 * a;
        double c2c = c23_1 + c23_2 * a + c23_3 * a * a;

        // we solve the first one
        double bb[2];
        if (!solve_quadratic_real(c1a, c1b, c1c, bb))
            continue;

        // and check the residuals of the other
        double res1 = c2a * bb[0] * bb[0] + c2b * bb[0] + c2c;
        double res2;

        // For data where X(3,:) = 0 there is only a single solution
        // In this case the second solution will be NaN
        if (std::isnan(bb[1]))
            res2 = std::numeric_limits<double>::max();
        else
            res2 = c2a * bb[1] * bb[1] + c2b * bb[1] + c2c;

        double b = (std::abs(res1) > std::abs(res2)) ? bb[1] : bb[0];

        Eigen::Matrix<double, 8, 1> p = N.col(0) * a + N.col(1) * b + N.col(2);

        // This gives us the first two rows of the camera matrix
        Eigen::Matrix<double,3,4> pose;        
        pose.row(0) << p(0), p(1), p(2), p(3);
        pose.row(1) << p(4), p(5), p(6), p(7);
                
        double scale = pose.block<1,3>(0,0).norm();
        pose.row(0) /= scale;
        pose.row(1) /= scale;
        pose.block<1,3>(2,0) = pose.block<1,3>(0,0).cross(pose.block<1,3>(1,0));

        // now we need to solve for focal length and t3
        Eigen::Matrix<double, 5, 2> A;
        Eigen::Matrix<double, 5, 1> B;

        for(size_t k = 0; k < 5; ++k) {
            // f * P(1:2,:) * X - x * (P(3,1:3)*X + t3) = 0
            // f * x'*P(1:2,:)*X - |x|^2 * (P(3,1:3)*X + t3) = 0
            // (x'*P(1:2,:)*X  -|x|^2) * [f t3] = |x|^2 * P(3,1:3)*X
            // (x'*P(1:2,:)*X/|x|  -|x|) * [f t3] = |x| * P(3,1:3)*X

            double norm_x = points2D[k].norm();
            Eigen::Vector3d RX = pose.block<3,3>(0,0) * points3D[k];
            
            A(k,0) = points2D[k].dot(RX.topRows<2>() + pose.block<2,1>(0,3)) / norm_x;
            A(k,1) = -norm_x;
            B(k) = norm_x * RX(2);            
        }

        //Eigen::Vector2d ft3 = A.colPivHouseholderQr().solve(B);
        Eigen::Vector2d ft3 = (A.transpose() * A).inverse() * (A.transpose() * B);
        
        double focal = ft3(0);
        pose(2,3) = ft3(1);
        
        if(focal < 0) {
            focal *= -1.0;
            pose.row(0) *= -1.0;
            pose.row(1) *= -1.0;
        }

        // To avoid numerical troubles we project to rotation here. This adds some slight runtime overhead but we might 
        // avoid some headaches later.
        Eigen::JacobiSVD<Eigen::Matrix3d> svd(pose.leftCols<3>(), Eigen::ComputeFullU | Eigen::ComputeFullV);
        pose.leftCols<3>() = svd.matrixU() * svd.matrixV().transpose();

        (*poses)[i] = pose;
        (*focal_lengths)[i] = focal;
    }
    return n_roots;
}

}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    
    if (nrhs != 1) {
        mexErrMsgTxt("One input argument required.");
    }
    if (nlhs > 3) {
        mexErrMsgTxt("Too many output arguments.");
    }
    double const * const matrixAux= static_cast<double const *>(mxGetData(prhs[0]));
    const mwSize *sizeInputMatrix= mxGetDimensions(prhs[0]);
    double*  xX= (double*)malloc(sizeInputMatrix[0] *sizeInputMatrix[1]* sizeof(double));
    double *pindR, *pindt, *pindf;
    const int size0 = sizeInputMatrix[0];
    const int size1 = sizeInputMatrix[1];
    for (int k = 0; k < size0; k++)
    {
        int kOffset = k*size1; // this saves re-computation time
        for (int i = 0; i < size1; i++)
        {
            int iOffset = i*size0;
            xX[i  + kOffset] = matrixAux[i  + kOffset];
        }
    }
    
    std::vector<Eigen::Matrix<double, 3, 4>> poses;
    std::vector<double> focal_lengths;
    std::vector<double> dist_params;
    
    Eigen::Vector2d x1,x2,x3,x4,x5;
    Eigen::Vector3d X1,X2,X3,X4,X5;
    x1 << xX[0],xX[1];
    x2 << xX[3], xX[4];
    x3 << xX[6],xX[7];
    x4 << xX[9],xX[10];
    x5 << xX[12],xX[13];
    X1 << xX[15],xX[16],xX[17];
    X2 << xX[18],xX[19],xX[20];
    X3 << xX[21],xX[22],xX[23];
    X4 << xX[24],xX[25],xX[26];
    X5 << xX[27],xX[28],xX[29];
    
    std::vector<Eigen::Vector2d> x;
    x.push_back(x1); x.push_back(x2); x.push_back(x3); x.push_back(x4); x.push_back(x5);
    std::vector<Eigen::Vector3d> X;
    X.push_back(X1); X.push_back(X2); X.push_back(X3); X.push_back(X4); X.push_back(X5);
//     
  //  for (int j = 0; j < 2; j++) {
   //     mexPrintf("%f  \n", X3(j));
    //}
//     int noOfSol = p4pfr_plnar(x, X, &poses, &focal_lengths, &dist_params);
    int noOfSol = colmap::p5pf(x, X, &poses, &focal_lengths,false);
    plhs[0] = mxCreateDoubleMatrix(3, noOfSol*4, mxREAL);
    plhs[1] = mxCreateDoubleMatrix(1, noOfSol, mxREAL);
    pindR = mxGetPr(plhs[0]);
    pindf = mxGetPr(plhs[1]);
    
    
    for(int i =0; i<noOfSol; i++) {
        
        for (int k = 0; k < 3; k++) {
            for (int j = 0; j < 4; j++) {
                pindR[(k)+3*j+12*i] = poses[i](k,j);
            }
//             mexPrintf(" \n ");
        }
    }
    for(int i =0; i<noOfSol; i++) {
        pindf[i] = focal_lengths[i];        
    }
}