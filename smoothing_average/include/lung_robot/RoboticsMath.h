#pragma once

#define _USE_MATH_DEFINES
#include <cmath>
#include <Eigen/Dense>
#include <vector>

namespace RoboticsMath
{
  // define custom data types
  using Matrix6d = Eigen::Matrix<double, 6, 6>;
  using Vector5d = Eigen::Matrix<double, 5, 1>;
  using Vector6d = Eigen::Matrix<double, 6, 1>;
  using Vector7d = Eigen::Matrix<double, 7, 1>;
  using Vector8d = Eigen::Matrix<double, 8, 1>;

  // template <typename Derived>
  // Eigen::MatrixXd Inverse(Eigen::DenseBase<Derived> const& g) {
  //     Eigen::MatrixXd ginv = Eigen::MatrixXd::Zero(4,4);
  //     ginv.block(0,0,3,3) = g.block(0,0,3,3).transpose(); //R^T
  //     ginv.block(0,3,3,1) = -g.block(0,0,3,3).transpose()*g.block(0,3,3,1); //-R^T*p
  //     ginv(3,3) = 1.0;
  //     return ginv;
  // }

  double deg2rad(double degrees) {
    return degrees * 4.0 * atan(1.0) / 180.0;
  }
  double sgn(double x) {
    double s = (x > 0) - (x < 0);
    return s;
  }
  double trace(Eigen::Matrix<double,3,3> R) {
    double tr = R(0,0) + R(1,1) + R(2,2);
    return tr;
  }
  double vectornorm(Eigen::VectorXd v) {
          return sqrt(v.transpose()*v);
  }
  Eigen::Matrix3d orthonormalize(Eigen::Matrix3d R) {
    Eigen::Matrix3d R_ortho;
    R_ortho.fill(0);
    // Normalize the first column:
    R_ortho.col(0) = R.col(0) / vectornorm(R.col(0));

    // Orthogonalize & normalize second column:
    R_ortho.col(1) = R.col(1);
    double c = (R_ortho.col(1).transpose()*R_ortho.col(0));
    c = c / (R_ortho.col(0).transpose()*R_ortho.col(0));
    R_ortho.col(1) = R_ortho.col(1) - c*R_ortho.col(0);
    R_ortho.col(1) = R_ortho.col(1) / vectornorm(R_ortho.col(1));

    // Orthogonalize & normalize third column:
    R_ortho.col(2) = R.col(2);
    double d = (R_ortho.col(2).transpose()*R_ortho.col(0));
    d = d / (R_ortho.col(0).transpose()*R_ortho.col(0));
    R_ortho.col(2) = R_ortho.col(2) - d*R_ortho.col(0);
    double e = (R_ortho.col(2).transpose()*R_ortho.col(1));
    e = e / (R_ortho.col(1).transpose()*R_ortho.col(1));
    R_ortho.col(2) = R_ortho.col(2) - e*R_ortho.col(1);
    R_ortho.col(2) = R_ortho.col(2) / vectornorm(R_ortho.col(2));
    return R_ortho;
  }

  Eigen::Matrix3d hat3(Eigen::Vector3d v) {
      Eigen::Matrix3d H = Eigen::Matrix<double, 3, 3>::Zero();
      H(0, 1) = -1 * v(2);
      H(0, 2) = v(1);
      H(1, 0) = v(2);
      H(1, 2) = -1 * v(0);
      H(2, 0) = -1 * v(1);
      H(2, 1) = v(0);

      return H;
  }
  Eigen::Matrix4d hat6(const Eigen::Matrix<double,6,1> x) {
  //function calculates the hat operator for a 6x1 matrix (result is 4x4 double)
      Eigen::Matrix4d cross_product;
      cross_product << 0.0, -x(5), x(4), x(0),
                      x(5), 0.0, -x(3), x(1),
                      -x(4), x(3), 0.0, x(2),
                      0.0, 0.0, 0.0, 0.0;
      return cross_product;
  }

  Eigen::Matrix<double,3,1> hom2cart(Eigen::Matrix<double,4,1> ph) {
    Eigen::Matrix<double,3,1> p;
    p << ph(0), ph(1), ph(2);
    return p;
  }
  Eigen::Matrix<double,4,1> cart2hom(Eigen::Matrix<double,3,1> p) {
    Eigen::Matrix<double,4,1> ph;
    ph << p(0), p(1), p(2), 1.0;
    return ph;
  }
  Eigen::Matrix3d quat2rotm(const Eigen::Matrix<double,4,1> x) {
    //function calculates the rotation from a 4x1 quaternion (result is 3x3 double)
    Eigen::Matrix3d rot_mat;
    double qw= x(0); double qx= x(1); double qy= x(2); double qz= x(3);
    rot_mat <<     1.0-2.0*(qy*qy+qz*qz), 2.0*(qx*qy-qz*qw), 2.0*(qx*qz+qy*qw),
                    2.0*(qx*qy+qz*qw), 1.0-2.0*(qx*qx+qz*qz), 2.0*(qy*qz-qx*qw),
                    2.0*(qx*qz-qy*qw), 2.0*(qy*qz+qx*qw), 1.0-2.0*(qx*qx+qy*qy);
    return rot_mat;
  }
  Eigen::Vector4d rotm2quat(Eigen::Matrix3d R) {
      R = orthonormalize(R);
      Eigen::Vector4d Q;
      Q.fill(0);

      double trace = R(0, 0) + R(1, 1) + R(2, 2);
      if (trace > 0)
      {
              double s = 2*sqrt(trace + 1.0);
              Q(0) = 0.25*s; //w eqn
              Q(1) = (R(2, 1) - R(1, 2)) / s; //x eqn
              Q(2) = (R(0, 2) - R(2, 0)) / s; //y eqn
              Q(3) = (R(1, 0) - R(0, 1)) / s; //z eqn
      }
      else
      {
              if (R(0, 0)>R(1, 1) && R(0, 0)>R(2, 2))
              {
                      double s = 2*sqrt(1.0 + R(0, 0) - R(1, 1) - R(2, 2));
                      Q(0) = (R(2, 1) - R(1, 2))/s; //w eqn
                      Q(1) = 0.25*s; //x eqn
                      Q(2) = (R(0, 1) + R(1, 0))/s; //y eqn
                      Q(3) = (R(0, 2) + R(2, 0))/s; //z eqn
              }
              else if (R(1, 1)>R(2, 2))
              {
                      double s = 2*sqrt(1.0 + R(1, 1) - R(0, 0) - R(2, 2));
                      Q(0) = (R(0, 2) - R(2, 0))/s; //w eqn
                      Q(1) = (R(0, 1) + R(1, 0))/s; //x eqn
                      Q(2) = 0.25*s; //y eqn
                      Q(3) = (R(1, 2) + R(2, 1))/s; //z eqn
              }
              else
              {
                      double s = 2*sqrt(1.0 + R(2, 2) - R(0, 0) - R(1, 1));
                      Q(0) = (R(1, 0) - R(0, 1))/s; //w eqn
                      Q(1) = (R(0, 2) + R(2, 0))/s; //x eqn
                      Q(2) = (R(1, 2) + R(2, 1))/s; //y eqn
                      Q(3) = 0.25*s; //z eqn
              }
      }

      return Q;
  }
  Eigen::Matrix4d quat2tform(const Eigen::Matrix<double,4,1> x) {
    //function calculates the rotation from a 4x1 quaternion (result is 3x3 double)
    Eigen::Matrix4d transform_mat;
    double qw= x(0); double qx= x(1); double qy= x(2); double qz= x(3);
    transform_mat <<     1.0-2.0*(qy*qy+qz*qz), 2.0*(qx*qy-qz*qw), 2.0*(qx*qz+qy*qw), 0.0,
                         2.0*(qx*qy+qz*qw), 1.0-2.0*(qx*qx+qz*qz), 2.0*(qy*qz-qx*qw), 0.0,
                         2.0*(qx*qz-qy*qw), 2.0*(qy*qz+qx*qw), 1.0-2.0*(qx*qx+qy*qy), 0.0,
                         0.0,                0.0,                0.0,                1.0;
    return transform_mat;
  }
  double rotationDiff(Eigen::Matrix3d R1, Eigen::Matrix3d R2) {
    Eigen::Matrix3d R_diff;
    R_diff = R1.transpose()*R2;

    double num = R_diff.trace() - 1;
    double den = 2.0;
    double theta_rad = acos(num/den);

    double theta_deg = theta_rad*180.0/M_PI;
    return theta_deg;
  }
  // double quatChordalDiff(Eigen::Matrix<double,4,1> q1, Eigen::Matrix<double,4,1> q2);
  Eigen::Matrix<double,3,3> rotd(Eigen::Matrix<double,3,1> ax, double angle_deg) {
    // Note: angle assumed to be degrees as input
    Eigen::Matrix<double,3,3> R_ax = hat3(ax);
    double ang = angle_deg*M_PI/180.0;

    double ctheta = cos(ang);
    double stheta = sin(ang);
    double vtheta = 1 - ctheta;

    double u1 = ax(0);
    double u2 = ax(1);
    double u3 = ax(2);

    Eigen::Matrix<double,3,3> R;
    R(0,0) = u1*u1*vtheta + ctheta;
    R(0,1) = u1*u2*vtheta - u3*stheta;
    R(0,2) = u1*u2*vtheta + u2*stheta;

    R(1,0) = u1*u2*vtheta + u3*stheta;
    R(1,1) = u2*u2*vtheta + ctheta;
    R(1,2) = u2*u3*vtheta - u1*stheta;

    R(2,0) = u1*u3*vtheta - u2*stheta;
    R(2,1) = u2*u3*vtheta + u1*stheta;
    R(2,2) = u3*u3*vtheta + ctheta;

    return R;

  }
  Eigen::Matrix<double,4,1> rotm2axang(Eigen::Matrix<double,3,3> R) {
    Eigen::AngleAxisd angle_axis(R);
    Eigen::Matrix<double,3,1> axis = angle_axis.axis();

    Eigen::Matrix<double,4,1> axang;
    axang << axis(0), axis(1), axis(2), angle_axis.angle();
    return axang;
  }

  Eigen::Matrix4d assembleTransformation(Eigen::Matrix3d Rot, Eigen::Vector3d Trans) {
    Rot = orthonormalize(Rot);
    Eigen::Matrix4d T;
    T.fill(0);
    T.topLeftCorner(3, 3) = Rot;
    T.topRightCorner(3, 1) = Trans;
    T(3, 3) = 1;
    return T;
  }
  Eigen::Matrix4d inverseTransform(Eigen::Matrix4d T) {
    Eigen::Matrix4d Tinv;
    Tinv.fill(0);
    Tinv.topLeftCorner(3, 3) = T.topLeftCorner(3, 3).transpose();
    Tinv.topRightCorner(3, 1) = -1 * T.topLeftCorner(3, 3).transpose()*T.topRightCorner(3, 1);
    Tinv(3, 3) = 1.0;
    return Tinv;
  }
  Vector7d collapseTransform(Eigen::Matrix4d T) {
    Eigen::Matrix<double, 7, 1> x;
    x.fill(0);
    x.head(3) = T.topRightCorner(3, 1);
    x.tail(4) = rotm2quat(T.topLeftCorner(3, 3));
    return x;
  }

  void getCofactor(double A[6][6], double temp[6][6], int p, int q, int n) {
    // Function to get cofactor of A[p][q] in temp[][]
    int i = 0, j = 0;

    // Looping for each element of the matrix
    for (int row = 0; row < n; row++)
    {
            for (int col = 0; col < n; col++)
            {
                    //  Copying into temporary matrix only those element
                    //  which are not in given row and column
                    if (row != p && col != q)
                    {
                            temp[i][j++] = A[row][col];

                            // Row is filled, so increase row index and
                            // reset col index
                            if (j == n - 1)
                            {
                                    j = 0;
                                    i++;
                            }
                    }
            }
    }
  }
  double determinant6D(double A[6][6], int n) {
    // Recursive function for finding determinant of 6x6 matrix.
    double D = 0; // Initialize result

                              //  Base case : if matrix contains single element
    if (n == 1)
            return A[0][0];

    double temp[6][6]; // To store cofactors

    int sign = 1;  // To store sign multiplier

                               // Iterate for each element of first row
    for (int f = 0; f < n; f++)
    {
            // Getting Cofactor of A[0][f]
            getCofactor(A, temp, 0, f, n);
            D += sign * A[0][f] * determinant6D(temp, n - 1);

            // terms are to be added with alternate sign
            sign = -sign;
    }

    return D;
  }
  void adjoint(double A[6][6], double adj[6][6]) {
      // Function to get adjoint of A[N][N] in adj[N][N].

      // temp is used to store cofactors of A[][]
      int sign = 1;
      double temp[6][6];

      for (int i = 0; i<6; i++)
      {
              for (int j = 0; j<6; j++)
              {
                      // Get cofactor of A[i][j]
                      getCofactor(A, temp, i, j, 6);

                      // sign of adj[j][i] positive if sum of row
                      // and column indexes is even.
                      sign = ((i + j) % 2 == 0) ? 1 : -1;

                      // Interchanging rows and columns to get the
                      // transpose of the cofactor matrix
                      adj[j][i] = (sign)*(determinant6D(temp, 6 - 1));
              }
      }
  }
  void inverse(double A[6][6], double inverse[6][6]) {
      // Function to calculate and store inverse, returns false if
      // matrix is singular

      // Find determinant of A[][]
      double det = determinant6D(A, 6);
      if (det != 0) //matrix is not singular
      {
              // Find adjoint
              double adj[6][6];
              adjoint(A, adj);

              // Find Inverse using formula "inverse(A) = adj(A)/det(A)"
              for (int i = 0; i<6; i++)
                      for (int j = 0; j<6; j++)
                              inverse[i][j] = adj[i][j] / double(det);
      }
  }

  Matrix6d Adjoint_p_q(Eigen::Vector3d p, Eigen::Vector4d q) {
      Eigen::Matrix3d R = quat2rotm(q);
      Eigen::Matrix3d phat = hat3(p);
      Eigen::Matrix<double, 6, 6> Ad = Eigen::Matrix<double, 6, 6>::Zero();
      Ad.topLeftCorner<3, 3>() = R;
      Ad.bottomRightCorner<3, 3>() = R;
      Ad.topRightCorner<3, 3>() = phat*R;

      return Ad;
  }
  Eigen::Vector4d slerp(Eigen::Vector4d qa, Eigen::Vector4d qb, double t) {
      Eigen::Vector4d qm;
      qm.fill(0);

      double cosHalfTheta = qa.transpose()*qb;
      if (std::abs(cosHalfTheta) >= 1.0)
      {
              qm = qa;
              return qm;
      }

      double halfTheta = acos(cosHalfTheta);
      double sinHalfTheta = sqrt(1.0 - cosHalfTheta*cosHalfTheta);

      if (std::abs(sinHalfTheta)<0.001)
      {
              qm(0) = 0.5*qa(0) + 0.5*qb(0);
              qm(1) = 0.5*qa(1) + 0.5*qb(1);
              qm(2) = 0.5*qa(2) + 0.5*qb(2);
              qm(3) = 0.5*qa(3) + 0.5*qb(3);
              return qm;
      }

      double ratioA = sin((1 - t)*halfTheta) / sinHalfTheta;
      double ratioB = sin(t*halfTheta) / sinHalfTheta;

      qm(0) = ratioA*qa(0) + ratioB*qb(0);
      qm(1) = ratioA*qa(1) + ratioB*qb(1);
      qm(2) = ratioA*qa(2) + ratioB*qb(2);
      qm(3) = ratioA*qa(3) + ratioB*qb(3);
      return qm;
  }
  Eigen::Matrix<double, 4, Eigen::Dynamic> quatInterp(Eigen::Matrix<double,4,Eigen::Dynamic> refQuat, Eigen::VectorXd refArcLengths, Eigen::VectorXd interpArcLengths) {
      int count = 0;
      int N = interpArcLengths.size();
      Eigen::MatrixXd quatInterpolated(4, N);
      quatInterpolated.fill(0);

      for (int i = 0; i<N; i++)
      {
              // setLinSpaced sometimes returns very small negative number instead of 0
              if (interpArcLengths(i) < 1E-10)
              {
                interpArcLengths(i) = 0.0;
              }
              if (interpArcLengths(i) < refArcLengths(count + 1)) {
                      count = count + 1;
              }
              double L = refArcLengths(count) - refArcLengths(count + 1);
              double t = (refArcLengths(count) - interpArcLengths(i)) / L;
              quatInterpolated.col(i) = slerp(refQuat.col(count), refQuat.col(count + 1), t);
      }

      return quatInterpolated;
  }
}
