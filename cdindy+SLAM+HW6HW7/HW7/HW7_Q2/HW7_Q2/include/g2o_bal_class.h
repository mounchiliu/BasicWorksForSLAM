#include <Eigen/Core>
#include "g2o/core/base_vertex.h"
#include "g2o/core/base_binary_edge.h"

#include <sophus/se3.h>

#include "rotation.h"
#include "projection.h"

typedef Eigen::Matrix< double, 9, 1 > Vector9d;

//Camera para
class CameraBAL
{
public:
  Sophus::SE3 mSE3;
  double mf;
  double mk1;
  double mk2;
public:
  CameraBAL(){}
  CameraBAL(Vector9d cam){
    Eigen::Matrix<double,6,1> se3;
    se3.head<3>() = cam.block<3,1>(3,0);//Sophus: t, R //g2o: R, t;
    se3.tail<3>() = cam.head<3>();

    mSE3 = Sophus::SE3::exp(se3);//group
    mf = cam[6];//f
    mk1 = cam[7];//k1
    mk2 = cam[8];//k2
  }
};


class VertexCameraBAL : public g2o::BaseVertex<9,CameraBAL>//dimension & estimate type
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
    VertexCameraBAL() {}

    virtual bool read ( std::istream& /*is*/ )
    {
        return false;
    }

    virtual bool write ( std::ostream& /*os*/ ) const
    {
        return false;
    }

    virtual void setToOriginImpl() {}

    virtual void oplusImpl ( const double* update )
    {
      Eigen::Matrix<double,6,1> t_up;//T -> R + t
      t_up<<update[3],update[4],update[5],update[0],update[1],update[2];//R, t

      _estimate.mSE3 = Sophus::SE3::exp(t_up)*_estimate.mSE3;

      _estimate.mf += update[6];//f
      _estimate.mk1 += update[7];//k1
      _estimate.mk2 += update[8];//k2

    }

};


class VertexPointBAL : public g2o::BaseVertex<3, Eigen::Vector3d>
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
    VertexPointBAL() {}

    virtual bool read ( std::istream& /*is*/ )
    {
        return false;
    }

    virtual bool write ( std::ostream& /*os*/ ) const
    {
        return false;
    }

    virtual void setToOriginImpl() {}

    virtual void oplusImpl ( const double* update )
    {
        Eigen::Vector3d::ConstMapType v ( update );
        _estimate += v;
    }
};

class EdgeObservationBAL : public g2o::BaseBinaryEdge<2, Eigen::Vector2d,VertexCameraBAL, VertexPointBAL>
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
    EdgeObservationBAL() {}

    virtual bool read ( std::istream& /*is*/ )
    {
        return false;
    }

    virtual bool write ( std::ostream& /*os*/ ) const
    {
        return false;
    }

    virtual void computeError() override   // The virtual function comes from the Edge base class. Must define if you use edge.
    {
      //Vertex 0 -> cam pose
      const VertexCameraBAL* cam = static_cast<const VertexCameraBAL*> (vertex(0));
      //Vertex 1 -> 3D point
      const VertexPointBAL* point = static_cast<const VertexPointBAL*> (vertex(1));

      Eigen::Vector3d Pc = cam->estimate().mSE3 * point->estimate(); //camera coordinate

      //Z=1
      double xn = -Pc[0]/Pc[2];
      double yn = -Pc[1]/Pc[2];

      const double& f = cam->estimate().mf;
      const double& k1 = cam->estimate().mk1;
      const double& k2 = cam->estimate().mk2;

      double r2 = xn*xn + yn*yn;
      //undistorted
      double dist = 1.0 + r2*(k1+k2*r2);

      _error[0] = f*dist*xn - _measurement[0];
      _error[1] = f*dist*yn - _measurement[1];

    }



    virtual void linearizeOplus() override
    {
      const VertexCameraBAL* cam = static_cast<const VertexCameraBAL*>(vertex(0));
      const VertexPointBAL* point = static_cast<const VertexPointBAL*>(vertex(1));

      Eigen::Vector3d Pc = cam->estimate().mSE3 * point->estimate(); // position
      double xc = Pc[0];
      double yc = Pc[1];
      double zc = Pc[2];
      double xn = -xc/zc;
      double yn = -yc/zc;
      const double& f = cam->estimate().mf;
      const double& k1 = cam->estimate().mk1;
      const double& k2 = cam->estimate().mk2;
      double r2 = xn*xn + yn*yn;
      double rp = 1.0 + k1*r2 + k2*r2*r2;

      Eigen::Matrix<double, 2, 6> de_dkesi;
      Eigen::Matrix<double, 2, 3> de_dPc;
      Eigen::Matrix<double, 3, 6> dPc_dkesi = Eigen::Matrix<double, 3, 6>::Zero();
      Eigen::Matrix<double, 2, 3> de_dPw;
      Eigen::Matrix<double, 3, 3> dPc_dPw;
      Eigen::Vector2d de_df, de_dk1, de_dk2;

      double zc_2 = zc * zc;
      double zc_3 = zc_2 * zc;
      double rp2 = k1 + 2*k2*r2;


      de_dPc(0, 0) = f*rp/zc + 2*f*xc*xc*rp2/zc_3;
      de_dPc(0, 1) = 2*f*xc*yc*rp2/zc_3;
      de_dPc(0, 2) = -f*xc*rp/zc_2 - 2*f*xc*r2*rp2/zc_2;
      de_dPc(1, 0) = 2*f*xc*yc*rp2/zc_3;
      de_dPc(1, 1) = f*rp/zc + 2*f*yc*yc*rp2/zc_3;
      de_dPc(1, 2) = -f*yc*rp/zc_2 - 2*f*yc*r2*rp2/zc_2;


      dPc_dkesi(0, 3) = 1;
      dPc_dkesi(0, 1) = zc;
      dPc_dkesi(0, 2) = -yc;
      dPc_dkesi(1, 4) = 1;
      dPc_dkesi(1, 0) = -zc;
      dPc_dkesi(1, 2) = xc;
      dPc_dkesi(2, 5) = 1;
      dPc_dkesi(2, 0) = yc;
      dPc_dkesi(2, 1) = -xc;

      de_dkesi = de_dPc * dPc_dkesi;

      de_df(0, 0) = xc*rp/zc;
      de_df(1, 0) = yc*rp/zc;
      de_dk1(0, 0) = f*xc*r2/zc;
      de_dk1(1, 0) = f*yc*r2/zc;
      de_dk2(0, 0) = f*xc*r2*r2/zc;
      de_dk2(1, 0) = f*yc*r2*r2/zc;

      // de_dPw
      dPc_dPw = cam->estimate().mSE3.rotation_matrix();
      de_dPw = de_dPc * dPc_dPw;

      // jacobian
      _jacobianOplusXi.block<2,6>(0,0) = -de_dkesi;
      _jacobianOplusXi.block<2,1>(0,6) = -de_df;
      _jacobianOplusXi.block<2,1>(0,7) = -de_dk1;
      _jacobianOplusXi.block<2,1>(0,8) = -de_dk2;


      _jacobianOplusXj = -de_dPw;


    }
};
