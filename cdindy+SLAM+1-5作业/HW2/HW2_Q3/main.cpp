#include <iostream>
using namespace std;

#include <Eigen/Core>
// Eigen 几何模块
#include <Eigen/Geometry>

int main()
{
  //Quaternion
  Eigen::Quaterniond q1(0.55,0.3,0.2,0.2);
  Eigen::Quaterniond q2(-0.1,0.3,-0.7,0.2);

  //normalized
  q1 = q1.normalized();
  q2 = q2.normalized();

  //vector t
  Eigen::Vector3d t1(0.7,1.1,0.2);
  Eigen::Vector3d t2(-0.1,0.4,0.8);

  //coordinates in the first coordinate system
  Eigen::Vector3d p1(0.5,-0.1,0.2);

  //Transform
  //Method 1
  Eigen::Matrix3d R1;//3x3 rotation matrix
  R1 = q1.toRotationMatrix();
  Eigen::Matrix3d R2;//3x3 rotation matrix
  R2 = q2.toRotationMatrix();

  //x_c = Rx_w + t
  Eigen::Vector3d p_w = R1.inverse()*(p1-t1);
  Eigen::Vector3d p_2 = R2*p_w + t2;
  cout<<"Method 1 p_2 = "<<endl<<p_2 << endl<<endl;

  //Method 2
  //v_rotated = q1*x_w
  //v_1 = q1*x_w + t
  Eigen::Vector3d p_w_2 = q1.inverse() * (p1-t1);
  Eigen::Vector3d p_2_2 = q2*p_w_2+t2;
  cout<<"Method 1 p_2_2 = "<<endl<<p_2_2 << endl<<endl;




  return 0;
}

