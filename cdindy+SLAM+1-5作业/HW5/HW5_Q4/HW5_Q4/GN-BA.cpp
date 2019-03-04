//
// Created by xiang on 12/21/17.
//

#include <Eigen/Core>
#include <Eigen/Dense>

using namespace Eigen;

#include <vector>
#include <fstream>
#include <iostream>
#include <iomanip>

#include "sophus/se3.h"

using namespace std;

typedef vector<Vector3d, Eigen::aligned_allocator<Vector3d>> VecVector3d;
typedef vector<Vector2d, Eigen::aligned_allocator<Vector3d>> VecVector2d;
typedef Matrix<double, 6, 1> Vector6d;

string p3d_file = "../HW5_Q4/p3d.txt";
string p2d_file = "../HW5_Q4/p2d.txt";

int main(int argc, char **argv) {

    VecVector2d p2d;
    VecVector3d p3d;
    Matrix3d K;
    double fx = 520.9, fy = 521.0, cx = 325.1, cy = 249.7;
    K << fx, 0, cx, 0, fy, cy, 0, 0, 1;

    // load points in to p3d and p2d 

    // START YOUR CODE HERE
    ifstream File_2d,File_3d;
    File_2d.open(p2d_file);
    File_3d.open(p3d_file);

    while(!File_2d.eof() && !File_3d.eof())
    {
      Vector2d p2d_v;
      File_2d>>p2d_v(0)>>p2d_v(1);
      p2d.push_back(p2d_v);

      Vector3d p3d_v;
      File_3d>>p3d_v(0)>>p3d_v(1)>>p3d_v(2);
      p3d.push_back(p3d_v);
    }

    File_2d.close();
    File_3d.close();
    // END YOUR CODE HERE
    assert(p3d.size() == p2d.size());

    int iterations = 100;
    double cost = 0, lastCost = 0;
    int nPoints = p3d.size();
    cout << "points: " << nPoints << endl;

    Sophus::SE3 T_esti; // estimated pose//Initial value == Identity Matix I



    for (int iter = 0; iter < iterations; iter++) {

        Matrix<double, 6, 6> H = Matrix<double, 6, 6>::Zero();
        Vector6d b = Vector6d::Zero();

        cost = 0;

        // compute cost
        for (int i = 0; i < nPoints; i++) {
            // compute cost for p3d[I] and p2d[I]
            // START YOUR CODE HERE

            //Updata R, t
            Matrix3d R = T_esti.matrix().block<3,3>(0,0);//3x3, start from 0,0
            Vector3d t = T_esti.matrix().block<3,1>(0,3);
            //camera coordinate
            Vector3d Pc =  R*p3d[i]+t;
            //pixel coordinate
            Vector3d u = 1/Pc(2) * K * Pc;

            Vector2d e;

            //form of 2D (2x1)
            e <<(p2d[i](0)-u(0)),(p2d[i](1)-u(1));


	    // END YOUR CODE HERE

	    // compute jacobian
            Matrix<double, 2, 6> J;
            // START YOUR CODE HERE
            J<<-fx/Pc(2),0,(fx*Pc(0))/(Pc(2)*Pc(2)),(fx*Pc(0)*Pc(1))/(Pc(2)*Pc(2)),-fx-(fx*Pc(0)*Pc(0))/(Pc(2)*Pc(2)),(fx*Pc(1))/Pc(2),
                0,-fy/Pc(2),(fy*Pc(1))/(Pc(2)*Pc(2)),fy+(fy*Pc(1)*Pc(1))/(Pc(2)*Pc(2)),-fy*Pc(0)*Pc(1)/(Pc(2)*Pc(2)),-fy*Pc(0)/Pc(2);


	    // END YOUR CODE HERE

            H += J.transpose() * J;
            b += -J.transpose() * e;

            cost += (e.transpose()*e)(0);

        }

	// solve dx 
        Vector6d dx;

        // START YOUR CODE HERE 
        dx = H.ldlt().solve(b);

        // END YOUR CODE HERE

/*        if (isnan(dx[0])) {
            cout << "result is nan!" << endl;
            break;
        }*/
        if (iter > 0 && cost >= lastCost) {
            // cost increase, update is not good
            cout << "cost: " << cost << ", last cost: " << lastCost << endl;
            break;
        }

        // update your estimation
        // START YOUR CODE HERE 
        //dx当作扰动
        T_esti =  Sophus::SE3::exp(dx)*T_esti;

        // END YOUR CODE HERE
        
        lastCost = cost;

        cout << "iteration " << iter << " cost=" << cout.precision(12) << cost << endl;
    }

    cout << "estimated pose: \n" << T_esti.matrix() << endl;


    return 0;
}
