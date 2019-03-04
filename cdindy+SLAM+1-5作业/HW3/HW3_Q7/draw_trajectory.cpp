#include <sophus/se3.h>
#include <string>
#include <iostream>
#include <fstream>

// need pangolin for plotting trajectory
#include <pangolin/pangolin.h>

using namespace std;

// path to trajectory file
string ground_file = "../HW3_Q7/groundtruth.txt";
string estimated_file = "../HW3_Q7/estimated.txt";

// function for plotting trajectory, don't edit this code
// start point is red and end point is blue
void DrawTrajectory(vector<Sophus::SE3, Eigen::aligned_allocator<Sophus::SE3>>,vector<Sophus::SE3, Eigen::aligned_allocator<Sophus::SE3>>);
void ShowRMSE(vector<Sophus::SE3, Eigen::aligned_allocator<Sophus::SE3>>, vector<Sophus::SE3, Eigen::aligned_allocator<Sophus::SE3>>);

int main() {

    vector<Sophus::SE3, Eigen::aligned_allocator<Sophus::SE3>> poses_g, poses_e;


    /// implement pose reading code
    // start your code here (5~10 lines)

    //Open 1st file and read data
    ifstream File;
    File.open(ground_file);
    double data[8];

    while(!File.eof())
    {
      //read one line (8 elements)
      for(int j=0; j<8; j++)
        File>>data[j];
      //vector t
      Eigen::Vector3d t(data[1], data[2], data[3]);
      //Quaternion //For Eigen, q = s + ai + bj + ck = [s,a,b,c]
      Eigen::Quaterniond q(data[4],data[5],data[6],data[7]);
      //SE3
      Sophus::SE3 pose_SE3(q,t);
      poses_g.push_back(pose_SE3);
    }
    File.close();



    //Open 2nd file and read data
    File.open(estimated_file);

    while(!File.eof())
    {
      //read one line (8 elements)
      for(int j=0; j<8; j++)
        File>>data[j];
      //vector t
      Eigen::Vector3d t(data[1], data[2], data[3]);
      //Quaternion //For Eigen, q = s + ai + bj + ck = [s,a,b,c]
      Eigen::Quaterniond q(data[4],data[5],data[6],data[7]);
      //SE3
      Sophus::SE3 pose_SE3(q,t);
      poses_e.push_back(pose_SE3);
    }
    File.close();
    // end your code here

    //compute RMSE
    ShowRMSE(poses_g,poses_e);

    // draw trajectory in pangolin
    DrawTrajectory(poses_g,poses_e);



    return 0;
}

/*******************************************************************************************/
void DrawTrajectory(vector<Sophus::SE3, Eigen::aligned_allocator<Sophus::SE3>> poses1,vector<Sophus::SE3, Eigen::aligned_allocator<Sophus::SE3>> poses2) {
    if (poses1.empty()||poses2.empty()) {
        cerr << "Trajectory is empty!" << endl;
        return;
    }

    // create pangolin window and plot the trajectory
    pangolin::CreateWindowAndBind("Trajectory Viewer", 1024, 768);
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    pangolin::OpenGlRenderState s_cam(
            pangolin::ProjectionMatrix(1024, 768, 500, 500, 512, 389, 0.1, 1000),
            pangolin::ModelViewLookAt(0, -0.1, -1.8, 0, 0, 0, 0.0, -1.0, 0.0)
    );

    pangolin::View &d_cam = pangolin::CreateDisplay()
            .SetBounds(0.0, 1.0, pangolin::Attach::Pix(175), 1.0, -1024.0f / 768.0f)
            .SetHandler(new pangolin::Handler3D(s_cam));


    while (pangolin::ShouldQuit() == false) {
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        d_cam.Activate(s_cam);
        glClearColor(1.0f, 1.0f, 1.0f, 1.0f);

        //Draw ground truth
        glLineWidth(2);
        for (size_t i = 0; i < poses1.size() - 1; i++) {
            glColor3f(1 - (float) i / poses1.size(), 0.0f, (float) i / poses1.size());
            glBegin(GL_LINES);
            auto p1 = poses1[i], p2 = poses1[i + 1];
            glVertex3d(p1.translation()[0], p1.translation()[1], p1.translation()[2]);
            glVertex3d(p2.translation()[0], p2.translation()[1], p2.translation()[2]);
            //两点连线
            glEnd();
        }


        //Draw estimated
        glLineWidth(2);
        for (size_t i = 0; i < poses2.size() - 1; i++) {
            glColor3f(1 - (float) i / poses2.size(), 5.0f, (float) i / poses2.size());
            glBegin(GL_LINES);
            auto p1 = poses2[i], p2 = poses2[i + 1];
            glVertex3d(p1.translation()[0], p1.translation()[1], p1.translation()[2]);
            glVertex3d(p2.translation()[0], p2.translation()[1], p2.translation()[2]);
            glEnd();
        }
        pangolin::FinishFrame();
        usleep(5000);   // sleep 5 ms
    }

}


void ShowRMSE(vector<Sophus::SE3, Eigen::aligned_allocator<Sophus::SE3>> poses_g, vector<Sophus::SE3, Eigen::aligned_allocator<Sophus::SE3>> poses_e){

  double Sum_ei_d = 0;
  size_t size = min(poses_g.size(), poses_e.size());

  for (size_t i = 0; i < size ; i++) {

    Sophus::SE3 Tg_inv = poses_g[i].inverse();


    Sophus::SE3 Ei_SE3;
    Ei_SE3 = Tg_inv*poses_e[i];
    //Result 6x1 vector
    typedef Eigen::Matrix<double,6,1> Vector6d;
    Vector6d ei_v= Ei_SE3.log();


    //||ei||2
    Eigen::MatrixXd ei_2 = ei_v.transpose()*ei_v;

    //Sum
    Sum_ei_d = Sum_ei_d + ei_2(0,0);

  }

  //n = poses.size
  double dRMSE = sqrt((1/double(size)) * Sum_ei_d );

  cout << "RMSE = " << dRMSE << endl;



}

