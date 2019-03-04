#include <opencv2/opencv.hpp>
#include <string>
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Dense>

using namespace std;
using namespace cv;

// this program shows how to use optical flow

string file_1 = "../HW6_Q4/left.png";  // first image
string file_2 = "../HW6_Q4/right.png";  // second image
string file_3 = "../HW6_Q4/disparity.png";//disparoty image

// TODO implement this funciton
/**
 * single level optical flow
 * @param [in] img1 the first image
 * @param [in] img2 the second image
 * @param [in] kp1 keypoints in img1
 * @param [in|out] kp2 keypoints in img2, if empty, use initial guess in kp1
 * @param [out] success true if a keypoint is tracked successfully
 * @param [in] inverse use inverse formulation?
 */
void OpticalFlowSingleLevel(
        const Mat &img1,
        const Mat &img2,
        const vector<KeyPoint> &kp1,
        vector<KeyPoint> &kp2,
        vector<bool> &success,
        bool inverse = false
);

// TODO implement this funciton
/**
 * multi level optical flow, scale of pyramid is set to 2 by default
 * the image pyramid will be create inside the function
 * @param [in] img1 the first pyramid
 * @param [in] img2 the second pyramid
 * @param [in] kp1 keypoints in img1
 * @param [out] kp2 keypoints in img2
 * @param [out] success true if a keypoint is tracked successfully
 * @param [in] inverse set true to enable inverse formulation
 */
void OpticalFlowMultiLevel(
        const Mat &img1,
        const Mat &img2,
        const vector<KeyPoint> &kp1,
        vector<KeyPoint> &kp2,
        vector<bool> &success,
        bool inverse = false
);

/**
 * get a gray scale value from reference image (bi-linear interpolated)
 * @param img
 * @param x
 * @param y
 * @return
 */
inline float GetPixelValue(const cv::Mat &img, float x, float y) {
    uchar *data = &img.data[int(y) * img.step + int(x)];
    float xx = x - floor(x);
    float yy = y - floor(y);
    return float(
            (1 - xx) * (1 - yy) * data[0] +
            xx * (1 - yy) * data[1] +
            (1 - xx) * yy * data[img.step] +
            xx * yy * data[img.step + 1]
    );
}


int main(int argc, char **argv) {

    // images, note they are CV_8UC1, not CV_8UC3
    Mat img1 = imread(file_1, 0);
    Mat img2 = imread(file_2, 0);
    Mat mDisparity = imread(file_3,0);

    // key points, using GFTT here.
    vector<KeyPoint> kp1;
    Ptr<GFTTDetector> detector = GFTTDetector::create(1000, 0.01, 20); // maximum 500 keypoints
    detector->detect(img1, kp1);

    // now lets track these key points in the second image
    // first use single level LK in the validation picture

    // then test multi-level LK
    vector<KeyPoint> kp2_multi;
    vector<bool> success_multi;
    OpticalFlowMultiLevel(img1, img2, kp1, kp2_multi, success_multi);


    Mat img_left;
    cv::cvtColor(img1, img_left, CV_GRAY2BGR);
    for (int i = 0; i < kp2_multi.size(); i++) {
        if (success_multi[i]) {
            int iDisparity = GetPixelValue(mDisparity,kp1[i].pt.x,kp1[i].pt.y);
            int iMeasuredDisparity = kp1[i].pt.x - kp2_multi[i].pt.x;
            int diff = iMeasuredDisparity - iDisparity;
            if(diff<=5 && diff>=-5)
              cv::putText(img_left, to_string(diff), kp2_multi[i].pt,cv::FONT_HERSHEY_COMPLEX,0.5, cv::Scalar(0, 250, 0));
            else if(diff<=10 && diff>=-10)
              cv::putText(img_left, to_string(diff), kp2_multi[i].pt,cv::FONT_HERSHEY_COMPLEX,0.5, cv::Scalar(250, 0, 0));
            else
              cv::putText(img_left, to_string(diff), kp2_multi[i].pt,cv::FONT_HERSHEY_COMPLEX,0.5, cv::Scalar(0, 0, 250));
            cv::circle(img_left, kp2_multi[i].pt, 2, cv::Scalar(0, 250, 0), 2);
            cv::line(img_left, kp1[i].pt, kp2_multi[i].pt, cv::Scalar(0, 250, 0));
        }

    }


    cv::imshow("tracked multi level", img_left);
    cv::waitKey(0);

    return 0;
}

void OpticalFlowSingleLevel(
        const Mat &img1,
        const Mat &img2,
        const vector<KeyPoint> &kp1,
        vector<KeyPoint> &kp2,
        vector<bool> &success,
        bool inverse
) {

    // parameters
    int half_patch_size = 4;
    int iterations = 10;
    bool have_initial = !kp2.empty();

    for (size_t i = 0; i < kp1.size(); i++) {
        auto kp = kp1[i];
        double dx = 0, dy = 0; // dx,dy need to be estimated
        if (have_initial) {
            dx = kp2[i].pt.x - kp.pt.x;
            dy = kp2[i].pt.y - kp.pt.y;
        }

        double cost = 0, lastCost = 0;
        bool succ = true; // indicate if this point succeeded

        // Gauss-Newton iterations
        for (int iter = 0; iter < iterations; iter++) {
            Eigen::Matrix2d H = Eigen::Matrix2d::Zero();
            Eigen::Vector2d b = Eigen::Vector2d::Zero();
            cost = 0;

            if (kp.pt.x + dx <= half_patch_size || kp.pt.x + dx >= img1.cols - half_patch_size ||
                kp.pt.y + dy <= half_patch_size || kp.pt.y + dy >= img1.rows - half_patch_size) {   // go outside
                succ = false;
                break;
            }

            // compute cost and jacobian
            for (int x = -half_patch_size; x < half_patch_size; x++)
                for (int y = -half_patch_size; y < half_patch_size; y++) {

                    // TODO START YOUR CODE HERE (~8 lines)

                    double I1,I2;
                    Eigen::Vector2d GradientI;

                    I1 = GetPixelValue(img1,kp.pt.x+x,kp.pt.y+y);
                    I2 = GetPixelValue(img2,kp.pt.x+x+dx,kp.pt.y+y+dy);
                    double error = I2-I1;
                    Eigen::Vector2d J;  // Jacobian

                    if (inverse == false) {
                        // Forward Jacobian
                      GradientI(0) = GetPixelValue(img2,kp.pt.x+x+dx,kp.pt.y+y+dy)
                                     - GetPixelValue(img2,kp.pt.x+x+dx-1,kp.pt.y+y+dy);//x
                      GradientI(1) = GetPixelValue(img2,kp.pt.x+x+dx,kp.pt.y+y+dy)
                                     - GetPixelValue(img2,kp.pt.x+x+dx,kp.pt.y+y+dy-1);//y

                      J = GradientI;


                    } else {
                        // Inverse Jacobian
                        // NOTE this J does not change when dx, dy is updated, so we can store it and only compute error
                      GradientI(0) = 0.5*(GetPixelValue(img1,kp.pt.x+x+1,kp.pt.y+y)
                          - GetPixelValue(img1,kp.pt.x+x-1,kp.pt.y+y));//x
                      GradientI(1) = 0.5*(GetPixelValue(img1,kp.pt.x+x,kp.pt.y+y+1)
                          - GetPixelValue(img1,kp.pt.x+x,kp.pt.y+y-1));//y

                      J = GradientI;



                    }

                    // compute H, b and set cost;
                    H += J*J.transpose();
                    b += -J*error;
                    cost += error*error;
                    // TODO END YOUR CODE HERE
                }

            // compute update
            // TODO START YOUR CODE HERE (~1 lines)
            Eigen::Vector2d update;
            update = H.ldlt().solve(b);
;
            // TODO END YOUR CODE HERE

/*            if (isnan(update[0])) {
                // sometimes occurred when we have a black or white patch and H is irreversible
                cout << "update is nan" << endl;
                succ = false;
                break;
            }*/
            if (iter > 0 && cost > lastCost) {
                cout << "cost increased: " << cost << ", " << lastCost << endl;
                break;
            }

            // update dx, dy
            dx += update[0];
            dy += update[1];
            lastCost = cost;
            succ = true;
        }

        success.push_back(succ);

        // set kp2
        if (have_initial) {
            kp2[i].pt = kp.pt + Point2f(dx, dy);
        } else {
            KeyPoint tracked = kp;
            tracked.pt += cv::Point2f(dx, dy);
            kp2.push_back(tracked);
        }
    }
}

void OpticalFlowMultiLevel(
        const Mat &img1,
        const Mat &img2,
        const vector<KeyPoint> &kp1,
        vector<KeyPoint> &kp2,
        vector<bool> &success,
        bool inverse) {

    // parameters
    int pyramids = 4;
    double pyramid_scale = 0.5;
    double scales[] = {1.0, 0.5, 0.25, 0.125};

    // create pyramids
    vector<Mat> pyr1, pyr2; // image pyramids
    // TODO START YOUR CODE HERE (~8 lines)
    for (int i = 0; i < pyramids; i++) {

      Mat temp1, temp2;
      if(i==0){
        pyr1.push_back(img1);
        pyr2.push_back(img2);
      }
      else{

        pyrDown(pyr1[pyr1.size()-1],temp1);
        pyrDown(pyr2[pyr2.size()-1],temp2);
        pyr1.push_back(temp1);
        pyr2.push_back(temp2);

      }
    }
    // TODO END YOUR CODE HERE

    // coarse-to-fine LK tracking in pyramids
    // TODO START YOUR CODE HERE
    for (int i = pyramids-1; i >= 0; i--) {
      KeyPoint kp_temp;
      vector<KeyPoint> kp1_temp;

      for(int j=0; j<kp1.size();j++)
      {
        kp_temp.pt.x = kp1[j].pt.x*scales[i];
        kp_temp.pt.y = kp1[j].pt.y*scales[i];
        kp1_temp.push_back(kp_temp);

        if(i!=pyramids-1)
        {
          if(j==10)
          {
            cout<<kp2[j].pt.x<<endl;
          }
          kp2[j].pt.x = kp2[j].pt.x*(1/pyramid_scale);
          kp2[j].pt.y = kp2[j].pt.y*(1/pyramid_scale);

          if(j==10)
            cout<<kp2[j].pt.x<<endl;
        }

      }


      OpticalFlowSingleLevel(pyr1[i],pyr2[i],kp1_temp,kp2,success,inverse);
    }



    // TODO END YOUR CODE HERE
    // don't forget to set the results into kp2
}