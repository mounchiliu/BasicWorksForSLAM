
#include <iostream>
using namespace std;

#include <ctime>
// Eigen 部分
#include <Eigen/Core>
// 稠密矩阵的代数运算（逆，特征值等）
#include <Eigen/Dense>

//size
//set size 100x100
#define MATRIX_SIZE 100

int main()
{
  //dynamic matrix
  //Eigen的固定大小矩阵支持到50,所以需要用动态大小的矩阵
  Eigen::MatrixXd matrix_A(MATRIX_SIZE,MATRIX_SIZE);
  //set random matrix
  matrix_A = Eigen::MatrixXd::Random(MATRIX_SIZE,MATRIX_SIZE);


  //matrix b
  Eigen::MatrixXd v_b;
  v_b = Eigen::MatrixXd::Random(MATRIX_SIZE,1);

  cout<<"A = "<<endl<<matrix_A<<endl;
  cout<<"b = "<<endl<<v_b<<endl;

  //timing
  clock_t time_stt = clock();

  //Inverse
  Eigen::MatrixXd x = matrix_A.inverse()*v_b;
  //print result and processing time
  //Result
  cout<<"x = "<<endl<<x<<endl;
  //Time
  cout<<"Time used in normal inverse is "<<1000*(clock()-time_stt)/(double)CLOCKS_PER_SEC<<"ms"<<endl<<endl;



  //QR
  time_stt = clock();
  x = matrix_A.colPivHouseholderQr().solve(v_b);
  //print result and processing time
  //Result
  cout<<"x = "<<endl<<x<<endl;
  //Time
  cout<<"Time used in QR decomposition is "<<1000*(clock()-time_stt)/(double)CLOCKS_PER_SEC<<"ms"<<endl<<endl;



  //Cholesky

  //only for positive definite matrix
  Eigen::MatrixXd T = Eigen::MatrixXd::Random(MATRIX_SIZE,1);
  bool bok = true;
  //con1. A=AT
  if(matrix_A != matrix_A.transpose())
  {
    bok = false;
    cout<<"Can't use LLT!"<<endl;
  }

  else
  {
    //con2. X!=0, XTAX>0
    Eigen::MatrixXd result = T.transpose()*matrix_A*T;
    if(result(0,0)<=0)
    {
      bok = false;
       cout<<"Can't use LLT!"<<endl;
    }
  }

  if(bok==true)
  {
    time_stt = clock();
    x = matrix_A.llt().solve(v_b);
    //print result and processing time
    //Result
    cout<<"x = "<<endl<<x<<endl;
    //Time
    cout<<"Time used in Cholesky(LLT) decomposition is "<<1000*(clock()-time_stt)/(double)CLOCKS_PER_SEC<<"ms"<<endl;

  }


  return 0;
}

