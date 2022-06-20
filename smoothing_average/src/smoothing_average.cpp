
#include<iostream>
#include<string>
#include<vector>
#include<random>
#include<algorithm>
#include<cmath>
#include<time.h>
#include<sstream>
#include<Eigen/Dense>
#include<unsupported/Eigen/MatrixFunctions>

// Include special ros libraries
#include "ros/ros.h"

// Message Libraries
#include "std_msgs/String.h"
#include "std_msgs/Int32.h"
#include "std_msgs/Float32.h"
#include "std_msgs/Bool.h"
#include "geometry_msgs/Point.h"
#include "geometry_msgs/Pose.h"
#include <random>

using Eigen::MatrixXd;








const int window = 10;
Eigen::Matrix<double, 4, window> window_matrix = Eigen::MatrixXd::Zero(4, window); // define with zeros: full window matrix (4xwindow)

Eigen::Matrix<double, 4, 1> window_row_averages = Eigen::MatrixXd::Zero(4,1); // define with zeros: avg across window rows (4x1)







void subscriber_callback(const geometry_msgs::Pose::ConstPtr& msg)
{
	Eigen::Matrix<double,4,1> window_vector = Eigen::MatrixXd::Zero(4,1); // make a temp vector 4x1 zeros
  window_vector(0) = msg->position.x;				// populate temp vector with msg.position data (taken from new_rand_mat)
  window_vector(1) = msg->position.y;
  window_vector(2) = msg->position.z;
  window_vector(3) = msg->orientation.x;

  for (int ii=0; ii<window-1; ii++) // from 0 to window-1
  {
  	window_matrix.col(ii) << window_matrix.col(ii+1); // shift all window columns leftward, adding the temp vector of new values to the rightmost column
  }
  window_matrix.col(window-1) << window_vector;

  window_row_averages = window_matrix.rowwise().mean(); // recalculate the average across each row, rewrite window_row_averages 4x1
  std::cout<<"\n"<<window_row_averages.transpose()<<std::endl; 
  //std::cout<<window_matrix<<std::endl;
}







int main(int argc, char *argv[])
{
  
  ros::init(argc, argv, "smoothing_average");
  ros::NodeHandle nh;

  double controller_rate = 40.0;        // Hz
  double dt_loop = 1.0/controller_rate; // period of each loop
  ros::Rate loop_rate(controller_rate);

  ros::spinOnce();

  int count = 0;



// ----------------------------------------------------------------------------------------------

// PROBLEM 1 (matrix with 0's, 1's, and rand numbers)

  int rows=4, cols=10;
  double mean=0.5, stddev=20.0;
  std::default_random_engine generator;
  std::normal_distribution<double> distribution(mean, stddev);


  MatrixXd A = MatrixXd::Zero(4, 10); 
  A.col(3) << 1,1,1,1;
  A.col(4) << 1,1,1,1;
  A.col(5) << 1,1,1,1;
  A.col(6) << 1,1,1,1;
  A.col(7) << 1,1,1,1;
  
  Eigen::Matrix<double,4,2> new_rand_mat = Eigen::MatrixXd::Zero(4,2);
  
  for (int ii=0; ii<4; ii++)
  {
  	double val1 = distribution(generator);
  	double val2 = distribution(generator);
  	new_rand_mat(ii,0)=val1;
  	new_rand_mat(ii,1)=val2;
  }

  A.block(0,8,4,2) = new_rand_mat; // BLOCK FCTN: top left row col, height x width
  std::cout<<A<<std::endl;


// PROBLEM 2 (calculating row averages)
  Eigen::Matrix<double, 4,1> row_avgs = Eigen::MatrixXd::Zero(4,1);
  row_avgs = A.rowwise().mean();
  std::cout<<"\n"<<row_avgs.transpose()<<std::endl;


// PROBLEM 3 (shifting values leftward)
  std::cout<<"\n"<<A<<std::endl;

  for (int ii=0; ii<10; ii++)
  {
  	A.col(0) << A.col(1);
  	A.col(1) << A.col(2);
  	A.col(2) << A.col(3);
  	A.col(3) << A.col(4);
  	A.col(4) << A.col(5);
  	A.col(5) << A.col(6);
  	A.col(6) << A.col(7);
  	A.col(7) << A.col(8);
  	A.col(8) << A.col(9);
  	A.col(9) << 0,0,0,0;

  	std::cout<<"\n"<<A<<std::endl;
  }
// -------------------------------------------------------------------------------------------------



// PUBLISHERS
	ros::Publisher random_vector_pub = nh.advertise<geometry_msgs::Pose>("vector", 1); // random_vector_pub "msg" to topic "vector" (new rand values)
	geometry_msgs::Pose msg;


	ros::Publisher window_row_averages_pub = nh.advertise<geometry_msgs::Pose>("averages",1); // window_row_averages_pub "avg_msg" to topic "averages" (avg row values)
	geometry_msgs::Pose avg_msg;

// SUBSCRIBERS
	ros::Subscriber sub = nh.subscribe("vector", 1, subscriber_callback); // sub to topic "vector" (new rand values)



// -------------------------------------------------------------------------------------

  while (ros::ok())
  {
    count++;

    Eigen::Matrix<double,4,1> new_rand_mat = Eigen::MatrixXd::Zero(4,1); // create a vector to fill with new random numbers (4x1 zeros)
    for (int ii=0; ii<4; ii++)
	  {
	  	double val1 = distribution(generator);  // get a random value
	  	new_rand_mat(ii,0)=val1;                // put it in new_rand_mat
	  }
	  msg.position.x = new_rand_mat(0);     // put those random numbers from new_rand_mat into "msg"
	  msg.position.y = new_rand_mat(1);
	  msg.position.z = new_rand_mat(2);
	  msg.orientation.x = new_rand_mat(3);

	  avg_msg.position.x = window_row_averages(0);   // put the calculated window_row_averages into "avg_msg"
	  avg_msg.position.y = window_row_averages(1);
	  avg_msg.position.z = window_row_averages(2);
	  avg_msg.orientation.x = window_row_averages(3);

	  random_vector_pub.publish(msg);     // publish msg
	  window_row_averages_pub.publish(avg_msg);  // publish avg_msg

    ros::spinOnce();
    if( !loop_rate.sleep() ) std::cout<<"Loop Rate was not satisfied.!" <<std::endl;
  }

  return 0;


}
