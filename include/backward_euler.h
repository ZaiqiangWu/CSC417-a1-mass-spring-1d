#include <Eigen/Dense>

//Input:
//  q - generalized coordiantes for the mass-spring system
//  qdot - generalized velocity for the mass spring system
//  dt - the time step in seconds
//  mass - the mass
//  force(q, qdot) - a function that computes the force acting on the mass as a function. This takes q and qdot as parameters.
//  stiffness(q, qdot) - a function that computes the stiffness (negative second derivative of the potential energy). This takes q and qdot as parameters.  
//Output:
//  q - set q to the updated generalized coordinate using Backward Euler time integration
//  qdot - set qdot to the updated generalized velocity using Backward Euler time integration

template<typename FORCE, typename STIFFNESS> 
inline void backward_euler(Eigen::VectorXd &q, Eigen::VectorXd &qdot, double dt, double mass,  FORCE &force, STIFFNESS &stiffness) {


    Eigen::MatrixXd stiff;
    stiffness(stiff, q, qdot);
    double a = stiff(0,0)/mass;
    Eigen::MatrixXd A;
    A.resize(2,2);
    A<< 1.0, -dt*a, -dt, 1.0;
    Eigen::VectorXd Q_old(2);
    //Q_old.resize(2,3);
    Q_old(0)=qdot(0);


    Q_old(1)=q(0);

    auto Q = A.inverse()*Q_old;
    qdot=Q.row(0);
    q=Q.row(1);





}