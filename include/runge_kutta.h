//Input:
//  q - generalized coordiantes for the mass-spring system
//  qdot - generalized velocity for the mass spring system
//  dt - the time step in seconds
//  mass - the mass
//  force(q, qdot) - a function that computes the force acting on the mass as a function. This takes q and qdot as parameters.
//Output:
//  q - set q to the updated generalized coordinate using Runge-Kutta time integration
//  qdot - set qdot to the updated generalized velocity using Runge-Kutta time integration

template<typename FORCE>
inline Eigen::VectorXd compute_diff(Eigen::VectorXd &Q, double mass,  FORCE &force)
{
    Eigen::VectorXd f;
    Eigen::VectorXd q=Q.row(1);
    Eigen::VectorXd qdot=Q.row(0);
    force(f, q, qdot);
    Eigen::VectorXd Qdot(2);
    Qdot.row(0)=f/mass;
    Qdot.row(1)=qdot;
    return Qdot;
}

template<typename FORCE> 
inline void runge_kutta(Eigen::VectorXd &q, Eigen::VectorXd &qdot, double dt, double mass,  FORCE &force) {

    Eigen::VectorXd Q(2);
    Q.row(0)=qdot;
    Q.row(1)=q;
    Eigen::VectorXd Q_old = Q;
    Eigen::VectorXd Qdot = compute_diff(Q,mass,force);
    Eigen::VectorXd k1=Qdot;
    Q=Q_old+0.5*dt*k1;
    Qdot = compute_diff(Q,mass,force);
    Eigen::VectorXd k2 = Qdot;
    Q=Q_old+0.5*dt*k2;
    Qdot = compute_diff(Q,mass,force);
    Eigen::VectorXd k3 = Qdot;
    Q=Q_old+dt*k3;
    Qdot = compute_diff(Q,mass,force);
    Eigen::VectorXd k4 = Qdot;
    Eigen::VectorXd Q_new = Q_old+1.0/6.0*dt*(k1+2*k2+2*k3+k4);
    qdot=Q_new.row(0);
    q=Q_new.row(1);


    
}