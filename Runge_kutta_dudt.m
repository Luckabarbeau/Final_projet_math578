function Y= Runge_kutta_dudt(Y,t,dt,A)
% Y the solution at time T dt the delta of time X the Operator that
% modified in function of dY/dt in function of time( in heat equation that is zero) A is the
% differential operator apply to Y to get dY/dt

k1=dt*A*Y;%+X(t);
k2=dt*A*(Y+k1*0.5);%+X(t+dt/2);
k3=dt*A*(Y+k2*0.5);%+X(t/+dt2);
k4=dt*A*(Y+k3);%+X(t+dt);
Y=1/6*(k1+2*k2+2*k3+k4);

end
