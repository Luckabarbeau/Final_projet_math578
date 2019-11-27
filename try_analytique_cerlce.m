rin=0.05
rout=0.2
Tin=1.3
Tout=1

nb=100
X=linspace(rin,rout,nb)
for i=1:nb
    T(i)=Tout+log(X(i)/rout)/log(rin/rout)*(Tin-Tout);
end

    plot(X,T,'.')