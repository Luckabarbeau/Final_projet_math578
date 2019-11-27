
x=linspace(0,1,100000000);
tic

for i=1:length(x)
    y(i)=sin(log(x(i)+1));
end
toc
parfor i=1:length(x)
    y(i)=sin(log(x(i)+1));
end
toc