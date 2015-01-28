tic;
nSim = 10000;

x = zeros(1,nSim);
parfor i=1:nSim
    a = rand(1000);
    x(i) = prod(sum(a));
    if mod(i,100)==0
        disp(i);
    end
end
toc;
