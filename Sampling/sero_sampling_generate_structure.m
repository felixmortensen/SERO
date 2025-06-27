function [S, D, T, V] = sero_sampling_generate_structure(z, iters)
%generate line phantom

if nargin < 1
    z = 50;
    iters = 1;
end

S = zeros(z,iters);
D = zeros(z,iters);
V = zeros(z,iters);
T = zeros(z,iters);

n_obj = randi(7,[iters,1])+3;


for i = 1:iters
    
    w = zeros(z, n_obj(i));
    p = randi(z, [n_obj(i),1]);
    r = abs(randn([n_obj(i),1])*10);
    
    s = rand(n_obj(i),1)*2.8+0.2;
    d = rand(n_obj(i),1)*2+0.3;
    v = rand(n_obj(i),1)*0.5;
    t = rand(n_obj(i),1)*2+0.5;
    
    st = round(p-r);
    st(st<1) = 1;
    
    en = round(p+r);
    en(en>z) = z;
    
    for j = 1:n_obj(i)
        
        if j == 1 % background
            w(:,j) = rand(size(w(:,j)))*0.2;
        else
            w(st(j):en(j),j) = rand(1)*0.5+0.5;
        end
        
        w(:,j) = smooth(w(:,j), ceil(sqrt(rand(1)*10)));
    end
    
    w = w./sum(w,2);
    w(isnan(w)) = 0;

    S(:,i) = w*s;
    D(:,i) = w*d;
    V(:,i) = w*v;
    T(:,i) = w*t;
    
end

if nargout == 1
    tmp(:,1,:) = S;
    tmp(:,2,:) = D;
    tmp(:,3,:) = T;
    tmp(:,4,:) = V;
    S = tmp;
end








