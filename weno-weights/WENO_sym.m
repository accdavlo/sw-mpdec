% WENO5 Uses 5 averages values U(-2:2) and
% 3 candidate stencils S_0 =(x_{0:2}), S_1=(x_{-1:1}), S_2=(x_{-2:0})
% 3 reconstructions v_r, with primitives V_r(x)=\int_{-inf}^x v_r(x) dx
% v_r(x) has degree 2 so V_r(x) has degree 3 and it must interpolate 

clear all
clc

syms x y

J=5;       %number of cell averages (odd)
K=(J+1)/2; %number of stencils and stencil lenght

% Choose the quadrature points in which we compute the reconstruction
%x_quads=[-1/2*sqrt(sym(3/5)), 0 ,1/2*sqrt(sym(3/5))];
%x_quads=[-1/2/sqrt(3),1/2/sqrt(3)];
%x_quads=[0];
%x_quads=[-sym(1/2),sym(1/2)];
x_quads=[-0.5*sqrt(sym(3/7)+sym(2/7)*sqrt(sym(6/5))),-0.5*sqrt(sym(3/7)-sym(2/7)*sqrt(sym(6/5))),...
    0.5*sqrt(sym(3/7)-sym(2/7)*sqrt(sym(6/5))),0.5*sqrt(sym(3/7)+sym(2/7)*sqrt(sym(6/5)))];



xj=[-K+1:K-1]; % index of cell averages

x_sten=[-K+1-1/2:1/2]; %first stencil interface points

% Computing lagrangian interpolation 
for k=1:K+1
    phiBase(k)=sym(1);
    for j=1:K+1
        if k~=j
            phiBase(k)=phiBase(k)*(x-x_sten(j))/(x_sten(k)-x_sten(j));
        end
    end
end
%phiBase={-(x-1)*(x-2)*(x-3)/6,x*(x-2)*(x-3)/2,-x*(x-1)*(x-3)/2,x*(x-1)*(x-2)/6};


% Auxiliary matrix to assemble the linear system
M=zeros(K+1,K);
for i=1:K+1
    for j=1:K
        if i>j
            M(i,j)=1;
        end
    end
end
    

% Compute the coefficients of the low order polynomials
for k=1:K %stencil
    MM = zeros(K+1,J);
    MM(:,k:k+K-1) = M;
    for i = 1:J % cell average
        c{k,i}= 0;
        for r=1:K+1
            z=diff(subs(phiBase(r),x,y-k+1));
            c{k,i}= simplify(c{k,i}+z*MM(r,i));
        end
        
    end
end


% Compute smoothness indicators
beta=sym(zeros(K,J,J));

for m=1:K % stencil
    for j=m:m+K-1 %coeffient U_j
        for l=m:m+K-1 %coefficient U_l
            for dd=1:K-1 % sum on derivatives
                beta(m,l,j)=beta(m,l,j)+...
                    int(diff(c{m,j},dd)*diff(c{m,l},dd),y,-1/2,1/2);
            end
        end 
    end    
end


fprintf("!------------------------------------!\n")
fprintf("!         Smoothness Indicators      !\n")
fprintf("!------------------------------------!\n")
for m=1:K
    fprintf("beta%1d = ",m)
    for j=m:m+K-1
        for l=m:m+K-1
            fprintf("+(%s*Q(%d)*Q(%d))",char(beta(m,l,j)),j-K,l-K)
        end
    end
    fprintf("\n")
end

fprintf("\n")


% Compute the coefficients of the high order polynomials

M2=zeros(J+1,J);
for i=1:J+1
    for j=1:J
        if i>j
            M2(i,j)=1;
        end
    end
end
x_sten=[-K+0.5:K-0.5];
for k=1:J+1
    phiBig(k)=sym(1);
    for j=1:J+1
        if k~=j
            phiBig(k)=phiBig(k)*(y-x_sten(j))/(x_sten(k)-x_sten(j));
        end
    end
end

hoCoefficients= simplify(diff(phiBig)*M2);


%for every quadrature point substitute it into the coefficient formula and
%solve the least square method
for x_quad = x_quads

    for i= 1:size(c,1)
        for j= 1:size(c,2)
            c_quad(i,j)=subs(c{i,j},y,x_quad);
        end
    end
    disp("x_quad")
    disp(x_quad)
    disp(latex(x_quad))
%   
    disp("c_quad")
    disp(c_quad)
    disp(latex(c_quad))
    disp(double(c_quad))

    rhs=subs(hoCoefficients,y,x_quad);
    Amat=inv((c_quad*c_quad.'));
    d=simplify(Amat*(c_quad*rhs.'));
    disp("d")
    disp(d)
    disp(latex(d))
    disp(double(d))
    fprintf("Residual error")

    residu=double(d'*c_quad-rhs);
    disp(residu)

    fprintf("!------------------------------------!\n")
    fprintf("! node %2.14f             !\n",x_quad)
    fprintf("!------------------------------------!\n")
    
    for k=1:K
        fprintf("gamma%1d =%2.16f\n",k,double(d(k)));
    end
%   
    
    for k=1:K
        fprintf("W%1d = ",k)
        for j=1:J
           fprintf("%+2.16f * Q(%d)",double(c_quad(k,j)),j-K)
        end
        fprintf("\n")
    end
end