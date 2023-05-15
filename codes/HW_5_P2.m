X = [-1, -1, 0, 0, 1, 1;...
    1, -1, 1, -1, 1, -1];

LV = [1, 3, 4, 3;...
    2, 2, 6, 6;...
    3, 4, 3, 5];

BE = [1, 2, 4, 6, 5, 3;...
    2, 4, 6, 5, 3, 1;...
    1, 2, 3, 4, 4, 1];
HH = [1, 0, 0, 0, 0, 0];
nh = size(BE, 2);
EtaG = [5, 6];
nel = 4;
nod = 6;
ke = 1;
ng=length(EtaG); II=eye(nod);
centrd = zeros([2, 4]);
K=zeros(nod,nod); F=zeros(nod,1);
for ii=1:nel
    xe = X(:,LV(:,ii));
    centrd(:,ii) = mean(xe,2);
    lge=LV(:,ii);
    dN=[xe(2,2)-xe(2,3),xe(2,3)-xe(2,1),xe(2,1)-xe(2,2);...
    xe(1,3)-xe(1,2),xe(1,1)-xe(1,3),xe(1,2)-xe(1,1)];
    Ae2=dN(2,3)*dN(1,2)-dN(1,3)*dN(2,2);
    % local shape functions
    N1 = @(x, y) (1/Ae2)*(dN(1,1).*(x - xe(1,2))+ dN(2,1).*(y - xe(2,2)));
    N2 = @(x, y) (1/Ae2)*(dN(1,2).*(x - xe(1,3))+ dN(2,2).*(y - xe(2,3)));
    N3 = @(x, y) (1/Ae2)*(dN(1,3).*(x - xe(1,1))+ dN(2,3).*(y - xe(2,1)));

    if ii==2
        N_test = N3;
        gradN_test = (1/Ae2)*[dN(1,3), dN(2,3)];
        N_test(-0.5, -0.1)
    end
    
    % mass term
    N11 = @(x, y) N1(x,y).*N1(x,y); N12 = @(x, y) N1(x,y).*N2(x,y); N13 = @(x, y) N1(x,y).*N3(x,y);
    N21 = @(x, y) N2(x,y).*N1(x,y); N22 = @(x, y) N2(x,y).*N2(x,y); N23 = @(x, y) N2(x,y).*N3(x,y);
    N31 = @(x, y) N3(x,y).*N1(x,y); N32 = @(x, y) N3(x,y).*N2(x,y); N33 = @(x, y) N3(x,y).*N3(x,y);

    N_mat = [TriIntegral(N11,xe(1,:),xe(2,:)), TriIntegral(N12,xe(1,:),xe(2,:)), TriIntegral(N13,xe(1,:),xe(2,:));...
        TriIntegral(N21,xe(1,:),xe(2,:)), TriIntegral(N22,xe(1,:),xe(2,:)), TriIntegral(N23,xe(1,:),xe(2,:));...
        TriIntegral(N31,xe(1,:),xe(2,:)), TriIntegral(N32,xe(1,:),xe(2,:)), TriIntegral(N33,xe(1,:),xe(2,:))];
    
    N_p = @(x, y) N1(x,y).*N1(x,y);
%     TriIntegral(N_p,xe(1,:),xe(2,:))

    lp1 = @(x, y) (x +y).*N1(x,y);
    lp2 = @(x, y) (x +y).*N2(x,y);
    lp3 = @(x, y) (x +y).*N3(x,y);
    Fe = [TriIntegral(lp1,xe(1,:),xe(2,:)); TriIntegral(lp2,xe(1,:),xe(2,:)); TriIntegral(lp3,xe(1,:),xe(2,:))];
    dN=dN/Ae2;
    Ke=Ae2/2*ke*dN'*dN + N_mat;
    K(lge,lge) = K(lge,lge) + Ke;
    F(lge) = F(lge) + Fe;
    Ke, Fe
end

for ig=1:ng
    K(EtaG(ig),:)=II(EtaG(ig),:);
    F(EtaG(ig))=X(2,EtaG(ig));
end

for ied=1:nh
    if HH(1,ied)~=0
        lged=BE(1:2,ied);
        xed(:,1:2)=X(:,lged);
        Led=norm(xed(:,1)-xed(:,2));
        pts = LV(:,BE(3,ied));
        xe = X(:,pts);
        dN=[xe(2,2)-xe(2,3),xe(2,3)-xe(2,1),xe(2,1)-xe(2,2);...
        xe(1,3)-xe(1,2),xe(1,1)-xe(1,3),xe(1,2)-xe(1,1)];
        Ae2=dN(2,3)*dN(1,2)-dN(1,3)*dN(2,2);
        N1 = @(y) (1/Ae2)*(dN(2,1).*(y - xe(2,2))).*(y.^2 - 1);
        N2 = @(y) (1/Ae2)*(dN(2,2).*(y - xe(2,3))).*(y.^2 - 1);
        q1 = integral(N1, -1, 1);
        q2 = integral(N2, -1, 1);
        q1, q2
        F(lged)=F(lged)+[q1;q2];
    end
end
%% solve algebraic system
U=K\F;
%% plot
trisurf(LV',X(1,:),X(2,:),U)
colorbar
xlabel('x1'), ylabel('x2'), zlabel('Solution')
x = linspace(-1,1, 50);
y = linspace(-1,1, 50);
[X1,Y1] = meshgrid(x,y);
Z = X1+Y1;

for ii=1:50
    for jj=1:50
        Z(ii,jj) = pllt(N3,X1(ii,jj),Y1(ii,jj),X(1,LV(:,4)),X(2,LV(:,4)));
    end
end

figure(2)
mesh(X1,Y1,Z)
colorbar

% for ii=1:nel
%     xe = X(:,LV(:,ii));
%     centrd(:,ii) = mean(xe,2);
%     lge=LV(:,ii);
%     dN=[xe(2,2)-xe(2,3),xe(2,3)-xe(2,1),xe(2,1)-xe(2,2);...
%     xe(1,3)-xe(1,2),xe(1,1)-xe(1,3),xe(1,2)-xe(1,1)];
%     Ae2=dN(2,3)*dN(1,2)-dN(1,3)*dN(2,2);
%     N1 = @(x, y) (1/Ae2)*(dN(1,1).*(x - xe(1,2))+ dN(2,1).*(y - xe(2,2)));
%     N2 = @(x, y) (1/Ae2)*(dN(1,2).*(x - xe(1,3))+ dN(2,2).*(y - xe(2,3)));
%     N3 = @(x, y) (1/Ae2)*(dN(1,3).*(x - xe(1,1))+ dN(2,3).*(y - xe(2,1)));
% 
%     for j=1:50
%         for k=1:50
%             Z(j,k) = Z(j,k)+ computeSol(xe(1,:), xe(2,:), U(lge), N1, N2, N3,X1(j,k),Y1(j,k));
%         end
%     end
% end
% 
% figure(3)
% mesh(X1,Y1,Z)
% colorbar

%%
% Reference for TriIntegral:
% https://www.mathworks.com/matlabcentral/answers/430420-integrate-a-function-over-a-triangle-area
function I = TriIntegral(f, Tx, Ty)
% I = TriIntegral(f, Tx, Ty)
% 2D integration of f on a triangle
% INPUTS:
%   - f is the vectorized function handle that when calling f(x,y) returns
%       function value at (x,y), x and y are column vectors
%   - Tx,Ty are is two vectors of length 3, coordinates of the triangle
% OUTPUT
%   I: integral of f in T
T = [Tx(:), Ty(:)];
I = integral2(@(s,t) fw(f,s,t,T),0,1,0,1);
A = det(T(2:3,:)-T(1,:));
I = I*abs(A);
end % TriIntegral
%%
function y = fw(f, s, t, T)
sz = size(s);
w1 = (1-s); % Bug fix
w2 = s.*t;
w3 = 1-w1-w2;
P = [w1(:),w2(:),w3(:)] * T;
y = feval(f,P(:,1),P(:,2));
y = s(:).*y(:);
y = reshape(y,sz);
end

function z = pllt(N,xq,yq,xv,yv)
    if inpolygon(xq,yq,xv,yv)==1 
        z = N(xq, yq);
    else
        z = 0;
    end
end

function z = computeSol(xv, yv, Unod, N1, N2, N3,xq,yq)
    if inpolygon(xq,yq,xv,yv)==1 
        z = Unod(1)*pllt(N1,xq,yq,xv,yv) + Unod(2)*pllt(N2,xq,yq,xv,yv) + Unod(3)*pllt(N3,xq,yq,xv,yv);
    else
        z = 0;
    end
end
