function [X, LV, BE, BN] = CP3Mesh(Hmax, eccentricity)
% Given Hmax, eccentricity, construct a mesh
% eccentricity is larger than 0 and less or equal than 1, and is the ratio between the large and small axis of an ellipse
% If eccentricity is equal to 0, then no central hole is considered.
% X: Coordinates - X(i,a) is the i-th coordinate of node a
% LV: Connectivity - LV(a,e) node number of the a-th node of element e
% BE: Boundary elements - BE(a,e), a=1,2, node numbers, BE(3,e), line
% BN: Boundary nodes - BN(1,a), node number in column a, BN(2,a), line 

%% Set the geometry
model = createpde;
R1 = [3,4,-1,1,1,-1,-1,-1,1,1]';
if(eccentricity~=0)
    bigR=0.25;
    smallR=bigR*min(abs(eccentricity),1);
    C1 = [4,0,0,bigR,smallR,0,0,0,0,0]';
    gm = [R1,C1];
    sf = 'R1-C1';
    ns = char('R1','C1');
else
    gm = [R1];
    sf = 'R1';
    ns = char('R1');
end
ns = ns';
g = decsg(gm,sf,ns);
geometryFromEdges(model,g);
%pdegplot(model,'EdgeLabels','on');
axis equal
xlim([-1.1,1.1])
generateMesh(model,'Hmax',Hmax);
figure(1)
pdeplot(model);
%% Extract coordinates and connectivity
X=model.Mesh.Nodes;
LV=model.Mesh.Elements;

%% Extract edges on outer boundary of square
% All edges in the mesh
LE=[LV(1,:), LV(2,:), LV(3,:);
    LV(2,:), LV(3,:), LV(1,:);
    LV(4,:), LV(5,:), LV(6,:)];

% Edges that appear only once 
sortedLE=sort(LE);              % Sort vertices in each edge
sortedLE=sortrows(sortedLE')';    % Dictionary sort of edges

i=1;
while(i<size(sortedLE,2))
    if(sortedLE(1,i)==sortedLE(1,i+1) & sortedLE(2,i)==sortedLE(2,i+1))
        sortedLE(:,i)=0;
        sortedLE(:,i+1)=0;
        i=i+1;
    end
    i=i+1;
end

% Keep only all edges on all boundaries
LE=sortedLE(:,find(sortedLE(1,:)~=0));

% x=1  is line 1
% y=1  is line 2
% x=-1 is line 3
% y=-1 is line 4
% ellipse is line 5

% Boundary elements on the external boundaries
i1=find( X(1,LE(1,:))==1  & X(1,LE(2,:))==1);
i3=find( X(1,LE(1,:))==-1 & X(1,LE(2,:))==-1);
i2=find( X(2,LE(1,:))==1  & X(2,LE(2,:))==1);
i4=find( X(2,LE(1,:))==-1 & X(2,LE(2,:))==-1);

% Remaining elements should be the circle
allindices=[1:size(LE,2)];
allindices([i1, i2, i3, i4])=0;
i5=find(allindices~=0);

% Move middle nodes on curved boundaries to the mid point of the edge
X(:,LE(3,i5))=(X(:,LE(1,i5))+X(:,LE(2,i5)))/2;

%% Form BE matrix
BE = [LE(:,i1), LE(:,i2), LE(:,i3), LE(:,i4), LE(:,i5);
    1*ones(1,size(i1,2)), 2*ones(1,size(i2,2)), 3*ones(1,size(i3,2)), 4*ones(1,size(i4,2)), 5*ones(1,size(i5,2))];

%% Boundary nodes per line - BN(1,:)=node numbers, BN(2,:)=line it belongs to
BN=[BE(1,:), BE(2,:), BE(3,:); BE(4,:), BE(4,:), BE(4,:)];
BN=unique(BN','rows')';
