clear 
clc
% load the input information.
load nnodes.txt;
load coord.txt;
load concen.txt;
load fixity.txt;
load nele.txt;
load ends.txt;
load Area.txt;
load Izz.txt;
load Iyy.txt;
load J.txt;
load E.txt;
load v.txt;
load beta_ang.txt;

% Generate the node id. 
% Each line represents the 6 dofs for one node.
node_id=zeros(nnodes,6);
for i=1:nnodes
    for j=1:6
        node_id(i,j)=(i-1)*6+j;
    end
end

% Generate the member id.
% Each line represents the 12 DOFs for one element.
mem_id=zeros(nele,12);
for i=1:nele
    mem_id(i,1:6)=node_id(ends(i,1),1:6);
    mem_id(i,7:12)=node_id(ends(i,2),1:6);   
end

% Generate the displacement boundary conditions.
% Each line in fixity represents the disp BC on 6 DOFs of one node.
% 0 represents that the DOF is fixed. 
% NaN represnts that the DOF is free.
fixity_tran=fixity';
D=fixity_tran(:);
free_dof=find(isnan(D));
fixed_dof=find(D==0);             


% Generate the force boundary conditions.
% Each line in concen represents the force BC on 6 DOFs of one node.
% 0 represents that the DOF does not has external force. 
% The number represnts that the force values on the DOF.
concen_tran=concen';
P_total=concen_tran(:);               
P_free=P_total(free_dof);

% Generate the xaxis.
% xaxis will be used in the rotation matrix.
% The rotation matrix is to transform the stiffness matrix 
% from global to local coordinate system.
for i=1:nele
    xaxis(i,:)=(coord(ends(i,2),:)-coord(ends(i,1),:));
    Length(i) = norm(xaxis(i,:));
    xaxisunit(i,:)=xaxis(i,:)/Length(i) ;
end

% Generate the local stiffness matrix for each element.
k_stack_local=zeros(nele,12,12);
gamma_stack_local=zeros(nele,12,12);
for i=1:nele
    k_stack_local(i,:,:)=estiff(Area(i),Izz(i),Iyy(i),J(i),E(i),v(i),Length(i));
    gamma_stack_local(i,:,:)=etran(beta_ang(i),xaxisunit(i,1:3));
end

 % Transform the local stiffness matrix to global stiffness matrix.
k_stack_global=zeros(nele,12,12);
for i=1:nele
    k_stack_global(i,:,:)=(squeeze(gamma_stack_local(i,:,:))')...
                          *squeeze(k_stack_local(i,:,:))...
                          *squeeze(gamma_stack_local(i,:,:));
end
  
% Assemble the global stiffness matrix together. 
ndof=6*nnodes;
k_total=zeros(ndof,ndof);
for i=1:nele
    k_total(mem_id(i,1:6),mem_id(i,1:6))=...
        k_total(mem_id(i,1:6),mem_id(i,1:6))+squeeze(k_stack_global(i,1:6,1:6));   
    k_total(mem_id(i,1:6),mem_id(i,7:12))=...
        k_total(mem_id(i,1:6),mem_id(i,7:12))+squeeze(k_stack_global(i,1:6,7:12));
    k_total(mem_id(i,7:12),mem_id(i,1:6))=...
        k_total(mem_id(i,7:12),mem_id(i,1:6))+squeeze(k_stack_global(i,7:12,1:6));
    k_total(mem_id(i,7:12),mem_id(i,7:12))=...
        k_total(mem_id(i,7:12),mem_id(i,7:12))+squeeze(k_stack_global(i,7:12,7:12));
end

% Extract the K_ff. 
K_ff=k_total(free_dof,free_dof);  
K_sf=k_total(fixed_dof,free_dof);  


% C = rcond(A) returns an estimate for the reciprocal condition of A in 1-norm. 
% If A is well conditioned, rcond(A) is near 1.0. 
% If A is badly conditioned, rcond(A) is near 0.
if abs(det(K_ff))<(10^(-15))
    AFLAG=0;
    display('Unstable Structure.');
    display('Maybe you should fix more DOFs .');
else
    AFLAG=1;
    display('The structure is stable.');
end

if abs(rcond(K_ff))<(10^(-15))
    BFLAG=0;
    display('The matrix is badly conditioned. ');
    display('The results may have a large error. ');
    display('Some element may have a way bigger stiffness than others.');
    display('Maybe you need to check the units of element properties.');
else
    BFLAG=1;
    display('The matrix is well-conditioned.');
end

defl=K_ff\P_free;

% Put the deflection to its appropriate location.
% Include the fixed DOF.
defl_vector_total=zeros(ndof,1);
defl_vector_total(free_dof)=defl;

% Each line in DEFL' represents the 6 DOFs at each node.
DEFL=reshape(defl_vector_total,6,nnodes);
DEFL=DEFL';
display(DEFL)

% Generate the reaction.
react=K_sf*defl;

react_vector_total=zeros(ndof,1);
react_vector_total(fixed_dof)=react;
REACT=reshape(react_vector_total,6,nnodes);
REACT=REACT';
display(REACT)
                    
% Calculate the element force.
% Transform the element from from global to local coordiante system.
ELE_FOR=zeros(nele,12);
for i=1:nele
    element_disp_global=zeros(1,12);
    element_disp_global(1:6)=DEFL(ends(i,1),1:6);
    element_disp_global(7:12)=DEFL(ends(i,2),1:6);
    element_disp_global=element_disp_global';
    element_force_local=squeeze(k_stack_local(i,:,:))...
                        *squeeze(gamma_stack_local(i,:,:))*element_disp_global;
    ELE_FOR(i,:)=element_force_local(:);
end
display(ELE_FOR)