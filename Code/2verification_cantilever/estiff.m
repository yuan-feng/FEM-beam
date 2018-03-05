function [stiffmatrix] = estiff(A,Izz,Iyy,J,E,v,Length)
   	% A      -- Area
    % Izz    -- Second moment of area over axis zz 
    % Iyy    -- Second moment of area over axis zz
    % J      -- torsion constant
    % E      -- Elastic modulus
    % v      -- Poisson's ratio
    % Length -- Element Length

	% write 4 (2 by 2 ) small matrix first and then put them together. 
	KK1=diag([A/Length, 12*Izz/Length^3, 12*Iyy/Length^3, J/(2*(1+v)*Length), 4*Iyy/Length, 4*Izz/Length]);
	KK1(2,6)=6*Izz/Length^2;
	KK1(6,2)=KK1(2,6);
	KK1(3,5)=-6*Iyy/Length^2;
	KK1(5,3)=KK1(3,5);

	KK2=diag([-A/Length, -12*Izz/Length^3, -12*Iyy/Length^3, -J/(2*(1+v)*Length), 2*Iyy/Length, 2*Izz/Length]);
	KK2(2,6)=6*Izz/Length^2;
	KK2(6,2)=-KK2(2,6);
	KK2(3,5)=-6*Iyy/Length^2;
	KK2(5,3)=6*Iyy/Length^2;

	KK3=diag([-A/Length, -12*Izz/Length^3, -12*Iyy/Length^3, -J/(2*(1+v)*Length), 2*Iyy/Length, 2*Izz/Length]);
	KK3(2,6)=-6*Izz/Length^2;
	KK3(6,2)=6*Izz/Length^2;
	KK3(3,5)=6*Iyy/Length^2;
	KK3(5,3)=-6*Iyy/Length^2;

	KK4=diag([A/Length, 12*Izz/Length^3, 12*Iyy/Length^3, J/(2*(1+v)*Length), 4*Iyy/Length, 4*Izz/Length]);
	KK4(2,6)=-6*Izz/Length^2;
	KK4(6,2)=KK4(2,6);
	KK4(3,5)=6*Iyy/Length^2;
	KK4(5,3)=KK4(3,5);

	% Put them together.
	stiffmatrix=E*[KK1 KK2;KK3 KK4];

end

