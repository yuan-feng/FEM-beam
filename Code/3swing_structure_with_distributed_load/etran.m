function [gamma] = etran(beta_ang,xaxis)
	% Make the transformation matrix by two steps.
	% Step 1: take xaxis by matrix StepOne. 
	% Step 2: rotate xaxis by matrix StepTwo.
	
	% Initialize the matrix StepOne and StepTwo.
	StepOne=zeros(3);
	StepTwo=zeros(3);

	% StepOne starts!!!
	% 	Deal with the special case when xaxis is parallel to y direction. 
	if (xaxis(1)==0 & xaxis(3)==0)
		StepOne(1,:)=xaxis;
		StepOne(2,:)=[-1,0,0];
		StepOne(3,:)=cross(StepOne(1,:),StepOne(2,:));

	% 	Deal with the normal case
	else
		StepOne(1,:)=xaxis;
		% Change to unit vector
		StepOne(1,:)=StepOne(1,:)/norm(StepOne(1,:));
		yaxis=[0 1 0];
		StepOne(3,:)=cross(StepOne(1,:),yaxis) ;
		StepOne(3,:)=StepOne(3,:)/norm(StepOne(3,:));
		StepOne(2,:)=cross(StepOne(3,:),StepOne(1,:));
		StepOne(2,:)=StepOne(2,:)/norm(StepOne(2,:));
	end
	% % StepTwo starts!!!
	StepTwo(1,1)=1;
	StepTwo(2,2)=cos(beta_ang);
	StepTwo(3,3)=cos(beta_ang);
	StepTwo(2,3)=sin(beta_ang);
	StepTwo(3,2)=-sin(beta_ang);

	% The final gamma:
	gammaPart=StepTwo*StepOne;

	% To 12 by 12
	gamma=zeros(12);
	for i=1:4
		gamma(3*i-2:3*i,3*i-2:3*i)=gammaPart;
	end

end

