function om = Drag(v,Cd1,Cd2)
om = zeros(size(v));
om(v>0)=Cd1*v(v>0).^2 + Cd2./(v(v>0).^2);