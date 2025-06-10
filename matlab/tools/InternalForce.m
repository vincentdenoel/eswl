function [Axial, Shear, Bendi]=InternalForce(NDOF,X,Ke,CORRES,NElem,ELEMDOF)
XX=zeros(NDOF,1); XX(CORRES)=X;

Axial = zeros(NElem,2);
Shear = zeros(NElem,2);
Bendi = zeros(NElem,2);

for iel=1:NElem
    xloc = XX(ELEMDOF(:,iel));
    Fint = Ke(:,:,iel)*xloc;
    Axial(iel,1)= Fint(1);  Axial(iel,2)=-Fint(4);
    Shear(iel,1)= Fint(2);  Shear(iel,2)=-Fint(5);
    Bendi(iel,1)=-Fint(3);  Bendi(iel,2)= Fint(6);
end