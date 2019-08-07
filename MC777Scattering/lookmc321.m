% lookmc321.m
clear

fid = fopen('mc321_.out');

for i=1:5
    s = fgetl(fid);
end

A = fscanf(fid,'%f %f %f',[5 inf])';
nm  = A(:,1);
Fsph = A(:,2);
Fcyl = A(:,3);
Fpla = A(:,4);
Fobl = A(:,5);

fclose(fid)


figure('Name','Fluence  Rates Across Symmetrical Surfaces');clf
semilogy(nm,Fsph,'ro-','linewidth',2)
hold on
semilogy(nm,Fcyl,'gd-','linewidth',2)
semilogy(nm,Fpla,'bs-','linewidth',2)
semilogy(nm,Fobl,'ys-','linewidth',2)

set(gca,'fontsize',18)
xlabel('Range, r[cm]')
ylabel('Fluence rate [W/cm^2]')
legend('F_{sph}','F_{cyl}','F_{pla}','F_{obl}')
axis([0 6 1e-5 1e3])
