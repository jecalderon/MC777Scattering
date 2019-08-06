% lookmc321.m
clear

fid = fopen('mc321_.out');

for i=1:4
    s = fgetl(fid);
end

A = fscanf(fid,'%f %f %f',[4 inf])';
nm  = A(:,1);
Fsph = A(:,2);
Fcyl = A(:,3);
Fpla = A(:,4);

fclose(fid)


figure(1);clf
semilogy(nm,Fsph,'ro-','linewidth',2)
hold on
semilogy(nm,Fcyl,'gd-','linewidth',2)
semilogy(nm,Fpla,'bs-','linewidth',2)

set(gca,'fontsize',18)
xlabel('Range, \r[cm]')
ylabel('Fluence rate [W/cm^2]')
legend('F_{sph}','F_{cyl}','F_{pla}',1)
axis([0 3.1 1e-4 1e3])
