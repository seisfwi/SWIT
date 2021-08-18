function cmap=cmap_rwb
col(1,:)=[1 0 0];
col(2,:)=[1 1 1];
col(3,:)=[0 0 0];

N=20;
for i=1:3;
    cmap1(1:N,i)=linspace(col(1,i),col(2,i),N);
    cmap2(1:N,i)=linspace(col(2,i),col(3,i),N);
end
cmap=[cmap1;cmap2(2:N,:)];
