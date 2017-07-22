function gaussmixm_cart
% create a CART tree for gaussmixm
% 5 parameters
% (1) d = distance between points on x axis [0, Inf]
% (2,3) xd, yd = mean of gaussian normalized by d [0,Inf]
% (4) rv = correlation coefficient [-1, 1]
% (5) vx = 1-vy = variance in x direction (0, 1)
pth=mfilename('fullpath');
nx=100000;
invprec=100; % precision is 1/iprec
rhvi={[0.5 1 2 5 20] 'd'; ...
    [0 0.1 0.3 0.5 1 2 5 20] 'x/d'; ...
    [0 0.1 0.3 0.5 1 2 5 20] 'y/d';  ...
    [-1 -0.9 -0.7 -0.5 0 0.5 1] 'rho';  ...
    [0.05 0.1 0.2 0.3 0.4 0.5 0.8 0.95] 'vx'};
nq=size(rhvi,1); % number of parameters
rhvs=cell(nq,1);
np=zeros(1,nq);
% create axis labels for each parameter
for i=1:nq
    np(i)=length(rhvi{i}); % number of values for parameter i
    rhvs{i}=cell(np(i),1);
    for j=1:np(i)
        rhvs{i}{j}=sprintf('%.2g',rhvi{i}(j));
    end
end
ntr=prod(np);
cin=zeros(ntr,nq); % input to tree
cout=zeros(ntr,1); % output from tree
ic=zeros(1,nq);
for i=1:ntr
    [ic(1),ic(2),ic(3),ic(4),ic(5)]=ind2sub(np,i);
    for j=1:nq
        cin(i,j)=rhvi{j}(ic(j));
    end
    d=cin(i,1); % distance between origins
    vx=cin(i,5);
    vy=1-vx;
    vxy=cin(i,4)*sqrt(vx*vy);
    m=d*cin(i,2:3);
    v=[vx vxy; vxy vy];
    x=randvec(nx,m,v);
    rhos=corr(sqrt([(x(:,1)+0.5*d).^2+x(:,2).^2 (x(:,1)-0.5*d).^2+x(:,2).^2]));
    cout(i)=rhos(2);
end
cart=classregtree(cin,round(cout*invprec)/invprec); % round to required precision
save(pth,'cart');
%% evaluate the tree at the training points
rhe=eval(cart,cin); % now do the inverse
[maxerr,maxeri]=max(abs(rhe(:)-cout(:)));
[id,ixd,iyd,irv,ivx]=ind2sub(np,maxeri);
figure(1);
plot(sort(abs(rhe(:)-cout(:))));
axis([1 ntr 0 maxerr]);
xlabel('Cumulative count');
ylabel('Absolute error');
fprintf('Max error: %.2f (true=%.2f est=%.2f) @\n    d(%d)=%.2g, x/d(%d)=%.2g, y/d(%d)=%.2g, rho(%d)=%.2f, vx(%d)=%.2f, vy=%.2f\n', ...
    rhe(maxeri)-cout(maxeri), cout(maxeri),rhe(maxeri),id,rhvi{1}(id),ixd,rhvi{2}(ixd),iyd,rhvi{3}(iyd),irv,rhvi{4}(irv),ivx,rhvi{5}(ivx),1-rhvi{5}(ivx));
%% plot a series of graphs
pli=[0 -1 -2 1 5]; % 0 for graph series, -1 for x, -2 fo y
[pld,plj]=sort(pli);
rhe=reshape(rhe,np);
cout=reshape(cout,np);
% create the axis tick labels
xtick=1:np(plj(1));
xtl=rhvs{plj(1)};
ytick=1:np(plj(2));
ytl=rhvs{plj(2)};
for i=1:np(plj(3))
    j=[pli;pli];
    j(:,plj(3))=[i;i];
    j(:,plj(1:2))=[ones(1,2);np(plj(1:2))];
    rhp=permute(rhe(j(1,1):j(2,1),j(1,2):j(2,2),j(1,3):j(2,3),j(1,4):j(2,4),j(1,5):j(2,5)),plj);
    rht=permute(cout(j(1,1):j(2,1),j(1,2):j(2,2),j(1,3):j(2,3),j(1,4):j(2,4),j(1,5):j(2,5)),plj);
    maxerp=max(abs(rhp(:)-rht(:)));
    figure(100+i);

    for isp=1:2
        subplot(1,2,isp);
        if isp==1
            imagesc(rht);
            title(sprintf('True (%s=%.2g,%s=%.2g,%s=%.2g)', ...
                rhvi{plj(3),2},rhvi{plj(3),1}(i),rhvi{plj(4),2},rhvi{plj(4),1}(pld(4)),rhvi{plj(5),2},rhvi{plj(5),1}(pld(5))));
        else
            imagesc(rhp);
            title(sprintf('Est (Max error=%.2f)',maxerp));
        end
        v_colormap('v_bipliny');
        set(gca,'clim',[-1 1],'xtick',xtick,'xticklabel',xtl,'ytick',ytick,'yticklabel',ytl);
        axis('xy');
        xlabel(rhvi{plj(2),2});
        ylabel(rhvi{plj(1),2});
        %         title(sprintf('d=%.2g, rv=%.2g, vx=%.2g, vy=%.2g',ds(id),rvs(irv),vxs(ivx),1-vxs(ivx)));
    end
end
tilefigs
