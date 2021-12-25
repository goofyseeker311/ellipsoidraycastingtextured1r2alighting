output_precision = 32; format long; pkg load statistics; pkg load matgeom;

function [k,n,v,u] = raysphereintersection(sphere, ray)
  spherecount = size(sphere,1); raycount = size(ray,1); kmatrix = zeros(raycount,2);
  for m=1:spherecount
    sw = (sphere(m,5:7).^2);
    aa = sum( sw.*ray(:,4:6).^2 ,2);
    bb = 2*sum( sw.*ray(:,4:6).*(ray(:,1:3)-sphere(m,1:3)) ,2);
    cc = sum( sw.*((ray(:,1:3)-sphere(m,1:3)).^2) ,2) - sphere(m,4).^2;
    dd = sqrt(bb.^2 - 4*aa.*cc);
    kv1 = (-bb + dd) ./ (2*aa); k1 = [((((abs(imag(kv1))==0).*kv1)>0).*kv1) ones(raycount,1)*m];
    kv2 = (-bb - dd) ./ (2*aa); k2 = [((((abs(imag(kv2))==0).*kv2)>0).*kv2) ones(raycount,1)*m];
    k1g0 = (k1(:,1)>0); k2g0 = (k2(:,1)>0); kmle0 = (kmatrix(:,1)<=0); kmg0 = !kmle0; k1gkm = (k1(:,1)>kmatrix(:,1)); k2gkm = (k2(:,1)>kmatrix(:,1));
    kmatrix = (kmle0.*k1g0.*k1).+kmg0.*(k1g0.*(((!k1gkm).*k1).+(k1gkm.*kmatrix)).+(!k1g0).*kmatrix);
    kmatrix = (kmle0.*k2g0.*k2).+kmg0.*(k2g0.*(((!k2gkm).*k2).+(k2gkm.*kmatrix)).+(!k2g0).*kmatrix);
  endfor
  nmatrix = kmatrix(:,1).*ray(:,4:6)+ray(:,1:3); uvi=find(kmatrix(:,2)>0); uvn=kmatrix(uvi,2); uvp=nmatrix(uvi,:);
  sn=sphere(uvn,1:3);uvv=uvp.-sn;uvd=sqrt(sum(uvv.^2,2))*ones(1,3); vmatrix=zeros(raycount,3);vmatrix(uvi,:)=uvv./uvd;
  umatrix=zeros(raycount,2); ratio=(1/360); ratio2=(1/180); uvicount=size(uvi,1);
  umatrix(uvi,1)=0.5+ratio.*(acosd(dot(vmatrix(uvi,1:2),ones(uvicount,1)*[1 0],2)./sqrt(sum(vmatrix(uvi,1:2).^2,2))).*((-1).^(vmatrix(uvi,2)<0)));
  umatrix(uvi,2)=1-(0.5+ratio2.*asind(vmatrix(uvi,3)));
  k = kmatrix; n = nmatrix; v=vmatrix; u = umatrix;
endfunction

wimagefiles={"physical-world-map.jpg","Fondvulcain.jpg","world116cE.jpg", "galaxy2.jpg", "skygalaxy-hdri-114148912.jpg", "voie-lactee-hd-eso-l3.jpg", "bush_texture.jpg", "fir_tree_PNG5661.png", "pine-tree-texture-3.jpg", "tree_texture.jpg"};
wimagescount=size(wimagefiles,2); wimages={}; wimagesizes=zeros(wimagescount,3);
for u=1:wimagescount; wimages{u}=double(imread(wimagefiles{u}))/255; wimagesizes(u,:)=size(wimages{u}); endfor;
lightvalue=1; lightvalue2=0.04; minlightdist=1.1; pointlightscount = 5;
pointlights = [100*(rand(pointlightscount,3)-0.5) randi(255,pointlightscount,3)/255 -2*ones(pointlightscount,1)];
pointlightscount+=1; pointlights=[pointlights; [10000 10000 10000 1 1 1 0]];
nspheres = 100; scenespheres = [100*(rand(nspheres,3)-0.5) 6*rand(nspheres,1)+4 2*rand(nspheres,3)+0.5 randi(255,nspheres,3)/255 randi(wimagescount+1,nspheres,1)-1];
scenespheres = [scenespheres; [pointlights(:,1:3) ones(pointlightscount,1) ones(pointlightscount,3) pointlights(:,4:6) zeros(pointlightscount,1)]]; nspheres+=pointlightscount;
cameraorigin = [0 0 0]; cv=[cameraorigin 0 0 0]; cameravector = cv + createLine3d(pi/2,0); cameraupvector = cv + createLine3d(0,0); camerarightvector = cv + createLine3d(pi/2,-pi/2);
spherebufferwidth = 1*3840; spherebufferheight = 1*1080; numbytes = 4;
raycount = spherebufferwidth*spherebufferheight; horizontalstep = 360/(spherebufferwidth-1); verticalstep = 180/(spherebufferheight-1);
camerazbuffersphere = zeros(spherebufferwidth, spherebufferheight); camerargbabuffersphere = zeros(spherebufferwidth, spherebufferheight,numbytes);
ratio = pi/180; thetasteps = flip(-180:horizontalstep:180,2)*ratio; phisteps = (0:verticalstep:180)*ratio; rays = zeros(raycount,6); raycounter = 1;
for theta=thetasteps; for phi=phisteps; rays(raycounter,:) = createLine3d(phi,theta); raycounter++; endfor; endfor; trays = (ones(raycount,1)*cv) .+ rays;

[k, kn, kv, ku] = raysphereintersection(scenespheres(:,1:7), trays); kd=k(:,1);kb=k(:,2); m = max(kd);
kc=zeros(raycount,3);kci=find(kb);kbi=kb(kci);kc(kci,:)=scenespheres(kbi,8:10); kt=scenespheres(kbi,11); kx=kn(kci,1);ky=kn(kci,2);kz=kn(kci,3);kcicount=size(kci,1);knc=kn(kci,:);
kc1=reshape(kc(:,1),spherebufferheight,spherebufferwidth); kc2=reshape(kc(:,2),spherebufferheight,spherebufferwidth); kc3=reshape(kc(:,3),spherebufferheight,spherebufferwidth);
kuv=ku(kci,:); kut=ones(raycount,3); for m=1:wimagescount; kti=find(kt==m); kcti=kci(kti); spx=round(kuv(kti,2).*(wimagesizes(m,1)-1))+1; spy=round(kuv(kti,1).*(wimagesizes(m,2)-1))+1;
wimage1=wimages{m}(:,:,1);wimage2=wimages{m}(:,:,2);wimage3=wimages{m}(:,:,3); flati=wimagesizes(m,1)*(spy-1)+spx; wmst=wimagesizes(m,1)*wimagesizes(m,2);
kut(kcti,1)=reshape(wimage1,wmst,1)(flati);kut(kcti,2)=reshape(wimage2,wmst,1)(flati);kut(kcti,3)=reshape(wimage3,wmst,1)(flati); endfor;
kut1=reshape(kut(:,1),spherebufferheight,spherebufferwidth);kut2=reshape(kut(:,2),spherebufferheight,spherebufferwidth);kut3=reshape(kut(:,3),spherebufferheight,spherebufferwidth);
renderim=zeros(spherebufferheight,spherebufferwidth,3); renderim(:,:,1) = kc1.*kut1; renderim(:,:,2) = kc2.*kut2; renderim(:,:,3) = kc3.*kut3;
figure(2); clf; imshow(renderim, "xdata", [-180 180], "ydata", [-90 90]); axis on; figure(3); hist(nonzeros(kd),30); imwrite(renderim,"testCCCA0.png");

[X,Y,Z] = sphere(); figure(1); clf; view(3); axis([-65 65 -65 65 -65 65]); hold on;
plot3(cameravector(1),cameravector(2),cameravector(3),'ro'); plot3(cameravector(1)+10*cameravector(4),cameravector(2)+10*cameravector(5),cameravector(3)+10*cameravector(6),'bx');
plot3([cameravector(1) cameravector(1)+10*cameravector(4)],[cameravector(2) cameravector(2)+10*cameravector(5)],[cameravector(3) cameravector(3)+10*cameravector(6)],'g');
plot3([cameravector(1) cameraupvector(1)+30*cameraupvector(4)],[cameravector(2) cameraupvector(2)+30*cameraupvector(5)],[cameravector(3) cameraupvector(3)+30*cameraupvector(6)],'b');
plot3([cameravector(1) camerarightvector(1)+20*camerarightvector(4)],[cameravector(2) camerarightvector(2)+20*camerarightvector(5)],[cameravector(3) camerarightvector(3)+20*camerarightvector(6)],'r');
plot3(cameraupvector(1)+30*cameraupvector(4),cameraupvector(2)+30*cameraupvector(5),cameraupvector(3)+30*cameraupvector(6),'kx');
plot3(camerarightvector(1)+20*camerarightvector(4),camerarightvector(2)+20*camerarightvector(5),camerarightvector(3)+20*camerarightvector(6),'bx');
for u=1:nspheres; a=scenespheres(u,4);ax=scenespheres(u,5);ay=scenespheres(u,6);az=scenespheres(u,7);
x=scenespheres(u,1);y=scenespheres(u,2);z=scenespheres(u,3); surf((a/ax)*X+x,(a/ay)*Y+y,(a/az)*Z+z); endfor;
plot3(kx,ky,kz,'y.'); plot3(pointlights(:,1),pointlights(:,2),pointlights(:,3),'rx'); hold off;

kl=zeros(raycount,3);
for v=1:pointlightscount;
  printf("v/plcount: %i/%i\n",v,pointlightscount);
  pla=(ones(kcicount,1)*pointlights(v,1:3)); pllc=pointlights(v,7);
  plv=pla.-knc; pld=sqrt(sum(plv.^2,2)); pld3=pld*ones(1,3); plnv=plv./pld3; lv=[knc.+(minlightdist.*plnv) plnv];
  [ck, ckn, ckv] = raysphereintersection(scenespheres(:,1:7), lv); ckl=(ck(:,2)==(nspheres-pointlightscount+v));
  cki=find(ckl); ckir=find(!ckl); ckd=ck(cki,1); ckdr=ckd.^pllc; ckicount=size(cki,1);ckircount=size(ckir,1);
  ckni=ckn(cki,:);cknir=ckn(ckir,:); kcni=kci(cki);kcnir=kci(ckir); pld3n=pld3(cki,:);
  plc=ones(ckicount,1)*pointlights(v,4:6); kl(kcni,:).+=(lightvalue.*ckdr.*plc); plcr=ones(ckircount,1)*pointlights(v,4:6); kl(kcnir,:)+=(lightvalue2.*plcr);
  figure(1); hold on; plot3(ckni(:,1),ckni(:,2),ckni(:,3),'r.'); hold off;
endfor;
for v=1:pointlightscount; ktz=find(kb==(nspheres-pointlightscount+v)); kl(ktz,:)=1; endfor;
kl1=reshape(kl(:,1),spherebufferheight,spherebufferwidth); kl2=reshape(kl(:,2),spherebufferheight,spherebufferwidth); kl3=reshape(kl(:,3),spherebufferheight,spherebufferwidth);
renderim2=zeros(spherebufferheight,spherebufferwidth,3); renderim2(:,:,1) = kc1.*kl1.*kut1; renderim2(:,:,2) = kc2.*kl2.*kut2; renderim2(:,:,3) = kc3.*kl3.*kut3;
figure(4); clf; view(2); hold on; imshow(renderim2, "xdata", [-180 180], "ydata", [-90 90]); hold off; axis on; imwrite(renderim2,"testCCCB0.png");

saveas(1,"raster0.png");
saveas(2,"caster0.png");
saveas(3,"hist0.png");
saveas(4,"caster0a.png");
