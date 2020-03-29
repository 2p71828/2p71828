%% Read habenula ROI mask and QSM data
%% Then look into anterior-half and posterior-half of ROI
% ==> AP contrast
%% Then we may look into 45 degree-adjacent region and compare chi mean in
% (1) anterior half, (2) posterior half, (3) adjacent region.


%% Read mask and QSM, for middle slice display
rt1 = 'W:\seungkyun_lab_public\collaboration\MtSinai_JWKim\from JWKim\to_skku';
fn1 = dir([rt1,filesep,'**\ha03_LR.mat']);      % 48 = 21 + 3*(10-1)
rt2 = rt1;
fn2 = dir([rt2,filesep,'**\qsm.mat']);      % check: 48x1 struct
tic;
h = waitbar(0);
for j=1:length(fn2)             % length(fn1) = length(fn2)
    waitbar(j/length(fn2),h);
    fn_m = [fn1(j).folder,filesep,fn1(j).name]; % mask file name
    fn_q = [fn2(j).folder,filesep,fn2(j).name]; % qsm file name
    load(fn_m,'ha03_L','ha03_R');
    load(fn_q,'q');
    ha03_L = logical(ha03_L);
    ha03_R = logical(ha03_R);
    [T,i1,i2,i3]=trim3d(ha03_L);
    [S,j1,j2,j3]=trim3d(ha03_R);
    Li(j,:)=[i1,i2,i3];             % left mask trim index
    Ri(j,:)=[j1,j2,j3];             % right mask trim index
    L2mid = ceil( (i2(1)+i2(2))/2 );    % if (1,2,3)==>2, (1,2,3,4)==>3.
                                        % we can make this mid number part
                                        % of the latter half.
    LP = ha03_L;    % Posterior is the earlier one in i2 direction
    LP(:,L2mid:end,:) = 0;
    LA = ha03_L;
    LA(:,1:L2mid-1,:) = 0;
    R2mid = ceil( (j2(1)+j2(2))/2 );
    RP = ha03_R;
    RP(:,R2mid:end,:) = 0;
    RA = ha03_R;
    RA(:,1:R2mid-1,:) = 0;
    dat(j,:) = [mean(q(LP)),mean(q(LA)),mean(q(RP)),mean(q(RA)),...
                mean(q(LP|RP)),mean(q(LA|RA)),...    
                mean(q(ha03_L)),mean(q(ha03_R))];   % main data
    vat(j,:) = [sum(LP(:)),sum(LA(:)),sum(RP(:)),sum(RA(:)),...
                    sum(LP(:)+RP(:)),sum(LA(:)+RA(:)),...
                    sum(ha03_L(:)),sum(ha03_R(:))];     % voxel count
    
    smin = min(i3(1),j3(1));
    smax = max(i3(2),j3(2));
    s = round((smin+smax)/2);           % middle slice
    qat(:,:,j) = q(:,:,s);              % qsm data
    mat(:,:,j) = ha03_L(:,:,s) + ha03_R(:,:,s); % mask data
                
    st = fn2(j).folder;         % folder label (chararcters after last \)
    loc = find(st=='\');
    lab = st(loc(end)+1:end);
    lat{j,1} = j;               % label data
    lat{j,2} = lab;
end
close(h);
toc; % 224 sec

%% Save for now
outname = 'data.mat'
save(outname,'dat','vat','Li','Ri','lat','qat','mat');
disp([outname,' < dat, vat, Li, Ri, lat, qat, mat']);

%% Show for individual scans
[T,i1,i2,i3]=trim3d(mat);
r1 = [i1(1)-5:i1(2)+5]; % display range
r2 = [i2(1)-5:i2(2)+5];
qdis = qat(r1,r2,:);
mdis = mat(r1,r2,:);
j = 1
a = qdis(:,:,j);
figure;imagesc(a,[-0.15,0.15]);axis image;colormap gray;hold on;
title([num2str(j),'.',lat{j,2}]);
B = bwboundaries_out(mdis(:,:,j),8);
for k=1:length(B)
    plot(B{k}(:,2),B{k}(:,1),'r');
end
    
% SliceBrowser(pdis,-0.15,0.15);colormap gray;

%% collecting all in one figure (subplot-based)
%{
ah = subpl_pos
6
8
48
[1:48]
colormap gray
set(ah,'fontsize',8)
set(ah,'xtick',[],'ytick',[])
set(gcf,'posi',[-1800,10,1300,950]); % -1800 if using 2nd screen on left.
%}





