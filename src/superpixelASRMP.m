%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Similarity Ratio with Mahalanobis distance (SRMP)  %
% with adaptive scheme                               %
%                                                    %
% Emre AKYILMAZ                                      %
%                                                    % 
% This work is a part of thesis titled as            % 
% "SIMILARITY RATIO BASED ALGORITHMS TO GENERATE SAR % 
% SUPERPIXELS"                                       % 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ww, imm, iwm, edges, alphaMat] = superpixelASRMP(fileName, S, numberOfIterations)
im = uint8(imread(fileName)); 
N1 = size(im,1); 
N2 = size(im,2);
im = double(im);
parameterValue = 1.0; %0.5;

%% Initialize variables
nc = ceil(N1/S)*ceil(N2/S);
sp_list=zeros(4, nc); % SP centers: center value (mean intensity later), column, row, number of pixels, Eigen vector
lbl=-ones(size(im)); % image of superpixel labels
d=inf(size(im)); % distance of each pixel to the closest sp

%% Assign cluster centers and form grid structure
% counter for superpixels
kk=1;

% initialize first row's center pixel
r=floor(S/2);

% create the first grid structure
% for all superpixel rows
for ri=1:ceil(N1/S)
    % initialize column center
    c=floor(S/2);
    % for all superpixel columns
    for ci=1:ceil(N2/S)
        sp_list(1:3, kk)=[squeeze(im(r,c));c;r]; % set center value, row, column
        c=c+S;
        kk=kk+1;
        if(c > size(im,2))
            break;
        end
    end
    r=r+S;
    if(r > size(im,1))
        break;
    end
end

% initially all superpixels have only one member
sp_list(4,:) = 1;

for k=1:nc
    % the search window limits
    imin=floor(max(sp_list(3,k)-S/2,1)+0.5); imax=floor(min(sp_list(3,k)+S/2,N1)+0.5);
    jmin=floor(max(sp_list(2,k)-S/2,1)+0.5); jmax=floor(min(sp_list(2,k)+S/2,N2)+0.5);
    lbl(imin:imax, jmin:jmax) = k;
    intensityValues(k) = mean(mean(im(imin:imax, jmin:jmax)));
    sp_list(1,k) = mean(mean(im(imin:imax, jmin:jmax)));
    sp_list(4,k) = (imax-imin)*(jmax-jmin);
end

std_superpixels = std(intensityValues);
mean_superpixels = mean(intensityValues);

%% Main Routine
for n=1:numberOfIterations
    % for each sp
    for k=1:nc
        % the search window limits
        imin=floor(max(sp_list(3,k)-S,1)+0.5); imax=floor(min(sp_list(3,k)+S,N1)+0.5);
        jmin=floor(max(sp_list(2,k)-S,1)+0.5); jmax=floor(min(sp_list(2,k)+S,N2)+0.5);
        
        % search window image and indices
        subim = im(imin:imax, jmin:jmax);
        [x,y] = meshgrid(jmin:jmax, imin:imax);
        alphaMat = ones(size(subim));
        
        % organize search image as rows
        sz = size(subim);
        
        % sp center value
        M1M = double(sp_list(1,k));
        
        % number of members of the first400/ set for radiometric terms
        m1 = (sp_list(4,k));
        
        if((sp_list(2,k) == 0 && sp_list(3,k) == 0) || (sp_list(4,k) == 0))
            continue;
        end
        
        % the image of radiometric distances
        distanceSubim = zeros(sz);
        
        % calculate radiometric distances
        for yy = 1:sz(1)
            for xx = 1:sz(2)
                % find 3x3 pixel neighborhood
                index = 1;
                neighbors = [];
                if(yy-1 > 0 && xx-1 > 0)
                    neighbors(index) = subim(yy-1,xx-1);  % Upper left.  r = row, c = column.
                    index = index + 1;
                end
                if(yy-1 > 0)
                    neighbors(index) = subim(yy-1,xx);    % Upper middle.  r = row, c = column.
                    index = index + 1;
                end
                if(yy-1 > 0 && xx+1 < sz(2))
                    neighbors(index) = subim(yy-1,xx+1);  % Upper right.  r = row, c = column.
                    index = index + 1;
                end
                if(xx-1 > 0)
                    neighbors(index) = subim(yy,xx-1);    % left.  r = row, c = column.
                    index = index + 1;
                end
                if(xx+1 < sz(2))
                    neighbors(index) = subim(yy,xx+1);    % right. r = row, c = column.
                    index = index + 1;
                end
                if(yy+1 < sz(1) && xx+1 < sz(2))
                    neighbors(index) = subim(yy+1,xx+1);  % Lowerleft.  r = row, c = column.
                    index = index + 1;
                end
                if(yy+1 < sz(1))
                    neighbors(index) = subim(yy+1,xx);    % lower middle.  r = row, c = column.
                    index = index + 1;
                end
                if(yy+1 < sz(1) && xx-1 > 0)
                    neighbors(index) = subim(yy+1,xx-1);  % Lower left.  r = row, c = column.
                    index = index + 1;
                end
                neighbors(index) = subim(yy,xx);
                
                % mean of the second set
                M2M = double(mean(neighbors));
                
                % number of members of the first set for radiometric term
                m2 = size(neighbors,2);
                
                % total number of elements
                m3=m1+m2;
                
                if(M1M == 0)
                    M1M = M1M + 0.00001;
                end
                
                if(M2M == 0)
                    M2M = M2M + 0.00001;
                end
                
                % radiometric distance term
                R = m3*log((m1*M1M+m2*M2M)/m3)-(m1*log(M1M)+m2*log(M2M));
                distanceSubim(yy,xx) = R;       

                diff = abs(sp_list(1,k)-sp_list(1, lbl(yy,xx)));
                val1 = mean_superpixels - std_superpixels;
                val2 = mean_superpixels + std_superpixels;
                alphaMat(yy, xx) = 1/(1 + exp(parameterValue*(diff - val1))) + 1/(1 + exp(parameterValue*(-diff + val2)));
            end
        end % end of radiometric distances
        
        % coordinate list of the search window
        coords = [x(:) y(:)]; 

        % coordinate list of the current sp
        %[y2 x2] = find(previouslbl==k);
        [y2 x2] = find(lbl==k);
        coordsSP = [x2(:) y2(:)];
        
        if(size(coordsSP,1) == 0)
            continue;
        end

        % Mahal. dist. of all search window elements to the current sp 
        ds2 = exp((-1)*mahal(coords(:,:), coordsSP(:,:)));
        ds2 = reshape(ds2, [size(subim,1) size(subim,2)]);
        dc2 = distanceSubim;
        
        % Cost Function for Mahalanobis        
        D = dc2 - alphaMat.*ds2;
        
        subd =  d(imin:imax, jmin:jmax);
        subl =  lbl(imin:imax, jmin:jmax);
        updateMask = D < subd;
        subd(updateMask) = D(updateMask);
        subl(updateMask) = k;
        
        d(imin:imax, jmin:jmax) = subd;
        lbl(imin:imax, jmin:jmax) = subl;
    end % for each sp
    
    sp_list(:) = 0;
    for r = 1:N1
        for c = 1:N2
            tmp = [im(r,c); c; r; 1];
            sp_list(:, lbl(r,c)) = sp_list(:, lbl(r,c)) + tmp;
        end
    end
    
    for k = 1:nc
        sp_list(1:3,k) = sp_list(1:3,k)/sp_list(4,k);
    end
end

C1=(sortrows(sp_list'))';

%% Morphological Cleanup
% se=strel('disk', 2);
% mask=zeros(size(lbl));
% for p=1:nc
%     b=lbl==p;
%     mask=mask | b-imopen(b,se);
% end
% 
% [idx,idx]=bwdist(~mask);
% lbl(mask)=lbl(idx(mask));

for k=1:nc
    w(lbl==k)=sp_list(1,k);
end
w=reshape(w,[N1 N2]);
ww=uint8(w);
superpixelResult = ww;

%% Combine Formed Superpixel Result with original Image
h=[1 0 -1;2 0 -2;1 0 -1];
gx=filter2(h,lbl);
gy=filter2(h',lbl);
maskx=(gx.^2+gy.^2)>0;
maskx=bwmorph(maskx,'thin',Inf);
maskx(1,:)=0;maskx(end,:)=0;
maskx(:,1)=0;maskx(:,end)=0;
mask=255*(maskx);
edges = mask;
mm=im+mask;
imm=uint8(mm);
wm=double(w)+mask;
iwm=uint8(wm);

