%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Similarity Ratio with Euclidean distance (SREP)    %
%                                                    %
% Emre AKYILMAZ                                      %
%                                                    % 
% This work is a part of thesis titled as            % 
% "SIMILARITY RATIO BASED ALGORITHMS TO GENERATE SAR % 
% SUPERPIXELS"                                       % 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ww, imm, iwm, edges] = superpixelSREP(fileName, S, numberOfIterations)
im = uint8(imread(fileName));
N1 = size(im,1); 
N2 = size(im,2);
im = double(im);

nc = (N1/S)*(N2/S);
C=zeros(4,nc);
l=-ones(size(im));
d=inf(size(im));
kk=1;
r=S;

for ri=1:(N1/S)
    c=S;
    for ci=1:(N2/S)
        C(1:3,kk)=[squeeze(im(r,c));c;r];
        c=c+S;
        kk=kk+1;
    end
    r=r+S;
end

C(4,:) = 1;

for n=1:numberOfIterations
    for k=1:nc
        imin=max(C(3,k)-S,1); imax=min(C(3,k)+S,N1);
        jmin=max(C(2,k)-S,1); jmax=min(C(2,k)+S,N2);
        subim = im(imin:imax, jmin:jmax);
        [x,y] = meshgrid(jmin:jmax, imin:imax);
        x = x-C(2,k);
        y = y-C(3,k);
        
        sz = size(subim);
        subim_org = reshape(subim,sz(1)*sz(2),1);
        
        M1M = double(C(1,k));
        m1 = 1;
        
        distanceSubim = zeros(sz);
        for yy = 1:sz(1)
            for xx = 1:sz(2)
                %Find 3x3 pixel neighborhood
                index = 1;
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
                 
                %Calculate distance
                M2M = double(mean(neighbors));
                m2 = size(neighbors,2);
                
                m3=m1+m2;           
                
                if(M1M == 0)
                    M1M = M1M + 0.00001;
                end
                
                if(M2M == 0)
                    M2M = M2M + 0.00001;
                end          
                
                R = m3*log((m1*M1M+m2*M2M)/m3)-(m1*log(M1M)+m2*log(M2M));
                %R2 = abs(1-exp(R));
                                
                distanceSubim(yy,xx) = R;
            end
        end
        
        dc2 = distanceSubim;        
        ds2 = sqrt(x.^2 + y.^2);
        %D = dc2 + 0.5*ds2/S;
        D = dc2 + (0.005)*ds2;
   
        subd =  d(imin:imax, jmin:jmax);
        subl =  l(imin:imax, jmin:jmax);
        updateMask = D < subd;
        subd(updateMask) = D(updateMask);
        subl(updateMask) = k;
        
        d(imin:imax, jmin:jmax) = subd;
        l(imin:imax, jmin:jmax) = subl;
    end

    C(:) = 0;
    for r = 1:N1
        for c = 1:N2
            tmp = [im(r,c); c; r; 1];
            C(:, l(r,c)) = C(:, l(r,c)) + tmp;
        end
    end
    
    for k = 1:nc
        C(1:3,k) = round(C(1:3,k)/C(4,k));
    end
end

t1 = size(C,2);
for i=1:t1
    C(1,i) = i;
end

C1=(sortrows(C'))';

%%%%%%%%Morfolojik temizlik
% se=strel('disk', 2);
% mask=zeros(size(l));
% for p=1:nc
%     b=l==p;
%     mask=mask | b-imopen(b,se);
% end
% 
% [idx,idx]=bwdist(~mask);
% 
% l(mask)=l(idx(mask));

for k=1:nc
    w(l==k)=C(1,k);
end
w=reshape(w,[N1 N2]);
ww=uint8(w);


%%%%%%%%%%%%%%%%%%%%Combine Result with original Image%%%%%%%
h=[1 0 -1;2 0 -2;1 0 -1];
gx=filter2(h,l);
gy=filter2(h',l);
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