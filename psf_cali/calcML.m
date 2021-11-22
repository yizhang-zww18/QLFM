function MLARRAY = calcML(fml, k, x1MLspace, x2MLspace, x1space, x2space)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x1length = length(x1space);%�����x��������
x2length = length(x2space);%�����y��������
x1MLdist = length(x1MLspace);%����΢͸����x��������
x2MLdist = length(x2MLspace);%����΢͸����y��������
x1center = find(x1space==0);%�������
x2center = find(x2space==0);%�������
x1centerALL = [  (x1center: -x1MLdist:1)  (x1center + x1MLdist: x1MLdist :x1length)];%��ÿ��΢͸����Ϊһ�������أ�x��
x1centerALL = sort(x1centerALL);%sort�����Ծ����ÿһ�зֱ������������
x2centerALL = [  (x2center: -x2MLdist:1)  (x2center + x2MLdist: x2MLdist :x2length)];%��ÿ��΢͸����Ϊһ�������أ�y��
x2centerALL = sort(x2centerALL);

zeroline = zeros(1, length(x2space) );
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
patternML = zeros( length(x1MLspace), length(x2MLspace) );%��һ΢͸������λ����
patternMLcp = zeros( length(x1MLspace), length(x2MLspace) );
for a=1:length(x1MLspace),
    for b=1:length(x2MLspace),        
        x1 = x1MLspace(a);
        x2 = x2MLspace(b);
        xL2norm = x1^2 + x2^2;
        

        patternML(a,b) = exp(-i*k/(2*fml)*xL2norm);   %���׹�ʽ4
        patternMLcp(a,b) = exp(-0.05*i*k/(2*fml)*xL2norm);  
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% [xx,yy] = meshgrid(1-round(length(x1MLspace)+1)/2:round(length(x1MLspace)+1)/2-1,1-round(length(x2MLspace)+1)/2:round(length(x2MLspace)+1)/2-1);
% mask = (xx.^2+yy.^2)<=((length(x1MLspace)-1)/2)^2;
% 
% patternML = patternML.*mask;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MLspace = zeros( length(x1space), length(x2space) );
MLcenters = MLspace;
for a=1:length(x1centerALL),
    for b=1:length(x2centerALL),
        MLcenters( x1centerALL(a), x2centerALL(b)) = 1;%������lensletƽ����չ�ɶ��������lenslets
    end
end
MLARRAY = conv2(MLcenters, patternML, 'same');%΢͸�����е���λ����
% MLARRAYcp = conv2(MLcenters, patternMLcp, 'same');
% 
% MLARRAYcpANG = angle(MLARRAYcp);
% MLARRAYcpANG = MLARRAYcpANG - min(min(MLARRAYcpANG)) + 0.0;
% MLARRAYcpANGnorm = MLARRAYcpANG/max(max(MLARRAYcpANG));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%