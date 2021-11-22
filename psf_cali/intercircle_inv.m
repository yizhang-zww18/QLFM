function calcu_phase = intercircle_inv(calcu_dephase,maxIte)
[rr,cc,~] = size(calcu_dephase);
ddx = 1/rr;
ddy = 1/cc;
dfx=squeeze(calcu_dephase(:,:,1));
dfy=squeeze(calcu_dephase(:,:,2));
ra = (rr-1)/2;
[X,Y]=meshgrid([-ra:ra],[-ra:ra]);
dfx( (X.^2+Y.^2)>(ra)^2 )=0;
dfy( (X.^2+Y.^2)>(ra)^2 )=0;
mask = X.^2+Y.^2<=(ra^2);

calcu_phase=zeros(size(calcu_dephase,1),size(calcu_dephase,2));
centerX = (rr+1)/2;
centerY = (rr+1)/2;

Nnum=rr;
for u=1:Nnum
    for v=1:Nnum
        if(mask(u,v)==0)
            continue;
        end
        record = 0;
        count = 0;
        for i = 1:maxIte
            step = zeros(1,u+v-2);
            sr = v-1;
            sd = u-1;
            cs=randperm(length(step));
            step(cs(1:sr)) = 1;
            r2 = 0;
            rr = Nnum;
            cc = Nnum;
            for ids = 1:length(step)
                if(step(ids)==1)
                    r2 = r2-dfy(rr,cc)*ddy;
                    cc = cc-1;
                else
                    r2 = r2-dfx(rr,cc)*ddx;
                    rr = rr-1;
                end
                
            end
            record = record+r2;
        end
        
        calcu_phase(u,v) = record/(maxIte);
        disp(['record = ',num2str(record/maxIte)]);
        disp(['u = ',num2str(u),'; v = ',num2str(v)]);
        
    end
end
%f=f-f(7,7);
calcu_phase( (X.^2+Y.^2)>ra^2 )=0;
end

