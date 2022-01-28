% The I-function generalizes the Fox H-function implemented in [R1].
%[R1] Abdelaziz Soulimani, Mustapha Benjillali, Hatim Chergui, Daniel B. da Costa,%Multihop Weibull-fading communications: Performance analysis framework and applications,
%Journal of the Franklin Institute, Volume 358, Issue 15, 2021, Pages 8012-8044,
function out = I_Func(an, An, AAn,ap, Ap, AAp, bm, Bm , BBm, bq, Bq, BBq, z,v)
%vc = 0,1,2... to ensure convergence
F = @(s)(GammaProd(bm,Bm,BBm,s).* GammaProd(1-an,-An,AAn,s).* z.^-s )./ (GammaProd(1-bq,-Bq,BBq,s).* GammaProd(ap,Ap,AAp,s));
%% Contour preparation:
% epsilon = 10^1.2;
% Sups = min((1-an)./An); Infs = max(-bm./Bm);
% if(isempty(Sups) && isempty(Infs))
% WPx=1;
% elseif(isempty(Sups) && ~isempty(Infs))
% WPx = Infs +epsilon;
% elseif(~isempty(Sups) && isempty(Infs))
% WPx = Sups -epsilon;
% else
% WPx = (Sups + Infs)/1.35;% s between Sups and Infs
% end
% %% integration:
% infity = 10;
% out =real( (1/(2*1i*pi))*integral(F,WPx-1i*infity, WPx+1i*infity));
% return
%% Parameters:
p = length([An Ap]);
q = length([Bm Bq]);%length([An Ap]);%
alphaFox = sum(An)-sum(Ap)+sum(Bm)-sum(Bq);
mu = sum([Bm Bq]) - sum([An Ap]);
betaFox = prod([An Ap].^-[An Ap])* prod([Bm Bq].^-[Bm Bq]);
delta = sum([bm bq]) - sum([an ap]) + (p-q)/2 ;
Tol = 10^-5;
%% Conditions per contour:
% Contour L_(c+i*infinity):
condition01= alphaFox>0 && abs(angle(z))<pi*alphaFox/2;
condition02=alphaFox==0&&(delta*mu+real(delta))<-1 && angle(z)==0;
condition0 = condition01 || condition02;
% contour L_(-infinity)
condition11 = (mu>0)&& z~=0;
condition12 = (mu==0) && abs(z)<betaFox && abs(z)>0;
condition13=(mu==0)&&abs(z)==betaFox &&real(delta)<-1;
condition1 = condition11||condition12 || condition13;
% contour L_(+infinity)
condition21 = (mu<0)&& z~=0;
condition22 = (mu==0) && abs(z)>betaFox;
condition2 = condition21 || condition22;
%% Contour preparation:
epsilon = 10^1.2;
Sups = (min((v+1-an)./An)); Infs = max(-(bm+v)./Bm);
if(isempty(Sups) && isempty(Infs))
    WPx=1;
elseif(isempty(Sups) && ~isempty(Infs))
    WPx = Infs +epsilon;
elseif(~isempty(Sups) && isempty(Infs))
    WPx = Sups -epsilon;
else
    WPx = (Sups+ Infs)/2 %WPx = (Sups+ Infs)/3.1; for OP cap  light WPx = (Sups+ Infs)/1.741  % mod thick cap WPx = (Sups+ Infs)/1.9741
    % BER light WPx = (Sups+ Infs)/12 % OP mod. WPx =0.1250
end
WayPoints = [WPx-1i*epsilon WPx+1i*epsilon];
  %WPx=0.1;
%% integration:
if(condition0 || (~condition1 && ~condition2))
    infity = 10;
    out = real((1/(2*1i*pi))*integral(F,WPx-1i*infity,WPx+1i*infity,'AbsTol' , Tol , 'RelTol' , Tol));
    return
end
if(condition1)
    infity = 100;
    if(~isempty(min(-bm./Bm)))
        infity = infity - min(-bm./Bm);
    end
    out = real((1/(2*1i*pi))*integral(F,-infity,-infity, 'Waypoints',WayPoints));
    return
end
if(condition2)
    infity = 100;
    if(~isempty(max((1-an)./An)))
        infity = infity + max((1-an)./An);
    end
    %   Tol = 10^-5;
    out = real((1/(2*1i*pi))*integral(F,infity,infity,'Waypoints',WayPoints, 'AbsTol' ,Tol, 'RelTol' ,Tol));
end
return
    function output = GammaProd(p,x,y,X)
        [pp, XX] = meshgrid(p,X);
        xx = meshgrid(x,X);
        if (isempty(p))
            output = ones(size(X));
        else
            output = reshape(prod(double(gammaz(pp+xx.*XX).^y),2),size(X));
        end
    end
end
