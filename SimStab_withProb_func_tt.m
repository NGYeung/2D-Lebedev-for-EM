function [R,R0,nextloop]=SimStab_withProb_func_tt(D,D0,src,NoiseMethod)
% I added stabilizatoin for M0 here, whcih it should not need in theory.


%% Construct
    nt=size(D,3);%timingST.nt;
    n=src.nsrc;
    %ht=timingST.ht;
    %construct matricees
   

statOut=[];
thld=1e-11; % threshold of stability, conditioning 1e-11
nextloop = 0;

if(NoiseMethod==1)    
    unstable=true;
    mu=1.00-0.01;
    %reg = 0;
    while unstable
        mu=mu+0.01; % increase multiplication factor of zero data
        D_tr=D;
        D0_tr=D0;

        D_tr(:,:,1)=mu*D(:,:,1);
        D0_tr(:,:,1)=mu*D0(:,:,1);

        
        UTU=mutu_vec(D_tr, n, nt);
        
        
        
        
        %UTU=0.5*(UTU+UTU');
        
        
        %UTU = UTU + reg.*eye(size(UTU,1)); 
        
        msigmat = eig(UTU,'vector');
 
        COND = min(msigmat)/msigmat(end);

        if(isempty(find( (msigmat./msigmat(end))<thld,1) )) % chek conditioning
        
        
            unstable=false;
        end
        
        
        if mu >= 5% just to make sure it doesn't take forever %changed for debugging!!!
            unstable=false;
            nextloop = 1;
        end
        
    end
    
   % disp(['mu is ' num2str(mu)])
    statOut=mu;
    D=D_tr;
    D0=D0_tr;
   


    M=mutu_vec(D,  src.nsrc, nt);
    M0=mutu_vec(D0,  src.nsrc, nt);


    %remember to unquote!!!!
    %M=0.5*(M+M');
    %M0=0.5*(M0+M0');
   
    mm=n*nt/2;

    Lorth_bg = mblockchol( M, n, nt/2);
    R=Lorth_bg';



    Lorth_bg = mblockchol( M0, n, nt/2);
    R0=Lorth_bg';

end