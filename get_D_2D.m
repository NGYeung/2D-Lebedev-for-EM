function [DaH,DHa,DbH,DHb,DEH,DHE,WH,WE,WH_sq,WE_sq,WH_sq_inv,WE_sq_inv]=get_D_2D(sizes)
%function [D,W,Dx2pol,Dy2pol,Dpol2y,Dpol2x]=get_D_2D(Input)
%
% This function generates the discetized spatial differntiation operator
% matrix D for a 2D field problem for 2D Lebedev in H polarization
%
% Inputs:       
% sizes     - structure
%
% Outputs:
% Da2b      - Differentiation matrix from a to b 
%             Example: Dx2pol differentiation from the x fielt to the
%             z/polarized field
%
%
%
% Date: 17-2-2022



delta=sizes.delta;

Nx=sizes.Nx;
Ny=sizes.Ny;

% homogenious stepsize
% Dirichlet boundary condition of tangetial and normal E at boundary

% this right now does not make any sense!!!!
% It should be deltahat_x_vec =delta*ones(Nx+1,1); and the 
% differentiation matrices need to be defined cleanly!!!!

    delta_x_vec    =delta*ones(Nx-1,1);
    deltahat_x_vec =delta*ones(Nx,1);

    delta_y_vec    =delta*ones(Ny-1,1);
    deltahat_y_vec =delta*ones(Ny,1);

    delta_x_vec_tilde=[delta/2;delta_x_vec;delta/2];
    delta_y_vec_tilde=[delta/2;delta_y_vec;delta/2];
    
% Primary and dual differentiation matrix on (00) cluster
    Xhat=spdiags([-1./deltahat_x_vec(2:Nx) 1./deltahat_x_vec(1:Nx-1)],[-1:0],Nx,Nx-1);
    X=spdiags([-1./delta_x_vec 1./delta_x_vec],[0:1],Nx-1,Nx);

    Yhat=spdiags([-1./deltahat_y_vec(2:Ny) 1./deltahat_y_vec(1:Ny-1)],[-1:0],Ny,Ny-1);
    Y=spdiags([-1./delta_y_vec 1./delta_y_vec],[0:1],Ny-1,Ny);

    tX=spdiags([-1./delta_x_vec_tilde(2:(Nx+1)) 1./delta_x_vec_tilde(1:Nx)],[-1:0],Nx+1,Nx);
	tY=spdiags([-1./delta_y_vec_tilde(2:(Ny+1)) 1./delta_y_vec_tilde(1:Ny)],[-1:0],Ny+1,Ny);
    
	tXhat=spdiags([-1./deltahat_x_vec 1./deltahat_x_vec],[0:1],Nx,Nx+1);
	tYhat=spdiags([-1./deltahat_y_vec 1./deltahat_y_vec],[0:1],Ny,Ny+1);
    
    
    

    


%     Wx=diag(kron( spdiags(deltahat_y_vec,0,Ny,Ny) , spdiags(delta_x_vec   ,0,Nx-1  ,Nx-1  ) ));
%     Wz=diag(kron( spdiags(delta_y_vec   ,0,Ny-1  ,Ny-1  ) , spdiags(delta_x_vec   ,0,Nx-1  ,Nx-1  ) ));
%     Wy=diag(kron( spdiags(delta_y_vec   ,0,Ny-1  ,Ny-1  ) , spdiags(deltahat_x_vec,0,Nx,Nx) ));

%     D=[ sparse(Nx*(Ny-1),Nx*(Ny-1)) -Dpol2x            sparse(Nx*(Ny-1),Ny*(Nx-1)) ; ...
%         -Dx2pol                     sparse(Nx*Ny,Nx*Ny) Dy2pol                      ; ...
%         sparse(Ny*(Nx-1),Nx*(Ny-1)) Dpol2y              sparse(Ny*(Nx-1),Ny*(Nx-1))];
%stepsize Matrix
%     W=spdiags([Wx;Wz;Wy],0,Nk,Nk);

    Dx2pol=kron(Y,speye(Nx));
    Dy2pol=kron(speye(Ny),X);
    Dpol2y=kron(speye(Ny),Xhat);
    Dpol2x=kron(Yhat,speye(Nx));
% Primary and dual differentiation matrix on (11) cluster
tDx2pol=kron(tYhat,speye(Nx+1));
tDy2pol=kron(speye(Ny+1),tXhat);
tDpol2y=kron(speye(Ny+1),tX);
tDpol2x=kron(tY,speye(Nx+1));

DaH=[-Dpol2x                        sparse(sizes.hz1p,sizes.ey2p);...
      sparse(sizes.hz2p,sizes.ex1p) tDpol2y];
DHa=[-Dx2pol                        sparse(sizes.ex1p,sizes.hz2p);...
      sparse(sizes.ey2p,sizes.hz1p) tDy2pol];
DbH=[ sparse(sizes.hz1p,sizes.ex2p) Dpol2y ;...
     -tDpol2x                       sparse(sizes.hz2p,sizes.ey1p)];

DHb=[sparse(sizes.ex2p,sizes.hz1p)  -tDx2pol;...
     Dy2pol                          sparse(sizes.ey1p,sizes.hz2p)];

DEH=[DaH  DbH];
DHE=[DHa ; DHb];


Wx          =spdiags(delta_x_vec       ,0,Nx-1,Nx-1);
Wxhat       =spdiags(deltahat_x_vec    ,0,Nx  ,Nx);
Wxtilde     =spdiags(delta_x_vec_tilde ,0,Nx+1,Nx+1);

Wy          =spdiags(delta_y_vec       ,0,Ny-1,Ny-1);
Wyhat       =spdiags(deltahat_y_vec    ,0,Ny  ,Ny);
Wytilde     =spdiags(delta_y_vec_tilde ,0,Ny+1,Ny+1);
    

WH=[kron(Wyhat,Wxhat)            , sparse(sizes.hz1p,sizes.hz2p);...
    sparse(sizes.hz2p,sizes.hz1p),kron(Wytilde,Wxtilde)];

WH_sq=spdiags(sqrt(diag(WH)),0,sizes.NH,sizes.NH);
WH_sq_inv=spdiags(1./sqrt(diag(WH)),0,sizes.NH,sizes.NH);
%diagonal of WE
ee=[diag(kron(Wy,Wxhat));...
    diag(kron(Wytilde,Wxhat));...
    diag(kron(Wyhat,Wxtilde));...
    diag(kron(Wyhat,Wx))];

WE=spdiags(ee,0,sizes.NE,sizes.NE);
WE_sq=spdiags(sqrt(ee),0,sizes.NE,sizes.NE);
WE_sq_inv=spdiags(1./sqrt(ee),0,sizes.NE,sizes.NE);

end
