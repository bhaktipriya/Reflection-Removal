function [I_t I_r ] = patch_gmm(I_in, h, w,c,dx,dy,ch)
%Let's construct the Ghosting Kernel
k_mat=construct_kernel(h,w,dx,dy,c); 
%Let's also construct an Identity matrix
I_mat=speye(h*w,h*w); 
%Let's define A as in paper. This maps an image to its 
%ghosted version
A=[I_mat k_mat]; 
%lambda=1/signma^2, sigma =10^-3
lambda=1e6; 
% patch size for patch-GMM
psize=8;
p_el=psize^2;
n_patches=(h-psize+1)*(w-psize+1);
mask=merge_two_patches(ones(p_el, n_patches),ones(psize^2, n_patches), h, w, psize);
%initialize
[I_t_i I_r_i]=grad_irls(I_in,dx,dy,c,h,w,ch);
% Setup for GMM prior
load GSModel_8x8_200_2M_noDC_zeromean.mat
excludeList=[];

% Create patches from the two layers.
est_t=im2patches(I_t_i, psize);
est_r=im2patches(I_r_i, psize);

niter=2;
beta=200;

for i = 1 : niter
  fprintf('Optimizing itn no %d ...\n', i);
    %beta->inf so we have to restrict patches Pix to be equal to auxiliary
    %variables

  % Merge the patches with bounded least squares
  %beta tends to infinity so update z such that z=pi_t
  sum_piX_zi=merge_two_patches(est_t,est_r,h,w,psize);
  z=beta*sum_piX_zi+lambda*A'*I_in(:); 

  % Non-neg. optimization by L-BFGSB
  % opts.factr  Tolerance setting (see this source code for more info)
  % opts.pgtol  Another tolerance setting, relating to norm(gradient,Inf)
  % opts.m      Number of limited-memory vectors to use in the algorithm
  opts=struct('factr',1e4,'pgtol',1e-8,'m',50);
  opts.printEvery=50;
  
  %7b
  sum_zi_2=norm(est_t(:))^2+norm(est_r(:))^2;
   %7a +7b 7b is expanded
    fcn = @(x)( lambda * norm(A*x - I_in(:))^2 + beta*( sum(x.*mask.*x - 2 * x.* sum_piX_zi(:)) + sum_zi_2));
  f_handle = @(x)(lambda * A'*(A*x) + beta*(mask.*x));
  grad = @(x)(2*(f_handle(x) - z));
  fun     = @(x)fminunc_wrapper( x, fcn, grad); 
   % The minimization problem that is solves is:
  % min_x  fun  (x)     subject to   l <= x <= u
  l=zeros(numel(sum_piX_zi),1);
  u=ones(numel(sum_piX_zi),1);
  % l and u are zeros and ones because the we want our optimising image to
  % be in double range.
  [out,~, info]=lbfgsb(fun,l,u,opts);
  % the above stmnt optimizes fun wrt out(which is X)
  % in other words fun has been optimized wrt I_t and I_r
  out=reshape(out,h, w,2);
  I_t=out(:,:,1); 
  I_r=out(:,:,2); 

  % restore the images I_t and I_r using patch based priors to get est_tand
  % est_r
  % Restore patches using the prior
  est_t=im2patches(I_t, psize);
  est_r=im2patches(I_r, psize);
  noiseSD=(1/beta)^0.5;
  [est_t t_cost]=aprxMAPGMM(est_t,psize,noiseSD,[h w],GS,excludeList);
  [est_r r_cost]=aprxMAPGMM(est_r,psize,noiseSD,[h w],GS,excludeList);
    
  %update for half quadratic functions, beta tends to infinity
  beta=beta*2;

end

end

