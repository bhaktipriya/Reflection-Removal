function out = merge_two_patches(est_t, est_r, h, w, psize)
  % Merge patches and concat to form X=[T R]
  t_merge=merge_patches(est_t,h,w,psize);
  %disp(t_merge);
  r_merge=merge_patches(est_r,h,w,psize);
  %disp(r_merge);
  out=[t_merge(:);r_merge(:)];

