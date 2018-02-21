module updates_baypass

use mt_stream
use mcmc_utils

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!calcul de la log_like_beta (a partir de la freq latente qui peut etre negative!)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function logl_ij_beta(Y_obs,N_obs,A_ij,LOG_CN)
 integer, intent(in) :: Y_obs,N_obs
 real (kind=8), intent(in) :: A_ij,LOG_CN
 real (kind=8) :: logl_ij_beta,eps=1e-16

 logl_ij_beta=0. !les cas impossible (a=0 et Y=N on ete ecarte dans la chaine lors de l'update des a)

 if(A_ij> eps .and. (1.-A_ij)>eps) then
   if(Y_obs==0) logl_ij_beta=(N_obs+0.)*log(1.-A_ij)
   if(N_obs==Y_obs) logl_ij_beta=(N_obs+0.)*log(A_ij)
   if(N_obs/=Y_obs .and. Y_obs/=0) logl_ij_beta=LOG_CN + (Y_obs+0.)*log(A_ij) + (N_obs - Y_obs + 0.)*log(1.-A_ij)
 end if

end function logl_ij_beta

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!calcul de la deviance
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function deviance(Y_ij,N_ij,A_ij,LOG_CN)
 integer, intent(in), dimension(:,:) :: Y_ij,N_ij
 real (kind=8), intent(in), dimension(:,:) :: A_ij,LOG_CN
 real (kind=8) :: deviance
 integer :: tmp_i, tmp_j,tmp_npop,tmp_nmrk

 tmp_npop=size(A_ij,2) ; tmp_nmrk=size(A_ij,1) ; deviance=0.
 do tmp_i=1,tmp_nmrk
  do tmp_j=1, tmp_npop
    deviance=deviance + logl_ij_beta(Y_ij(tmp_i,tmp_j),N_ij(tmp_i,tmp_j),A_ij(tmp_i,tmp_j),LOG_CN(tmp_i,tmp_j))
  end do
 end do
 deviance=-2.*deviance
end function deviance

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!perason coef
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function rho_pearson(vect_1,vect_2)
implicit none
 real (kind=8), intent(in), dimension(:) :: vect_1,vect_2
 real (kind=8) :: rho_pearson
 real (kind=8) :: m1,s1,m2,s2,s12
 integer :: dim_vect
 
 dim_vect=size(vect_1)
 m1=sum(vect_1)/dim_vect ; m2=sum(vect_2)/dim_vect
 s1=sqrt(sum((vect_1-m1)**2)/dim_vect)
 s2=sqrt(sum((vect_2-m2)**2)/dim_vect)
 s12=sum((vect_1-m1)*(vect_2-m2))/dim_vect
 rho_pearson=s12/(s1*s2)
 
end function rho_pearson
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!update des y
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine update_y (this,y_cur,popsize,y_lect,n_lect,y_minmax,alpha,delta_y,accept,y_out) 
implicit none
!! On fait du MH
!real, intent(in),dimension(:,:) :: log_pr_init !Y_l*log(y/n) + (N_l-Y-l)*log(1-y/n) + log(Cyn) ! en inout pour pouvoir le modifier
type(mt_state), intent(inout) :: this
real (kind=8), intent(in) :: alpha
integer,intent(in) :: y_lect,n_lect,y_cur,popsize,delta_y,y_minmax(2)

integer, intent(out) :: y_out,accept

real (kind=8) :: dum_real,rnd,diff_log,tmp_logcur,tmp_lognew!,inv_bwd,inv_fwd
integer :: b_inf,b_sup,ymax,ymin

tmp_logcur=0. ; tmp_lognew=0. 
ymin=y_minmax(1) ; ymax=y_minmax(2)

!vrais_cur
dum_real=(y_cur+0.)/(popsize+0.)
if(y_cur/=0 .and. y_cur/=popsize)  tmp_logcur=log_binomial_coef(popsize,y_cur) + &
                                              y_lect*log(dum_real) + (n_lect-y_lect)*log(1.-dum_real)
!autrement c'est egal a zero NB:vu l'update et les cdt initiales on ne peut pas avoir de cas pathos (ex y_cur=0 et y_lect>0)
!proposal
! y_out=y_cur
! b_inf=max(ymin,y_cur-delta_y) ; b_sup=min(ymax,y_cur+delta_y)
! do while(y_out==y_cur)
!  y_out=min(b_sup,b_inf+floor(grnd()*((b_sup - b_inf) + 1.))) !  y_out=b_inf+nint(rnd*(b_sup - b_inf))
!  !dans quelques tres rare cas le cast de fllor se fait mal: par exemple floor(2.99999999) se retrouve à 3 au lieu de 2 (pareil pour int)
!  !ce pb n'est pas resolu par du passage a real*8 d'ou la securite de min
! end do
!!1/p(q,q') et 1/p(q',q): forward et backward prop probability
! inv_fwd=b_sup - b_inf  !+ 1.: on ne met pas +1 car la valeur courante est exclue dans la proposal on a donc b_sup-b_inf+1 - 1 valeurs
! inv_bwd=min(ymax,y_out+delta_y)-max(ymin,y_out-delta_y)  ! + 1.:on ne met pas +1 car la valeur courante est exclue dans la proposal on a donc b_sup-b_inf+1 - 1 valeurs
! diff_log=log(inv_fwd) - log(inv_bwd) !prior uniforme

!vrais_prop
!dum_real=(y_out+0.)/(popsize+0.)
!if(y_out/=0 .and. y_out/=popsize)  tmp_lognew=log_binomial_coef(popsize,y_out) + &
!                                              y_lect*log(dum_real) + (n_lect-y_lect)*log(1.-dum_real)
! diff_log=diff_log+(y_out-y_cur)*(log(alpha) - log(1.-alpha)) + tmp_lognew -tmp_logcur

!!!!!!!!!!!!!!!!!

!proposal miroir
 y_out=y_cur
 b_inf=max(ymin,y_cur-delta_y) ; b_sup=min(ymax,y_cur+delta_y)
 rnd=y_cur+ (grnd(this)*(2.*delta_y + 1.))  - delta_y
  if(rnd<b_sup+1. .and. rnd>=b_inf) then 
    y_out=int(rnd)
  else
   if(rnd>b_sup) y_out=ceiling(2*b_sup - rnd + 1.)
   if(rnd<b_inf) y_out=int(2*b_inf - rnd)
  end if

if(y_out/=y_cur) then
!vrais_prop
 dum_real=(y_out+0.)/(popsize+0.)
 if(y_out/=0 .and. y_out/=popsize)  tmp_lognew=log_binomial_coef(popsize,y_out) + &
                                               y_lect*log(dum_real) + (n_lect-y_lect)*log(1.-dum_real)
  diff_log=(y_out-y_cur)*(log(alpha) - log(1.-alpha)) + tmp_lognew -tmp_logcur

  rnd=grnd(this)
  if(log(rnd)>diff_log) then
    accept=0 ; y_out=y_cur
  else
    accept=1
  end if
else
  accept=0 
end if
end subroutine update_y

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!update des alpha
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine update_alpha (this,Y_val,N_val,a_cur,a_std_cur,pi_val,cov_val,mu_up,var_up,delta_a,accept,a_out,a_std_out) !delta_a=taille de l'intervalle autour de la proposal de a
implicit none

!on fait une marche aleatoire pure
type(mt_state), intent(inout) :: this
integer, intent(in) :: Y_val, N_val 
real (kind=8), intent(in) :: a_cur,a_std_cur,pi_val,cov_val,mu_up,var_up, delta_a

integer, intent(out) :: accept 
real (kind=8), intent(out) :: a_out,a_std_out

real (kind=8) :: eps=1e-8, b_inf, b_sup,diff_log

accept=1
!on echantillonne la proposal
b_inf=a_std_cur - delta_a/2 ; b_sup=a_std_cur + delta_a/2
a_std_out=b_inf+grnd(this)*(b_sup - b_inf)
a_out=pi_val + cov_val + a_std_out*sqrt(pi_val*(1.0-pi_val)) 
!!!!!!
!calcul de log_new et log_cur
!!!!!!
diff_log= ((a_std_cur-mu_up)**2 - (a_std_out-mu_up)**2) / (2.0*var_up) 

!!!!CAS 1
 if((Y_val .ne. N_val) .and. Y_val>0) then
  if(a_out<0. .or. a_out>1.0) then !on rejette car pas possible
     accept=0 ; a_out=a_cur ; a_std_out=a_std_cur
  else 
   diff_log=diff_log + Y_val*(log(a_out)-log(a_cur)) +  (N_val-Y_val)*(log(1.0-a_out)-log(1.0-a_cur))
   if(log(grnd(this))>diff_log) then
      accept=0 ; a_out=a_cur ; a_std_out=a_std_cur
   end if
  end if
 end if

!!!!CAS 2
 if(Y_val .eq. N_val) then
  if(a_out<eps) then !on rejette car pas possible
     accept=0 ; a_out=a_cur ; a_std_out=a_std_cur
  else !normalement a_cur ne peut jamais etre egal à 0 dans cette configuration (cf au dessus)
    diff_log=diff_log + Y_val*(log(min(a_out,1.0))-log(min(a_cur,1.0)))
      if(log(grnd(this))>diff_log) then
          accept=0 ; a_out=a_cur ; a_std_out=a_std_cur
      end if
   end if
 end if

!!!!CAS 3
 if(Y_val==0) then
  if(a_out>(1.0-eps)) then !on rejette car pas possible
     accept=0 ; a_out=a_cur ; a_std_out=a_std_cur
  else !normalement a_cur ne peut jamais etre egal à 1 dans cette configuration (cf au dessus)
    diff_log=diff_log + N_val*(log(1.0-max(a_out,0.0))-log(1.0-max(a_cur,0.0)))
    if(log(grnd(this))>diff_log) then
      accept=0 ; a_out=a_cur ; a_std_out=a_std_cur
    end if
   end if
  end if
  
end subroutine update_alpha

!!!!

subroutine update_alpha_vect (this,Y_val,N_val,a_cur,a_std_cur,pi_val,cov_vals,delta_a,lambda_mat,c_matrix,accept,a_out,a_std_out) !delta_a=taille de l'intervalle autour de la proposal de a
implicit none

!on fait une marche aleatoire pure
type(mt_state), intent(inout) :: this
integer, intent(in),dimension(:) :: Y_val, N_val 
real (kind=8), intent(in),dimension(:,:) :: lambda_mat,c_matrix
real (kind=8), intent(in),dimension(:) :: a_cur,a_std_cur,cov_vals
real (kind=8), intent(in) :: pi_val,delta_a

integer, intent(out) :: accept 
real (kind=8), intent(out),dimension(:) :: a_out,a_std_out

real (kind=8) :: eps=1e-8, diff_log
integer :: tmp_npop,tmp_i

tmp_npop=size(a_cur)
accept=1

!on echantillonne la proposal
do tmp_i=1,tmp_npop
 a_out(tmp_i)=random_normal(this)*delta_a
end do
a_out=a_cur + matmul(c_matrix,a_out)
a_std_out=(a_out-pi_val-cov_vals)/sqrt(pi_val*(1.-pi_val))

!!!!!!
!calcul de log_new et log_cur
!!!!!!
diff_log= -0.5*(dot_product(a_std_out,matmul(lambda_mat,a_std_out)) - dot_product(a_std_cur,matmul(lambda_mat,a_std_cur)) )

tmp_i=0
do while (accept==1 .and. tmp_i<tmp_npop)
 tmp_i=tmp_i+1
!!!!CAS 1
 if((Y_val(tmp_i) .ne. N_val(tmp_i)) .and. Y_val(tmp_i)>0) then
  if(a_out(tmp_i)<0. .or. a_out(tmp_i)>1.0) then !on rejette car pas possible
     accept=0 
  else 
   diff_log=diff_log + Y_val(tmp_i)*(log(a_out(tmp_i))-log(a_cur(tmp_i))) +  (N_val(tmp_i)-Y_val(tmp_i))*(log(1.0-a_out(tmp_i))-log(1.0-a_cur(tmp_i)))
  end if
 end if
!!!!CAS 2
 if(Y_val(tmp_i) .eq. N_val(tmp_i)) then
  if(a_out(tmp_i)<eps) then !on rejette car pas possible
     accept=0 
  else !normalement a_cur ne peut jamais etre egal à 0 dans cette configuration (cf au dessus)
    diff_log=diff_log + Y_val(tmp_i)*(log(min(a_out(tmp_i),1.0))-log(min(a_cur(tmp_i),1.0)))
  end if
 end if

!!!!CAS 3
 if(Y_val(tmp_i)==0) then
  if(a_out(tmp_i)>(1.0-eps)) then !on rejette car pas possible
     accept=0 
  else !normalement a_cur ne peut jamais etre egal à 1 dans cette configuration (cf au dessus)
    diff_log=diff_log + N_val(tmp_i)*(log(1.0-max(a_out(tmp_i),0.0))-log(1.0-max(a_cur(tmp_i),0.0)))
  end if
 end if
end do

if(accept==0) then
  a_out=a_cur ; a_std_out=a_std_cur
else
  if(log(grnd(this))>diff_log) then
    accept=0 ; a_out=a_cur ; a_std_out=a_std_cur
  end if
end if

  
end subroutine update_alpha_vect 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!update des P_i
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine update_pi (this,p_cur,vect_a_cur,cov_vals,lambda_mat,beta_params,delta_p,accept,p_out,dum_real_pop) 
implicit none

!integer, intent(in) :: i_cur !locus considere
type(mt_state), intent(inout) :: this
real (kind=8), intent(in), dimension(:,:) :: lambda_mat,vect_a_cur
real (kind=8), intent(in), dimension(:) :: cov_vals
real (kind=8), intent(in) :: delta_p , p_cur, beta_params(2)

integer, intent(out) :: accept !0 si accepte, 1 sinon
real (kind=8), intent(out) :: p_out
real (kind=8), intent(out), dimension(:) :: dum_real_pop !on stocke vect_a_std_out

real (kind=8) :: eps=1e-4,b_inf, b_sup, inv_fwd, inv_bwd, diff_log,coef

 accept=1
 b_inf=max(eps,p_cur-delta_p/2) ; b_sup=min(1.-eps,p_cur+delta_p/2)
 p_out=b_inf+grnd(this)*(b_sup - b_inf)

!1/p(q,q') et 1/p(q',q): forward et backward prop probability
 inv_fwd=b_sup - b_inf
 inv_bwd=min(1.-eps,p_out+delta_p/2) - max(eps,p_out-delta_p/2)

coef=size(dum_real_pop)/2. + 1.
diff_log=(beta_params(1) - coef)*(log(p_out) - log(p_cur)) + &
         (beta_params(2) - coef)*(log(1.0-p_out) - log(1.0-p_cur)) + &
          log(inv_fwd) - log(inv_bwd)


 dum_real_pop=(vect_a_cur(:,1)-p_out-cov_vals)/sqrt(p_out*(1.-p_out))
 diff_log=diff_log - 0.5*(dot_product(dum_real_pop,matmul(lambda_mat,dum_real_pop)) - &
                     dot_product(vect_a_cur(:,2),matmul(lambda_mat,vect_a_cur(:,2))))

 if(log(grnd(this))>diff_log) then
   accept=0 ;  p_out=p_cur ; dum_real_pop=vect_a_cur(:,2)
 end if   

 end subroutine update_pi
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!update de pi en slice
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function slc_pi(this,p_cur,vect_a_cur,lambda_mat,beta_params) 
implicit none

!integer, intent(in) :: i_cur !locus considere
type(mt_state), intent(inout) :: this
real (kind=8), intent(in), dimension(:,:) :: lambda_mat,vect_a_cur
real (kind=8), intent(in) :: p_cur, beta_params(2)

real (kind=8) :: slc_pi
real (kind=8), dimension(:) :: dum_real_pop(size(vect_a_cur,1)) !on stocke vect_a_std_out

real (kind=8) :: eps=1e-4,b_inf, b_sup, F_cur,F_out,F_rnd,p_out,coef

  coef=size(dum_real_pop)/2. + 1.
  F_cur=(beta_params(1) - coef)*log(p_cur) + (beta_params(2) - coef)*log(1.0-p_cur) &
         - 0.5*dot_product(vect_a_cur(:,2),matmul(lambda_mat,vect_a_cur(:,2)))
  F_rnd=log(grnd(this)) + F_cur
  b_inf=eps ; b_sup=1.-eps

  p_out=b_inf+grnd(this)*(b_sup - b_inf)
  dum_real_pop=(vect_a_cur(:,1)-p_out)/sqrt(p_out*(1.-p_out))
  F_out=(beta_params(1) - coef)*log(p_out) + (beta_params(2) - coef)*log(1.0-p_out) &
         - 0.5*dot_product(dum_real_pop,matmul(lambda_mat,dum_real_pop))
  do while (F_out<F_rnd)
   if (F_out<F_cur) then
    b_inf=p_out
   else
    b_sup=p_out
   end if
   p_out=b_inf+grnd(this)*(b_sup - b_inf)
   dum_real_pop=(vect_a_cur(:,1)-p_out)/sqrt(p_out*(1.-p_out))
   F_out=(beta_params(1) - 1.5)*log(p_out) + (beta_params(2) - 1.5)*log(1.0-p_out) &
         - 0.5*dot_product(dum_real_pop,matmul(lambda_mat,dum_real_pop))
  end do
  
  slc_pi=p_out
end function slc_pi

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!     UPDATE DE MU BATA_PI  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine update_beta_pi_params(this,beta_params,tmp_nmrk,params_save,delta_mu,delta_phi,accept_mu,accept_phi,val_out)
implicit none
!val_out(5): a_out,b_out,gamma_log(phi_out),gamma_log(phi_out*mu_out),gamma_log(phi_out*(1-mu_out))
type(mt_state), intent(inout) :: this
real  (kind=8), intent(in) :: delta_mu,delta_phi, beta_params(2),params_save(5)
integer, intent(in) :: tmp_nmrk 

integer, intent(out) :: accept_mu,accept_phi !0 si accepte, 1 sinon
real  (kind=8), intent(out) :: val_out(5)

real  (kind=8) :: b_inf, b_sup, mu_out,inv_fwd, inv_bwd,mu_cur,phi_cur,phi_out,diff_log,tmp_val_out(5,2)

 phi_cur=sum(beta_params) ; mu_cur=beta_params(1)/phi_cur

!!!!!!!!!!!!!
!update de mu
!!!!!!!!!!!!!

 b_inf=max(1e-3,mu_cur-delta_mu/2) ; b_sup=min(1.-1e-3,mu_cur+delta_mu/2)
 mu_out=b_inf+grnd(this)*(b_sup - b_inf)
!1/p(q,q') et 1/p(q',q): forward et backward prop probability
 inv_fwd=b_sup - b_inf
 inv_bwd=min(1.-1e-3,mu_out+delta_mu/2)-max(1e-3,mu_out-delta_mu/2)

 tmp_val_out(3,1)=params_save(3)
 tmp_val_out(4,1)=gamma_log(mu_out*phi_cur) ; tmp_val_out(5,1)=gamma_log((1.-mu_out)*phi_cur)

!prior uniforme sur les pi
 diff_log=log(inv_fwd) - log(inv_bwd) + tmp_nmrk*(sum(params_save(4:5)-tmp_val_out(4:5,1))) + &
          phi_cur*(mu_out-mu_cur)*(params_save(1)-params_save(2))

 if(log(grnd(this))>diff_log) then
   accept_mu=0 ; mu_out=mu_cur ; tmp_val_out(4:5,1)=params_save(4:5)
 else
   accept_mu=1 
 end if

!!!!!!!!!!!!!
!update de phi
!!!!!!!!!!!!!

  phi_out=exp(random_normal(this)*delta_phi + log(phi_cur))
  tmp_val_out(4,2)=gamma_log(mu_out*phi_out) ; tmp_val_out(5,2)=gamma_log((1.-mu_out)*phi_out) ; tmp_val_out(3,2)=gamma_log(phi_out)
  
  diff_log=log(phi_out) - log(phi_cur) + tmp_nmrk*( sum(tmp_val_out(4:5,1) - tmp_val_out(4:5,2)) + &
                  tmp_val_out(3,2) - tmp_val_out(3,1)) + &
                  (phi_out-phi_cur)*(mu_out*params_save(1) + (1.-mu_out)*params_save(2)-1.)

 if(log(grnd(this))>diff_log) then
   accept_phi=0 ; phi_out=phi_cur ; val_out(3:5)=tmp_val_out(3:5,1)
 else
   accept_phi=1 ; val_out(3:5)=tmp_val_out(3:5,2)
 end if

 val_out(1)=mu_out*phi_out ; val_out(2)=(1.-mu_out)*phi_out

end subroutine update_beta_pi_params


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!UPDATE DES DELTA_i!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine up_delta(this,delta_cur,beta_cur,tmp_npop,P,alpha_vect, chol_lambda_in,phi_j,cov_val_in,pi_val,delta_out,cov_val_out)
implicit none

type(mt_state), intent(inout) :: this
real (kind=8), intent(in), dimension(:,:) :: chol_lambda_in
real (kind=8), intent(in), dimension(:)   :: alpha_vect,phi_j,cov_val_in
real (kind=8), intent(in) :: beta_cur,P,pi_val
integer , intent(in) :: tmp_npop,delta_cur

integer ,intent(out) :: delta_out
real (kind=8), intent(out), dimension(:) :: cov_val_out

real (kind=8) :: ddot_alpha(tmp_npop),tilde_phi(tmp_npop),f_0, f_1, prob_ber

 cov_val_out = cov_val_in ! toujours pour delta_0
 if(delta_cur==1) cov_val_out=cov_val_in - beta_cur*phi_j 
 ddot_alpha=(alpha_vect - pi_val - cov_val_out)/sqrt(pi_val*(1.-pi_val))
 ddot_alpha=matmul(chol_lambda_in,ddot_alpha)

 f_0 = log(1.-P) - 0.5*sum(ddot_alpha**2)
 tilde_phi=matmul(chol_lambda_in,phi_j)/sqrt(pi_val*(1.-pi_val))
 f_1 = log(P) - 0.5*sum((ddot_alpha - beta_cur*tilde_phi)**2)
 prob_ber = 1./(1. + exp(f_0 - f_1))
 delta_out = 0
 if(grnd(this) < prob_ber) delta_out = 1
 cov_val_out=cov_val_out + delta_out*beta_cur*phi_j

end subroutine up_delta

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!UPDATE DES DELTA_i (Ising prior)    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine up_delta_ising(this,i_cur,delta_vect,ising_b,beta_cur,tmp_npop,tmp_nmrk,P,alpha_vect, chol_lambda_in,phi_j,cov_val_in,pi_val,delta_out,cov_val_out)
!si ising_b=0. up_delta_ising <=> up_delta
implicit none

type(mt_state), intent(inout) :: this
real (kind=8), intent(in), dimension(:,:) :: chol_lambda_in
real (kind=8), intent(in), dimension(:)   :: alpha_vect,phi_j,cov_val_in
real (kind=8), intent(in) :: beta_cur,P,pi_val,ising_b
integer , intent(in), dimension(:) :: delta_vect
integer , intent(in)               :: tmp_npop,i_cur,tmp_nmrk

integer ,intent(out) :: delta_out
real (kind=8), intent(out), dimension(:) :: cov_val_out

real (kind=8) :: ddot_alpha(tmp_npop),tilde_phi(tmp_npop),f_0, f_1, prob_ber,delta_cur
integer, dimension(2):: cpt_ising=(/-1,1/),pen_ising=(/-1,1/),val_cpt_ising=(/0,0/)

 cov_val_out = cov_val_in ! toujours pour delta_0
 delta_cur = delta_vect(i_cur)
 if(delta_cur==1) cov_val_out=cov_val_in - beta_cur*phi_j 
 ddot_alpha=(alpha_vect - pi_val - cov_val_out)/sqrt(pi_val*(1.-pi_val))
 ddot_alpha=matmul(chol_lambda_in,ddot_alpha)
!ising stuff
 if(i_cur<tmp_nmrk) val_cpt_ising=cpt_ising*pen_ising(delta_vect(i_cur+1)+1)
 if(i_cur>1) val_cpt_ising=val_cpt_ising + cpt_ising*pen_ising(delta_vect(i_cur-1)+1)

 f_0 = log(1.-P) - 0.5*sum(ddot_alpha**2) + ising_b*val_cpt_ising(1)
 tilde_phi=matmul(chol_lambda_in,phi_j)/sqrt(pi_val*(1.-pi_val))
 f_1 = log(P) - 0.5*sum((ddot_alpha - beta_cur*tilde_phi)**2) + ising_b*val_cpt_ising(2)
 prob_ber = 1./(1. + exp(f_0 - f_1))
 delta_out = 0
 if(grnd(this) < prob_ber) delta_out = 1
 cov_val_out=cov_val_out + delta_out*beta_cur*phi_j

end subroutine up_delta_ising


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!UPDATE DES Beta_i!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine up_beta_i(this,beta_cur,tmp_npop,alpha_vect, chol_lambda_in,delta_i,phi_j,cov_val_in,pi_val,tau_beta,beta_out,cov_val_out)
implicit none

type(mt_state), intent(inout) :: this
real (kind=8), intent(in), dimension(:,:) :: chol_lambda_in !en fait inv_c_mat
real (kind=8), intent(in), dimension(:)   :: alpha_vect,phi_j,cov_val_in
real (kind=8), intent(in) :: tau_beta,beta_cur,pi_val
integer , intent(in) :: delta_i,tmp_npop

real (kind=8),intent(out) :: beta_out
real (kind=8), intent(out), dimension(:) :: cov_val_out

real (kind=8) :: ddot_alpha(tmp_npop),tilde_phi(tmp_npop)
real (kind=8) :: tmp_mu, tmp_tau

 cov_val_out=cov_val_in
 if(delta_i==0) then
  beta_out=random_normal(this)/sqrt(tau_beta)
 else
  tilde_phi=matmul(chol_lambda_in,phi_j)/sqrt(pi_val*(1.-pi_val))
  ddot_alpha=(alpha_vect - pi_val - cov_val_in + beta_cur*phi_j)/sqrt(pi_val*(1.-pi_val))
  ddot_alpha=matmul(chol_lambda_in,ddot_alpha)
  tmp_tau = tau_beta + sum(tilde_phi**2)
  tmp_mu  = sum(phi_j*ddot_alpha)/tmp_tau
  beta_out=(random_normal(this)/sqrt(tmp_tau)) + tmp_mu
  cov_val_out=cov_val_in + (beta_out-beta_cur)*phi_j
 end if

end subroutine up_beta_i

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!UPDATE DES Beta_i!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine up_beta_i_coop(this,beta_cur,tmp_npop,alpha_vect, chol_lambda_in,delta_i,phi_j,cov_val_in,pi_val,min_beta,max_beta,delta_b,beta_out,cov_val_out,accept)
implicit none

type(mt_state), intent(inout) :: this
real (kind=8), intent(in), dimension(:,:) :: chol_lambda_in
real (kind=8), intent(in), dimension(:)   :: alpha_vect,phi_j,cov_val_in
real (kind=8), intent(in) :: beta_cur,pi_val,min_beta, max_beta, delta_b
integer , intent(in) :: delta_i,tmp_npop

real (kind=8),intent(out) :: beta_out
real (kind=8), intent(out), dimension(:) :: cov_val_out
integer , intent(out) :: accept

real (kind=8) :: ddot_alpha(tmp_npop),tilde_phi(tmp_npop)
real (kind=8) :: diff_log,b_sup,b_inf,inv_fwd,inv_bwd 

 cov_val_out=cov_val_in
 accept=1
 if(delta_i==1) then
  b_inf=max(min_beta,beta_cur-delta_b) ; b_sup=min(max_beta,beta_cur+delta_b)
  beta_out=b_inf+grnd(this)*(b_sup - b_inf)
  inv_fwd=b_sup - b_inf
  inv_bwd=min(max_beta,beta_out+delta_b)-max(min_beta,beta_out-delta_b)
  tilde_phi=matmul(chol_lambda_in,phi_j)/sqrt(pi_val*(1.-pi_val))
  ddot_alpha=(alpha_vect - pi_val - cov_val_in + beta_cur*phi_j)/sqrt(pi_val*(1.-pi_val))
  ddot_alpha=matmul(chol_lambda_in,ddot_alpha)
  diff_log=log(inv_fwd) - log(inv_bwd) - 0.5*( &
            sum((ddot_alpha - beta_out*tilde_phi)**2) - sum((ddot_alpha - beta_cur*tilde_phi)**2))
  if(log(grnd(this))>diff_log) then 
   beta_out=beta_cur ; accept=0
  end if
  cov_val_out=cov_val_in + (beta_out-beta_cur)*phi_j 
 else 
  beta_out=min_beta+grnd(this)*(max_beta - min_beta)
 end if

end subroutine up_beta_i_coop

!~ subroutine up_beta_i_coop(beta_cur,aux_model,tmp_npop,alpha_vect, chol_lambda_in,delta_i,phi_j,cov_val_in,pi_val,min_beta,max_beta,delta_b,beta_out,cov_val_out,accept)
!~ implicit none
!~ 
!~ real (kind=8), intent(in), dimension(:,:) :: chol_lambda_in
!~ real (kind=8), intent(in), dimension(:)   :: alpha_vect,phi_j,cov_val_in
!~ real (kind=8), intent(in) :: beta_cur,pi_val,min_beta, max_beta, delta_b
!~ integer , intent(in) :: delta_i,tmp_npop
!~ logical , intent(in) :: aux_model
!~ 
!~ real (kind=8),intent(out) :: beta_out
!~ real (kind=8), intent(out), dimension(:) :: cov_val_out
!~ integer , intent(out) :: accept
!~ 
!~ real (kind=8) :: ddot_alpha(tmp_npop),tilde_phi(tmp_npop)
!~ real (kind=8) :: diff_log,b_sup,b_inf,inv_fwd,inv_bwd
!~ integer :: tmp_i
!~ 
!~ cov_val_out=cov_val_in ; accept=1
!~ ! beta_out=min(min_beta,max(min_beta,beta_cur + random_normal()*sd_prop))
!~ if(delta_i==1) then
!~   b_inf=max(min_beta,beta_cur-delta_b) ; b_sup=min(max_beta,beta_cur+delta_b)
!~   beta_out=b_inf+grnd()*(b_sup - b_inf)
!~   inv_fwd=b_sup - b_inf
!~   inv_bwd=min(max_beta,beta_out+delta_b)-max(min_beta,beta_out-delta_b)
!~   tilde_phi=matmul(chol_lambda_in,phi_j)/sqrt(pi_val*(1.-pi_val))
!~   ddot_alpha=(alpha_vect - pi_val - cov_val_in + beta_cur*phi_j)/sqrt(pi_val*(1.-pi_val))
!~   ddot_alpha=matmul(chol_lambda_in,ddot_alpha)
!~   diff_log=log(inv_fwd) - log(inv_bwd) - 0.5*( &
!~             sum((ddot_alpha - beta_out*tilde_phi)**2) - sum((ddot_alpha - beta_cur*tilde_phi)**2))
!~   if(aux_model) then !prior triangulaire
!~    inv_fwd=beta_out/beta_cur ! inv_fwd recycle
!~    inv_bwd=max_beta/min_beta ! inv_bwd recycle
!~    if(inv_fwd>0.) then
!~     diff_log=diff_log + log(inv_fwd)
!~    else
!~     if(beta_out<0.) then
!~      diff_log=diff_log + log(inv_fwd*inv_bwd)
!~     else
!~      diff_log=diff_log + log(inv_fwd/inv_bwd)
!~     end if
!~    end if
!~   end if
!~   if(log(grnd())>diff_log) then
!~    beta_out=beta_cur ; accept=0
!~   end if
!~   cov_val_out=cov_val_in + (beta_out-beta_cur)*phi_j
!~ else ! forcement aux_model=>echantillonne dans prior triangulaire via un slice initilaise sur valeur absolue de beta_cur
!~   beta_out=abs(beta_cur)
!~   if(grnd()<max_beta/(max_beta-min_beta)) then !triangle positif
!~     b_sup=max_beta
!~     do tmp_i=1,100
!~      inv_fwd=grnd()*beta_out/(2.*b_sup)
!~      b_inf=2*b_sup*inv_fwd
!~      beta_out=b_inf + grnd()*(b_sup - b_inf)
!~     end do
!~   else !triangle negatif
!~     b_sup=abs(min_beta)
!~     do tmp_i=1,100
!~      inv_fwd=grnd()*beta_out/(2.*b_sup)
!~      b_inf=2*b_sup*inv_fwd
!~      beta_out=b_inf + grnd()*(b_sup - b_inf)
!~     end do
!~     beta_out=-1.*beta_out
!~   end if
!~ end if
!~ 
!~ end subroutine up_beta_i_coop

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!     UPDATE DE P  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

function up_P (this,Delta_i,beta_aP,beta_bP)
implicit none

type(mt_state), intent(inout) :: this
integer, intent(in), dimension(:) :: Delta_i
real (kind=8), intent(in) :: beta_aP,beta_bP
real (kind=8) :: up_P,sum_delta, tmp_I

 tmp_I=size(Delta_i) + 0. ; sum_delta=sum(Delta_i) + 0.
 up_P=random_beta(beta_aP + sum_delta,beta_bP + tmp_I - sum_delta,this)

 end function up_P

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!UPDATE DES tau_beta!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

function up_tau_beta(this,vect_beta_i,delta_i,k_gamma,l_gamma)
implicit none

type(mt_state), intent(inout) :: this
real (kind=8), intent(in), dimension(:) :: vect_beta_i
real (kind=8), intent(in) :: k_gamma,l_gamma
integer , intent(in), dimension(:) :: delta_i

real (kind=8) :: up_tau_beta

real (kind=8) :: k_gam_full, l_gam_full

  k_gam_full  = k_gamma + sum(delta_i)/2.
  l_gam_full  = l_gamma + sum( delta_i*vect_beta_i*vect_beta_i)/2.
  up_tau_beta = random_gamma1(k_gam_full,this) / l_gam_full
  if(isnan(up_tau_beta)) up_tau_beta=100.
!  up_tau_beta  = max(0.0001,min(10000.,random_gamma1(k_gam_full) / l_gam_full))

end function up_tau_beta

end module updates_baypass
