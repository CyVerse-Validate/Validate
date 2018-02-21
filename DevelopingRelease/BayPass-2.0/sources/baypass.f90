!17/01/2015: passage à la command line
!02/02/2015: introduction de la possibilité d'avoir des missing data dans le cas de données Poolseq seulement
!24/07/2015: passage à bypass avec modifs marginales (i) BF en dB, ii) resolution des Inf dans BPis, iii) retrait impression summary_omega si -omegafile, iv) ajout delta0yij option, v) si aux on active automatique covmcmc,)
!8/12/2015: corrections fuites memoires (cur_stream=0 qd pas dans boucle parallele, initialisation de dum_char_array la premiere fois), correction lecture opt_reg, opt_pheno et opt_aux avec initialisation propre, modification affichage ecran (anglais, matrice plus affichée, initialisation stream plus affichee, heure..) 

include 'mcmc_mat_utils.f90'
include 'updates_baypass.f90'
include 'M_kracken.f90' ! ligne de commande

program baypass
use mt_stream
use omp_lib
use mcmc_utils !fonctions diverses pour les summary statistics
use updates_baypass !update
use M_kracken

implicit none

integer:: err, i_thin,npop, nmrk, pop, mrk, tst, tst_pi,tst_pij,tst_b,iter, pilot,&
          val_systime(2),tmp_jhms(4),seed, nvaleurs, thin,burn_in,npilot,pilot_length,y_out,tst_y,&
          dum_int,dum_int2,ddl_wish,rho,tmp_y_rep,npheno,pheno,ninter_beta,delta0_yij,&
          iflen,igot,ibegin(2),iterm(2),ilen,ier,nstream,cur_stream
integer, allocatable :: Y_OBS(:,:), N_OBS(:,:) ,& !observations
                        PoolSize(:),Y_READS(:,:),N_READS(:,:),delta_y(:,:),Yminmax(:,:,:),&
                        INITS_Delta(:,:),&
                        omega_index_subpop(:,:),dum_int_pop(:),& !tab de dim npop,npop-1 dans lesquels sont stockés les indices relevant pour les extractions de chque pop (cf alpha_up_muvar)
                        dum_read_line(:)
real (kind=8), allocatable :: INITS_Pij(:,:,:),mean_pij(:,:,:),delta_pij(:,:),acc_pij(:,:),& !frequences alleliques
                              INITS_Pi(:),delta_pi(:),acc_pi(:),mean_pi(:,:),& !frequences ancestrales
                              PHENO_VAL(:,:),INITS_Beta_i(:,:),INITS_tau_beta(:),INIT_Pdelta(:),&
                              mean_beta_i(:,:,:),mean_tau_beta(:,:),mean_delta_i(:,:),mean_p_delta(:,:),sum_covariates(:,:),&
                              delta_betai(:,:),acc_betai(:,:),&!au cas ou up_betai_coop
                              INITS_LDA(:,:),mean_lda(:,:,:),omega_mat(:,:),mean_omega(:,:,:),mean_xtx(:,:),& 
                              acc_y(:,:),CPO(:,:),log_cn(:,:),mean_yij(:,:,:),muvar_cur(:,:),&
                              dum_real_pop(:),dum_real_pop2(:),dum_real_pheno(:),c_mat(:,:),dum_vect_grid(:,:),&
                              reg_mat(:,:),inv_c_mat(:,:),mean_rho_coef(:,:,:),mean_isbf(:,:,:),grid_is(:,:,:),dum_pheno_std(:,:),&
                              array_real_pop(:,:,:)
real (kind=8):: delta0_pij,delta0_pi,acc_inf=0.25,acc_sup=0.4,rate_adjust=1.25,& !pour etudes pilotes
             !  beta_pi=0.7,lda_a_prior=3.,lda_b_prior=4.,& !parmatres hyperpriors (a rentrer)
                pi_beta_params(2),acc_beta_params(2),mean_beta_params(2,2),pi_beta_update_save(5),& !Sum(log(Pi)) ; Sum(log(1-pi)) ; gamma_log(Phi) ; gamma_log(Phi*mu) ; gamma_log(Phi*(1-mu)) avec Phi=(a_pi+b_pi) et mu=a_pi/Phi
                P_beta_params(2),&
                delta_pi_beta_mu=0.05,delta_pi_beta_phi=1,&
                bpval,kld,bf,tmp_mean,tmp_min,tmp_max,dum_real,dum_real2,dum_real5(5),tmp_a,tmp_b,&
                alpha_up_muvar_param(2),tmp_t_rep,tmp_t,tmp_t_obs,tmp_vij,k_tau_prior=1.,l_tau_prior=10.,tau_beta0=20.,&
                min_beta=-0.3 , max_beta=0.3 ,ising_b
real (kind=8) :: dum_real16,dum_real16_2, mean_deviance,lpml !passage a kind=8 car pb avc certains compulateur
!attention tres important d'initialiser ici les opt et dum_char_array
logical  :: up_params_beta, poolseq , up_alpha_coop, opt_pheno=.false. ,& !up_slc_pi=.false.
            opt_aux=.false. , estim_tau_beta=.false. , out_pilot , estim_omega,opt_reg=.false.,up_betai_coop,opt_ising,scale_cov
logical, allocatable :: missing_data(:,:)

character(len=255) :: Geno_file,Pheno_file,omega_file,poolsize_file,dum_char='',dum_char_array(2)='',&
                      sum_yij_pij_file,sum_pij_file,sum_pi_file,sum_omega_file,sum_mat_omega_file,&
                      sum_betaparams_file,sum_betai_file,sum_Pdelta_file,sum_taufile,&
                      sum_beta_i_reg_file,sum_DIC_file,sum_phenostd_file,out_prefix

type(mt_state),allocatable :: mts(:) !mts(0:NSTREAM-1)

call system_clock(val_systime(1),tst, tst_pi)
call void_read_input()
!~ call sgrnd(seed) !seed du mersenne twister

print *,''
print *,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
print *,'        PILOT RUNS START'
print *,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
print *,''
!~ call timestamp()

pilot=1 
 call void_initialisation()
 print *,'' ; print *,'###########INITS VALUES#############' ; print *,''
 call void_printval()

do while (pilot<=npilot)
 acc_pij=0. ; acc_pi=0.  ; acc_beta_params=0.
 if(opt_pheno .and. up_betai_coop) acc_betai=0.
 if(poolseq) acc_y=0.
 print *,'PILOT RUN: ',pilot
 do iter=1,pilot_length
  if(mod(iter,pilot_length/10)==0) print *,'  iteration=',iter
  call void_mcmc_iter()
 end do
 call void_adjust_proposals()
 print *,'' ; print *,'########### NODE VALUES #############' ; print *,''
 call void_printval()
 pilot=pilot+1
 call void_initialisation()
 print *,''
end do

print *,''
print *,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
print *,'        BURN-IN: START'
print *,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
print *,''

if(opt_aux) then
 if(estim_tau_beta) then
  allocate(dum_real_pheno(npheno)) ;  dum_real_pheno=0.
 else
  out_pilot=.true.
 end if
end if

do iter=1,burn_in
 if(mod(iter,burn_in/10)==0) print *,'  iteration=',iter
 call void_mcmc_iter()
 if(opt_aux .and. estim_tau_beta) dum_real_pheno=dum_real_pheno + INITS_tau_beta
end do

if(opt_aux .and. estim_tau_beta) then
 INITS_tau_beta=(dum_real_pheno/burn_in)*(P_beta_params(1)/sum(P_beta_params))
 out_pilot=.true.
 print *,'' ; print *,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!' ; print *,'     2nd BURN-IN: START'
 print *,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!' ; print *,''
 write(*,'(A,1x,100(f12.6,1x))') 'Estimated Precision(s): ',INITS_tau_beta
 print *,''
 do iter=1,burn_in
  if(mod(iter,burn_in/10)==0) print *,'  iteration=',iter
  call void_mcmc_iter()
 end do
end if

call void_printval()

print *,''
print *,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
print *,'        MCMC CHAIN STARTED'
print *,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
print *,''
print *,''

 allocate(CPO(nmrk,npop)) ; CPO=0. ; mean_deviance=0.
 acc_pij=0. ; acc_pi=0. ; acc_beta_params=0.
 allocate(mean_pij(nmrk,npop,4),mean_pi(nmrk,2),mean_xtx(nmrk,2),mean_lda(npop,npop,2),mean_omega(npop,npop,2))
 mean_lda=0. ; mean_omega=0. ; mean_pij=0.; mean_pi=0.; mean_beta_params=0. ; mean_xtx=0.
 if(poolseq) then
  allocate(mean_yij(nmrk,npop,2))
  mean_yij=0. ; acc_y=0.
 end if

 if(opt_pheno) then
  allocate(mean_beta_i(nmrk,npheno,2),mean_tau_beta(npheno,2),mean_delta_i(nmrk,npheno),mean_p_delta(npheno,2))
  mean_beta_i=0. ; mean_tau_beta=0. ; mean_delta_i=0. ; mean_p_delta=0.
 end if

 do iter=1,nvaleurs
  if(mod(iter,nvaleurs/100)==0) print *,'  iteration=',iter
  do i_thin=1,thin
   call void_mcmc_iter()
  end do
!  write(7,'(100(f14.8,1x))') (/(INITS_LDA(dum_int,dum_int),dum_int=1,npop)/)
  call void_update_summary()
 end do

 call void_print_summary()

 do pop=0,nstream-1
    call delete(mts(pop))
 enddo 
call system_clock(val_systime(2),tst, tst_pi)
tmp_jhms=elapse_time(val_systime(1),val_systime(2),tst)
write(*,'(A,4(1x,i2,A))') 'Analysis took:',tmp_jhms(1),' days',tmp_jhms(2), ' h',tmp_jhms(3),' min',tmp_jhms(4),' sec'

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!! calcul des differentes reg_mat (il y en a npop)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine compute_reg_mat() !delta_a=taille de l'intervalle autour de la proposal de a
 do dum_int=1,npop
  reg_mat(dum_int,:)=matmul(omega_mat(dum_int,omega_index_subpop(dum_int,:)),inv(omega_mat(omega_index_subpop(dum_int,:),omega_index_subpop(dum_int,:))))
 end do
end subroutine compute_reg_mat

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   UPDATE DE MATRICE LDA !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

function up_lda_mat()
!NB ici p_wish=J=npop (ca simplifie un peu les calculs mais pas complique a changer au cas ou)
 integer (kind=8) :: tmp_i,tmp_j
 real (kind=8) :: up_lda_mat(npop,npop),mat_wish(npop,npop) 

  mat_wish=0.
  if(opt_pheno) then
   do tmp_i=1,nmrk
    mat_wish=mat_wish + txx(INITS_Pij(tmp_i,:,2) - sum_covariates(tmp_i,:)/(sqrt(INITS_Pi(tmp_i)*(1.-INITS_Pi(tmp_i)))))
   end do
  else
   do tmp_i=1,nmrk
    mat_wish=mat_wish + txx(INITS_Pij(tmp_i,:,2))
   end do
  end if
  do tmp_j=1,npop
    mat_wish(tmp_j,tmp_j)=mat_wish(tmp_j,tmp_j) + rho
  end do
  mat_wish=inv(mat_wish)
  call wishart_sample ( npop, ddl_wish, mat_wish, up_lda_mat, mts(cur_stream) )

 end function up_lda_mat


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!LECTURE DES INPUT/ALLOCATION 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine void_read_input()
!~  call kracken('cmd',&
!~    &'-accinf 0.25 -accsup 0.4 -adjrate 1.25 -auxmodel no -auxPbetaprior 0.02 1.98 -betapiprior 1.0 1.0 &
!~      &-burnin 5000 -covmcmc no -d0pi 0.5 -d0pij 0.05 -efile -estpibetapar no -esttaubeta no -gfile -help no &
!~      &-isingbeta 0.0 -maxbeta 0.3 -minbeta -0.3 -nbetagrid 201 -npilot 20 -npop -1 -nval 1000 &
!~      &-omegafile -outprefix -pilotlength 1000 -poolsizefile -pppval no -rho 1 -scalecov no -seed 5001 -tauB0 -1. -thin 25 -upalphaalt no ') 
 call kracken('cmd','-accinf 0.25 -accsup 0.4 -adjrate 1.25 -auxmodel no -auxPbetaprior 0.02 1.98 -betapiprior 1.0 1.0 -burnin 5000 -covmcmc no -d0pi 0.5 -d0pij 0.05 -d0yij 1 -efile -setpibetapar no -esttaubeta no -gfile -help no -isingbeta 0.0 -maxbeta 0.3 -minbeta -0.3 -nbetagrid 201 -npilot 20 -npop -1 -nthreads 1 -nval 1000 -omegafile -outprefix -pilotlength 1000 -poolsizefile -rho 1 -scalecov no -seed 5001 -tauB0 -1. -thin 25 -upalphaalt no')      

 call retrev("cmd_help",dum_char,iflen,ier)
 if(iflen==0) then !on imprime l'aide
  call void_print_help()
  stop
 end if

 npop = iget("cmd_npop")
 if(npop<0) then
  write(*,'(A)') 'ERROR: Please provide the number of populations'
  stop
 end if

!recuperation de nmrk
 call retrev("cmd_gfile",dum_char,iflen,ier)
 if(iflen==0) then
  write(*,'(A)') 'ERROR: Please provide a genotyping file'
  stop
 else !on recupere le nombre de marker et on verifie que ca s'ouvre bien
  Geno_file=dum_char(:iflen)
  open(1,file=Geno_file,status='old') 
  nmrk=0
  do 
   read(1,*,end=1,iostat=err) 
   nmrk=nmrk+1
  end do
  1 continue
  close(1)
 end if

 seed=iget("cmd_seed") ; ninter_beta=iget("cmd_nbetagrid")
 npilot = iget("cmd_npilot") ; burn_in =  iget("cmd_burnin") ; thin=iget("cmd_thin")
 nvaleurs=iget("cmd_nval") ; pilot_length = iget("cmd_pilotlength")
 delta0_yij=iget("cmd_d0yij") ; delta0_pij=rget("cmd_d0pij") ; delta0_pi=rget("cmd_d0pi") ; rate_adjust=rget("cmd_adjrate")
 acc_inf=rget("cmd_accinf") ; acc_sup=rget("cmd_accsup")
 min_beta=rget("cmd_minbeta") ; max_beta=rget("cmd_maxbeta")

 write(*,'(A,i9)') ' No of Markers                      = ',nmrk ;   write(*,'(A,i9)') ,' No of Populations                  = ',npop
 write(*,'(A,A)') ' Genotype File name                 =    ',trim(Geno_file)
 print *,''
 print *,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
 print *,'        Model Specifications'
 print *,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
 print *,''

!poolseq?
 call retrev("cmd_poolsizefile",dum_char,iflen,ier) 
 if(iflen>0) then
  poolseq=.true. ; poolsize_file=dum_char(:iflen)
  allocate(PoolSize(npop))
  poolsize_file=dum_char(:iflen)
  open(1,file=poolsize_file,status='old')
  read(1,*) PoolSize(1:npop)
  write(*,'(A,1000(i5,1x))') ' -> POOLSEQ MODE with pool sizes: ',PoolSize
  close(1)
  else
   poolseq=.false.
  end if

!beta prior sur les Pi?
 call retrev("cmd_setpibetapar",dum_char,iflen,ier)
 up_params_beta=.true. ; if(iflen==0) up_params_beta=.false.
 print *,''
 call retrev("cmd_betapiprior",dum_char,iflen,ier)
 call delim(dum_char,dum_char_array(1:2),2,igot,ibegin,iterm,ilen,' ,:')
 call string_to_real(dum_char_array(1),pi_beta_params(1),ier)
 call string_to_real(dum_char_array(2),pi_beta_params(2),ier)

!!Type de modele possible: i)   sans covariate (avec ou sans estim_omega)  =>opt_pheno=opt_reg=opt_aux=false et estim_omega=true ou false
!                          ii)  IS covariate (avec ou sans estim_omega)     =>opt_reg=true ; opt_pheno=opt_aux=false   et estim_omega=true ou false
!                          iii) MCMC covariate (forcement sans estim_omega) => STD model (opt_pheno=true et opt_aux=false) OU AUX model (opt_pheno=true et opt_aux=true);  et estim_omega=false,opt_ppp=.false
!                               iii.1) prior uniforme sur beta avec ou sans aux
!                               iii.2) prior gaussienne avec ou sans aux

!on recupere l'ensemble des infos et on gere les incompatibilites'

 call retrev("cmd_omegafile",dum_char,iflen,ier)  
 if(iflen>0) then
  omega_file=dum_char(:iflen)
  estim_omega=.false. 
 else 
  estim_omega=.true.
  rho=iget("cmd_rho") ; ddl_wish=rho + nmrk 
 end if
 
 call retrev("cmd_auxPbetaprior",dum_char,iflen,ier)
 call delim(dum_char,dum_char_array(1:2),2,igot,ibegin,iterm,ilen,' ,:')
 call string_to_real(dum_char_array(1),P_beta_params(1),ier)
 call string_to_real(dum_char_array(2),P_beta_params(2),ier)

 call retrev("cmd_esttaubeta",dum_char,iflen,ier)
 estim_tau_beta=.false. ; if(iflen==0) estim_tau_beta=.true.
 tau_beta0=rget("cmd_tauB0") 
 if(tau_beta0<0) then
  up_betai_coop=.true.
 else
  up_betai_coop=.false.
 end if

!!!!!!!!!!!!!!!!!!!!
!!!les modeles:
!!!!!!!!!!!!!!!!!!!!

 call retrev("cmd_efile",dum_char,iflen,ier) 
 !ouvrir et recuperer npheno 
 if(iflen>0) then
  opt_reg=.true. ! par defaut modele IS
  Pheno_file=dum_char(:iflen) !recuperation de npheno
  open(1,file=Pheno_file,status='old') 
  npheno=0
  do 
   read(1,*,end=2,iostat=err) 
   npheno=npheno+1
  end do
  2 continue
  close(1)
  call retrev("cmd_scalecov",dum_char,iflen,ier)
  scale_cov=.false. ; if(iflen==0) scale_cov=.true.
 end if

!is ou mcmc mode?
 call retrev("cmd_covmcmc",dum_char,iflen,ier)
 if(iflen==0 .or. opt_aux) then
  opt_pheno=.true.
  if(opt_reg) then !on avait bien un fichier efile mais on veut desactiver IS
   opt_reg=.false.
  else
    print *,'The covmcmc mode requires a covariate file (-efile option)'
    stop
  end if
 end if

 call retrev("cmd_auxmodel",dum_char,iflen,ier)
 if(iflen==0) then
  opt_aux=.true. ;  opt_pheno=.true.
  if(opt_reg) then !on avait bien un fichier efile mais on veut desactiver IS
   opt_reg=.false.
  else
    print *,'The covmcmc mode requires a covariate file (-efile option)'
    stop
  end if
 end if
 ising_b=rget("cmd_isingbeta") 

 if(opt_pheno) then
  if(estim_omega) then
    print *,'ERROR: MCMC STD or AUX covariate mode is activated but no omega matrix file was provided'
    print *,' Please provide an Omega matrix'
    print *,' AND/OR switch to the IS covariate mode (i.e. remove -covmcmc or -aux option)'
    stop
   end if
 end if

!recuperation du prefixe output
 call retrev("cmd_outprefix",dum_char,iflen,ier)
 if(iflen==0) then
  out_prefix=''
 else 
   out_prefix=trim(dum_char(:iflen))//'_'
 end if 
 sum_yij_pij_file=trim(out_prefix)//'summary_yij_pij.out' ;  sum_pij_file=trim(out_prefix)//'summary_pij.out'
 sum_pi_file=trim(out_prefix)//'summary_pi_xtx.out'  ; sum_omega_file=trim(out_prefix)//'summary_lda_omega.out'
 sum_mat_omega_file=trim(out_prefix)//'mat_omega.out' ; sum_betaparams_file=trim(out_prefix)//'summary_beta_params.out'
 sum_betai_file=trim(out_prefix)//'summary_betai.out' ; sum_Pdelta_file=trim(out_prefix)//'summary_Pdelta.out'
 sum_taufile=trim(out_prefix)//'summary_tau.out' ; sum_beta_i_reg_file=trim(out_prefix)//'summary_betai_reg.out'
 sum_DIC_file=trim(out_prefix)//'DIC.out' ; sum_phenostd_file=trim(out_prefix)//'covariate.std'

!recuperation des threads et initialisation du RNG
 nstream=iget("cmd_nthreads")
 call omp_set_num_threads ( nstream )
 call set_mt19937
 allocate(mts(0:(nstream-1)))
 call new(mts(0))
 call init(mts(0),seed)  ! init by scalar
 do pop=1,nstream-1
    call create_stream(mts(0),mts(pop),pop)
 enddo 

!!!!!!!!!!
!Synthese:
!!!!!!!!!!

!modele basique
 if(.not. opt_pheno .and. .not. opt_reg) then
  print *,' -> Basic Mode (no covariate):'
  if(estim_omega) &
   write(*,'(A,i5)') '  (*) Omega matrix is estimated with prior inv(Omega) ~ Wish_npop((1/rho)Id_npop,rho)  where rho= ',rho
  if(up_params_beta) then 
   write(*,'(A)') '  (*) Parameters from the Pi prior dist. will be estimated'
  else
   write(*,'(A,f6.3,A,f6.3,A)') '  (*) Prior on Pi is Beta(',pi_beta_params(1),',',pi_beta_params(2),')'
  end if
   write(*,'(A)') '  (*) Outlier detection: XtX stat. for each SNP        '
 end if

!modele IS covariate
 if(opt_reg) then
   write(*,'(A,i3,A)') '  -> IS Covariate Mode with ',npheno, ' covariates (i.e. BF and Beta coef will be estimated using an Importance Sampling algorithm):'
  if(estim_omega) &
   write(*,'(A,i5)') '  (*) Omega matrix is estimated with prior inv(Omega) ~ Wish_npop((1/rho)Id_npop,rho)  where rho= ',rho
  if(up_params_beta) then 
   write(*,'(A)') '  (*) Parameters from the Pi prior dist. will be estimated'
  else
   write(*,'(A,f6.3,A,f6.3,A)') '  (*) Prior on Pi is Beta(',pi_beta_params(1),',',pi_beta_params(2),')'
  end if
   write(*,'(A)') '  (*) Outlier detection: XtX stat. for each SNP        '
   write(*,('(A,f8.4,A,f8.4,A)')) '  (*) Prior on Beta coef is U(',min_beta,',',max_beta,')'
   write(*,('(A,i4,A,f8.4,A,f8.4,A)')) '  (*) IS estimates are approximate on the grid with ',ninter_beta,' points uniformly distributed over the (',min_beta,',',max_beta,') interval' 
 end if

!modele MCMC covariate
 if(opt_pheno) then
  write(*,'(A,i3,A)') '  -> MCMC Covariate Mode with ',npheno, ' covariates (outlier detection and estimation of omega are inactivated)'
   if(up_betai_coop) then
    write(*,('(A,f8.4,A,f8.4,A)')) '  (*) Prior on Beta coef is U(',min_beta,',',max_beta,')'
   else
    write(*,('(A,f8.4,A,f8.4,A)')) '  (*) Prior on Beta coef is N(0,',1./tau_beta0,')'
   end if
   if(opt_aux) then
    write(*,('(A,f6.3,A,f6.3,A)')) '  (*) Auxiliary variable model with P~Beta(',P_beta_params(1),',',P_beta_params(2),') => BF will be calculated'
    if(ising_b>-0.0001 .and. ising_b<0.0001) then
     write(*,'(A)') '        => No spatial dependancy between markers (isingbeta=0.)'
     ising_b=0.
     opt_ising=.false.
    else 
     opt_ising=.true.
     write(*,('(A,f6.3)')) '        => Spatial dependancy between markers is modeled via an Ising model on the auxiliary variable with beta_ising=',ising_b
    end if
   else
    write(*,'(A)') '  (*) Full model (BPvalue will be computed)'
   end if
 end if

 print *,''
 print *,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
 print *,'           MCMC specifications'
 print *,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
 print *,''
 print *,' Nb. of sampled parameter values       = ',nvaleurs
 print *,' Thinning Rate                         = ',thin
 print *,' Burn in Period Length                 = ',burn_in
 print *,' Max Number of Pilot runs              = ',npilot
 print *,' Pilot run Length                      = ',pilot_length
 print *,''

 write(*,'(A,1x,f6.4)')  ' Init. deltas for Pi proposals          = ',delta0_pi
 call retrev("cmd_upalphaalt",dum_char,iflen,ier)
 up_alpha_coop=.true. ; if(iflen==0) up_alpha_coop=.false.
 if(up_alpha_coop) then
  write(*,'(A,1x,f6.4)') ' Coop Pij proposal with init. sig. for Pi. multivariate prop.          = ',delta0_pij 
 else
  write(*,'(A,1x,f6.4)') ' Alt. Pij proposal with init. sig. for Pij prop.                       = ',delta0_pij
 end if
 write(*,'(A,1x,2(f4.2,1x))') ' Pilot Run Adj. Factor Pilot            = ',rate_adjust
 write(*,'(A,1x,f4.2,A,f4.2)') ' Targeted Rej./Acc. rates              = ',acc_inf,'-',acc_sup
 write(*,'(A,1x,i6)')          ' R.N.G. seed                            = ',seed 

!!on rentre dans le dur

 if(opt_pheno) then
  allocate(PHENO_VAL(npop,npheno),dum_pheno_std(0:(nstream-1),npop))
  allocate(INITS_Delta(nmrk,npheno),INITS_Beta_i(nmrk,npheno),INITS_tau_beta(npheno),INIT_Pdelta(npheno))
  if(up_betai_coop) then
   allocate(acc_betai(nmrk,npheno),delta_betai(nmrk,npheno))
   delta_betai=0.05 ; acc_betai=0.
  end if
 end if

 if(opt_reg) then !forcement opt_pheno=TRUE donc pas besoin d'allouer
  allocate(PHENO_VAL(npop,npheno),dum_pheno_std(0:(nstream-1),npop))
  allocate(mean_rho_coef(nmrk,npheno,2),mean_isbf(nmrk,npheno,4),grid_is(0:(nstream-1),ninter_beta,npheno+1),dum_vect_grid(0:(nstream-1),ninter_beta))
  mean_rho_coef=0. !pearson
  mean_isbf=0. !1=m_bf ; 2=sd_bf ; 3=m_beta ; 4=sd_beta ! par importance sampling sur grille
  dum_real=(max_beta - min_beta)/(ninter_beta-1.)
  grid_is(:,1,1)=min_beta
  do dum_int=2,ninter_beta
   grid_is(:,dum_int,1)=grid_is(:,dum_int-1,1) + dum_real
  end do
 end if

 if(opt_pheno .or. opt_reg) then
  print *,''
  print *,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
  print *,'        READING COVARIATE DATA'
  print *,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
  print *,''
  open(1,file=Pheno_file,status='old') 
  if(scale_cov) open(10,file=sum_phenostd_file,status='unknown') 
  do pheno=1,npheno 
   read(1,*,iostat=err) PHENO_VAL(1:npop,pheno) 
   dum_real = sum(PHENO_VAL(:,pheno))/npop
   dum_real2 = sqrt((sum(PHENO_VAL(:,pheno)**2)/npop - dum_real**2)*(npop+0.)/(npop-1.))
   write(*,'(A,1x,i3,1x,A,f6.4,1x,A,f6.4,A)') 'Cov:',pheno,' (Mean =',dum_real,' SD =',dum_real2,')'
   if(scale_cov) then
    PHENO_VAL(:,pheno)=(PHENO_VAL(:,pheno)-dum_real )/dum_real2
    write(10,'(1000(f10.6,1x))') PHENO_VAL(:,pheno)
   end if
  end do
  close(1)
  if(scale_cov) close(10)
 end if

 print *,''
 print *,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
 print *,'        READING COUNT DATA'
 print *,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
 print *,''
 allocate(Y_OBS(nmrk,npop), N_OBS(nmrk,npop) ,missing_data(nmrk,npop),log_cn(nmrk,npop),dum_read_line(2*npop),dum_int_pop(npop))
 open(1,file=Geno_file,status='old') 
 do pop=1,npop
  dum_int_pop(pop)=2*pop-1 !on stocke la position
 end do
 do mrk=1,nmrk
  read(1,*,iostat=err) dum_read_line(1:(2*npop))
  Y_OBS(mrk,1:npop) = dum_read_line(dum_int_pop)
  N_OBS(mrk,1:npop) = Y_OBS(mrk,1:npop)  + dum_read_line(dum_int_pop+1) 
 end do
 close(1)
 write(*,'(A,1x,1000(i4,1x))') 'First SNP ref. Allele count ',Y_OBS(1,:) 
 write(*,'(A,1x,1000(i4,1x))') 'First SNP total Gene count  ',N_OBS(1,:) 
 write(*,'(A,1x,1000(i4,1x))') 'Last  SNP ref. Allele count ',Y_OBS(nmrk,:) 
 write(*,'(A,1x,1000(i4,1x))') 'Last  SNP total Gene count  ',N_OBS(nmrk,:)

!calcul des coef n=binomiaux utiles pour la suite (calcule de deviance...)
!a ce stade Y_obs peuvent etre des lectures ou des comptages
 missing_data=.false.
 do mrk=1,nmrk
  do pop=1,npop
   if(N_OBS(mrk,pop)<1) then
    missing_data(mrk,pop)=.true.
    log_cn(mrk,pop)=0.
   else 
    log_cn(mrk,pop)=log_binomial_coef(N_OBS(mrk,pop),Y_OBS(mrk,pop))
   end if
  end do
 end do

 if(poolseq) then !les Y_OBS et N_OBS deviennent des comptades et on introduit les lectures
  allocate(Y_READS(nmrk,npop),N_READS(nmrk,npop),delta_y(nmrk,npop),Yminmax(nmrk,npop,2),acc_y(nmrk,npop))
  Y_READS=Y_OBS ; N_READS=N_OBS ; delta_y=delta0_yij
  do pop=1,npop
   N_OBS(:,pop)=PoolSize(pop)
   do mrk=1,nmrk
    if(missing_data(mrk,pop)) then
      Yminmax(mrk,pop,1)=0 ; Yminmax(mrk,pop,2)=PoolSize(pop)
    else
     if(Y_READS(mrk,pop)/=0 .and. Y_READS(mrk,pop)/=N_READS(mrk,pop)) then
       Yminmax(mrk,pop,1)=1 ; Yminmax(mrk,pop,2)=PoolSize(pop)-1 !ymax=popsize-1 ; ymin=1 
     end if
     if(Y_READS(mrk,pop)==0) then
      Yminmax(mrk,pop,1)=0 ; Yminmax(mrk,pop,2)=PoolSize(pop)-1 !ymax=popsize-1 ; ymin=0
     end if
     if(Y_READS(mrk,pop)==N_READS(mrk,pop)) then
      Yminmax(mrk,pop,1)=1 ; Yminmax(mrk,pop,2)=PoolSize(pop)  ! ymax=popsize ; ymin=1
     end if
    end if
   end do
  end do
 end if

 allocate(INITS_Pij(nmrk,npop,2) , INITS_Pi(nmrk) , INITS_LDA(npop,npop),omega_mat(npop,npop),&
         sum_covariates(nmrk,npop),delta_pi(nmrk),acc_pi(nmrk),muvar_cur(npop,2),&
         dum_real_pop(npop),dum_real_pop2(npop),c_mat(npop,npop),inv_c_mat(npop,npop),&
         array_real_pop(0:(nstream-1),npop,2))
 if(up_alpha_coop) then
   allocate(delta_pij(nmrk,1),acc_pij(nmrk,1)) 
 else
  allocate(delta_pij(nmrk,npop),acc_pij(nmrk,npop),omega_index_subpop(npop,npop-1),reg_mat(npop,npop-1))
  dum_int_pop=(/(tst,tst=1,npop)/)
  do pop=1,npop
   omega_index_subpop(pop,:)=pack(dum_int_pop,dum_int_pop/=pop)
  end do
 end if

 if(.not. estim_omega) then 
  print *,''
  print *,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
  print *,'        READING OMEGA MATRIX'
  print *,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
  print *,''
  open(1,file=omega_file,status='old') 
  do pop=1,npop
   read(1,*,iostat=err)  omega_mat(pop,1:npop) 
  end do
  close(1)
  INITS_LDA=inv(omega_mat)
  c_mat = omega_mat ; call chol(c_mat) 
  inv_c_mat=inv(c_mat)
  if(.not. up_alpha_coop) call compute_reg_mat()
 end if

 delta_pi=delta0_pi ; delta_pij=delta0_pij 

end subroutine void_read_input

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!PRINT HELP PAGE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine void_print_help()

   print *,"Version 2.1"
   print *,''
   print *,'Usage: BayPass [options]'
   print *,''
   print *,'Options:'
   print *,' I)   General Options:'
   write(*,'(A)') '  -help                   Display the help page '
   write(*,'(A)') '  -npop            INT    Number of populations                                           (always required) '
   write(*,'(A)') '  -gfile           CHAR   Name of the Genotyping Data File                                (always required)'
   write(*,'(A)') '  -efile           CHAR   Name of the covariate file: activate covariate mode             (def="")'
   write(*,'(A)') '  -scalecov        CHAR   Scale covariates                                                (def="")'
   write(*,'(A)') '  -poolsizefile    CHAR   Name of the Pool Size file => activate PoolSeq mode             (def="")'
   write(*,'(A)') '  -outprefix       CHAR   Prefix used for the output files                                (def="")'   
   print *,''
   print *,' II)  Modeling Options:'
   write(*,'(A)') '  -omegafile       CHAR   Name of the omega matrix file => inactivate estim. of omega     (def="")'
   write(*,'(A)') '  -rho             INT    Rho parameter of the Wishart prior on omega                     (def=1)'   
   write(*,'(A)') '  -setpibetapar           Inactivate estimation of the Pi beta priors parameters            '
   write(*,'(A)') '  -betapiprior     FLOAT2 Pi Beta prior parameters (if -setpibetapar)                    (def=1.0 1.0)'   
   write(*,'(A)') '  -minbeta         FLOAT  Lower beta coef. for the grid                                   (def=-0.3) '
   write(*,'(A)') '  -maxbeta         FLOAT  Upper beta coef. for the grid                                   (def= 0.3) '
   print *,''
   print *,'  I.1)  IS covariate mode (default covariate mode):'
   write(*,'(A)') '    -nbetagrid       INT    Number of grid points (IS covariate mode)                     (def=201) '
   print *,''
   print *,'  I.2)  MCMC covariate mode:'
   write(*,'(A)') '    -covmcmc                Activate mcmc covariate mode (desactivate estim. of omega) '
   write(*,'(A)') '    -auxmodel               Activate Auxiliary variable mode to estimate BF        '
   write(*,'(A)') '    -isingbeta       FLOAT  Beta (so-called inverse temperature) of the Ising model       (def=0.0)       '   
   write(*,'(A)') '    -auxPbetaprior   FLOAT2 auxiliary P Beta prior parameters                             (def=0.02 1.98)                         '   
!   write(*,'(A)') '    -tauB0           FLOAT  Set Gauss. Prior on beta coef: b~N(0,1/tauB0)              '
!   write(*,'(A)') '  -esttaubeta        ' ! a virer
   print *,''
   print *,' III) MCMC Options:'
   write(*,'(A)') '  -nthreads        INT    Number of threads                                               (def=1) '
   write(*,'(A)') '  -nval            INT    Number of post-burnin and thinned samples to generate           (def=1000) '
   write(*,'(A)') '  -thin            INT    Size of the thinning (record one every thin post-burnin sample) (def=25) '
   write(*,'(A)') '  -burnin          INT    Burn-in length                                                  (def=5000) '
   write(*,'(A)') '  -npilot          INT    Number of pilot runs (to adjust proposal distributions)         (def=20) '
   write(*,'(A)') '  -pilotlength     INT    Pilot run length                                                (def=1000) '
   write(*,'(A)') '  -accinf          FLOAT  Lower target acceptance rate bound                              (def=0.25) '
   write(*,'(A)') '  -accsup          FLOAT  Upper target acceptance rate bound                              (def=0.40) '
   write(*,'(A)') '  -adjrate         FLOAT  Adjustement factor                                              (def=1.25) '
   write(*,'(A)') '  -d0pi            FLOAT  Initial delta for the pi all. freq. proposal                    (def=0.5)'
   write(*,'(A)') '  -upalphaalt             Alternative update of the pij  '
   write(*,'(A)') '  -d0pij           FLOAT  Initial delta for the pij all. freq. proposal (alt. update)     (def=0.05)'
   write(*,'(A)') '  -d0yij           INT    Initial delta for the yij all. count (PoolSeq mode)             (def=1)'
   write(*,'(A)') '  -seed            INT    Random Number Generator seed                                    (def=5001) '

end subroutine void_print_help

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!INITIALISATION
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine void_initialisation()
 !1) initialisation des Pij, Alpha_i et Pi_i
 INITS_Pi=0.0 !on initialise à la moyenne des A_IJ
 do mrk=1,nmrk
  do pop=1,npop
   if(poolseq) then !Cas du poolseq: initialisation des Y_OBS
    if(missing_data(mrk,pop)) then
     Y_OBS(mrk,pop)=PoolSize(pop)/2
    else
     dum_real=(Y_READS(mrk,pop)*PoolSize(pop))/N_READS(mrk,pop)
     Y_OBS(mrk,pop)=nint(dum_real)
     if(Y_READS(mrk,pop)>0 .and. Y_OBS(mrk,pop)==0) Y_OBS(mrk,pop)=1
     if(Y_READS(mrk,pop)<N_READS(mrk,pop) .and. Y_OBS(mrk,pop)==PoolSize(pop)) Y_OBS(mrk,pop)=PoolSize(pop)-1
    end if
   end if
  INITS_Pij(mrk,pop,1)=(Y_OBS(mrk,pop)+1.)/(N_OBS(mrk,pop)+2.)
  end do
  INITS_Pi(mrk)=min(0.999,max(0.001,sum(INITS_Pij(mrk,:,1))/npop))
  INITS_Pij(mrk,:,2)= (INITS_Pij(mrk,:,1) - INITS_Pi(mrk))/sqrt(INITS_Pi(mrk)*(1.-INITS_Pi(mrk)))
 end do

 if(estim_omega) then
 !2) initialisation de LDA: comme un phylogenie en etoile cad diagonale=1/cj
  INITS_LDA=0. ; omega_mat=0.
  do pop=1,npop
    do mrk=1,nmrk
     INITS_LDA(pop,pop)= INITS_LDA(pop,pop) + INITS_Pij(mrk,pop,2)**2 
    end do
    INITS_LDA(pop,pop) = (nmrk+0.)/INITS_LDA(pop,pop)
    omega_mat(pop,pop) = 1./INITS_LDA(pop,pop)
  end do  
  c_mat = omega_mat ; call chol(c_mat) 
  inv_c_mat=inv(c_mat)
  if(.not. up_alpha_coop) call compute_reg_mat()
 end if

!initialisation des parametres de la beta_pi
 if(up_params_beta) then
  tmp_a=sum(INITS_Pi)/(nmrk+0.) ; tmp_b=sum(INITS_Pi*INITS_Pi-tmp_a**2)/(nmrk-1.)
  dum_real=(tmp_a*(1-tmp_a))/tmp_b - 1.
  pi_beta_params(1)=max(0.001,tmp_a*dum_real)  ; pi_beta_params(2)=max(0.001,(1.-tmp_a)*dum_real)  !estimateur des moments
  tmp_b=sum(pi_beta_params)
  pi_beta_update_save(3)=gamma_log(tmp_b)
  pi_beta_update_save(4)=gamma_log(tmp_b*tmp_a) ; pi_beta_update_save(5)=gamma_log(tmp_b*(1-tmp_a))
 end if
 
 sum_covariates=0.
 if(opt_pheno) then
  INITS_Delta = 1 ; INITS_Beta_i=0. ; INIT_Pdelta=P_beta_params(1)/sum(P_beta_params) ; INITS_tau_beta = tau_beta0
!~   if(opt_aux) INITS_tau_beta=INITS_tau_beta*INIT_Pdelta !=>on suppose que la transformation est deja dans la valeur passée en commande
 end if
 
end subroutine void_initialisation

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!PRINT node values
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine void_printval()
 if(poolseq) then
  tmp_min=minval(Y_OBS+0.) ; tmp_max=maxval(Y_OBS+0.)
  write(*,'(A,(1x,f5.1,2x),A,(1x,f5.1,2x),A,(1x,f5.1))')   'Yij min         =',tmp_min,'Yij max         =',tmp_max,&
                                                           'Average Yij     =',sum(Y_OBS+0.)/(nmrk*npop)
 end if
 tmp_min=minval(INITS_Pij(:,:,1)) ; tmp_max=maxval(INITS_Pij(:,:,1))
 write(*,'(A,(1x,f12.5,2x),A,(1x,f12.5,2x),A,(1x,f12.5))') 'Pij min         =',tmp_min,'Pij max         =',tmp_max,&
                                                           'Average Pij     =',sum(INITS_Pij(:,:,1))/(nmrk*npop)
 tmp_min=minval(INITS_Pij(:,:,2)) ; tmp_max=maxval(INITS_Pij(:,:,2))
 write(*,'(A,(1x,f12.5,2x),A,(1x,f12.5,2x),A,(1x,f12.5))') 'Pij_til min     =',tmp_min,'Pij_til max     =',tmp_max,&
                                                           'Average Pij_til =',sum(INITS_Pij(:,:,2))/(nmrk*npop)
 tmp_min=minval(INITS_Pi) ; tmp_max=maxval(INITS_Pi)
 write(*,'(A,(1x,f12.5,2x),A,(1x,f12.5,2x),A,(1x,f12.5))') 'PI min          =',tmp_min,'PI max          =',tmp_max,&
                                                           'Average PI      =',sum(INITS_Pi)/nmrk

 if(up_params_beta) write(*,'(A,1x,2(f12.5,1x))')          'Pi_beta_params values ',pi_beta_params

!~  if(estim_omega) then
!~   print *,'OMEGA node values '
!~   do pop=1,npop
!~    write(*,'(1000(f12.8,1x))') ( omega_mat(pop,dum_int), dum_int=1,npop )
!~   end do
!~  end if

end subroutine void_printval

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!AJUSTEMENT DES PROPOSALS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine void_adjust_proposals()
 tst_pij=0 ; tst_pi=0 ; tst_y=0 ; tst_b=0

 if(poolseq) then  !ajustement des y_ij
  do mrk=1,nmrk
   do pop=1,npop
    if(acc_y(mrk,pop)/pilot_length>acc_sup) then
     delta_y(mrk,pop)=nint(min(PoolSize(pop)+0.,delta_y(mrk,pop)+1.)) ; tst_y=tst_y+1
    end if
    if(acc_y(mrk,pop)/pilot_length<acc_inf) then
     delta_y(mrk,pop)=nint(max(1.,delta_y(mrk,pop)-1.)) ; tst_y=tst_y+1
    end if
   end do
  end do
 write(*,'(A,(1x,i8,2x),A,(1x,i8,2x),A)') '  Acceptance Rate not achieved for ',tst_y,' out of ',nmrk*npop ,' y'
 end if

!ajustement des alpha proposals
 dum_int=npop
 if(up_alpha_coop) dum_int=1
 do mrk=1,nmrk
  do pop=1,dum_int
   if(acc_pij(mrk,pop)/pilot_length>acc_sup) then
     delta_pij(mrk,pop)=delta_pij(mrk,pop)*rate_adjust ; tst_pij=tst_pij+1
   end if
   if(acc_pij(mrk,pop)/pilot_length<acc_inf) then
    delta_pij(mrk,pop)=delta_pij(mrk,pop)/rate_adjust ; tst_pij=tst_pij+1
   end if
  end do
 end do
 if(up_alpha_coop) then
  write(*,'(A,(1x,i8,2x),A,(1x,i8,2x),A)') '  Acceptance Rate not achieved for ',tst_pij,' out of ',nmrk ,' p_ij vectors'
 else
  write(*,'(A,(1x,i8,2x),A,(1x,i8,2x),A)') '  Acceptance Rate not achieved for ',tst_pij,' out of ',nmrk*npop ,' p_ij'
 end if

!ajustement des pi proposals (le cas echeant)
! if(.not. up_slc_pi) then
  do mrk=1,nmrk
   if(acc_pi(mrk)/pilot_length>acc_sup) then
    delta_pi(mrk)=delta_pi(mrk)*rate_adjust ; tst_pi=tst_pi+1
   end if    
   if(acc_pi(mrk)/pilot_length<acc_inf) then
    delta_pi(mrk)=delta_pi(mrk)/rate_adjust ; tst_pi=tst_pi+1
   end if  
  end do
  write(*,'(A,(1x,i8,2x),A,(1x,i8,2x),A)') '  Acceptance Rate not achieved for ',tst_pi,' out of ',nmrk ,' pi'
! end if

 if(opt_pheno .and. up_betai_coop) then
  do mrk=1,nmrk
   do pheno=1,npheno
    if(acc_betai(mrk,pheno)/pilot_length>acc_sup) then
     delta_betai(mrk,pheno)=delta_betai(mrk,pheno)*rate_adjust ; tst_b=tst_b+1
    end if    
    if(acc_betai(mrk,pheno)/pilot_length<acc_inf) then
     delta_betai(mrk,pheno)=delta_betai(mrk,pheno)/rate_adjust ; tst_b=tst_b+1
    end if  
   end do
  end do
  write(*,'(A,(1x,i8,2x),A,(1x,i8,2x),A)') '  Acceptance Rate not achieved for ',tst_b,' out of ',nmrk*npheno ,' betai'
 end if


 if(poolseq) write(*,'(A,(1x,f7.5,2x),A,(1x,f6.1))') '      Mean Acceptance Rate Y  = ',sum(acc_y)/(pilot_length*nmrk*npop),&
                                             ' mean delta_y= ',sum(delta_y+0.)/(nmrk*npop)
 if(up_alpha_coop) then
  write(*,'(A,(1x,f7.5,2x),A,(1x,f7.5))') '      Mean Acceptance Rate Pij  = ',sum(acc_pij)/(pilot_length*nmrk),&
                                         ' mean delta_p= ',sum(delta_pij)/(nmrk)
 else
  write(*,'(A,(1x,f7.5,2x),A,(1x,f7.5))') '      Mean Acceptance Rate Pij  = ',sum(acc_pij)/(pilot_length*nmrk*npop),&
                                         ' mean delta_p= ',sum(delta_pij)/(nmrk*npop)
 end if
! if(.not. up_slc_pi) write(*,'(A,(1x,f7.5,2x),A,(1x,f7.5))') '      Mean Acceptance Rate Pi   = ',sum(acc_pi)/(pilot_length*nmrk),&
!                                         ' mean delta_pi= ',sum(delta_pi)/nmrk
 write(*,'(A,(1x,f7.5,2x),A,(1x,f7.5))') '      Mean Acceptance Rate Pi   = ',sum(acc_pi)/(pilot_length*nmrk),&
                                         ' mean delta_pi= ',sum(delta_pi)/nmrk

 if(opt_pheno .and. up_betai_coop) then
   write(*,'(A,(1x,f7.5,2x),A,(1x,f7.5))') '      Mean Acceptance Rate Betai   = ',sum(acc_betai)/(pilot_length*nmrk*npheno),&
                                          ' mean delta_betai= ',sum(delta_betai)/(nmrk*npheno)
 end if


 if(up_params_beta) then !ajustement des beta_pi_mu et beta_pi_phi
  if(acc_beta_params(1)/pilot_length>acc_sup) then
     delta_pi_beta_mu=delta_pi_beta_mu*rate_adjust
  end if    
  if(acc_beta_params(1)/pilot_length<acc_inf) then
     delta_pi_beta_mu=delta_pi_beta_mu/rate_adjust
  end if 
 write(*,'(A,(1x,f7.5,2x),A,(1x,f7.5))') '      Mean Acceptance Rate b_mu = ',acc_beta_params(1)/pilot_length,&
                                         '      delta  = ',delta_pi_beta_mu
  if(acc_beta_params(2)/pilot_length>acc_sup) then
    delta_pi_beta_phi=delta_pi_beta_phi*rate_adjust
  end if    
  if(acc_beta_params(2)/pilot_length<acc_inf) then
    delta_pi_beta_phi=delta_pi_beta_phi/rate_adjust
  end if   
 write(*,'(A,(1x,f7.5,2x),A,(1x,f7.5))') '      Mean Acceptance Rate b_phi= ',acc_beta_params(2)/pilot_length,&
                                        '      delta  = ',delta_pi_beta_phi
 end if

end subroutine void_adjust_proposals

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!MCMC_ITER
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine void_mcmc_iter()
 if(poolseq) then
 !$OMP PARALLEL private(pop,cur_stream,dum_int,y_out)
 !$OMP do schedule(guided) 
  do mrk=1,nmrk
   cur_stream=omp_get_thread_num()
   do pop=1,npop
!! ATTENTION PB dans fonction random_bionmial2: de toute facon pas necessaire (OK avec MH car y_lect=n_lect=0 si missing_data)
!~     if(missing_data(mrk,pop)) then
!~      acc_y(mrk,pop)=acc_y(mrk,pop)+1.
!~      if(INITS_Pij(mrk,pop,1)>0. .and. INITS_Pij(mrk,pop,1)<1.) then
!~       dum_int=random_binomial2(PoolSize(pop),INITS_Pij(mrk,pop,1),mts(cur_stream))
!~       Y_OBS(mrk,pop)=dum_int
!~       if(Y_OBS(mrk,pop)<0) print *,mrk,pop,Y_OBS(mrk,pop),missing_data(mrk,pop),INITS_Pij(mrk,pop,1),PoolSize(pop),dum_int 
!~      else
!~       if(INITS_Pij(mrk,pop,1)<1e-16) Y_OBS(mrk,pop)=0
!~       if(INITS_Pij(mrk,pop,1)>0.9999) Y_OBS(mrk,pop)=PoolSize(pop)    
!~      end if
!~     else
     call update_y(mts(cur_stream),Y_OBS(mrk,pop),PoolSize(pop),Y_READS(mrk,pop),N_READS(mrk,pop),&
                   Yminmax(mrk,pop,:),max(0.,min(1.,INITS_Pij(mrk,pop,1))),delta_y(mrk,pop),dum_int,y_out)
     acc_y(mrk,pop)=acc_y(mrk,pop)+dum_int ; Y_OBS(mrk,pop)=y_out
!~     end if
!~     if(Y_OBS(mrk,pop)<0) print *,mrk,pop,Y_OBS(mrk,pop),missing_data(mrk,pop),INITS_Pij(mrk,pop,1),PoolSize(pop),random_binomial2(PoolSize(pop),INITS_Pij(mrk,pop,1),mts(cur_stream)) 
   end do
  end do
  !$OMP END do   
  !$OMP END PARALLEL  
 end if

  if(up_alpha_coop) then
 !$OMP PARALLEL private(cur_stream,dum_int)
 !$OMP do schedule(guided) 
   do mrk=1,nmrk
    cur_stream=omp_get_thread_num()
!    print *,cur_stream,size(dum_real_pop2) => illustre probleme avec ifort et -openmp: siz(dum_real_pop2)=0!!!
    call update_alpha_vect(mts(cur_stream),Y_OBS(mrk,:),N_OBS(mrk,:),INITS_Pij(mrk,:,1),INITS_Pij(mrk,:,2),INITS_Pi(mrk),&
                           sum_covariates(mrk,:),delta_pij(mrk,1),INITS_LDA,c_mat,dum_int,array_real_pop(cur_stream,:,1),array_real_pop(cur_stream,:,2))
     INITS_Pij(mrk,:,1)=array_real_pop(cur_stream,:,1) ; INITS_Pij(mrk,:,2)=array_real_pop(cur_stream,:,2) ; acc_pij(mrk,1)=acc_pij(mrk,1)+dum_int 
    end do
  !$OMP END do   
  !$OMP END PARALLEL  
  else
   do pop=1,npop
    alpha_up_muvar_param(2)=omega_mat(pop,pop)-dot_product(reg_mat(pop,:),omega_mat(pop,omega_index_subpop(pop,:)))
 !$OMP PARALLEL private(cur_stream,dum_int,dum_real,dum_real2)
 !$OMP do schedule(guided) 
   do mrk=1,nmrk
    cur_stream=omp_get_thread_num()
    alpha_up_muvar_param(1)=dot_product(reg_mat(pop,:),INITS_Pij(mrk,omega_index_subpop(pop,:),2))
    call update_alpha (mts(cur_stream),Y_OBS(mrk,pop),N_OBS(mrk,pop),INITS_Pij(mrk,pop,1),INITS_Pij(mrk,pop,2),INITS_Pi(mrk),&
                       sum_covariates(mrk,pop),alpha_up_muvar_param(1),alpha_up_muvar_param(2),delta_pij(mrk,pop),dum_int,dum_real,dum_real2)
     INITS_Pij(mrk,pop,1)=dum_real ; INITS_Pij(mrk,pop,2)=dum_real2 ; acc_pij(mrk,pop)=acc_pij(mrk,pop)+dum_int 
    end do
  !$OMP END do   
  !$OMP END PARALLEL 
   end do
  end if

!  if(up_slc_pi) then
!   do mrk=1,nmrk
!    INITS_Pi(mrk)=slc_pi(INITS_Pi(mrk),INITS_Pij(mrk,:,:),INITS_LDA,pi_beta_params) 
!    INITS_Pij(mrk,:,2)= (INITS_Pij(mrk,:,1) - INITS_Pi(mrk))/sqrt(INITS_Pi(mrk)*(1-INITS_Pi(mrk)))
!   end do 
!  else
 !$OMP PARALLEL private(cur_stream,dum_int,dum_real)
 !$OMP do schedule(guided) 
   do mrk=1,nmrk
    cur_stream=omp_get_thread_num()
    call update_pi(mts(cur_stream),INITS_Pi(mrk),INITS_Pij(mrk,:,:),sum_covariates(mrk,:),INITS_LDA,pi_beta_params,delta_pi(mrk),dum_int,dum_real,array_real_pop(cur_stream,:,1))
    acc_pi(mrk)=acc_pi(mrk)+dum_int ; INITS_Pi(mrk)=dum_real ; INITS_Pij(mrk,:,2)= array_real_pop(cur_stream,:,1)
   end do 
  !$OMP END do   
  !$OMP END PARALLEL 

  if(up_params_beta) then
   cur_stream=0
   pi_beta_update_save(1:2)=0.
   do mrk=1,nmrk
     pi_beta_update_save(1)=pi_beta_update_save(1) + log(INITS_Pi(mrk))
     pi_beta_update_save(2)=pi_beta_update_save(2) + log(1.-INITS_Pi(mrk))
   end do
    call update_beta_pi_params(mts(cur_stream),pi_beta_params,nmrk,pi_beta_update_save,delta_pi_beta_mu,&
                                delta_pi_beta_phi,dum_int,dum_int2,dum_real5)
    acc_beta_params(1)=acc_beta_params(1) + dum_int ; acc_beta_params(2)=acc_beta_params(2) + dum_int2
    pi_beta_params=dum_real5(1:2) ; pi_beta_update_save(3:5)=dum_real5(3:5)
  end if

  if(estim_omega) then
   cur_stream=0
   INITS_LDA=up_lda_mat() ; omega_mat=inv(INITS_LDA)
   c_mat=omega_mat ; call chol(c_mat)
   inv_c_mat=inv(c_mat)
   if(.not. up_alpha_coop) call compute_reg_mat()
  end if

  if(opt_pheno) then
   do pheno=1,npheno
    if(up_betai_coop) then
 !$OMP PARALLEL private(cur_stream,dum_int,dum_real)
 !$OMP do schedule(guided) 
     do mrk=1,nmrk
      cur_stream=omp_get_thread_num()
      call up_beta_i_coop(mts(cur_stream),INITS_Beta_i(mrk,pheno),npop,INITS_Pij(mrk,:,1),inv_c_mat,INITS_Delta(mrk,pheno),PHENO_VAL(:,pheno),&
                          sum_covariates(mrk,:),INITS_Pi(mrk),min_beta,max_beta,delta_betai(mrk,pheno),dum_real,array_real_pop(cur_stream,:,1),dum_int) 
      INITS_Beta_i(mrk,pheno)=dum_real ; sum_covariates(mrk,:)=array_real_pop(cur_stream,:,1) ; acc_betai(mrk,pheno)=acc_betai(mrk,pheno) + dum_int
     end do
  !$OMP END do   
  !$OMP END PARALLEL 
    else
 !$OMP PARALLEL private(cur_stream,dum_real)
 !$OMP do schedule(guided) 
     do mrk=1,nmrk
      cur_stream=omp_get_thread_num()
      call up_beta_i(mts(cur_stream),INITS_Beta_i(mrk,pheno),npop,INITS_Pij(mrk,:,1),inv_c_mat,INITS_Delta(mrk,pheno),PHENO_VAL(:,pheno),&
                          sum_covariates(mrk,:),INITS_Pi(mrk),INITS_tau_beta(pheno),dum_real,array_real_pop(cur_stream,:,1)) 
      INITS_Beta_i(mrk,pheno)=dum_real ; sum_covariates(mrk,:)=array_real_pop(cur_stream,:,1)
     end do
  !$OMP END do   
  !$OMP END PARALLEL 
    end if
    if(opt_aux .and. out_pilot) then
     if(opt_ising) then
      cur_stream=0
      do mrk=1,nmrk
       call up_delta_ising(mts(cur_stream),mrk,INITS_Delta(:,pheno),ising_b,INITS_Beta_i(mrk,pheno),npop,nmrk,INIT_Pdelta(pheno),&
                          INITS_Pij(mrk,:,1),inv_c_mat,PHENO_VAL(:,pheno),sum_covariates(mrk,:),INITS_Pi(mrk),dum_int,dum_real_pop)  
       INITS_Delta(mrk,pheno)=dum_int ; sum_covariates(mrk,:)=dum_real_pop
      end do
     else
 !$OMP PARALLEL private(cur_stream,dum_int)
 !$OMP do schedule(guided) 
      do mrk=1,nmrk
       cur_stream=omp_get_thread_num()
       call up_delta(mts(cur_stream),INITS_Delta(mrk,pheno),INITS_Beta_i(mrk,pheno),npop,INIT_Pdelta(pheno),INITS_Pij(mrk,:,1),inv_c_mat,&
                     PHENO_VAL(:,pheno),sum_covariates(mrk,:),INITS_Pi(mrk),dum_int,array_real_pop(cur_stream,:,1))  
       INITS_Delta(mrk,pheno)=dum_int ; sum_covariates(mrk,:)=array_real_pop(cur_stream,:,1)
      end do
  !$OMP END do   
  !$OMP END PARALLEL 
     end if
    end if
   end do

    if(estim_tau_beta) then
     cur_stream=0
     do pheno=1,npheno
      INITS_tau_beta(pheno)=up_tau_beta(mts(cur_stream),INITS_Beta_i(:,pheno),INITS_Delta(:,pheno),k_tau_prior,l_tau_prior)
     end do
    end if
    
    if(opt_aux .and. out_pilot) then
     cur_stream=0 !pas la peine de paralleliser ca
     do pheno=1,npheno
      INIT_Pdelta(pheno)=up_P(mts(cur_stream),INITS_Delta(:,pheno),P_beta_params(1),P_beta_params(2))
     end do
    end if

  end if

end subroutine void_mcmc_iter

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!! CALCUL SUMMARY STATS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine void_update_summary()

mean_deviance=(mean_deviance*(iter-1) + deviance(Y_OBS,N_OBS,INITS_Pij(:,:,1),log_cn))/iter
!$OMP PARALLEL private(pop)
!$OMP do schedule(guided) 
do mrk=1,nmrk
 do pop=1,npop
  CPO(mrk,pop)=(CPO(mrk,pop)*(iter-1) + &
               exp(-1.*logl_ij_beta(Y_OBS(mrk,pop),N_OBS(mrk,pop),INITS_Pij(mrk,pop,1),log_cn(mrk,pop))))/iter
 end do
end do
!$OMP END do   
!$OMP END PARALLEL 

!$OMP PARALLEL private(cur_stream,pop,dum_real)
!$OMP do schedule(guided) 
do mrk=1,nmrk
 cur_stream=omp_get_thread_num()
 mean_pi(mrk,1)=(mean_pi(mrk,1)*(iter-1)+INITS_Pi(mrk))/iter
 mean_pi(mrk,2)=(mean_pi(mrk,2)*(iter-1)+(INITS_Pi(mrk))**2)/iter
 !pij standardisée et xtx stat
 array_real_pop(cur_stream,:,1) = matmul(inv_c_mat,INITS_Pij(mrk,:,2))
 dum_real=dot_product(array_real_pop(cur_stream,:,1),array_real_pop(cur_stream,:,1))
! dum_real=dot_product(INITS_Pij(mrk,:,2),matmul(INITS_LDA,INITS_Pij(mrk,:,2)))
 mean_xtx(mrk,1)=(mean_xtx(mrk,1)*(iter-1)+dum_real)/iter
 mean_xtx(mrk,2)=(mean_xtx(mrk,2)*(iter-1)+dum_real**2)/iter
 do pop=1,npop
  if(poolseq) then
    mean_yij(mrk,pop,1)=(mean_yij(mrk,pop,1)*(iter-1)+Y_OBS(mrk,pop))/iter
    mean_yij(mrk,pop,2)=(mean_yij(mrk,pop,2)*(iter-1)+(Y_OBS(mrk,pop)+0.)**2)/iter
  end if
  dum_real=max(0.,min(1.,INITS_Pij(mrk,pop,1)))
  mean_pij(mrk,pop,1)=(mean_pij(mrk,pop,1)*(iter-1)+dum_real)/iter
  mean_pij(mrk,pop,2)=(mean_pij(mrk,pop,2)*(iter-1)+dum_real**2)/iter
  mean_pij(mrk,pop,3)=(mean_pij(mrk,pop,3)*(iter-1)+array_real_pop(cur_stream,pop,1))/iter
  mean_pij(mrk,pop,4)=(mean_pij(mrk,pop,4)*(iter-1)+array_real_pop(cur_stream,pop,1)**2)/iter
 end do 
end do
!$OMP END do   
!$OMP END PARALLEL 

if(up_params_beta) then
 mean_beta_params(1,1)=(mean_beta_params(1,1)*(iter-1)+pi_beta_params(1))/iter
 mean_beta_params(1,2)=(mean_beta_params(1,2)*(iter-1)+(pi_beta_params(1))**2)/iter
 mean_beta_params(2,1)=(mean_beta_params(2,1)*(iter-1)+pi_beta_params(2))/iter
 mean_beta_params(2,2)=(mean_beta_params(2,2)*(iter-1)+(pi_beta_params(2))**2)/iter
end if

 mean_lda(:,:,1)=(mean_lda(:,:,1)*(iter-1)+ INITS_LDA)/iter
 mean_lda(:,:,2)=(mean_lda(:,:,2)*(iter-1)+ INITS_LDA**2)/iter
 mean_omega(:,:,1)=(mean_omega(:,:,1)*(iter-1)+ omega_mat)/iter
 mean_omega(:,:,2)=(mean_omega(:,:,2)*(iter-1)+ omega_mat**2)/iter

 if(opt_pheno) then
  if(opt_aux) then
   do pheno=1,npheno
    do mrk=1,nmrk
     dum_real=INITS_Beta_i(mrk,pheno)*INITS_Delta(mrk,pheno)
     mean_beta_i(mrk,pheno,1)=(mean_beta_i(mrk,pheno,1)*(iter-1)+dum_real)/iter
     mean_beta_i(mrk,pheno,2)=(mean_beta_i(mrk,pheno,2)*(iter-1)+dum_real**2)/iter
     mean_delta_i(mrk,pheno)=(mean_delta_i(mrk,pheno)*(iter-1)+INITS_Delta(mrk,pheno))/iter
    end do
   mean_p_delta(pheno,1)=(mean_p_delta(pheno,1)*(iter-1)+INIT_Pdelta(pheno))/iter
   mean_p_delta(pheno,2)=(mean_p_delta(pheno,2)*(iter-1)+INIT_Pdelta(pheno)**2)/iter
   end do
  else
   do pheno=1,npheno
    do mrk=1,nmrk
     mean_beta_i(mrk,pheno,1)=(mean_beta_i(mrk,pheno,1)*(iter-1)+INITS_Beta_i(mrk,pheno))/iter
     mean_beta_i(mrk,pheno,2)=(mean_beta_i(mrk,pheno,2)*(iter-1)+INITS_Beta_i(mrk,pheno)**2)/iter
    end do
    if(estim_tau_beta) then
     mean_tau_beta(pheno,1)=(mean_tau_beta(pheno,1)*(iter-1)+ INITS_tau_beta(pheno))/iter
     mean_tau_beta(pheno,2)=(mean_tau_beta(pheno,2)*(iter-1)+ INITS_tau_beta(pheno)**2)/iter
    end if 
   end do
  end if
 end if

 if(opt_reg) then
  do pheno=1,npheno
   dum_real_pop2 = matmul(inv_c_mat,PHENO_VAL(:,pheno))
!$OMP PARALLEL private(cur_stream,dum_int,dum_real)
!$OMP do schedule(guided) 
   do mrk=1,nmrk
    cur_stream=omp_get_thread_num()
    dum_pheno_std(cur_stream,:)=dum_real_pop2 / sqrt(INITS_Pi(mrk)*(1-INITS_Pi(mrk)))
    do dum_int=1,ninter_beta
     grid_is(cur_stream,dum_int,pheno+1)= -0.5*(grid_is(cur_stream,dum_int,1)**2)*(sum(dum_pheno_std(cur_stream,:)**2))
    end do
    array_real_pop(cur_stream,:,1) = matmul(inv_c_mat,INITS_Pij(mrk,:,2))
    dum_real= rho_pearson(dum_pheno_std(cur_stream,:),array_real_pop(cur_stream,:,1))
    mean_rho_coef(mrk,pheno,1)=(mean_rho_coef(mrk,pheno,1)*(iter-1)+ dum_real)/iter
    mean_rho_coef(mrk,pheno,2)=(mean_rho_coef(mrk,pheno,2)*(iter-1)+ dum_real**2)/iter
    !calcul des stats par importance
    do dum_int=1,ninter_beta
     dum_vect_grid(cur_stream,dum_int)=exp(grid_is(cur_stream,dum_int,1)*sum(array_real_pop(cur_stream,:,1)*dum_pheno_std(cur_stream,:))  + grid_is(cur_stream,dum_int,pheno+1))
    end do
    dum_real=sum(dum_vect_grid(cur_stream,2:ninter_beta) + dum_vect_grid(cur_stream,1:(ninter_beta-1)))/(2*ninter_beta-2.)  !pour les bf: car prior uniforme (tester si pas uniforme)
    mean_isbf(mrk,pheno,1)=(mean_isbf(mrk,pheno,1)*(iter-1)+ dum_real)/iter
    mean_isbf(mrk,pheno,2)=(mean_isbf(mrk,pheno,2)*(iter-1)+ dum_real**2)/iter
    dum_vect_grid(cur_stream,1:(ninter_beta-1))=dum_vect_grid(cur_stream,2:ninter_beta) + dum_vect_grid(cur_stream,1:(ninter_beta-1))
    dum_real=(sum( dum_vect_grid(cur_stream,1:(ninter_beta-1)) * 0.5 * (grid_is(cur_stream,2:ninter_beta,1) + grid_is(cur_stream,1:(ninter_beta-1),1)) ))/sum(dum_vect_grid(cur_stream,:))
    mean_isbf(mrk,pheno,3)=(mean_isbf(mrk,pheno,3)*(iter-1)+ dum_real)/iter
    mean_isbf(mrk,pheno,4)=(mean_isbf(mrk,pheno,4)*(iter-1)+ dum_real**2)/iter
   end do
!$OMP END do   
!$OMP END PARALLEL 
  end do
 end if

end subroutine void_update_summary

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!! IMPRIMER SUMMARY STATS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine void_print_summary()
 if(up_alpha_coop) then
  write(*,'(A,(1x,f7.5,2x),A,(1x,f7.5))') '      Mean Final Acceptance Rate Pij  = ',sum(acc_pij)/(thin*nvaleurs*nmrk),&
                                          ' mean delta_pij= ',sum(delta_pij)/(nmrk)
 else 
  write(*,'(A,(1x,f7.5,2x),A,(1x,f7.5))') '      Mean Final Acceptance Rate Pij  = ',sum(acc_pij)/(thin*nvaleurs*nmrk*npop),&
                                          ' mean delta_pij= ',sum(delta_pij)/(nmrk*npop)
 end if                                        
 write(*,'(A,(1x,f7.5,2x),A,(1x,f7.5))') '      Mean Final Acceptance Rate Pi   = ',sum(acc_pi)/(thin*nvaleurs*nmrk),&
                                        ' mean delta_pi= ',sum(delta_pi)/nmrk
 if(up_params_beta) then 
  write(*,'(A,(1x,f7.5,2x),A,(1x,f7.5))') '      Mean Final Acceptance Rate b_mu = ',acc_beta_params(1)/(thin*nvaleurs),&
                                         '      delta  = ',delta_pi_beta_mu
  write(*,'(A,(1x,f7.5,2x),A,(1x,f7.5))') '      Mean Final Acceptance Rate b_phi= ',acc_beta_params(2)/(thin*nvaleurs),&
                                         '      delta  = ',delta_pi_beta_phi
 end if
 
 !!!!!!!!!!!
 !!impression summary stats
 !!!!!!!!!!!


 if(estim_omega) then
  open(3,file=sum_omega_file,status='unknown') ;  write (3,*) 'I J M_lda_ij SD_lda_ij M_omega_ij SD_omega_ij'
  open(4,file=sum_mat_omega_file,status='unknown') 
  do pop=1,npop 
   write(4,'(1000(f12.6,1x))') ( mean_omega(pop,dum_int,1), dum_int=1,npop ) 
   do mrk=1,npop
    write(3,'(2(i5,1x),4(f15.8,1x))') pop,mrk,mean_lda(pop,mrk,1),sqrt(mean_lda(pop,mrk,2)-(mean_lda(pop,mrk,1))**2),&
                                              mean_omega(pop,mrk,1),sqrt(mean_omega(pop,mrk,2)-(mean_omega(pop,mrk,1))**2)
   end do
  end do
 end if
 close(3) ; close(4)

!Impression des données pop et locus
 if(poolseq) then
   open(1,file=sum_yij_pij_file,status='unknown')
   write (1,*) 'POP MRK M_Y SD_Y M_P SD_P M_Pstd SD_Pstd DELTA_Y ACC_Y DELTA_P ACC_P'
 else
   open(1,file=sum_pij_file,status='unknown')
   write (1,*) 'POP MRK M_P SD_P M_Pstd SD_Pstd DELTA_P ACC_P'
 end if
 
 open(2,file=sum_pi_file,status='unknown') 
 write (2,*) 'MRK M_P SD_P DELTA_P ACC_P M_XtX SD_XtX'
 dum_int=npop ; if(up_alpha_coop) dum_int=1
 do mrk=1,nmrk
   write(2,'(1i8,6(f12.8,1x))') mrk,mean_pi(mrk,1),sqrt(mean_pi(mrk,2)-(mean_pi(mrk,1))**2),&
                                delta_pi(mrk),acc_pi(mrk)/(nvaleurs*thin),&
                                mean_xtx(mrk,1),sqrt(mean_xtx(mrk,2)-(mean_xtx(mrk,1))**2)
  do pop=1,npop
   dum_int=pop ; if(up_alpha_coop) dum_int=1
   if(poolseq) then
     write(1,'(i4,1x,i8,2(f8.2,1x),6(f12.8,1x),1x,i5,1x,f6.4)') pop,mrk,mean_yij(mrk,pop,1),sqrt(mean_yij(mrk,pop,2)-(mean_yij(mrk,pop,1))**2),&
                                        mean_pij(mrk,pop,1),sqrt(mean_pij(mrk,pop,2)-(mean_pij(mrk,pop,1))**2),&
                                        mean_pij(mrk,pop,3),sqrt(mean_pij(mrk,pop,4)-(mean_pij(mrk,pop,3))**2),&
                                        delta_pij(mrk,dum_int),acc_pij(mrk,dum_int)/(nvaleurs*thin),&
                                        delta_y(mrk,pop),acc_y(mrk,pop)/(nvaleurs*thin)
   else
     write(1,'(i4,1x,i8,6(f12.8,1x))') pop,mrk,mean_pij(mrk,pop,1),sqrt(mean_pij(mrk,pop,2)-(mean_pij(mrk,pop,1))**2),&
                                        mean_pij(mrk,pop,3),sqrt(mean_pij(mrk,pop,4)-(mean_pij(mrk,pop,3))**2),&
                                        delta_pij(mrk,dum_int),acc_pij(mrk,dum_int)/(nvaleurs*thin)
   end if
  end do
 end do
  close(1) ; close(2)

  if(up_params_beta) then
    open(4,file=sum_betaparams_file,status='unknown') ;  write (4,*) 'PARAM Mean SD'
    write(4,'(A,1x,2(f12.6,1x))') 'a_beta_pi ',mean_beta_params(1,1),sqrt(mean_beta_params(1,2)-(mean_beta_params(1,1))**2)
    write(4,'(A,1x,2(f12.6,1x))') 'b_beta_pi ',mean_beta_params(2,1),sqrt(mean_beta_params(2,2)-(mean_beta_params(2,1))**2)
    close(4)
  end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!cas avec covariables
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 if(opt_pheno) then
  if(.not. estim_tau_beta) then 
   mean_tau_beta(:,1)=tau_beta0 ; mean_tau_beta(:,2)=tau_beta0**2 
  end if
  if(opt_aux) then !Bayes Factor
   open(25,file=sum_betai_file,status='unknown') !on ne met pas les accRate car pas de sens
   write (25,*) 'COVARIABLE MRK M_Beta SD_Beta M_Delta BF(dB) '
   open(26,file=sum_Pdelta_file,status='unknown')
   write (26,*) 'COVARIABLE M_P SD_P'!' M_Tau SD_tau'
   dum_real=P_beta_params(1)/sum(P_beta_params) ; dum_real=(1.-dum_real)/dum_real !inverse prior odd
!~    dum_real=(1.-mean_p_delta(pheno,1))/mean_p_delta(pheno,1) 
   do pheno=1,npheno
    write(26,'(i3,1x,4(f12.8,1x))') pheno,mean_p_delta(pheno,1),sqrt(mean_p_delta(pheno,2)-mean_p_delta(pheno,1)**2)!,mean_tau_beta(pheno,1),sqrt(mean_tau_beta(pheno,2)-mean_tau_beta(pheno,1)**2)
    do mrk=1,nmrk
     tmp_mean=mean_delta_i(mrk,pheno)/(1.-mean_delta_i(mrk,pheno)) !posterior odds
     if(mean_delta_i(mrk,pheno)<1./nvaleurs) tmp_mean=0.5/(nvaleurs-0.5)
     if(mean_delta_i(mrk,pheno)>(nvaleurs-1.)/nvaleurs) tmp_mean=2.*(nvaleurs-0.5)
     bf=log10(dum_real*tmp_mean)
     write(25,'(i3,1x,i7,1x,10(f12.8,1x))') pheno,mrk,mean_beta_i(mrk,pheno,1),&
                                          sqrt(mean_beta_i(mrk,pheno,2)-(mean_beta_i(mrk,pheno,1))**2),&
                                          mean_delta_i(mrk,pheno),10.*bf
    end do
   end do
   close(26) ; close(25)
  else !KLD et BPval
   if(up_betai_coop) then
    open(25,file=sum_betai_file,status='unknown')
    write (25,*) 'COVARIABLE MRK M_Beta SD_Beta AccRateB DeltaB eBPmc'
   else
    open(25,file=sum_betai_file,status='unknown')
    write (25,*) 'COVARIABLE MRK M_Beta SD_Beta KLD logBPval'
   end if
   if(estim_tau_beta) then
    open(26,file=sum_taufile,status='unknown')
    write (26,*) 'COVARIABLE M_Tau SD_tau'
   end if
   do pheno=1,npheno
    if(estim_tau_beta) then
      write(26,'(i3,1x,4(f12.8,1x))') pheno,mean_tau_beta(pheno,1),sqrt(mean_tau_beta(pheno,2)-mean_tau_beta(pheno,1)**2)
    end if
    do mrk=1,nmrk
     bpval=mean_beta_i(mrk,pheno,1)/sqrt(mean_beta_i(mrk,pheno,2)-(mean_beta_i(mrk,pheno,1))**2)
     bpval=min(8.1,max(bpval,-8.1))
     call normal_01_cdf(bpval,dum_real)
     bpval=-1.*log10(1.-2.*abs(0.5-dum_real))
     if(up_betai_coop) then
      write(25,'(i3,1x,i7,1x,10(f12.8,1x))') pheno,mrk,mean_beta_i(mrk,pheno,1),&
                                           acc_betai(mrk,pheno)/(nvaleurs*thin),delta_betai(mrk,pheno),&
                                           sqrt(mean_beta_i(mrk,pheno,2)-(mean_beta_i(mrk,pheno,1))**2),bpval
     else
      dum_real=(mean_beta_i(mrk,pheno,2)-(mean_beta_i(mrk,pheno,1))**2)*mean_tau_beta(pheno,1)
      kld=-1.*log(dum_real) + dum_real
      dum_real= (mean_beta_i(mrk,pheno,2)**2)*mean_tau_beta(pheno,1)
      kld=kld+dum_real-1. ; kld=0.5*kld
      write(25,'(i3,1x,i7,1x,8(f12.8,1x))') pheno,mrk,mean_beta_i(mrk,pheno,1),&
                                           sqrt(mean_beta_i(mrk,pheno,2)-(mean_beta_i(mrk,pheno,1))**2),&
                                           kld,bpval
     end if
    end do
   end do
   close(25)
   if(estim_tau_beta) close(26)
  end if
 end if

 if(opt_reg) then
  open(25,file=sum_beta_i_reg_file,status='unknown')
  write (25,*) 'COVARIABLE MRK M_Pearson SD_Pearson BF(dB) Beta_is SD_Beta_is eBPis' !pas de sens de mettre SD_BF
  do pheno=1,npheno
   do mrk=1,nmrk
    bpval=mean_isbf(mrk,pheno,3)/sqrt(mean_isbf(mrk,pheno,4)-(mean_isbf(mrk,pheno,3))**2)
    bpval=min(8.1,max(bpval,-8.1))
    call normal_01_cdf(bpval,dum_real)
    bpval=-1.*log10(1.-2.*abs(0.5-dum_real)) 
    write(25,'(i3,1x,i7,1x,8(f14.8,1x))') pheno,mrk,mean_rho_coef(mrk,pheno,1),&
                                       sqrt(mean_rho_coef(mrk,pheno,2)-(mean_rho_coef(mrk,pheno,1))**2),&
                                       10.*log10(mean_isbf(mrk,pheno,1)),&
                                       mean_isbf(mrk,pheno,3),sqrt(mean_isbf(mrk,pheno,4)-(mean_isbf(mrk,pheno,3))**2),bpval
   end do
  end do
  close(25)
 end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !impression deviance
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 lpml=0.
 do mrk=1,nmrk
  do pop=1,npop
   lpml=lpml - log(CPO(mrk,pop))
  end do
 end do
 open(1002,file=sum_DIC_file,status='unknown')
 write (1002,*) 'bar(D) pD DIC LPML'
 dum_real16=deviance(Y_OBS,N_OBS,mean_pij(:,:,1),log_cn)
 dum_real16_2=mean_deviance-dum_real16
 write(1002,'(4(f20.2,1x))') mean_deviance,dum_real16_2,dum_real16+2.*dum_real16_2,lpml
 close(1002)

end subroutine void_print_summary

end program baypass



