!!!!!!!!!!!!!!!!!!!!!!!!
!!Most of the functions and subroutines from this module were adapted from
!!Alan Miller Fortran module (random.f90) available in particular at http://jblevins.org/mirror/amiller/
!!!!!!!!!!!!!!!!!!!!!!!!

module mcmc_utils
use mt_stream
use mt_kind_defs
!~     integer, parameter :: defaultsd = 4357
!~     integer, parameter :: N = 624, N1 = N + 1
!~     integer, save, dimension(0:N-1) :: mt
!~     integer, save                   :: mti = N1
contains

function grnd(this)
implicit none
 type(mt_state), intent(inout) :: this
 real (kind=8) :: grnd
 grnd=genrand_double3(this)
 return
end function grnd

!____________________________________________________________________________
! A C-program for MT19937: Real number version
!   genrand() generates one pseudorandom real number (double)
! which is uniformly distributed on [0,1]-interval, for each
! call. sgenrand(seed) set initial values to the working area
! of 624 words. Before genrand(), sgenrand(seed) must be
! called once. (seed is any 32-bit integer except for 0).
! Integer generator is obtained by modifying two lines.
!   Coded by Takuji Nishimura, considering the suggestions by
! Topher Cooper and Marc Rieffel in July-Aug. 1997.
!
! This library is free software; you can redistribute it and/or
! modify it under the terms of the GNU Library General Public
! License as published by the Free Software Foundation; either
! version 2 of the License, or (at your option) any later
! version.
! This library is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
! See the GNU Library General Public License for more details.
! You should have received a copy of the GNU Library General
! Public License along with this library; if not, write to the
! Free Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
! 02111-1307  USA
!
! Copyright (C) 1997 Makoto Matsumoto and Takuji Nishimura.
! When you use this, send an email to: matumoto@math.keio.ac.jp
! with an appropriate reference to your work.
!
!***********************************************************************
! Fortran translation by Hiroshi Takano.  Jan. 13, 1999.
!
!   genrand()      -> double precision function grnd()
!   sgenrand(seed) -> subroutine sgrnd(seed)
!                     integer seed
!
! This program uses the following non-standard intrinsics.
!   ishft(i,n): If n>0, shifts bits in i by n positions to left.
!               If n<0, shifts bits in i by n positions to right.
!   iand (i,j): Performs logical AND on corresponding bits of i and j.
!   ior  (i,j): Performs inclusive OR on corresponding bits of i and j.
!   ieor (i,j): Performs exclusive OR on corresponding bits of i and j.
!
!***********************************************************************
! Fortran version rewritten as an F90 module and mt state saving and getting
! subroutines added by Richard Woloshyn. (rwww@triumf.ca). June 30, 1999

!Initialization subroutine
!~   subroutine sgrnd(seed)
!~     implicit none
!~ !
!~ !      setting initial seeds to mt[N] using
!~ !      the generator Line 25 of Table 1 in
!~ !      [KNUTH 1981, The Art of Computer Programming
!~ !         Vol. 2 (2nd Ed.), pp102]
!~ !
!~     integer, intent(in) :: seed
!~ 
!~     mt(0) = iand(seed,-1)
!~     do mti=1,N-1
!~       mt(mti) = iand(69069 * mt(mti-1),-1)
!~     enddo
!~ !
!~     return
!~   end subroutine sgrnd
!~ 
!~ !Random number generator
!~   real(8) function grnd()
!~     implicit integer(a-z)
!~ 
!~     integer, parameter :: M = 397, MATA  = -1727483681
!~     integer, parameter :: LMASK =  2147483647
!~     integer, parameter :: UMASK = -LMASK - 1
!~     integer, parameter :: TMASKB= -1658038656, TMASKC= -272236544
!~ 
!~     dimension mag01(0:1)
!~     data mag01/0, MATA/
!~     save mag01
!~ !                        mag01(x) = x * MATA for x=0,1
!~ 
!~     TSHFTU(y)=ishft(y,-11)
!~     TSHFTS(y)=ishft(y,7)
!~     TSHFTT(y)=ishft(y,15)
!~     TSHFTL(y)=ishft(y,-18)
!~ 
!~     if(mti.ge.N) then
!~ !                       generate N words at one time
!~       if(mti.eq.N+1) then
!~ !                            if sgrnd() has not been called,
!~         call sgrnd( defaultsd )
!~ !                              a default initial seed is used
!~       endif
!~ 
!~       do kk=0,N-M-1
!~           y=ior(iand(mt(kk),UMASK),iand(mt(kk+1),LMASK))
!~           mt(kk)=ieor(ieor(mt(kk+M),ishft(y,-1)),mag01(iand(y,1)))
!~       enddo
!~       do kk=N-M,N-2
!~           y=ior(iand(mt(kk),UMASK),iand(mt(kk+1),LMASK))
!~           mt(kk)=ieor(ieor(mt(kk+(M-N)),ishft(y,-1)),mag01(iand(y,1)))
!~       enddo
!~       y=ior(iand(mt(N-1),UMASK),iand(mt(0),LMASK))
!~       mt(N-1)=ieor(ieor(mt(M-1),ishft(y,-1)),mag01(iand(y,1)))
!~       mti = 0
!~     endif
!~ 
!~     y=mt(mti)
!~     mti = mti + 1 
!~     y=ieor(y,TSHFTU(y))
!~     y=ieor(y,iand(TSHFTS(y),TMASKB))
!~     y=ieor(y,iand(TSHFTT(y),TMASKC))
!~     y=ieor(y,TSHFTL(y))
!~ 
!~     if(y .lt. 0) then
!~       grnd=(dble(y)+2.0d0**32)/(2.0d0**32-1.0d0)
!~     else
!~       grnd=dble(y)/(2.0d0**32-1.0d0)
!~     endif
!~ 
!~     return
!~   end function grnd
!~ 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!FONCTION random_normal (avec appel au mersenne twister
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

FUNCTION random_normal(this)

implicit none
type(mt_state), intent(inout) :: this
REAL (kind=8) :: random_normal,&
         s = 0.449871, t = -0.386595, a = 0.19600, b = 0.25472,    &
         r1 = 0.27597, r2 = 0.27846, half=0.5, x, y, q ,u, v
!     Generate P = (u,v) uniform in rectangle enclosing acceptance region
DO
  u=grnd(this) ; v=grnd(this)
  v = 1.7156 * (v - half)
  x = u - s ; y = ABS(v) - t
  q = x**2 + y*(a*y - b*x)
  IF (q < r1) EXIT
  IF (q > r2) CYCLE
  IF (v**2 < -4.0*LOG(u)*u**2) EXIT
END DO
random_normal = v/u
RETURN

END FUNCTION random_normal

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!FONCTION random_gamma1 s>1 (avec appel au mersenne twister)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

FUNCTION random_gamma1(s,this) RESULT(fn_val)
! Generates a random gamma deviate for shape parameter s >= 1.
implicit none
type(mt_state), intent(inout) :: this
REAL (kind=8), INTENT(IN)    :: s
REAL (kind=8), SAVE  :: c, d
REAL (kind=8)                :: fn_val,u, v, x

  d = s - 1./3. ;  c = 1./SQRT(9.0*d)
! Start of main loop
DO
  DO
    x = random_normal(this) ;  v = (1.0 + c*x)**3
    IF (v > 0.) EXIT
  END DO
  u=grnd(this) 
  IF (u < 1.0 - 0.0331*x**4) THEN
    fn_val = d*v
    EXIT
  ELSE IF (LOG(u) < 0.5*x**2 + d*(1.0 - v + LOG(v))) THEN
    fn_val = d*v
    EXIT
  END IF
END DO

RETURN
END FUNCTION random_gamma1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!FONCTION random_gamma1 s<1 (avec appel au mersenne twister
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

FUNCTION random_gamma2(s,this) RESULT(fn_val)
implicit none
type(mt_state), intent(inout) :: this
REAL (kind=8), INTENT(IN)    :: s
REAL (kind=8), SAVE          :: a, p, c, uf, vr, d
REAL (kind=8)                :: fn_val,r, x, w,vsmall = TINY(1.0)


IF (s <= 0. .OR. s >= 1.) THEN
  WRITE(*, *) 'SHAPE PARAMETER VALUE OUTSIDE PERMITTED RANGE'
  STOP
END IF

  a = 1. - s ; p = a/(a + s*EXP(-a))
  IF (s < vsmall) THEN
    WRITE(*, *) 'SHAPE PARAMETER VALUE TOO SMALL'
    STOP
  END IF
  c = 1./s ;  uf = p*(vsmall/a)**s
  vr = 1. - vsmall ;  d = a*LOG(a)

DO
  r=grnd(this)
  IF (r >= vr) THEN
    CYCLE
  ELSE IF (r > p) THEN
    x = a - LOG((1. - r)/(1. - p)) ;  w = a*LOG(x)-d
  ELSE IF (r > uf) THEN
    x = a*(r/p)**c ;  w = x
  ELSE
    fn_val = 0.
    RETURN
  END IF

  r=grnd(this)
  IF (1.-r <= w .AND. r > 0.) THEN
    IF (r*(w + 1.) >= 1.) CYCLE
    IF (-LOG(r) <= w) CYCLE
  END IF
  EXIT
END DO

fn_val = x
RETURN

END FUNCTION random_gamma2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!Probabilite density function: gaussienne
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

function norm_pdf(x,mean,var)
  implicit none
  real, intent(in) :: x,mean,var
  real, parameter  :: pi = 3.141592653589793
  real :: norm_pdf

  norm_pdf = exp ( -0.5 * (x-mean) * (x-mean) / var ) / sqrt ( 2.0 * pi * var )

 ! return 
end function norm_pdf

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!Probabilite density function: log de la gaussienne simplifie (attention a utiliser dans unr apport pour justifier la siplification)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

function log_norm_pdf(x,mean,var)
  implicit none
  real (kind=8), intent(in) :: x,mean,var
  real (kind=8) :: log_norm_pdf, pi = 3.141592653589793

  log_norm_pdf = -0.5*log(var) - 0.5 * (x-mean) * (x-mean) / var -0.5*log( 2.0 * pi)

  return 
end function log_norm_pdf

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!Coefficient C(n,k) utilse pour la binomiale
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

function binomial_coef (n, k)
  implicit none
  integer, intent(in) :: n,k
  integer :: i,mn,mx
  real  (kind=8) :: binomial_coef !sinon probleme d'overflow!!

  mn = min ( k, n-k )
  if ( mn < 0 ) then
    binomial_coef = 0.
  else if ( mn == 0 ) then
    binomial_coef = 1.
  else
    mx = max ( k, n-k )
    binomial_coef = mx + 1.
    do i = 2, mn
      binomial_coef = ( binomial_coef * ( mx + i ) ) / i
    end do
  end if
end function binomial_coef

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!Log du Coefficient C(n,k) utilse pour la binomiale
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

function log_binomial_coef (n, k)
  implicit none
  integer, intent(in) :: n,k
  integer :: i,mn,mx
  real :: log_binomial_coef 

  mn = min ( k, n-k )
  if ( mn < 0 ) then
    print *,'calcul de logbinom impossible: k>n'
  else if ( mn == 0 ) then
    log_binomial_coef = 0.
  else
    mx = max ( k, n-k )
    log_binomial_coef = log(mx + 1.)
    do i = 2, mn
      log_binomial_coef = log_binomial_coef + log( mx + i + 0.) -log(i+0.)
    end do
  end if
end function log_binomial_coef


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!Coefficient BETA utilse pour la loi beta...
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function beta ( a, b )
  implicit none
  real (kind=8), intent(in) :: a,b
  real (kind=8) :: beta !,gamma_log

  if ( a <= 0.0D+00 .or. b <= 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'BETA - Fatal error!'
    write ( *, '(a)' ) '  Both A and B must be greater than 0.'
    stop
  end if
  beta = exp ( gamma_log ( a ) + gamma_log ( b ) - gamma_log ( a + b ) )

  return
end function

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!Coefficient GAMMA utilse pour la loi beta...
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function gamma_log ( x )
 implicit none
  real ( kind = 8 ), parameter, dimension ( 7 ) :: c = (/ &
    -1.910444077728D-03,8.4171387781295D-04,-5.952379913043012D-04, &
     7.93650793500350248D-04,-2.777777777777681622553D-03,&
     8.333333333333333331554247D-02,5.7083835261D-03 /)
  real    ( kind = 8 ), parameter :: d1 = -5.772156649015328605195174D-01,&
  d2 =  4.227843350984671393993777D-01,d4 =  1.791759469228055000094023D+00,&
  frtbig = 1.42D+09
  real    ( kind = 8 ), parameter, dimension ( 8 ) :: p1 = (/ &
    4.945235359296727046734888D+00, 2.018112620856775083915565D+02, &
    2.290838373831346393026739D+03, 1.131967205903380828685045D+04, &
    2.855724635671635335736389D+04, 3.848496228443793359990269D+04, &
    2.637748787624195437963534D+04, 7.225813979700288197698961D+03 /),&
   p2 = (/4.974607845568932035012064D+00, 5.424138599891070494101986D+02, &
    1.550693864978364947665077D+04, 1.847932904445632425417223D+05, &
    1.088204769468828767498470D+06, 3.338152967987029735917223D+06, &
    5.106661678927352456275255D+06, 3.074109054850539556250927D+06 /),&
   p4 = (/1.474502166059939948905062D+04, 2.426813369486704502836312D+06, &
    1.214755574045093227939592D+08, 2.663432449630976949898078D+09, &
    2.940378956634553899906876D+10, 1.702665737765398868392998D+11, &
    4.926125793377430887588120D+11, 5.606251856223951465078242D+11 /)
  real    ( kind = 8 ), parameter :: pnt68 = 0.6796875D+00
  real    ( kind = 8 ), parameter, dimension ( 8 ) :: q1 = (/ &
    6.748212550303777196073036D+01, 1.113332393857199323513008D+03, &
    7.738757056935398733233834D+03, 2.763987074403340708898585D+04, &
    5.499310206226157329794414D+04, 6.161122180066002127833352D+04, &
    3.635127591501940507276287D+04, 8.785536302431013170870835D+03 /)
  real    ( kind = 8 ), parameter, dimension ( 8 ) :: q2 = (/ &
    1.830328399370592604055942D+02, 7.765049321445005871323047D+03, &
    1.331903827966074194402448D+05, 1.136705821321969608938755D+06, &
    5.267964117437946917577538D+06, 1.346701454311101692290052D+07, &
    1.782736530353274213975932D+07, 9.533095591844353613395747D+06 /)
  real    ( kind = 8 ), parameter, dimension ( 8 ) :: q4 = (/ &
    2.690530175870899333379843D+03, 6.393885654300092398984238D+05, &
    4.135599930241388052042842D+07, 1.120872109616147941376570D+09, &
    1.488613728678813811542398D+10, 1.016803586272438228077304D+11, &
    3.417476345507377132798597D+11, 4.463158187419713286462081D+11 /)
  real    ( kind = 8 ), parameter :: sqrtpi = 0.9189385332046727417803297D+00
  real    ( kind = 8 ), parameter :: xbig = 4.08D+36
  integer ( kind = 4 ) :: i
  real    ( kind = 8 ) :: x,res,xden,xm1,xm2,xm4,xnum,xsq,corr,gamma_log

!  Return immediately if the argument is out of range.
  if ( x <= 0.0D+00 .or. xbig < x ) then
    gamma_log = huge ( gamma_log )
    return
  end if

  if ( x <= epsilon ( x ) ) then
    res = -log ( x )
  else if ( x <= 1.5D+00 ) then
    if ( x < pnt68 ) then
      corr = - log ( x ) ; xm1 = x
    else
      corr = 0.0D+00 ; xm1 = ( x - 0.5D+00 ) - 0.5D+00
    end if
    if ( x <= 0.5D+00 .or. pnt68 <= x ) then
      xden = 1.0D+00 ;  xnum = 0.0D+00
      do i = 1, 8
        xnum = xnum * xm1 + p1(i) ; xden = xden * xm1 + q1(i)
      end do
      res = corr + ( xm1 * ( d1 + xm1 * ( xnum / xden ) ) )
    else
      xm2 = ( x - 0.5D+00 ) - 0.5D+00
      xden = 1.0D+00 ; xnum = 0.0D+00
      do i = 1, 8
        xnum = xnum * xm2 + p2(i)
        xden = xden * xm2 + q2(i)
      end do
      res = corr + xm2 * ( d2 + xm2 * ( xnum / xden ) )
    end if
  else if ( x <= 4.0D+00 ) then
    xm2 = x - 2.0D+00
    xden = 1.0D+00 ; xnum = 0.0D+00
    do i = 1, 8
      xnum = xnum * xm2 + p2(i)
      xden = xden * xm2 + q2(i)
    end do
    res = xm2 * ( d2 + xm2 * ( xnum / xden ) )
  else if ( x <= 12.0D+00 ) then
    xm4 = x - 4.0D+00
    xden = - 1.0D+00 ; xnum = 0.0D+00
    do i = 1, 8
      xnum = xnum * xm4 + p4(i)
      xden = xden * xm4 + q4(i)
    end do
    res = d4 + xm4 * ( xnum / xden )
  else
    res = 0.0D+00
    if ( x <= frtbig ) then
      res = c(7) ; xsq = x * x
      do i = 1, 6
        res = res / xsq + c(i)
      end do
    end if
    res = res / x ; corr = log ( x )
    res = res + sqrtpi - 0.5D+00 * corr
    res = res + x * ( corr - 1.0D+00 )
  end if
  gamma_log = res
  return
end function gamma_log

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!CDF de la gaussienne.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine normal_01_cdf ( x, cdf )
  implicit none

  real ( kind = 8 ), parameter :: a1 = 0.398942280444D+00,&
           a2 = 0.399903438504D+00,a3 = 5.75885480458D+00,&
           a4 = 29.8213557808D+00,a5 = 2.62433121679D+00,&
           a6 = 48.6959930692D+00,a7 = 5.92885724438D+00,&
           b0 = 0.398942280385D+00,b1 = 3.8052D-08,&
           b2 = 1.00000615302D+00,b3 = 3.98064794D-04,&
           b4 = 1.98615381364D+00,b5 = 0.151679116635D+00,&
           b6 = 5.29330324926D+00,b7 = 4.8385912808D+00,&
           b8 = 15.1508972451D+00,b9 = 0.742380924027D+00,&
           b10 = 30.789933034D+00,b11 = 3.99019417011D+00
  real ( kind = 8 ) :: cdf,q,x,y

!  |X| <= 1.28.
  if ( abs ( x ) <= 1.28D+00 ) then
    y = 0.5D+00 * x * x
    q = 0.5D+00 - abs ( x ) * ( a1 - a2 * y / ( y + a3 - a4 / ( y + a5 &
      + a6 / ( y + a7 ) ) ) )
!  1.28 < |X| <= 12.7
  else if ( abs ( x ) <= 12.7D+00 ) then
    y = 0.5D+00 * x * x
    q = exp ( - y ) * b0 / ( abs ( x ) - b1 &
      + b2 / ( abs ( x ) + b3 + b4 / ( abs ( x ) - b5 &
      + b6 / ( abs ( x ) + b7 - b8 / ( abs ( x ) + b9 &
      + b10 / ( abs ( x ) + b11 ) ) ) ) ) )
!  12.7 < |X|
  else
    q = 0.0D+00
  end if
!  Take account of negative X.
  if ( x < 0.0D+00 ) then
    cdf = q
  else
    cdf = 1.0D+00 - q
  end if

  return
end subroutine normal_01_cdf

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!CDF de la gamma.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine gamma_cdf ( x, a, b, c, cdf )
 implicit none
  real ( kind = 8 ) :: a,b,c,cdf,p2,x,x2
 ! real    ( kind = 8 ) gamma_inc
  x2 = ( x - a ) / b ; p2 = c
  cdf = gamma_inc ( p2, x2 )
  return
end subroutine gamma_cdf

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!CDF de la Chi2.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine chi_square_cdf ( x, a, cdf )
  implicit none
  real ( kind = 8 ) :: a,a2,b2,c2,cdf,x,x2
  x2 = 0.5D+00 * x ; a2 = 0.0D+00
  b2 = 1.0D+00 ; c2 = 0.5D+00 * a
  call gamma_cdf ( x2, a2, b2, c2, cdf )
  return
end subroutine chi_square_cdf

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!Gamma Incomplete function
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

function gamma_inc ( p, x )
  implicit none

  real ( kind = 8 ) :: a,arg,b,c,cdf,gamma_inc,&
               p,pn1,pn2,pn3,pn4,pn5,pn6,rn,x
  !real    ( kind = 8 ) gamma_log
  real ( kind = 8 ), parameter :: exp_arg_min = -88.0D+00,&
                  overflow = 1.0D+37,plimit = 1000.0D+00,&
                  tol = 1.0D-07, xbig = 1.0D+08

  gamma_inc = 0.0D+00
  if ( p <= 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'GAMMA_INC - Fatal error!'
    write ( *, '(a)' ) '  Parameter P <= 0.'
    stop
  end if

  if ( x <= 0.0D+00 ) then
    gamma_inc = 0.0D+00
    return
  end if
!  Use a normal approximation if PLIMIT < P.
  if ( plimit < p ) then
    pn1 = 3.0D+00 * sqrt ( p ) * ( ( x / p ) ** ( 1.0D+00 / 3.0D+00 ) &
      + 1.0D+00 / ( 9.0D+00 * p ) - 1.0D+00 )
    call normal_01_cdf ( pn1, cdf )
    gamma_inc = cdf
    return
  end if
!  Is X extremely large compared to P?
  if ( xbig < x ) then
    gamma_inc = 1.0D+00
    return
  end if
!  Use Pearson's series expansion.
!  (P is not large enough to force overflow in the log of Gamma.
  if ( x <= 1.0D+00 .or. x < p ) then
    arg = p * log ( x ) - x - gamma_log ( p + 1.0D+00 )
    c = 1.0D+00
    gamma_inc = 1.0D+00
    a = p
    do
      a = a + 1.0D+00 ; c = c * x / a
      gamma_inc = gamma_inc + c
      if ( c <= tol ) then
        exit
      end if
    end do
    arg = arg + log ( gamma_inc )
    if ( exp_arg_min <= arg ) then
      gamma_inc = exp ( arg )
    else
      gamma_inc = 0.0D+00
    end if
  else
!  Use a continued fraction expansion.
    arg = p * log ( x ) - x - gamma_log ( p )
    a = 1.0D+00 - p ; b = a + x + 1.0D+00
    c = 0.0D+00 ; pn1 = 1.0D+00
    pn2 = x ; pn3 = x + 1.0D+00
    pn4 = x * b ; gamma_inc = pn3 / pn4
    do
      a = a + 1.0D+00 ; b = b + 2.0D+00
      c = c + 1.0D+00
      pn5 = b * pn3 - a * c * pn1
      pn6 = b * pn4 - a * c * pn2
      if ( 0.0D+00 < abs ( pn6 ) ) then
        rn = pn5 / pn6
        if ( abs ( gamma_inc - rn ) <= min ( tol, tol * rn ) ) then
          arg = arg + log ( gamma_inc )
          if ( exp_arg_min <= arg ) then
            gamma_inc = 1.0D+00 - exp ( arg )
          else
            gamma_inc = 1.0D+00
          end if
          return
        end if
        gamma_inc = rn
      end if
      pn1 = pn3 ; pn2 = pn4
      pn3 = pn5 ; pn4 = pn6
!  Rescale terms in continued fraction if terms are large.
      if ( overflow <= abs ( pn5 ) ) then
        pn1 = pn1 / overflow ; pn2 = pn2 / overflow
        pn3 = pn3 / overflow ; pn4 = pn4 / overflow
      end if
    end do
  end if
  return
end function gamma_inc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!Probabilite density function: BINOMIALE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

function binomial_pdf (x,a,b) !x=k;a=n,b=p
  implicit none
  real (kind=8), intent(in) :: b
  integer, intent(in) :: a,x
  real (kind=8) :: cnk , binomial_pdf

  if ( a < 1 ) then
    binomial_pdf = 0.0D+00
  else if ( x < 0 .or. a < x ) then
    binomial_pdf = 0.0D+00
  else if ( b == 0.0D+00 ) then
    if ( x == 0 ) then
    binomial_pdf = 1.0D+00
    else
    binomial_pdf = 0.0D+00
    end if
  else if ( b == 1.0D+00 ) then
    if ( x == a ) then
    binomial_pdf = 1.0D+00
    else
    binomial_pdf = 0.0D+00
    end if
  else
    
    cnk= binomial_coef(a , x)
    binomial_pdf = cnk * b**x * ( 1.0D+00 - b )**( a - x )

  end if

end function binomial_pdf

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!Probabilite density function: BINOMIALE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

function beta_pdf ( x, a, b)
  implicit none
  real (kind=8), intent(in) :: x,a,b
  real (kind=8) :: beta_pdf
  if ( x < 0.0D+00 .or. 1.0D+00 < x ) then
    beta_pdf = 0.0D+00
  else
    beta_pdf = x**( a - 1.0D+00 ) * ( 1.0D+00 - x )**( b - 1.0D+00 ) / beta ( a, b )
  end if
  return
end function

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!random beta
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

FUNCTION random_beta(aa, bb,this) RESULT(fn_val)
implicit none
type(mt_state), intent(inout) :: this
REAL (kind=8), INTENT(IN)    :: aa, bb
REAL (kind=8), PARAMETER  :: aln4 = 1.3862944
REAL (kind=8)    :: fn_val,a, b, g, r, s, x, y, z , vsmall = TINY(1.0), vlarge = HUGE(1.0)
REAL (kind=8), SAVE       :: d, f, h, t, c
LOGICAL, SAVE    :: swap

IF (aa <= 0. .OR. bb <= 0.) THEN
  WRITE(*, *) 'IMPERMISSIBLE SHAPE PARAMETER VALUE(S)'
  STOP
END IF

  a = aa ; b = bb
  swap = b > a
  IF (swap) THEN
    g = b ; b = a ; a = g
  END IF
  d = a/b ; f = a+b
  IF (b > 1.) THEN
    h = SQRT((2.*a*b - f)/(f - 2.))
    t = 1.
  ELSE
    h = b
    t = 1./(1. + (a/(vlarge*b))**b)
  END IF
  c = a+h

DO
  r=grnd(this) ; x=grnd(this) 
  s = r*r*x
  IF (r < vsmall .OR. s <= 0.) CYCLE
  IF (r < t) THEN
    x = LOG(r/(1. - r))/h
    y = d*EXP(x)
    z = c*x + f*LOG((1. + d)/(1. + y)) - aln4
    IF (s - 1. > z) THEN
      IF (s - s*z > 1.) CYCLE
      IF (LOG(s) > z) CYCLE
    END IF
    fn_val = y/(1. + y)
  ELSE
    IF (4.0*s > (1. + 1./d)**f) CYCLE
    fn_val = 1.
  END IF
  EXIT
END DO

IF (swap) fn_val = 1. - fn_val
RETURN
END FUNCTION random_beta

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!Calcul de la fonction digamma
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

function digamma ( x )
  implicit none
  real ( kind = 8 ), parameter :: c = 8.5D+00, d1 = -0.5772156649D+00,s = 0.00001D+00,&
                                  s3 = 0.08333333333D+00, s4 = 0.0083333333333D+00,s5 = 0.003968253968D+00

  real ( kind = 8 ) :: digamma,r,x,y

  if ( x <= 0.0D+00 ) then
    digamma = 0.0D+00
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'DIGAMMA - Fatal error!'
    write ( *, '(a)' ) '  X <= 0.'
    stop
  else if ( x <= s ) then
    digamma = d1 - 1.0D+00 / x
  else
    digamma = 0.0D+00
    y = x
    do while ( y < c )
      digamma = digamma - 1.0D+00 / y
      y = y + 1.0D+00
    end do
    r = 1.0D+00 / ( y * y )
    digamma = digamma + log ( y ) - 0.5D+00 / y &
      - r * ( s3 - r * ( s4 - r * s5 ) )
  end if
  return
end function digamma

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!Calcul temps ecoule entre 2 cpu
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

function elapse_time(t1,t2,cl_rate)
implicit none
integer, intent(in) :: t1,t2,cl_rate
integer :: elapse_time(4),tmp

tmp=(t2-t1)/cl_rate
elapse_time(1)=floor(tmp/86400.) ; tmp = tmp - elapse_time(1)*86400
elapse_time(2)=floor(tmp/3600.) ; tmp = tmp - elapse_time(2)*3600
elapse_time(3)=floor(tmp/60.) ; elapse_time(4) = tmp - elapse_time(3)*60
return

end function elapse_time

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!resolution de systeme lineraire
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine solve(a,rhs)
implicit none
! solves a*sol=rhs   (rhs = matrix)  
! sets rhs = sol; on output A = LU decomposition of A
! to get inverse set rhs = identity matrix
! optionally the log determinant is returned in det
  real (kind=8),intent(inout),dimension(:,:) :: a,rhs
  integer, save :: indx(8000)
  real (kind=8) :: d
  integer :: n,m,i

  n=size(a,1) ; m=size(rhs,2)
  call ludcmp(a,n,indx,d)
  do i=1,m
    call lubksb(a,n,indx,rhs(:,i))
  end do
end subroutine solve  

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Inversion de matrice
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

function inv(a)
implicit none
real (kind=8) , intent (in), dimension(:,:) :: a
real (kind=8), dimension(:,:) :: inv(size(a,1),size(a,1)),b(size(a,1),size(a,1)),aa(size(a,1),size(a,1))
integer :: n,i

n=size(a,1) ; b=0
do i=1,n
  b(i,i)=1
end do
aa=a
call solve(aa(1:n,1:n),b(1:n,1:n))
inv=b
end function inv

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Decomposition LU
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine ludcmp(a,n,indx,d)
!decomposition LU: a=matrice a decomposer, d=determinant
implicit none
real (kind=8),intent (inout), dimension(:,:) :: a
integer, intent (inout), dimension(:) :: indx
integer, intent (in) :: n

real (kind=8),intent (out) :: d
real (kind=8),dimension(:) :: vv(size(a,1))
integer :: i,j,k,imax!,nmax=8000
real (kind=8)  :: somme,dum,aamax,eps_t=1.0E-20

d=1.
do i=1,n
  aamax=0.
  do j=1,n
   if (abs(a(i,j))>aamax) aamax=abs(a(i,j))
  end do
  if (aamax==0.) stop 'Singular matrix.'
  vv(i)=1./aamax
end do

do j=1,n
 if (j>1) then
  do i=1,j-1
   somme=a(i,j)
   if (i>1) then
    do k=1,i-1
     somme=somme-A(I,K)*A(K,J)
    end do
   a(i,j)=somme
   end if
  end do
 end if
 aamax=0.
 do I=J,N
  somme=a(i,j)
  if (J>1) then
   do K=1,J-1
    somme=somme-A(I,K)*A(K,J)
   end do
   a(i,j)=somme
  end if
  dum=vv(i)*abs(somme)
  if (dum>=aamax) then
   imax=i ; aamax=dum
  end if
 end do
 if (J/=imax) then
  do k=1,n
   dum=a(imax,k) ; a(imax,k)=a(j,k) ; a(j,k)=dum
  end do
  d=-d ; vv(imax)=vv(j)
 end if
 indx(j)=imax
 if(j/=n) then
  if(a(j,j)==0.) a(j,j)=eps_t
  dum=1./a(j,j)
  do i=j+1,n
   a(i,j)=a(i,j)*dum
  end do
 end if
end do
if(a(n,n)==0.) a(n,n)=eps_t
return
end subroutine ludcmp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Resolution de Ax = b. (utilise avec ludcmp.)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine lubksb(a,n,indx,b)
implicit none
real (kind=8),intent(in),dimension(:,:) :: A
real (kind=8),intent (out), dimension(:) :: B
integer, intent(in), dimension (:) :: indx
integer, intent(in) :: n
integer :: ii,i,j,ll
real (kind=8) :: somme

ii=0
do i=1,n
 ll=indx(i) ; somme=b(ll) ; b(ll)=b(i)
 if (ii/=0) then
  do j=ii,i-1
   somme=somme-a(i,j)*b(j)
  end do
 else if (somme/=0.) then
  ii=i
 end if
 b(i)=somme
end do

do i=n,1,-1
 somme=b(i)
 if(i<n) then
  do j=i+1,n
   somme=somme-A(I,J)*B(J)
  end do
 end if
 b(i)=somme/a(i,i)
end do
return
end subroutine lubksb

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Decomposition Choleski: A = C*C'
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine chol(a)
implicit none
real (kind=8),intent (inout), dimension (:,:) :: a
real (kind=8) :: sum1,sum2
integer :: i,j,k,m,n

!a_sauv=a
n=size(a,1)
do k=1,n 
 sum1 = a(k,k) 
 do m=1,k-1 
  sum1 = sum1 - a(k,m)**2.0 
 end do
 a(k,k) = sqrt(sum1) 
 do i=k+1,n 
  sum2 = a(i,k) 
  do m=1,k-1 
   sum2 = sum2 - a(i,m)*a(k,m) 
  end do 
  a(i,k) = sum2/a(k,k) 
 end do
end do
 
do i=1,n
 do j=i+1,n
  a(i,j)=0.
 end do
end do

return
end subroutine chol


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Tirage dans une gaussienne multivariée
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

function rand_mvnorm(vect_mu,prec_mat,this)
!vect_mu=vecteur moyenne , prec_mat=matrice precision
implicit none
type(mt_state), intent(inout) :: this
real (kind=8), intent(in), dimension (:,:):: prec_mat
real (kind=8), intent(in), dimension (:):: vect_mu
real (kind=8) :: varcov_mat(size(prec_mat,1),size(prec_mat,1)),vect_z(size(prec_mat,1)),rand_mvnorm(size(prec_mat,1))
integer :: i,n

if(size(prec_mat,1)/=size(vect_mu)) stop 'rand_mvnorm: pb dimension mu et prec_mat'

n=size(prec_mat,1)
varcov_mat=inv(prec_mat)
call chol(varcov_mat) !decomposition de choleski de la matrice de variance covariance
do i=1,n
 vect_z(i)=random_normal(this)
end do

rand_mvnorm=vect_mu + matmul(varcov_mat,vect_z)

end function rand_mvnorm

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Calcul de la matrice carrée XX' à partir d'un vecteur X
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

function txx(vect)
  implicit none
  real (kind=8),dimension(:) :: vect
  real (kind=8),dimension(:,:) :: txx(size(vect),size(vect))
  integer :: tmp1,tmp2,dim_vect
 
  dim_vect=size(vect)
  do tmp1=1,dim_vect
   do tmp2=1,dim_vect
    txx(tmp1,tmp2)=vect(tmp1)*vect(tmp2)
   end do
  end do
end function txx

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!FONCTION random_chisq (avec appel au mersenne twister
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

FUNCTION random_chisq(df,this)

implicit none
type(mt_state), intent(inout) :: this
real (kind=8) :: random_chisq,df, arg

!  arg1 = 1.0E+00
  arg = df / 2.0E+00
  if(arg<1.0) then
   random_chisq = 2.0E+00 * random_gamma2(arg,this)
  else
   random_chisq = 2.0E+00 * random_gamma1(arg,this)
  end if

END FUNCTION random_chisq


subroutine bartlett_sample ( m, df, sigma, t,this )

!*****************************************************************************80
!
!! BARTLETT_SAMPLE samples the Bartlett distribution.
!
!  Discussion:
!
!    If the matrix T is sampled from the Bartlett distribution, then 
!    the matrix W = T' * T is a sample from the Wishart distribution.
!
!    This function requires functions from the PDFLIB and RNGLIB libraries.
!
!    The "initialize()" function from RNGLIB must be called before using
!    this function.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 August 2013
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Patrick Odell, Alan Feiveson,
!    A numerical procedure to generate a sample covariance matrix,
!    Journal of the American Statistical Association,
!    Volume 61, Number 313, March 1966, pages 199-203.
!
!    Stanley Sawyer,
!    Wishart Distributions and Inverse-Wishart Sampling,
!    Washington University,
!    30 April 2007, 12 pages.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the order of the matrix.
!
!    Input, integer ( kind = 4 ) DF, the number of degrees of freedom.
!    M <= DF.
!
!    Input, real ( kind = 8 ) SIGMA(M,M), the covariance matrix, which should be 
!    a symmetric positive definite matrix.
!
!    Output, real ( kind = 8 ) T(M,M), the sample matrix from 
!    the Bartlett distribution.
!
  implicit none
  type(mt_state), intent(inout) :: this
  integer ( kind = 4 ) m

  integer ( kind = 4 ) df
  integer ( kind = 4 ) flag;
  real ( kind = 8 ), allocatable :: r(:,:)
  real ( kind = 8 ) sigma(m,m)
  real ( kind = 8 ) t(m,m)
  real ( kind = 8 ), allocatable :: tu(:,:)

  if ( df < m ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'BARTLETT_SAMPLE - Fatal error!'
    write ( *, '(a)' ) '  DF = ', df, ' < M = ', m
    stop 1
  end if
!
!  Get the upper triangular Cholesky factor of SIGMA.
!
  allocate ( r(1:m,1:m) )
  call r8mat_cholesky_factor_upper ( m, sigma, r, flag )

  if ( flag /= 0 ) then
    write ( *, '(a)' ) ''
    write ( *, '(a)' ) 'BARTLETT_SAMPLE - Fatal error!'
    write ( *, '(a)' ) &
      '  Unexpected error return from R8MAT_CHOLESKY_FACTOR_UPPER.'
    write ( *, '(a,i4)' ) '  FLAG = ', flag
    stop 1
  end if
!
!  Sample the unit Bartlett distribution.
!
  allocate ( tu(1:m,1:m) )
  call bartlett_unit_sample ( m, df, tu, this )
!
!  Construct the matrix T = TU * R.
!
  t = matmul ( tu(1:m,1:m), r(1:m,1:m) )
!
!  Free memory.
!
  deallocate ( r )
  deallocate ( tu )

  return
end subroutine bartlett_sample

subroutine bartlett_unit_sample ( m, df, t, this )

!*****************************************************************************80
!
!! BARTLETT_UNIT_SAMPLE samples the unit Bartlett distribution.
!
!  Discussion:
!
!    If the matrix T is sampled from the unit Bartlett distribution, then 
!    the matrix W = T' * T is a sample from the unit Wishart distribution.
! 
!    This function requires functions from the PDFLIB and RNGLIB libraries.
!
!    The "initialize()" function from RNGLIB must be called before using
!    this function.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 October 2013
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Patrick Odell, Alan Feiveson,
!    A numerical procedure to generate a sample covariance matrix,
!    Journal of the American Statistical Association,
!    Volume 61, Number 313, March 1966, pages 199-203.
!
!    Stanley Sawyer,
!    Wishart Distributions and Inverse-Wishart Sampling,
!    Washington University,
!    30 April 2007, 12 pages.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the order of the matrix.
!
!    Input, integer ( kind = 4 ) DF, the number of degrees of freedom.
!    M <= DF.
!
!    Output, real ( kind = 8 ) T(M,M), the sample matrix from the 
!    unit Bartlett distribution.
!
  implicit none
  type(mt_state), intent(inout) :: this
  integer ( kind = 4 ) m

  integer ( kind = 4 ) df
  real ( kind = 8 ) df_chi
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
!  real ( kind = 8 ) r8_chi_sample
!  real ( kind = 8 ) r8_normal_01_sample
  real ( kind = 8 ) t(m,m)

  if ( df < m ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'BARTLETT_UNIT_SAMPLE - Fatal error!'
    write ( *, '(a)' ) '  DF = ', df, ' < M = ', m
    stop 1
  end if

  do i = 1, m

    t(i,1:i-1) = 0.0D+00

    df_chi = real ( df + 1 - i, kind = 8 )
    t(i,i) = sqrt ( random_chisq(df_chi,this) ) !sqrt ( r8_chi_sample ( df_chi ) )

    do j = i + 1, m
      t(i,j) = random_normal(this) !r8_normal_01_sample ( )
    end do

  end do

  return
end subroutine bartlett_unit_sample 

subroutine jacobi_eigenvalue ( n, a, it_max, v, d, it_num, rot_num )

!*****************************************************************************80
!
!! JACOBI_EIGENVALUE carries out the Jacobi eigenvalue iteration.
!
!  Discussion:
!
!    This function computes the eigenvalues and eigenvectors of a
!    real symmetric matrix, using Rutishauser's modfications of the classical
!    Jacobi rotation method with threshold pivoting. 
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 July 2013
!
!  Author:
!
!    FORTRAN90 version by John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input, real ( kind = 8 ) A(N,N), the matrix, which must be square, real,
!    and symmetric.
!
!    Input, integer ( kind = 4 ) IT_MAX, the maximum number of iterations.
!
!    Output, real ( kind = 8 ) V(N,N), the matrix of eigenvectors.
!
!    Output, real ( kind = 8 ) D(N), the eigenvalues, in descending order.
!
!    Output, integer ( kind = 4 ) IT_NUM, the total number of iterations.
!
!    Output, integer ( kind = 4 ) ROT_NUM, the total number of rotations.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n,n)
  real ( kind = 8 ) bw(n)
  real ( kind = 8 ) c
  real ( kind = 8 ) d(n)
  real ( kind = 8 ) g
  real ( kind = 8 ) gapq
  real ( kind = 8 ) h
  integer ( kind = 4 ) i
  integer ( kind = 4 ) it_max
  integer ( kind = 4 ) it_num
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  integer ( kind = 4 ) m
  integer ( kind = 4 ) p
  integer ( kind = 4 ) q
  integer ( kind = 4 ) rot_num
  real ( kind = 8 ) s
  real ( kind = 8 ) t
  real ( kind = 8 ) tau
  real ( kind = 8 ) term
  real ( kind = 8 ) termp
  real ( kind = 8 ) termq
  real ( kind = 8 ) theta
  real ( kind = 8 ) thresh
  real ( kind = 8 ) v(n,n)
  real ( kind = 8 ) w(n)
  real ( kind = 8 ) zw(n)

  call r8mat_identity ( n, v )

  call r8mat_diag_get_vector ( n, a, d )

  bw(1:n) = d(1:n)
  zw(1:n) = 0.0D+00
  it_num = 0
  rot_num = 0

  do while ( it_num < it_max )

    it_num = it_num + 1
!
!  The convergence threshold is based on the size of the elements in
!  the strict upper triangle of the matrix.
!
    thresh = 0.0D+00
    do j = 1, n
      do i = 1, j - 1
        thresh = thresh + a(i,j) ** 2
      end do
    end do

    thresh = sqrt ( thresh ) / real ( 4 * n, kind = 8 )

    if ( thresh == 0.0D+00 ) then
      exit 
    end if

    do p = 1, n
      do q = p + 1, n

        gapq = 10.0D+00 * abs ( a(p,q) )
        termp = gapq + abs ( d(p) )
        termq = gapq + abs ( d(q) )
!
!  Annihilate tiny offdiagonal elements.
!
        if ( 4 < it_num .and. &
             termp == abs ( d(p) ) .and. &
             termq == abs ( d(q) ) ) then

          a(p,q) = 0.0D+00
!
!  Otherwise, apply a rotation.
!
        else if ( thresh <= abs ( a(p,q) ) ) then

          h = d(q) - d(p)
          term = abs ( h ) + gapq

          if ( term == abs ( h ) ) then
            t = a(p,q) / h
          else
            theta = 0.5D+00 * h / a(p,q)
            t = 1.0D+00 / ( abs ( theta ) + sqrt ( 1.0D+00 + theta * theta ) )
            if ( theta < 0.0D+00 ) then 
              t = - t
            end if
          end if

          c = 1.0D+00 / sqrt ( 1.0D+00 + t * t )
          s = t * c
          tau = s / ( 1.0D+00 + c )
          h = t * a(p,q)
!
!  Accumulate corrections to diagonal elements.
!
          zw(p) = zw(p) - h                  
          zw(q) = zw(q) + h
          d(p) = d(p) - h
          d(q) = d(q) + h

          a(p,q) = 0.0D+00
!
!  Rotate, using information from the upper triangle of A only.
!
          do j = 1, p - 1
            g = a(j,p)
            h = a(j,q)
            a(j,p) = g - s * ( h + g * tau )
            a(j,q) = h + s * ( g - h * tau )
          end do

          do j = p + 1, q - 1
            g = a(p,j)
            h = a(j,q)
            a(p,j) = g - s * ( h + g * tau )
            a(j,q) = h + s * ( g - h * tau )
          end do

          do j = q + 1, n
            g = a(p,j)
            h = a(q,j)
            a(p,j) = g - s * ( h + g * tau )
            a(q,j) = h + s * ( g - h * tau )
          end do
!
!  Accumulate information in the eigenvector matrix.
!
          do j = 1, n
            g = v(j,p)
            h = v(j,q)
            v(j,p) = g - s * ( h + g * tau )
            v(j,q) = h + s * ( g - h * tau )
          end do

          rot_num = rot_num + 1

        end if

      end do
    end do

    bw(1:n) = bw(1:n) + zw(1:n)
    d(1:n) = bw(1:n)
    zw(1:n) = 0.0D+00

  end do
!
!  Restore upper triangle of input matrix.
!
  do j = 1, n
    do i = 1, j - 1
      a(i,j) = a(j,i)
    end do
  end do
!
!  Descending sort the eigenvalues and eigenvectors.
!
  do k = 1, n - 1

    m = k

    do l = k + 1, n
      if ( d(m) < d(l) ) then
        m = l
      end if
    end do

    if ( m /= k ) then

      t    = d(m)
      d(m) = d(k)
      d(k) = t

      w(1:n)   = v(1:n,m)
      v(1:n,m) = v(1:n,k)
      v(1:n,k) = w(1:n)

    end if

  end do

  return
end subroutine jacobi_eigenvalue

function r8_epsilon ( )

!*****************************************************************************80
!
!! R8_EPSILON returns the R8 roundoff unit.
!
!  Discussion:
!
!    The roundoff unit is a number R which is a power of 2 with the
!    property that, to the precision of the computer's arithmetic,
!      1 < 1 + R
!    but
!      1 = ( 1 + R / 2 )
!
!    FORTRAN90 provides the superior library routine
!
!      EPSILON ( X )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 September 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) R8_EPSILON, the round-off unit.
!
  implicit none

  real ( kind = 8 ) r8_epsilon

  r8_epsilon = 2.220446049250313D-016

  return
end function r8_epsilon

subroutine r8mat_cholesky_factor_upper ( n, a, c, flag )

!*****************************************************************************80
!
!! R8MAT_CHOLESKY_FACTOR_UPPER: upper Cholesky factor of a symmetric matrix.
!
!  Discussion:
!
!    The matrix must be symmetric and positive semidefinite.
!
!    For a positive semidefinite symmetric matrix A, the Cholesky factorization
!    is an upper triangular matrix R such that:
!
!      A = R * R'
!
!    The lower Cholesky factor is a lower triangular matrix L such that
!
!      A = L * L'
!
!    An R8MAT is an MxN array of R8's, stored by (I,J) -> [I+J*M].
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 August 2013
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of rows and columns of
!    the matrix A.
!
!    Input, real ( kind = 8 ) A(N,N), the N by N matrix.
!
!    Output, real ( kind = 8 ) C(N,N), the N by N upper triangular
!    Cholesky factor.
!
!    Output, integer ( kind = 4 ) FLAG:
!    0, no error occurred.
!    1, the matrix is not positive definite.
!    2, the matrix is not nonnegative definite.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n,n)
  real ( kind = 8 ) c(n,n)
  integer ( kind = 4 ) flag
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
 ! real ( kind = 8 ) r8_epsilon
  real ( kind = 8 ) sum2
 ! real ( kind = 8 ) tol

  flag = 0;

  c(1:n,1:n) = a(1:n,1:n);

  do j = 1, n

    c(j,1:j-1) = 0.0D+00

    do i = j, n

      sum2 = c(i,j) - dot_product ( c(1:j-1,j), c(1:j-1,i) )

      if ( i == j ) then
        if ( sum2 <= 0.0D+00 ) then
          flag = 1
          return
        else
          c(j,i) = sqrt ( sum2 )
        end if
      else
        if ( c(j,j) /= 0.0D+00 ) then
          c(j,i) = sum2 / c(j,j)
        else
          c(j,i) = 0.0D+00
        end if
      end if

    end do

  end do

  return
end subroutine r8mat_cholesky_factor_upper

subroutine r8mat_diag_get_vector ( n, a, v )

!*****************************************************************************80
!
!! R8MAT_DIAG_GET_VECTOR gets the value of the diagonal of an R8MAT.
!
!  Discussion:
!
!    An R8MAT is an MxN array of R8's, stored by (I,J) -> [I+J*M].
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 March 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of rows and columns of
!    the matrix.
!
!    Input, real ( kind = 8 ) A(N,N), the N by N matrix.
!
!    Output, real ( kind = 8 ) V(N), the diagonal entries
!    of the matrix.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n,n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) v(n)

  do i = 1, n
    v(i) = a(i,i)
  end do

  return
end subroutine r8mat_diag_get_vector

subroutine r8mat_diagonal ( n, diag, a )

!*****************************************************************************80
!
!! R8MAT_DIAGONAL returns a diagonal matrix as an R8MAT.
!
!  Discussion:
!
!    An R8MAT is an MxN array of R8's, stored by (I,J) -> [I+J*M].
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 July 2013
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of A.
!
!    Input, real ( kind = 8 ) DIAG(N), the diagonal entries.
!
!    Output, real ( kind = 8 ) A(N,N), the N by N diagonal matrix.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n,n)
  real ( kind = 8 ) diag(n)
  integer ( kind = 4 ) i

  a(1:n,1:n) = 0.0D+00

  do i = 1, n
    a(i,i) = diag(i)
  end do

  return
end subroutine r8mat_diagonal

subroutine r8mat_identity ( n, a )

!*****************************************************************************80
!
!! R8MAT_IDENTITY stores the identity matrix in an R8MAT.
!
!  Discussion:
!
!    An R8MAT is an MxN array of R8's, stored by (I,J) -> [I+J*M].
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of A.
!
!    Output, real ( kind = 8 ) A(N,N), the N by N identity matrix.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n,n)
  integer ( kind = 4 ) i

  a(1:n,1:n) = 0.0D+00

  do i = 1, n
    a(i,i) = 1.0D+00
  end do

  return
end subroutine r8mat_identity

function r8mat_norm_fro_affine ( m, n, a1, a2 )

!*****************************************************************************80
!
!! R8MAT_NORM_FRO_AFFINE returns the Frobenius norm of an R8MAT difference.
!
!  Discussion:
!
!    An R8MAT is an MxN array of R8's, stored by (I,J) -> [I+J*M].
!
!    The Frobenius norm is defined as
!
!      R8MAT_NORM_FRO = sqrt (
!        sum ( 1 <= I <= M ) sum ( 1 <= j <= N ) A(I,J) * A(I,J) )
!
!    The matrix Frobenius norm is not derived from a vector norm, but
!    is compatible with the vector L2 norm, so that:
!
!      r8vec_norm_l2 ( A * x ) <= r8mat_norm_fro ( A ) * r8vec_norm_l2 ( x ).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows.
!
!    Input, integer ( kind = 4 ) N, the number of columns.
!
!    Input, real ( kind = 8 ) A1(M,N), A2(M,N), the matrices for whose 
!    difference the Frobenius norm is desired.
!
!    Output, real ( kind = 8 ) R8MAT_NORM_FRO_AFFINE, the Frobenius 
!    norm of A1 - A2.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a1(m,n)
  real ( kind = 8 ) a2(m,n)
  real ( kind = 8 ) r8mat_norm_fro_affine

  r8mat_norm_fro_affine = sqrt ( sum ( ( a1(1:m,1:n) - a2(1:m,1:n) )**2 ) )

  return
end function r8mat_norm_fro_affine

subroutine r8mat_print ( m, n, a, title )

!*****************************************************************************80
!
!! R8MAT_PRINT prints an R8MAT.
!
!  Discussion:
!
!    An R8MAT is an MxN array of R8's, stored by (I,J) -> [I+J*M].
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows in A.
!
!    Input, integer ( kind = 4 ) N, the number of columns in A.
!
!    Input, real ( kind = 8 ) A(M,N), the matrix.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  character ( len = * ) title

  call r8mat_print_some ( m, n, a, 1, 1, m, n, title )

  return
end subroutine r8mat_print 

subroutine r8mat_print_some ( m, n, a, ilo, jlo, ihi, jhi, title )

!*****************************************************************************80
!
!! R8MAT_PRINT_SOME prints some of an R8MAT.
!
!  Discussion:
!
!    An R8MAT is an MxN array of R8's, stored by (I,J) -> [I+J*M].
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 September 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns.
!
!    Input, real ( kind = 8 ) A(M,N), an M by N matrix to be printed.
!
!    Input, integer ( kind = 4 ) ILO, JLO, the first row and column to print.
!
!    Input, integer ( kind = 4 ) IHI, JHI, the last row and column to print.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ), parameter :: incx = 5
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  character ( len = 14 ) ctemp(incx)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i2hi
  integer ( kind = 4 ) i2lo
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) inc
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) j2hi
  integer ( kind = 4 ) j2lo
  integer ( kind = 4 ) jhi
  integer ( kind = 4 ) jlo
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )

  if ( m <= 0 .or. n <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  (None)'
    return
  end if

  do j2lo = max ( jlo, 1 ), min ( jhi, n ), incx

    j2hi = j2lo + incx - 1
    j2hi = min ( j2hi, n )
    j2hi = min ( j2hi, jhi )

    inc = j2hi + 1 - j2lo

    write ( *, '(a)' ) ' '

    do j = j2lo, j2hi
      j2 = j + 1 - j2lo
      write ( ctemp(j2), '(i8,6x)' ) j
    end do

    write ( *, '(''  Col   '',5a14)' ) ctemp(1:inc)
    write ( *, '(a)' ) '  Row'
    write ( *, '(a)' ) ' '

    i2lo = max ( ilo, 1 )
    i2hi = min ( ihi, m )

    do i = i2lo, i2hi

      do j2 = 1, inc

        j = j2lo - 1 + j2

        if ( a(i,j) == real ( int ( a(i,j) ), kind = 8 ) ) then
          write ( ctemp(j2), '(f8.0,6x)' ) a(i,j)
        else
          write ( ctemp(j2), '(g14.6)' ) a(i,j)
        end if

      end do

      write ( *, '(i5,a,5a14)' ) i, ':', ( ctemp(j), j = 1, inc )

    end do

  end do

  return
end subroutine r8mat_print_some 

subroutine r8ut_inverse ( n, a )

!*****************************************************************************80
!
!! R8UT_INVERSE computes the inverse of an R8UT matrix.
!
!  Discussion:
!
!    The R8UT storage format is used for an M by N upper triangular 
!    matrix.  The format stores all M*N entries, even those which are zero.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    04 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Albert Nijenhuis, Herbert Wilf,
!    Combinatorial Algorithms,
!    Academic Press, 1978, second edition,
!    ISBN 0-12-519260-6
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input/output, real ( kind = 8 ) A(N,N).
!    On input, the upper triangular matrix to be inverted.
!    On output, the inverse of the upper triangular matrix.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
!
!  Check.
!
  do i = 1, n
    if ( a(i,i) == 0.0D+00 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R8UT_INVERSE - Fatal error!'
      write ( *, '(a)' ) '  Zero diagonal element.'
      stop 1
    end if
  end do

  do j = n, 1, -1

    do i = n, 1, -1

      if ( j < i ) then

        a(i,j) = 0.0D+00

      else if ( i == j ) then

        a(i,j) = 1.0D+00 / a(i,j)

      else if ( i < j ) then

        a(i,j) = - sum ( a(i,i+1:j) * a(i+1:j,j) ) / a(i,i)

      end if

    end do
  end do

  return
end subroutine r8ut_inverse 

subroutine r8vec_print ( n, a, title )

!*****************************************************************************80
!
!! R8VEC_PRINT prints an R8VEC.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 August 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of components of the vector.
!
!    Input, real ( kind = 8 ) A(N), the vector to be printed.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  integer ( kind = 4 ) i
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
  write ( *, '(a)' ) ' '

  do i = 1, n
    write ( *, '(2x,i8,a,1x,g16.8)' ) i, ':', a(i)
  end do

  return
end subroutine r8vec_print

subroutine timestamp ( )

!*****************************************************************************80
!
!! TIMESTAMP prints the current YMDHMS date as a time stamp.
!
!  Example:
!
!    31 May 2001   9:45:54.872 AM
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 May 2013
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    None
!
  implicit none

  character ( len = 8 ) ampm
  integer ( kind = 4 ) d
  integer ( kind = 4 ) h
  integer ( kind = 4 ) m
  integer ( kind = 4 ) mm
  character ( len = 9 ), parameter, dimension(12) :: month = (/ &
    'January  ', 'February ', 'March    ', 'April    ', &
    'May      ', 'June     ', 'July     ', 'August   ', &
    'September', 'October  ', 'November ', 'December ' /)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) s
  integer ( kind = 4 ) values(8)
  integer ( kind = 4 ) y

  call date_and_time ( values = values )

  y = values(1)
  m = values(2)
  d = values(3)
  h = values(5)
  n = values(6)
  s = values(7)
  mm = values(8)

  if ( h < 12 ) then
    ampm = 'AM'
  else if ( h == 12 ) then
    if ( n == 0 .and. s == 0 ) then
      ampm = 'Noon'
    else
      ampm = 'PM'
    end if
  else
    h = h - 12
    if ( h < 12 ) then
      ampm = 'PM'
    else if ( h == 12 ) then
      if ( n == 0 .and. s == 0 ) then
        ampm = 'Midnight'
      else
        ampm = 'AM'
      end if
    end if
  end if

  write ( *, '(i2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
    d, trim ( month(m) ), y, h, ':', n, ':', s, '.', mm, trim ( ampm )

  return
end subroutine timestamp

subroutine wishart_sample ( m, df, sigma, a, this )

!*****************************************************************************80
!
!! WISHART_SAMPLE samples the Wishart distribution.
!
!  Discussion:
!
!    This function requires functions from the PDFLIB and RNGLIB libraries.
!
!    The "initialize()" function from RNGLIB must be called before using
!    this function.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 August 2013
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Patrick Odell, Alan Feiveson,
!    A numerical procedure to generate a sample covariance matrix,
!    Journal of the American Statistical Association,
!    Volume 61, Number 313, March 1966, pages 199-203.
!
!    Stanley Sawyer,
!    Wishart Distributions and Inverse-Wishart Sampling,
!    Washington University,
!    30 April 2007, 12 pages.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the order of the matrix.
!
!    Input, integer ( kind = 4 ) DF, the number of degrees of freedom.
!    M <= DF.
!
!    Input, real ( kind = 8 ) SIGMA(M,M), the covariance matrix, which should be 
!    a symmetric positive definite matrix.
!
!    Output, real ( kind = 8 ) A(M,M), the sample matrix from 
!    the Wishart distribution.
!
  implicit none
  type(mt_state), intent(inout) :: this
  integer ( kind = 4 ) m

  real ( kind = 8 ) a(m,m)
  real ( kind = 8 ), allocatable :: au(:,:)
  real ( kind = 8 ), allocatable :: aur(:,:)
  integer ( kind = 4 ) df
  integer ( kind = 4 ) flag
  real ( kind = 8 ), allocatable :: r(:,:)
  real ( kind = 8 ) sigma(m,m)

  if ( df < m ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'WISHART_SAMPLE - Fatal error!'
    write ( *, '(a)' ) '  DF = ', df, ' < M = ', m
    stop 1
  end if
!
!  Get R, the upper triangular Cholesky factor of SIGMA.
!
  allocate ( r(1:m,1:m) )
  call r8mat_cholesky_factor_upper ( m, sigma, r, flag )

  if ( flag /= 0 ) then
    write ( *, '(a)' ) ''
    write ( *, '(a)' ) 'WISHART_SAMPLE - Fatal error!'
    write ( *, '(a)' ) &
      '  Unexpected error return from R8MAT_CHOLESKY_FACTOR_UPPER.'
    write ( *, '(a,i4)' ) '  FLAG = ', flag
    stop 1
  end if
!
!  Get AU, a sample from the unit Wishart distribution.
!
  allocate ( au(1:m,1:m) )
  call wishart_unit_sample ( m, df, au,this )
!
!  Construct the matrix A = R' * AU * R.
!
  allocate ( aur(1:m,1:m) )
  aur = matmul ( au(1:m,1:m), r(1:m,1:m) )

  a = matmul ( transpose ( r(1:m,1:m) ), aur(1:m,1:m) )
!
!  Free memory.
!
  deallocate ( au )
  deallocate ( aur )
  deallocate ( r )

  return
end subroutine wishart_sample 

subroutine wishart_sample_inverse ( m, df, sigma, a, this )

!*****************************************************************************80
!
!! WISHART_SAMPLE_INVERSE returns the inverse of a sample Wishart matrix.
!
!  Discussion:
!
!    This function requires functions from the PDFLIB and RNGLIB libraries.
!
!    The "initialize()" function from RNGLIB must be called before using
!    this function.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 October 2013
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Patrick Odell, Alan Feiveson,
!    A numerical procedure to generate a sample covariance matrix,
!    Journal of the American Statistical Association,
!    Volume 61, Number 313, March 1966, pages 199-203.
!
!    Stanley Sawyer,
!    Wishart Distributions and Inverse-Wishart Sampling,
!    Washington University,
!    30 April 2007, 12 pages.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the order of the matrix.
!
!    Input, integer ( kind = 4 ) DF, the number of degrees of freedom.
!    M <= DF.
!
!    Input, real ( kind = 8 ) SIGMA(M,M), the covariance matrix, which should be 
!    a symmetric positive definite matrix.
!
!    Output, real ( kind = 8 ) A(M,M), the inverse of a sample matrix from the 
!    Wishart distribution.
!
  implicit none
  type(mt_state), intent(inout) :: this
  integer ( kind = 4 ) m

  real ( kind = 8 ) a(m,m)
  integer ( kind = 4 ) df
  integer ( kind = 4 ) flag
  real ( kind = 8 ), allocatable :: r(:,:)
  real ( kind = 8 ), allocatable :: s(:,:)
  real ( kind = 8 ) sigma(m,m)
  real ( kind = 8 ), allocatable :: ua(:,:)

  if ( df < m ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'WISHART_SAMPLE_INVERSE - Fatal error!'
    write ( *, '(a)' ) '  DF = ', df, ' < M = ', m
    stop 1
  end if
!
!  Get R, the upper triangular Cholesky factor of SIGMA.
!
  allocate ( r(1:m,1:m) )
  call r8mat_cholesky_factor_upper ( m, sigma, r, flag )

  if ( flag /= 0 ) then
    write ( *, '(a)' ) ''
    write ( *, '(a)' ) 'WISHART_SAMPLE_INVERSE - Fatal error!'
    write ( *, '(a)' ) &
      '  Unexpected error return from R8MAT_CHOLESKY_FACTOR_UPPER.'
    write ( *, '(a,i4)' ) '  FLAG = ', flag
    stop 1
  end if
!
!  Get S, the inverse of R.
!
  allocate ( s(1:m,1:m) )
  s(1:m,1:m) = r(1:m,1:m)
  call r8ut_inverse ( m, s )
!
!  Get UA, the inverse of a sample from the unit Wishart distribution.
!
  allocate ( ua(1:m,1:m) )
  call wishart_unit_sample_inverse ( m, df, ua, this )
!
!  Construct the matrix A = S * UA * S'.
!
  a = matmul ( s, matmul ( ua, transpose ( s ) ) )
!
!  Free memory.
!
  deallocate ( r )
  deallocate ( s )
  deallocate ( ua )

  return
end subroutine wishart_sample_inverse
subroutine wishart_unit_sample ( m, df, a, this )

!*****************************************************************************80
!
!! WISHART_UNIT_SAMPLE samples the unit Wishart distribution.
!
!  Discussion:
!
!    This function requires functions from the PDFLIB and RNGLIB libraries.
!
!    The "initialize()" function from RNGLIB must be called before using
!    this function.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 October 2013
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Patrick Odell, Alan Feiveson,
!    A numerical procedure to generate a sample covariance matrix,
!    Journal of the American Statistical Association,
!    Volume 61, Number 313, March 1966, pages 199-203.
!
!    Stanley Sawyer,
!    Wishart Distributions and Inverse-Wishart Sampling,
!    Washington University,
!    30 April 2007, 12 pages.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the order of the matrix.
!
!    Input, integer ( kind = 4 ) DF, the number of degrees of freedom.
!    M <= DF.
!
!    Output, double A(M,M), the sample matrix from the 
!    unit Wishart distribution.
!
  implicit none
  type(mt_state), intent(inout) :: this
  integer ( kind = 4 ) m

  real ( kind = 8 ) a(m,m)
  real ( kind = 8 ), allocatable :: c(:,:)
  integer ( kind = 4 ) df
  real ( kind = 8 ) df_chi
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
!  real ( kind = 8 ) r8_chi_sample
!  real ( kind = 8 ) r8_normal_01_sample

  if ( df < m ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'WISHART_UNIT_SAMPLE - Fatal error!'
    write ( *, '(a)' ) '  DF = ', df, ' < M = ', m
    stop 1
  end if

  allocate ( c(1:m,1:m) )

  do i = 1, m

    do j = 1, i - 1
      c(i,j) = 0.0D+00
    end do

    df_chi = real ( df + 1 - i, kind = 8 )
    c(i,i) = sqrt ( random_chisq(df_chi, this)) !sqrt ( r8_chi_sample ( df_chi ) )

    do j = i + 1, m
      c(i,j) = random_normal(this) !r8_normal_01_sample ( )
    end do

  end do

  a = matmul ( transpose ( c ), c )
!
!  Free memory.
!
  deallocate ( c )

  return
end subroutine wishart_unit_sample
subroutine wishart_unit_sample_inverse ( m, df, a, this )

!*****************************************************************************80
!
!! WISHART_UNIT_SAMPLE_INVERSE inverts a unit Wishart sample matrix.
!
!  Discussion:
!
!    This function requires functions from the PDFLIB and RNGLIB libraries.
!
!    The "initialize()" function from RNGLIB must be called before using
!    this function.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 October 2013
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Patrick Odell, Alan Feiveson,
!    A numerical procedure to generate a sample covariance matrix,
!    Journal of the American Statistical Association,
!    Volume 61, Number 313, March 1966, pages 199-203.
!
!    Stanley Sawyer,
!    Wishart Distributions and Inverse-Wishart Sampling,
!    Washington University,
!    30 April 2007, 12 pages.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the order of the matrix.
!
!    Input, integer ( kind = 4 ) DF, the number of degrees of freedom.
!    M <= DF.
!
!    Output, real ( kind = 8 ) A(M,M), the inverse of a sample matrix from 
!    the unit Wishart distribution.
!
  implicit none
  type(mt_state), intent(inout) :: this
  integer ( kind = 4 ) m

  real ( kind = 8 ) a(m,m)
  real ( kind = 8 ), allocatable :: b(:,:)
  real ( kind = 8 ), allocatable :: c(:,:)
  integer ( kind = 4 ) df
  real ( kind = 8 ) df_chi
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
!  real ( kind = 8 ) r8_chi_sample
!  real ( kind = 8 ) r8_normal_01_sample

  if ( df < m ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'WISHART_UNIT_SAMPLE_INVERSE - Fatal error!'
    write ( *, '(a,i6,a,i6)' ) '  DF = ', df, ' < M = ', m
    stop 1
  end if
!
!  Compute C, an upper triangular matrix such that the
!  Wishart sample matrix is C' * C.
!
  allocate ( c(1:m,1:m) )

  do i = 1, m
    df_chi = real ( df - i + 1, kind = 8 )
    c(i,i) = sqrt ( random_chisq(df_chi,this) ) ! sqrt ( r8_chi_sample ( df_chi ) )
    do j = i + 1, m
      c(i,j) = random_normal(this) !r8_normal_01_sample ( )
    end do
  end do
!
!  Compute B, the inverse of C.
!
  allocate ( b(1:m,1:m) )
  b(1:m,1:m) = c(1:m,1:m)
  call r8ut_inverse ( m, b )
!
!  The inverse of the Wishart sample matrix C'*C is inv(C) * C'.
!
  a = matmul ( b, transpose ( b ) )
!
!  Free memory.
!
  deallocate ( b )
  deallocate ( c )

  return
end subroutine wishart_unit_sample_inverse 

end module mcmc_utils




