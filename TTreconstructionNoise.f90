!code that computes the noise for the reconstructed TT dipole
!Paper should be out soon(-ish).

program TTreconstructionNoise


  implicit none

  integer :: i,m,il,t,j,s, ilmin, il1, il2
  integer, parameter :: dl= KIND(1.d0)
  real(dl), allocatable :: powerspectra(:,:)
  real(dl), allocatable :: lensedps(:,:)
  real(dl), allocatable :: kszps(:)
  real(dl), allocatable :: NoiseSpectra(:,:)
  integer, allocatable :: ell(:)
  real(dl) :: NoiseInvPartial, NoiseInvSum, NoiseTotal, STN, errorbar
  integer :: lmax, lmin, ilsubmax
  real(dl)  :: a3j(0:20000)
  real(dl), allocatable  :: w3j(:,:,:)
  real(dl), parameter :: pi = 3.14159265359
  integer :: ell1, ell2, ell3
  integer :: cvon = 1 
  logical :: want_lensedTT, want_kSZ, want_exp
  integer :: lensed_on = 0, ksz_on = 0, forecast_on = 0
  real :: dummy
  want_lensedTT  = .false.
  want_kSZ = .false.
  !forecast using experimental noise
  want_exp = .false.

  if(want_kSZ) ksz_on = 1
  if(want_lensedTT) lensed_on = 1
  if(want_exp) then
     forecast_on = 1
     cvon = 0
     !lensed_on = 1
  endif
  
  !CAMB spectra
  !1 = TT, 2 = EE, 3 = BB, 4 = TE, 5 = PP, 6 = TP, 7 = EP
  !Noise spectra
  !1 = TT, 2 = EE, 3 = BB, 4 = TT(delensed), 5 = dd, 6 = PP
  lmax = 2000
  lmin = 2
  allocate(powerspectra(7,lmin:lmax))
  allocate(lensedps(4,lmin:lmax))
  allocate(NoiseSpectra(6,lmin:lmax))
  
  !kZS signal; should be added to the noise (see https://arxiv.org/pdf/1109.0553.pdf)
  allocate(kszps(lmin:lmax))
  open(unit=27,file="kszpowerspectrum.dat", status="old")

  do il1 = lmin, lmax
     read(27,*) m, kszps(il1)
  enddo
  close(27)
  
  allocate(ell(lmin:lmax))
  !get Cl's
  !note that the following file only includes unlensed spectra
  open(unit=28,file="Feb2017_highEllSpectra_lenspotentialCls.dat", status="old")
  !can't fully delens yet. If not delensed should use these spectra. 
  open(unit=29,file="Feb2017_highEllSpectra_lensedCls.dat", status="old")
  !specify file for experimental noise:
  !Stage4: 1 muK-arcmin, Stage4.5: .3 muK-arcmin, Stage5: 0.1 muK-arcmin
  !Satge5.5: 0.03 muK-arcmin, Stage6: 0.01 muK-arcmin, Stage6.5: 0.003 muK-arcmin
  !Stage7: 0.001 muK-arcmin
  !open(unit=30,file="Spectra/noise_spectra_with_kSZ_stage7.txt", status="old")
  !new CMBS4 noise, ell, TT, EE, phiphi
  open(unit=30,file="Spectra/CMBS4noise.txt", status="old")
  do il = 1, lmax
     if (il .eq. 1) then 
        read(28,*)
        read(29,*)
        !read(30,*)
        !write(*,*) ell(il), powerspectra(1:4,il)
     else
        read(28,*) ell(il), powerspectra(1:7,il)
        read(29,*) ell(il), lensedps(1:4,il)
        !read(30,*) dummy, NoiseSpectra(1:6,il)
        read(30,*) dummy, NoiseSpectra(1,il), NoiseSpectra(2,il), NoiseSpectra(6,il)
        !write(*,*) il,ell(il)
        !removing artifical noise
        !if(ell(il) .ge. 2999) then 
        !   NoiseSpectra(1,il) = NoiseSpectra(1,il) - 1000
           !NoiseSpectra(4,il) = lensedps(1,il)/fac(ell(il))
        !endif 
        !write(*,'(I5, 4E17.7)') ell(il), powerspectra(1:4,il)
        !write(*,*) ell(il), lensedPS(1,il), fac(ell(il))*NoiseSpectra(4,il)
     endif

  enddo
  close(28)
  close(29)
  close(30)
  STN = 0.d0
  open(unit=28,file="hatT2_nokSZ_CV_lmax2000.dat", status="replace")
  ilmin = 100
  do s = 2, 500
     ilmin = 500
     NoiseInvSum = 0.d0
     !ilsubmax = 1000+(s-1)*50
     !hardcode
     !ilsubmax = 5000
     !$OMP PARALLEL DO DEFAUlT(SHARED), &
     !$OMP PRIVATE(i,j,NoiseInvPartial,ell1,ell3,a3j) &
     !$OMP REDUCTION(+:NoiseInvSum)
     do il1 = lmin, lmax !l1 loop
        !ell1 = ell(il1)
        !mode you want to reconstruct 
        ell3 = ell(s) !dipole l = 1 !can change this to any mode...
        !write(*,*) ell(s)
        !get Wigner 3j for l2:
        call GetThreeJs(a3j((abs(ell3-il1)):(ell3+il1)),il1,ell3,0,0)
        !if(i .eq. ilmin .and. s .eq. 140) write(*,*) ilmin, a3j(ell(ilmin+1))

        !l2 loop. For a given L (ell3), triangle condition applies, i.e.
        !abs(l1-l3) < l2 < min(lmax, (l3+l1))
        do il2 = max(lmin,abs(il1-ell3)), min(lmax,il1+ell3)
!!$           NoiseInvPartial = ((2.*ell1+1.)*(2.*ell(j)+1.)/(4.*pi))*parity(ell1+ell(j)+ell3)* &
!!$                (prefactor(ell1,ell(j),ell3)*a3j(ell(j)))**2* &
!!$                (1./(((1.-lensed_on)*powerspectra(1,i)+(lensed_on)*(forecast_on-1)*lensedps(1,i)+ksz_on*kszps(i))/ &
!!$                fac(ell1)+(forecast_on)*(NoiseSpectra(1,j) + powerspectra(1,j))))* &
!!$                (powerspectra(5,j)/(ell(j)*(ell(j)+1.)*fac(ell(j))))**2/(powerspectra(5,j)/(ell(j)*(ell(j)+1.)*fac(ell(j)))+ &
!!$                (forecast_on)*NoiseSpectra(6,j))
           NoiseInvPartial = ((2.*il1+1.)*(2.*il2+1.)/(4.*pi))*parity(il1+il2+ell3)* &
                (prefactor(il1,il2,ell3)*a3j(il2))**2* &
                (1./((powerspectra(1,il1)+ksz_on*kszps(il1))/ &
                fac(il1)+(1-cvon)*(NoiseSpectra(1,il1))))* &
                (powerspectra(5,il2)/(fac(il2)*fac2(il2)))**2/(powerspectra(5,il2)/(fac2(il2)*fac2(il2))+ &
                (1-cvon)*NoiseSpectra(6,il2)/fac(il2))
           NoiseInvSum = NoiseInvSum + NoiseInvPartial     
           !write(*,'(I4,5E17.6)') j, NoiseInvSum, NoiseInvPartial, parity(ell(i)+ell(j)+1), a3j(j), powerspectra(5,j)!w3j(ell(i),ell(j),1)
        enddo
     enddo
     !$OMP END PARALLEL DO
     NoiseTotal = 1./NoiseInvSum
     STN = STN + (2.*ell3+1.)*powerspectra(1,s)**2/(ell3*(ell3+1.)*NoiseTotal/pi/2.0+powerspectra(1,s))**2
     errorbar = (ell3*(ell3+1.)*NoiseTotal/pi/2.0+powerspectra(1,s))/sqrt( (2.*ell3+1.))/sqrt(10.0)
     !write(*,*) "ellmax:",ell(ilsubmax),"TT reconstruction Noise:", NoiseTotal/pi
     write(*,'(I5,4E17.6)') ell3, powerspectra(1,s), ell3*(ell3+1.)*NoiseTotal/pi/2.0, sqrt(STN), errorbar
     !write(28,'(I5,E17.6)') ell(ilsubmax), NoiseTotal/pi
     
     write(28,'(I5,4E17.6)') ell3, powerspectra(1,s), ell3*(ell3+1.)*NoiseTotal/pi/2.0, sqrt(STN), errorbar
  enddo
  
  close(28)
  
  deallocate(powerspectra)
  deallocate(ell)
  deallocate(NoiseSpectra)
  deallocate(lensedps)
  deallocate(kszps)
   contains

     real function prefactor(l1,l2,l3)
       integer, intent(in) :: l1,l2,l3

       prefactor = 1./2.*(-1.0*l1*(l1+1.0)+1.0*l2*(l2+1.0)+1.0*l3*(l3+1.0))
     end function prefactor

     real function parity(n)
       integer, intent(in) :: n

       if (MOD(n,2)==0) then
          parity = 1.0
       else
          parity = 0.0
       endif
     end function parity

     real function fac(l)
       integer, intent(in) :: l
       fac = 1.d0*l*(l+1.)/(2.*pi)
     end function fac
     real function fac2(l)
       integer, intent(in) :: l
       fac2 =1.d0*l*(l+1.) 
     end function fac2
     real function Noise_l(sigmab,w0,ell)
       real(dl) :: sigmab
       real(dl) :: w0
       integer :: ell
       real(dl) :: arcmin
       real(dL) :: fac
       arcmin = 0.000290888
       fac = arcmin / 2.d0*sqrt(2.d0*log(2.d0))

       Noise_l = w0**2*arcmin**2*exp(ell**2*fac**2*sigmab**2)

     end function Noise_l

     subroutine GetThreeJs(thrcof,l2in,l3in,m2in,m3in)
       !Recursive evaluation of 3j symbols. Does minimal error checking on input
       !parameters.
       implicit none
       integer, parameter :: dl = KIND(1.d0)
       integer, intent(in) :: l2in,l3in, m2in,m3in
       real(dl), dimension(*) :: thrcof
       INTEGER, PARAMETER :: i8 = selected_int_kind(18)
       integer(i8) :: l2,l3,m2,m3
       integer(i8) :: l1, m1, l1min,l1max, lmatch, nfin, a1, a2

       real(dl) :: newfac, oldfac, sumfor, c1,c2,c1old, dv, denom, x, sum1,sumuni
       real(dl) :: x1,x2,x3, y,y1,y2,y3,sum2,sumbac, ratio,cnorm, sign1, thresh
       integer i,ier, index, nlim, sign2
       integer nfinp1,nfinp2,nfinp3, lstep, nstep2,n
       real(dl), parameter :: zero = 0._dl, one = 1._dl
       real(dl), parameter ::  tiny = 1.0d-30, srtiny=1.0d-15, huge = 1.d30,srhuge = 1.d15

       ! routine to generate set of 3j-coeffs (l1,l2,l3\\ m1,m2,m3)

       ! by recursion from l1min = max(abs(l2-l3),abs(m1)) 
       !                to l1max = l2+l3
       ! the resulting 3j-coeffs are stored as thrcof(l1-l1min+1)

       ! to achieve the numerical stability, the recursion will proceed
       ! simultaneously forwards and backwards, starting from l1min and l1max
       ! respectively.
       !
       ! lmatch is the l1-value at which forward and backward recursion are
       ! matched.
       !
       ! ndim is the length of the array thrcof
       !
       ! ier = -1 for all 3j vanish(l2-abs(m2)<0, l3-abs(m3)<0 or not integer)
       ! ier = -2 if possible 3j's exceed ndim
       ! ier >= 0 otherwise

       l2=l2in
       l3=l3in
       m2=m2in
       m3=m3in
       newfac = 0
       lmatch = 0
       m1 = -(m2+m3)

       ! check relative magnitude of l and m values
       ier = 0

       if (l2 < abs(m2) .or. l3 < m3) then
          ier = -1
          ! call MpiStop('error ier = -1')
          print*, 'error ier = -1'
          stop
          return
       end if

       ! limits for l1
       l1min = max(abs(l2-l3),abs(m1))
       l1max = l2+l3

       if (l1min >= l1max) then
          if (l1min/=l1max) then
             ier = -1

             !call MpiStop('error ier = -1')
             print*, 'error ier = -1' 
             stop
             return
          end if

          ! reached if l1 can take only one value, i.e.l1min=l1max
          thrcof(1) = (-1)**abs(l2+m2-l3+m3)/sqrt(real(l1min+l2+l3+1,dl))
          return

       end if

       nfin = l1max-l1min+1

       ! starting forward recursion from l1min taking nstep1 steps
       l1 = l1min
       thrcof(1) = srtiny
       sum1 = (2*l1 + 1)*tiny

       lstep = 1

30     lstep = lstep+1
       l1 = l1+1

       oldfac = newfac
       a1 = (l1+l2+l3+1)*(l1-l2+l3)*(l1+l2-l3)
       a2 = (l1+m1)*(l1-m1)*(-l1+l2+l3+1)
       newfac = sqrt(a2*real(a1,dl))
       if (l1 == 1) then
          !IF L1 = 1  (L1-1) HAS TO BE FACTORED OUT OF DV, HENCE
          c1 = -(2*l1-1)*l1*(m3-m2)/newfac
       else

          dv = -l2*(l2+1)*m1 + l3*(l3+1)*m1 + l1*(l1-1)*(m3-m2)
          denom = (l1-1)*newfac

          if (lstep > 2) c1old = abs(c1)
          c1 = -(2*l1-1)*dv/denom

       end if

       if (lstep<= 2) then

          ! if l1=l1min+1 the third term in the recursion eqn vanishes, hence
          x = srtiny*c1
          thrcof(2) = x
          sum1 = sum1+tiny*(2*l1+1)*c1*c1
          if(lstep==nfin) then
             sumuni=sum1
             go to 230
          end if
          goto 30

       end if

       c2 = -l1*oldfac/denom

       ! recursion to the next 3j-coeff x  
       x = c1*thrcof(lstep-1) + c2*thrcof(lstep-2)
       thrcof(lstep) = x
       sumfor = sum1
       sum1 = sum1 + (2*l1+1)*x*x
       if (lstep/=nfin) then

          ! see if last unnormalised 3j-coeff exceeds srhuge
          if (abs(x) >= srhuge) then

             ! REACHED IF LAST 3J-COEFFICIENT LARGER THAN SRHUGE
             ! SO THAT THE RECURSION SERIES THRCOF(1), ... , THRCOF(LSTEP)
             ! HAS TO BE RESCALED TO PREVENT OVERFLOW

             ier = ier+1
             do i = 1, lstep
                if (abs(thrcof(i)) < srtiny) thrcof(i)= zero
                thrcof(i) = thrcof(i)/srhuge
             end do

             sum1 = sum1/huge
             sumfor = sumfor/huge
             x = x/srhuge

          end if

          ! as long as abs(c1) is decreasing, the recursion proceeds towards
          ! increasing
          ! 3j-valuse and so is numerically stable. Once an increase of abs(c1) is 
          ! detected, the recursion direction is reversed.

          if (c1old > abs(c1)) goto 30

       end if !lstep/=nfin

       ! keep three 3j-coeffs around lmatch for comparison with backward recursion

       lmatch = l1-1
       x1 = x
       x2 = thrcof(lstep-1)
       x3 = thrcof(lstep-2)
       nstep2 = nfin-lstep+3

       ! --------------------------------------------------------------------------
       !
       ! starting backward recursion from l1max taking nstep2 stpes, so that
       ! forward and backward recursion overlap at 3 points 
       ! l1 = lmatch-1, lmatch, lmatch+1

       nfinp1 = nfin+1
       nfinp2 = nfin+2
       nfinp3 = nfin+3
       l1 = l1max
       thrcof(nfin) = srtiny
       sum2 = tiny*(2*l1+1)

       l1 = l1+2
       lstep=1

       do
          lstep = lstep + 1
          l1= l1-1

          oldfac = newfac
          a1 = (l1+l2+l3)*(l1-l2+l3-1)*(l1+l2-l3-1)
          a2 = (l1+m1-1)*(l1-m1-1)*(-l1+l2+l3+2)
          newfac = sqrt(a1*real(a2,dl))

          dv = -l2*(l2+1)*m1 + l3*(l3+1)*m1 +l1*(l1-1)*(m3-m2)

          denom = l1*newfac
          c1 = -(2*l1-1)*dv/denom
          if (lstep <= 2) then

             ! if l2=l2max+1, the third term in the recursion vanishes

             y = srtiny*c1
             thrcof(nfin-1) = y
             sumbac = sum2
             sum2 = sum2 + tiny*(2*l1-3)*c1*c1

             cycle

          end if

          c2 = -(l1-1)*oldfac/denom

          ! recursion to the next 3j-coeff y
          y = c1*thrcof(nfinp2-lstep)+c2*thrcof(nfinp3-lstep)

          if (lstep==nstep2) exit

          thrcof(nfinp1-lstep) = y
          sumbac = sum2
          sum2 = sum2+(2*l1-3)*y*y

          ! see if last unnormalised 3j-coeff exceeds srhuge
          if (abs(y) >= srhuge) then

             ! reached if 3j-coeff larger than srhuge so that the recursion series
             ! thrcof(nfin),..., thrcof(nfin-lstep+1) has to be rescaled to prevent
             ! overflow

             ier=ier+1
             do i = 1, lstep
                index=nfin-i+1
                if (abs(thrcof(index)) < srtiny) thrcof(index)=zero
                thrcof(index) = thrcof(index)/srhuge
             end do

             sum2=sum2/huge
             sumbac=sumbac/huge

          end if

       end do

       ! the forward recursion 3j-coeffs x1, x2, x3 are to be matched with the 
       ! corresponding backward recursion vals y1, y2, y3

       y3 = y
       y2 = thrcof(nfinp2-lstep)
       y1 = thrcof(nfinp3-lstep)

       ! determine now ratio such that yi=ratio*xi (i=1,2,3) holds with minimal
       ! error

       ratio = (x1*y1+x2*y2+x3*y3)/(x1*x1+x2*x2+x3*x3)
       nlim = nfin-nstep2+1

       if (abs(ratio) >= 1) then

          thrcof(1:nlim) = ratio*thrcof(1:nlim) 
          sumuni = ratio*ratio*sumfor + sumbac

       else

          nlim = nlim+1
          ratio = 1/ratio
          do n = nlim, nfin
             thrcof(n) = ratio*thrcof(n)
          end do
          sumuni = sumfor + ratio*ratio*sumbac

       end if
       ! normalise 3j-coeffs

230    cnorm = 1/sqrt(sumuni)

       ! sign convention for last 3j-coeff determines overall phase

       sign1 = sign(one,thrcof(nfin))
       sign2 = (-1)**(abs(l2+m2-l3+m3))
       if (sign1*sign2 <= 0) then
          cnorm = -cnorm
       end if
       if (abs(cnorm) >= one) then
          thrcof(1:nfin) = cnorm*thrcof(1:nfin)
          return
       end if

       thresh = tiny/abs(cnorm)

       do n = 1, nfin
          if (abs(thrcof(n)) < thresh) thrcof(n) = zero
          thrcof(n) = cnorm*thrcof(n)
       end do
       return 

     end subroutine GetThreeJs



   end program TTreconstructionNoise



  
