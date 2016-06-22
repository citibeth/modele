#include "rundeck_opts.h"

      subroutine obio_alkalinity_topaz(vrbos,kmax,i,j)

!@sum  online computation of alkalinity
!@auth Natassa Romanou

!-------------------------------------------------------------------------

      USE obio_dim
      USE obio_incom, only: rain_ratio,cpratio,sigma_Ca,d_Ca,
     .      npratio,uMtomgm3,cnratio,bn,zc,mgchltouMC
      USE obio_com, only: P_tend,p1d,pp2_1d,dp1d,A_tend,
     .      rhs,alk1d,caexp,kzc,co3_conc

#ifdef OBIO_ON_GARYocean
      USE MODEL_COM, only: nstep=> itime
      USE OCEAN, only: dxypo,lmm
#else
      USE hycom_scalars, only: nstep
      USE hycom_arrays, only: scp2
#endif

      implicit none

      integer nt,k,kmax,nchl1,nchl2,i,j
      real*8 J_PO4(kmax),pp,Jprod(kmax),Jprod_sum,Fc,zz,F_Ca(kmax+1),
     .       J_Ca(kmax),term,term1,term2,DOP,offterm
      logical vrbos
!--------------------------------------------------------------------------
!only compute tendency terms if total depth greater than conpensation depth
      if (p1d(kmax+1) .lt. p1d(kzc)) then
         A_tend(:) = 0.
         go to 100
      endif

!pp and therefore Jprod are defined at dp points (mid level)
!F_Ca is defined at pressure interfaces
!J_Ca is therefore defined mid-layer (dp points)

!compute sources/sinks of phosphate
!J_PO4 units uM/hr
      do k=1,kmax
!     J_PO4(k) =  P_tend(k,1)           !approximate by nitrate conc tendency
!                                       !NO3/PO4 ratio from Conkright et al, 1994
!                                       !expressed as nitrates
!                                       ! in notes npratio that multiplies tendency term cancells 
!                                       ! out with npratio that divides to get units of nitrate
!     J_PO4(k) =  rhs(k,1,5)+rhs(k,1,6)+rhs(k,1,7)+rhs(k,1,8)    !uptake of nitrate
!     J_PO4(k) =((rhs(k,1,5)+rhs(k,1,6)+rhs(k,1,7)+rhs(k,1,8))/bn
!    .         +  rhs(k,14,6)/mgchltouMC
!    .         +(rhs(k,2,5)+rhs(k,2,6)
!    .         + rhs(k,2,7)+rhs(k,2,8))/bn)/cnratio

      J_PO4(k) =  rhs(k,1,5)+rhs(k,1,6)+rhs(k,1,7)+rhs(k,1,8)
     .          +(rhs(k,5,13)+rhs(k,6,13)+rhs(k,7,13)+rhs(k,8,13))*bn
     .          + rhs(k,2,5)+rhs(k,2,6)+rhs(k,2,7)+rhs(k,2,8)

      term = -1.d0* J_PO4(k)            !uM,N/hr= mili-mol,N/m3/hr
      rhs(k,15,1) = term
      A_tend(k)= term 
      enddo

      call cadet_topaz(kmax,J_Ca,vrbos,i,j)

      caexp = J_Ca(kzc)        !mili-gC/m2/hr

!compute sources/sinks of CaCO3
      offterm= 0.d0
      do k=1,kmax
!        if (p1d(k) .le. p1d(kzc))  then
!            !formation of calcium carbonate above compensation depth
!            J_Ca(k) = -1.d0*rain_ratio*cpratio*(1.-sigma_Ca)*Jprod(k)   !mgC/m3/hr
!        else
!            !dissolution of calcium carbonate below compensation depth
!            J_Ca(k) = -1.d0* (F_Ca(k+1)-F_Ca(k)) / dp1d(k)   !mgC/m3/hr
!        endif
       term = 2.d0* J_Ca(k)/cnratio  ! mgC/m3/hr -> mili-mol,N/m3/hr  
                                     ! no need to multiply here by mol.weight
                                     ! because already in cnratio (see obio_init)
       rhs(k,15,5) = term
       A_tend(k) = A_tend(k) + term      

       offterm = offterm + term * dp1d(k)

      if(vrbos)
     .write(*,'(a,4i5,3e12.4)')'obio_alkalinity2:',
     .   nstep,i,j,k,term,offterm,rhs(k,15,5)


      !bottom boundary condition adjust bottom layer 
      if (k.eq.kmax) then
         rhs(kmax,15,5) = rhs(kmax,15,5) - offterm /dp1d(kmax)
         A_tend(kmax) = A_tend(kmax) - offterm / dp1d(kmax)
      endif

      if(vrbos)
     .write(*,'(a,4i5,3e12.4)')'obio_alkalinity3:',
     .   nstep,i,j,k,term,offterm,rhs(k,15,5)

      enddo

      !for consistency, keep term that goes into rhs table in uM/hr = mili-mol,N/m3/hr
      !convert A_tend terms into uE/kg/hr, the actual units of alkalinity 
      A_tend = A_tend /1024.5d0 *1.d3     ! mili-mol,N/m3/hr -> umol/m3/hr -> umol/kg/hr

 

 100  continue
      end subroutine obio_alkalinity_topaz

      subroutine cadet_topaz(kmax,J_Ca,vrbos,i,j)

      USE MODEL_COM, only: dtsrc, nstep=> itime
      USE obio_dim
      USE obio_incom, only: cnratio,mgchltouMC
      USE obio_com, only: co3_conc,temp1d,saln1d,obio_P,p1d
     .                   ,wsdet,ca_det_calc1d,Ca_tend

      implicit none
      
      integer k,kmax,i,j
      real*8 T,Salt,TK,PKSPC,wsink, gamma_cadet_calc
      real*8 co3_calc_sol,Omega_calc,J_cadet_calc
      real*8 lambda0,KEppley,P_star,P_min,CaN_calc_ratio,Omega_satmax
      real*8 N_cyanob,jgraz_cyanob_N,P_insitu,J_prod_cadet_calc
      real*8 J_Ca(kmax)
      logical vrbos

      do k=1,kmax
      !rate of constant dissolution of cadet_calc
      !wsink=1.d0
      wsink =  wsdet(k,1)    !wsdet for carbon????
      gamma_cadet_calc = wsink/1343.d0    ! in s-1
      gamma_cadet_calc = gamma_cadet_calc * 3600.d0     ! in hr-1

      if(vrbos)write(*,'(a,4i5,2e12.4)')'cadet1',
     .      nstep,i,j,k,wsink,gamma_cadet_calc

      T = temp1d(k)
      Salt = saln1d(k)
      TK = T + 273.15
      P_insitu = 0.1016*p1d(k)+1.013
      PKSPC = 171.9065 
     .      + 0.077993 * TK 
     .      - 2903.293 / TK 
     .      - 71.595* dlog10(TK) 
     .      - (-0.77712 + 2.8426d-3 *TK + 178.34/TK)*Salt**(1./2.)
     .      + 0.07711*salt 
     .      - 4.1249d-3 * salt**(3./2.) 
     .      - 0.02 
     .      - (48.76 - 0.5304*T) * (P_insitu - 1.013) / (191.46*TK) 
     .      + (1.d-3 *(11.76 -0.3692 *T))
     .      * (P_insitu -1.013)*(P_insitu-1.013)/(382.92*TK)
      co3_calc_sol = 10**(-PKSPC) / (2.937d-4 *max(5.d0, Salt))
      !??? co3_conc, how do I apply it over all depths?
      Omega_calc = co3_conc / co3_calc_sol

      !how do i compute ca_det_calc ????
      ! dissolution rate in units of mol kg-1 s-1 calculated in every layer, 
      ! and has to be multiplied by the thickness to get a layer flux.
      J_cadet_calc = 
     .     gamma_cadet_calc*max(0.d0,1.d0-Omega_calc)* Ca_det_calc1d(k) 

      if(vrbos)write(*,'(a,4i5,8e12.4)')'cadet2',
     .      nstep,i,j,k,T,Salt,P_insitu,PKSPC,co3_calc_sol,Omega_calc
     .     ,Ca_det_calc1d(k),J_cadet_calc

      !grazing rate constant at 0C
      lambda0 = 0.19/86400. *3600.   !hr^-1  
      !Temp. coeff for growth
      KEppley = 0.063   ! C^-1 (Eppley 1972) 
      !nitrogen in small phytoplankton (=here cyanobacteria)
      N_cyanob = obio_P(k,7) *mgchltouMC/cnratio  !uM,N/m3
      N_cyanob = 1.d-3*N_cyanob*1024.5d0   !mol,N/kg
      !pivot phyto. conc. for grazing allometry
      P_star = 1.9d-6 * 16/106   !mol,N/kg
      !min. phyto. conc. threshold for grazing
      P_min = 1.d-10   !mol,N/kg (added for stability)
      jgraz_cyanob_N = min(1.d0/dtsrc, 
     .                     lambda0 * exp(KEppley*T)
     .      *(N_cyanob * N_cyanob /
     .       (P_star *(N_cyanob + P_min)))) * N_cyanob
      !Calcite CaCO3 to nitrogen uptake ratio
      CaN_calc_ratio =  0.005 * 106/16    !mol,Ca/mol,N *** tunable to give right caexp
      !maximum saturation state
      Omega_satmax = 10.  !dimensionless; to limit possible extreme values
      J_prod_cadet_calc = jgraz_cyanob_N * CaN_calc_ratio 
     .                  * exp(-0.0539*T) 
     .                  * min(Omega_satmax,max(0.d0,Omega_calc - 1.d0))
      

      J_Ca(k) = J_cadet_calc - J_prod_cadet_calc
      Ca_tend(k) = J_Ca(k)

!change units in J_Ca (mol,Ca/kg/hr -> mili-mol,N/m3/hr)
      J_Ca(k) = 1.d3 * J_Ca(k) /CaN_calc_ratio * 1024.5d0


!     A_tend() = 2.d0 * (J_cadet_calc - J_prod_cadet_calc)

      enddo
  
      end subroutine cadet_topaz
