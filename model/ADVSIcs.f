#include "rundeck_opts.h"

      SUBROUTINE ADVSI
!@sum  ADVSI advects sea ice
!@+    Currently set up to advect ice on AGCM grid (i.e. usidt/vsidt are
!@+    on the AGCM grid, and RSI/MSI/HSI etc. are unchanged)
!@+    At some point this will change (USIDT/VSIDT on ice grid, and RSI
!@+    etc. will need to be interpolated back and forth).
!@auth Gary Russell/Gavin Schmidt
!@auth rewrite for cubed sphere: M. Kelley
c NOTE: CURRENTLY ASSUMING THAT THERE IS NO TRANSPORT OF ICE TO/FROM
c EQUATORIAL CUBE FACES.  WILL UPGRADE AS NEEDED.
      USE CONSTANT, only : byshi,lhm,grav
      USE MODEL_COM, only : im,focean,p,ptop,kocean,dts=>dtsrc
      USE DOMAIN_DECOMP_ATM, only : grid, GET, HALO_UPDATE
      USE GEOM, only : axyp,byaxyp,
     &     dlxsina,dlysina, ull2ucs,vll2ucs, ull2vcs,vll2vcs
      USE ICEDYN_COM, only : rsix,rsiy,rsisave,foa,byfoa
      USE SEAICE, only : ace1i,xsi,
     &     get_snow_ice_layer,set_snow_ice_layer,relayer,relayer_12
      USE SEAICE_COM, only : rsi,msi,snowi,hsi,ssi,lmi
#ifdef TRACERS_WATER
     *     ,trsi,ntm
#endif
      USE FLUXES, only : gtemp,apress,msicnv,fwsim,uisurf,visurf
#ifdef TRACERS_WATER
     *     ,gtracer
#endif
      USE DIAG_COM, only : oa,aij=>aij_loc,
     &     IJ_MUSI,IJ_MVSI,IJ_HUSI,IJ_HVSI,IJ_SUSI,IJ_SVSI
#ifdef TRACERS_WATER
      USE TRDIAG_COM, only : taijn=>taijn_loc,tij_tusi,tij_tvsi
#endif
      USE ICEDYN, only : grid_icdyn,usi,vsi
      USE ICEDYN_COM, only : i2a_uc,i2a_vc,UVLLATUC,UVLLATVC,CONNECT
      use cs2ll_utils, only : ll2csint_lij
      IMPLICIT NONE
!@var NTRICE max. number of tracers to be advected (mass/heat/salt+)
#ifndef TRACERS_WATER
      INTEGER, PARAMETER :: NTRICE=3*(LMI+2)
#else
      INTEGER, PARAMETER :: NTRICE=(3+NTM)*(LMI+2)
      INTEGER ITR
      REAL*8 TRSNOW(NTM,2), TRICE(NTM,LMI)
#endif
      INTEGER I,J,L,K
      REAL*8 DMHSI,ASI,YRSI,XRSI,FRSI,SICE,COUR,FAO,CNEW,
     &     ull,vll
C****
C**** USI, VSI  latlon-oriented U,V components of
C****           sea ice velocity (m/s)
C****
C**** FAW    flux of surface water area (m^2) = USIDT*DYP
C**** FASI   flux of sea ice area (m^2) = USIDT*DYP*RSIedge
C**** FMSI   flux of sea ice mass (kg) or heat (J) or salt (kg)
C****
      REAL*8, DIMENSION(ntrice,grid%I_STRT_HALO:grid%I_STOP_HALO) ::
     &     FMSI, FMSI_jm1
      REAL*8, DIMENSION(ntrice) :: AMSI, FMSI_im1
      REAL*8, DIMENSION(grid%I_STRT_HALO:grid%I_STOP_HALO) ::
     &     FASI, FXSI, FYSI, FAW,
     &     FASI_jm1, FXSI_jm1, FYSI_jm1, FAW_jm1
      REAL*8 ::
     &     FASI_im1, FXSI_im1, FYSI_im1, FAW_im1

!@var MHS mass/heat/salt content of sea ice
      REAL*8 MHS(NTRICE,grid%I_STRT_HALO:grid%I_STOP_HALO,
     &                  grid%J_STRT_HALO:grid%J_STOP_HALO)
      REAL*8 :: MICE(LMI),SICEg(LMI),HICE(LMI),TSIL(LMI),MSI1,FMSI2
      REAL*8 :: HSNOW(2),SNOWL(2),TSNW(2)

      REAL*8, DIMENSION(grid%I_STRT_HALO:grid%I_STOP_HALO,
     &                  grid%J_STRT_HALO:grid%J_STOP_HALO) ::
     &     UDYDT,VDXDT

      REAL*8, DIMENSION(2,grid_icdyn%im_world,
     &         grid_icdyn%J_STRT_HALO:grid_icdyn%J_STOP_HALO) ::
     &     uvll

      INTEGER I_0,I_1,J_0,J_1, I_0Y,I_1Y
!@var coastfac: A proportionality factor to compute the component
!@+   of advective velocity which limits ice buildup along
!@+   coastlines.  (At some gridcells, negative feedbacks on
!@+   ice production are not able to assert themselves when
!@+   sea ice does not reside on the ocean grid.)
      real*8 :: coastfac

C**** Get grid parameters
      CALL GET(grid, I_STRT=I_0, I_STOP=I_1, J_STRT=J_0, J_STOP=J_1)

      i_0y = max(1 ,i_0-1)
      i_1y = min(im,i_1+1)

      SNOWL=0.; HSNOW=0.; HICE=0.; SICEg=0.; TSNW=0.; TSIL=0.; MICE=0

      DO J=J_0,J_1
      DO I=I_0,I_1
C**** Reduce ice concentration gradients if ice amounts decreased
        IF (RSI(I,J).gt.1d-4) THEN
          IF (RSISAVE(I,J).gt.RSI(I,J)) THEN ! reduce gradients
            FRSI=(RSISAVE(I,J)-RSI(I,J))/RSISAVE(I,J)
            RSIX(I,J)=RSIX(I,J)*(1.-FRSI)
            RSIY(I,J)=RSIY(I,J)*(1.-FRSI)
          END IF
          IF(RSI(I,J)-RSIX(I,J).lt.0.)  RSIX(I,J) =    RSI(I,J)
          IF(RSI(I,J)+RSIX(I,J).lt.0.)  RSIX(I,J) =   -RSI(I,J)
          IF(RSI(I,J)-RSIX(I,J).gt.1d0) RSIX(I,J) =    RSI(I,J)-1d0
          IF(RSI(I,J)+RSIX(I,J).gt.1d0) RSIX(I,J) =1d0-RSI(I,J)
          IF(RSI(I,J)-RSIY(I,J).lt.0.)  RSIY(I,J) =    RSI(I,J)
          IF(RSI(I,J)+RSIY(I,J).lt.0.)  RSIY(I,J) =   -RSI(I,J)
          IF(RSI(I,J)-RSIY(I,J).gt.1d0) RSIY(I,J) =    RSI(I,J)-1d0
          IF(RSI(I,J)+RSIY(I,J).gt.1d0) RSIY(I,J) =1d0-RSI(I,J)
        ELSE
          RSIX(I,J) = 0.  ; RSIY(I,J) = 0.
        END IF
C**** update RSISAVE for diagnostics
        RSISAVE(I,J)=RSI(I,J)
C**** set up local MHS array to contain all advected quantities
C**** Currently this is on atmospheric grid

        call get_snow_ice_layer(
     &       SNOWI(I,J),MSI(I,J),HSI(:,I,J),SSI(:,I,J),
#ifdef TRACERS_WATER
     &       TRSI(:,:,I,J),TRSNOW,TRICE,
#endif 
     &       SNOWL,HSNOW,HICE,SICEg,TSNW,TSIL,MICE,.false.)

C-- MASS: ICE(LMI)
        MHS(1      :LMI    ,I,J) = MICE
C-- MASS: SNOW(2)
        MHS(LMI+1  :LMI+2  ,I,J) = SNOWL
C-- HEAT: ICE(LMI)
        MHS(3+LMI  :2+2*LMI,I,J) = HICE
C-- HEAT: SNOW(2)
        MHS(3+2*LMI:4+2*LMI,I,J) = HSNOW
C-- SALT: ICE(LMI)
        MHS(5+2*LMI:4+3*LMI,I,J) = SICEg
C-- SALT: SNOW (2, always =0)
        MHS(5+3*LMI:6+3*LMI,I,J) = 0.

#ifdef TRACERS_WATER
C**** add tracers to advected arrays
        DO ITR=1,NTM
          MHS(1+(3+ITR-1)*(LMI+2)    :LMI+(3+ITR-1)*(LMI+2),I,J)
     &        = TRICE(ITR,:)
          MHS(1+LMI+(3+ITR-1)*(LMI+2):(3+ITR)*(LMI+2)      ,I,J)
     &        = TRSNOW(ITR,:)
        ENDDO
#endif

      ENDDO ! i
      ENDDO ! j


      CALL HALO_UPDATE(grid, MSI)

C****
C**** Interpolate to obtain latlon-oriented ice velocities at
C**** cell edges, and transform to CS orientation.  ll2csint_lij
C**** fills any halo cells in its outputs.
C****
      do j=grid_icdyn%j_strt,grid_icdyn%j_stop
        do i=1,grid_icdyn%im_world
          uvll(1,i,j) = usi(i,j)*dts
          uvll(2,i,j) = vsi(i,j)*dts
        enddo
      enddo
      call ll2csint_lij(grid_icdyn,i2a_uc,uvll,uvllatuc,
     &     is_ll_vector=.true.)
      call ll2csint_lij(grid_icdyn,i2a_vc,uvll,uvllatvc,
     &     is_ll_vector=.true.)

      coastfac =
     &          1d-3 ! convert kg/m2 ice mass to ice thickness
     &         *1d-1 ! 10 cm/s speed for 1 m thickness difference over ~100 km
     &         *(real(im,kind=8)/90d0) ! scale by gridlength

      do j=j_0-1,j_1
      do i=i_0y,i_1y
        ull = uvllatvc(1,i,j+1)
        vll = uvllatvc(2,i,j+1)
        vdxdt(i,j) = (ull*ull2vcs(i,j+1)+vll*vll2vcs(i,j+1))
        if(connect(i,j)*connect(i,j+1).eq.0) then
          vdxdt(i,j)=0.
        elseif(connect(i,j)+connect(i,j+1).lt.30) then
          ! 
          if(focean(i,j).lt.focean(i,j+1)) then
            vdxdt(i,j) = vdxdt(i,j)
     &           + dts*min(+10d0,max(0d0,msi(i,j)-msi(i,j+1))*coastfac)
          else
            vdxdt(i,j) = vdxdt(i,j)
     &           + dts*max(-10d0,min(0d0,msi(i,j)-msi(i,j+1))*coastfac)
          endif
        endif
        vdxdt(i,j) = dlxsina(i,j+1)*vdxdt(i,j)
      enddo
      enddo
      do j=j_0,j_1
      do i=i_0-1,i_1
        ull = uvllatuc(1,i+1,j)
        vll = uvllatuc(2,i+1,j)
        udydt(i,j) = (ull*ull2ucs(i+1,j)+vll*vll2ucs(i+1,j))
        if(connect(i,j)*connect(i+1,j).eq.0) then
          udydt(i,j)=0.
        elseif(connect(i,j)+connect(i+1,j).lt.30) then
          if(focean(i,j).lt.focean(i+1,j)) then
            udydt(i,j) = udydt(i,j)
     &           + dts*min(+10d0,max(0d0,msi(i,j)-msi(i+1,j))*coastfac)
          else
            udydt(i,j) = udydt(i,j)
     &           + dts*max(-10d0,min(0d0,msi(i,j)-msi(i+1,j))*coastfac)
          endif
        endif
        udydt(i,j) = dlysina(i+1,j)*udydt(i,j)
      enddo
      enddo

c for now, no transport across cube edges
      if(i_0 .eq.  1) udydt(0 ,:) = 0.
      if(i_1 .eq. im) udydt(im,:) = 0.
      if(j_0 .eq.  1) vdxdt(:, 0) = 0.
      if(j_1 .eq. im) vdxdt(:,im) = 0.

C****
C**** Update halos of transported quantities
C****
      CALL HALO_UPDATE(grid, RSI)
      CALL HALO_UPDATE(grid, RSIX)
      CALL HALO_UPDATE(grid, RSIY)
      CALL HALO_UPDATE(grid, MHS, jdim=3)

C****
C**** Transport in the Y direction
C****
      j = j_0-1
      do i=i_0y,i_1y ! updating i-halo cells for subsequent x-sweep
        faw(i) = vdxdt(i,j)  ! should compute vdxdt here
        if(faw(i).le.0.) then
c**** sea ice velocity is southward at grid box edge
          cour = faw(i)*byaxyp(i,j+1)
          fao = faw(i)*focean(i,j+1)
          fasi(i)=fao*(rsi(i,j+1)-(1d0+cour)*rsiy(i,j+1))
          fxsi(i)=fao*rsix(i,j+1)
          fysi(i)=faw(i)*(cour*fao*rsiy(i,j+1)-3d0*fasi(i))
          fmsi(1:ntrice,i) = fasi(i)*mhs(1:ntrice,i,j+1)
        else
c**** sea ice velocity is northward at grid box edge
          cour = faw(i)*byaxyp(i,j)
          fao = faw(i)*focean(i,j)
          fasi(i)=fao*(rsi(i,j)+(1d0-cour)*rsiy(i,j))
          fxsi(i)=fao*rsix(i,j)
          fysi(i)=faw(i)*(cour*fao*rsiy(i,j)-3d0*fasi(i))
          fmsi(1:ntrice,i) = fasi(i)*mhs(1:ntrice,i,j)
        end if
        faw_jm1(i)  = faw(i)
        fasi_jm1(i) = fasi(i)
        fxsi_jm1(i) = fxsi(i)
        fysi_jm1(i) = fysi(i)
        fmsi_jm1(1:ntrice,i) = fmsi(1:ntrice,i)
      enddo

      do j=j_0,j_1
      do i=i_0y,i_1y ! updating i-halo cells for subsequent x-sweep

c _jm1 qtys are already zero. no need to re-zero.
        if(vdxdt(i,j-1).eq.0. .and. vdxdt(i,j).eq.0.) cycle

        faw(i) = vdxdt(i,j)  ! should compute vdxdt here
        if(faw(i).le.0.) then
c**** sea ice velocity is southward at grid box edge
          cour = faw(i)*byaxyp(i,j+1)
          fao = faw(i)*focean(i,j+1)
          fasi(i)=fao*(rsi(i,j+1)-(1d0+cour)*rsiy(i,j+1))
          fxsi(i)=fao*rsix(i,j+1)
          fysi(i)=faw(i)*(cour*fao*rsiy(i,j+1)-3d0*fasi(i))
          fmsi(1:ntrice,i) = fasi(i)*mhs(1:ntrice,i,j+1)
        else
c**** sea ice velocity is northward at grid box edge
          cour = faw(i)*byaxyp(i,j)
          fao = faw(i)*focean(i,j)
          fasi(i)=fao*(rsi(i,j)+(1d0-cour)*rsiy(i,j))
          fxsi(i)=fao*rsix(i,j)
          fysi(i)=faw(i)*(cour*fao*rsiy(i,j)-3d0*fasi(i))
          fmsi(1:ntrice,i) = fasi(i)*mhs(1:ntrice,i,j)
        end if
c accumulate transports
        aij(i,j,ij_mvsi)=aij(i,j,ij_mvsi)+sum(fmsi(1:lmi+2,i))
        aij(i,j,ij_hvsi)=aij(i,j,ij_hvsi)+sum(fmsi(3+lmi:4+2*lmi,i))
        aij(i,j,ij_svsi)=aij(i,j,ij_svsi)+sum(fmsi(5+2*lmi:6+3*lmi,i))
c#ifdef TRACERS_WATER
c        do itr=1,ntm
c           taijn(i,j,ticij_tvsi,itr)=taijn(i,j,ticij_tvsi,itr)+
c     &          sum(fmsi(1+(3+itr-1)*(lmi+2):(3+itr)*(lmi+2),i))
c         end do
c#endif

        if(faw_jm1(i).gt.0. .or. faw(i).lt.0.) then
! when there is inflow, use general-case formulas
          asi = rsi(i,j)*foa(i,j)
          do k=1,ntrice
            amsi(k) = asi*mhs(k,i,j) + (fmsi_jm1(k,i)-fmsi(k,i))
          enddo
          asi = asi + (fasi_jm1(i)-fasi(i))
          if(asi.le.foa(i,j)) then
            yrsi = (rsiy(i,j)*axyp(i,j)*foa(i,j)+
     &           (fysi_jm1(i)-fysi(i)) + 3d0*((faw_jm1(i)+
     &           faw(i))*asi-axyp(i,j)*(fasi_jm1(i)+fasi(i))))
     &           / (axyp(i,j) + (faw_jm1(i)-faw(i)))
            rsi(i,j)  = asi*byfoa(i,j)
            rsiy(i,j) = yrsi*byfoa(i,j)
            rsix(i,j) = rsix(i,j) + (fxsi_jm1(i)-fxsi(i))*byfoa(i,j)
            if(asi.gt.0) mhs(1:ntrice,i,j) = amsi(1:ntrice)/asi
          else

c**** sea ice crunches into itself and completely covers grid box
c**** move excess ice to lower layers, allow snow to pile up

            rsi(i,j)   = 1d0
            rsix(i,j)  = 0.
            rsiy(i,j)  = 0.
            do k=1,ntrice/(lmi+2)
              !**** k=1: mass, k=2: heat, k=3: salt
              mhs(1+(lmi+2)*(k-1),i,j) = amsi(1+(lmi+2)*(k-1)) / asi ! 4 ice layers
              mhs(2+(lmi+2)*(k-1),i,j) = amsi(2+(lmi+2)*(k-1)) / asi
              dmhsi = sum(amsi(1+(lmi+2)*(k-1):4+(lmi+2)*(k-1)))
     &             *(byfoa(i,j) -1d0/ asi)
              mhs(3+(lmi+2)*(k-1),i,j) = amsi(3+(lmi+2)*(k-1)) / asi +
     &             xsi(3)*dmhsi
              mhs(4+(lmi+2)*(k-1),i,j) = amsi(4+(lmi+2)*(k-1)) / asi +
     &             xsi(4)*dmhsi
              mhs(5+(lmi+2)*(k-1),i,j) =
     &             amsi(5+(lmi+2)*(k-1)) * byfoa(i,j) ! 2 snow layers
              mhs(6+(lmi+2)*(k-1),i,j) =
     &             amsi(6+(lmi+2)*(k-1)) * byfoa(i,j)
            enddo

          endif
        else
! when there is only outflow, use simpler formulas.
! why is mhs not updated here.
          rsi(i,j)  =  rsi(i,j) + (fasi_jm1(i)-fasi(i))*byfoa(i,j)
          cnew = 1d0+(faw_jm1(i)*focean(i,j-1)
     &               -faw(i)*focean(i,j))*byfoa(i,j)
          rsix(i,j) = rsix(i,j)*cnew
          rsiy(i,j) = rsiy(i,j)*cnew**2
        endif

        faw_jm1(i)  = faw(i)
        fasi_jm1(i) = fasi(i)
        fxsi_jm1(i) = fxsi(i)
        fysi_jm1(i) = fysi(i)
        fmsi_jm1(1:ntrice,i) = fmsi(1:ntrice,i)

C**** Limit RSIX and RSIY so that sea ice is positive at the edges.
        rsi(i,j) = max(0d0,rsi(i,j))
        if(rsi(i,j)-rsix(i,j).lt.0.)  rsix(i,j) =    rsi(i,j)
        if(rsi(i,j)+rsix(i,j).lt.0.)  rsix(i,j) =   -rsi(i,j)
        if(rsi(i,j)-rsix(i,j).gt.1d0) rsix(i,j) =    rsi(i,j)-1d0
        if(rsi(i,j)+rsix(i,j).gt.1d0) rsix(i,j) =1d0-rsi(i,j)
        if(rsi(i,j)-rsiy(i,j).lt.0.)  rsiy(i,j) =    rsi(i,j)
        if(rsi(i,j)+rsiy(i,j).lt.0.)  rsiy(i,j) =   -rsi(i,j)
        if(rsi(i,j)-rsiy(i,j).gt.1d0) rsiy(i,j) =    rsi(i,j)-1d0
        if(rsi(i,j)+rsiy(i,j).gt.1d0) rsiy(i,j) =1d0-rsi(i,j)

      enddo ! i
      enddo ! j

C****
C**** Transport in the X direction
C****

      do j=j_0,j_1

      i=i_0-1
      faw(i) = udydt(i,j)  ! should compute udydt here
      if(faw(i).le.0.) then
c**** sea ice velocity is westward at grid box edge
        cour = faw(i)*byaxyp(i+1,j)
        fao = faw(i)*focean(i+1,j)
        fasi(i)=fao*(rsi(i+1,j)-(1d0+cour)*rsix(i+1,j))
        fxsi(i)=faw(i)*(cour*fao*rsix(i+1,j)-3d0*fasi(i))
        fysi(i)=fao*rsiy(i+1,j)
        fmsi(1:ntrice,i) = fasi(i)*mhs(1:ntrice,i+1,j)
      else
c**** sea ice velocity is eastward at grid box edge
        cour = faw(i)*byaxyp(i,j)
        fao = faw(i)*focean(i,j)
        fasi(i)=fao*(rsi(i,j)+(1d0-cour)*rsix(i,j))
        fxsi(i)=faw(i)*(cour*fao*rsix(i,j)-3d0*fasi(i))
        fysi(i)=fao*rsiy(i,j)
        fmsi(1:ntrice,i) = fasi(i)*mhs(1:ntrice,i,j)
      endif
      faw_im1 = faw(i)
      fasi_im1 = fasi(i)
      fxsi_im1 = fxsi(i)
      fysi_im1 = fysi(i)
      fmsi_im1(1:ntrice) = fmsi(1:ntrice,i)

      do i=i_0,i_1

c _im1 qtys are already zero. no need to re-zero.
        if(udydt(i-1,j).eq.0. .and. udydt(i,j).eq.0.) cycle

        faw(i) = udydt(i,j) ! should compute udydt here
        if(faw(i).le.0.) then
c**** sea ice velocity is westward at grid box edge
          cour = faw(i)*byaxyp(i+1,j)
          fao = faw(i)*focean(i+1,j)
          fasi(i)=fao*(rsi(i+1,j)-(1d0+cour)*rsix(i+1,j))
          fxsi(i)=faw(i)*(cour*fao*rsix(i+1,j)-3d0*fasi(i))
          fysi(i)=fao*rsiy(i+1,j)
          fmsi(1:ntrice,i) = fasi(i)*mhs(1:ntrice,i+1,j)
        else
c**** sea ice velocity is eastward at grid box edge
          cour = faw(i)*byaxyp(i,j)
          fao = faw(i)*focean(i,j)
          fasi(i)=fao*(rsi(i,j)+(1d0-cour)*rsix(i,j))
          fxsi(i)=faw(i)*(cour*fao*rsix(i,j)-3d0*fasi(i))
          fysi(i)=fao*rsiy(i,j)
          fmsi(1:ntrice,i) = fasi(i)*mhs(1:ntrice,i,j)
        endif
c accumulate transports
        aij(i,j,ij_musi)=aij(i,j,ij_musi)+sum(fmsi(1:lmi+2,i))
        aij(i,j,ij_husi)=aij(i,j,ij_husi)+sum(fmsi(3+lmi:4+2*lmi,i))
        aij(i,j,ij_susi)=aij(i,j,ij_susi)+sum(fmsi(5+2*lmi:6+3*lmi,i))
c#ifdef TRACERS_WATER
c        do itr=1,ntm
c          taijn(i,j,ticij_tusi,itr)=taijn(i,j,ticij_tusi,itr)+
c     &         sum(fmsi(1+(3+itr-1)*(lmi+2):(3+itr)*(lmi+2),i))
c        enddo
c#endif

        if(faw_im1.gt.0. .or. faw(i).lt.0.) then
! when there is inflow, use general-case formulas
          asi = rsi(i,j)*foa(i,j)
          do k=1,ntrice
            amsi(k) = asi*mhs(k,i,j) + (fmsi_im1(k)-fmsi(k,i))
          enddo
          asi = asi + (fasi_im1-fasi(i))
          if(asi.le.foa(i,j)) then
            xrsi = (rsix(i,j)*axyp(i,j)*foa(i,j)+
     &        (fxsi_im1-fxsi(i)) + 3d0*((faw_im1+faw(i))*asi-
     &           axyp(i,j)*(fasi_im1+fasi(i))))
     &           / (axyp(i,j) + (faw_im1-faw(i)))
            rsi(i,j)  = asi*byfoa(i,j)
            rsix(i,j) = xrsi*byfoa(i,j)
            rsiy(i,j) = rsiy(i,j) + (fysi_im1-fysi(i))*byfoa(i,j)
            if (asi.gt.0) mhs(1:ntrice,i,j) = amsi(1:ntrice)/asi
          else

c**** sea ice crunches into itself and completely covers grid box
c**** move excess ice to lower layers, allow snow to pile up
            rsi(i,j)   = 1d0
            rsix(i,j)  = 0.
            rsiy(i,j)  = 0.
            do k=1,ntrice/(lmi+2)
              ! k=1: mass, k=2: heat, k=3: salt
              mhs(1+(lmi+2)*(k-1),i,j) = amsi(1+(lmi+2)*(k-1)) / asi ! 4 ice layers
              mhs(2+(lmi+2)*(k-1),i,j) = amsi(2+(lmi+2)*(k-1)) / asi
              dmhsi = sum(amsi(1+(lmi+2)*(k-1):4+(lmi+2)*(k-1)))
     &             *(byfoa(i,j) -1d0/ asi)
              mhs(3+(lmi+2)*(k-1),i,j) = amsi(3+(lmi+2)*(k-1)) / asi +
     &             xsi(3)*dmhsi
              mhs(4+(lmi+2)*(k-1),i,j) = amsi(4+(lmi+2)*(k-1)) / asi +
     &             xsi(4)*dmhsi
              mhs(5+(lmi+2)*(k-1),i,j) =
     &             amsi(5+(lmi+2)*(k-1)) * byfoa(i,j) ! 2 snow layers
              mhs(6+(lmi+2)*(k-1),i,j) =
     &             amsi(6+(lmi+2)*(k-1)) * byfoa(i,j)
            enddo

          endif
        else
! when there is only outflow, use simpler formulas
! why is mhs not updated here.
          rsi(i,j)  =  rsi(i,j) + (fasi_im1-fasi(i))*byfoa(i,j)
          cnew = 1d0+(faw_im1*focean(i-1,j)
     &               -faw(i )*focean(i  ,j))*byfoa(i,j)
          rsix(i,j) = rsix(i,j)*cnew**2
          rsiy(i,j) = rsiy(i,j)*cnew
        endif
        faw_im1  = faw(i)
        fasi_im1 = fasi(i)
        fxsi_im1 = fxsi(i)
        fysi_im1 = fysi(i)
        fmsi_im1(1:ntrice) = fmsi(1:ntrice,i)
        
C**** Limit RSIX and RSIY so that sea ice is positive at the edges.
c why is it necessary to do this after the advection?
        rsi(i,j) = max(0d0,rsi(i,j))
        if(rsi(i,j)-rsix(i,j).lt.0.)  rsix(i,j) =    rsi(i,j)
        if(rsi(i,j)+rsix(i,j).lt.0.)  rsix(i,j) =   -rsi(i,j)
        if(rsi(i,j)-rsix(i,j).gt.1d0) rsix(i,j) =    rsi(i,j)-1d0
        if(rsi(i,j)+rsix(i,j).gt.1d0) rsix(i,j) =1d0-rsi(i,j)
        if(rsi(i,j)-rsiy(i,j).lt.0.)  rsiy(i,j) =    rsi(i,j)
        if(rsi(i,j)+rsiy(i,j).lt.0.)  rsiy(i,j) =   -rsi(i,j)
        if(rsi(i,j)-rsiy(i,j).gt.1d0) rsiy(i,j) =    rsi(i,j)-1d0
        if(rsi(i,j)+rsiy(i,j).gt.1d0) rsiy(i,j) =1d0-rsi(i,j)

      enddo ! i
      enddo ! j

      IF (KOCEAN.ge.1) THEN ! full ocean calculation, adjust sea ice
C**** set global variables from local array
C**** Currently on atmospheric grid, so no interpolation necessary
        DO J=J_0,J_1
          DO I=I_0,I_1
            IF (FOCEAN(I,J).gt.0) THEN

C**** Fresh water sea ice mass convergence (needed for qflux model)
            MSICNV(I,J) = RSI(I,J)*(SUM(MHS(1:LMI+2,I,J))
     &           -SUM(MHS(5+2*LMI:4+3*LMI,I,J)))
     &                  - RSISAVE(I,J)*(ACE1I+SNOWI(I,J)
     &           +MSI(I,J)-SUM(SSI(1:LMI,I,J)))

C-- MASS: ICE(LMI)
            MICE  = MHS(1      :LMI    ,I,J)
C-- MASS: SNOW(2)
            SNOWL = MHS(LMI+1  :LMI+2  ,I,J)
C-- HEAT: ICE(LMI)
            HICE  = MHS(3+LMI  :2+2*LMI,I,J)
C-- HEAT: SNOW(2)
            HSNOW = MHS(3+2*LMI:4+2*LMI,I,J)
C-- SALT: ICE(LMI)
            SICEg = MHS(5+2*LMI:4+3*LMI,I,J)
C-- SALT: SNOW (2, always =0)
c         MHS(5+3*LMI:6+3*LMI,I,J) = 0.

#ifdef TRACERS_WATER
            DO ITR=1,NTM
              TRICE(ITR,:) = 
     &            MHS(1+(3+ITR-1)*(LMI+2)    :LMI+(3+ITR-1)*(LMI+2),I,J)
              TRSNOW(ITR,:) =
     &            MHS(1+LMI+(3+ITR-1)*(LMI+2):(3+ITR)*(LMI+2)      ,I,J)
            ENDDO
#endif

C-- Mass fluxes required to keep first layer ice = ACE1I
            FMSI2=MICE(1)+MICE(2)-ACE1I
C**** relayer lower ice layers
            call relayer(FMSI2,MICE,HICE,SICEg
#ifdef TRACERS_WATER
     &       ,TRICE
#endif
     &       )

C**** relayer upper two layers
            call relayer_12(HSNOW,HICE,SICEg,MICE,SNOWL
#ifdef TRACERS_WATER
     &       ,TRSNOW,TRICE
#endif 
     &       )

C**** reconstitute upper snow and ice layers
            call set_snow_ice_layer(HSNOW,HICE,SICEg,MICE,SNOWL,
#ifdef TRACERS_WATER
     &       TRSNOW,TRICE,TRSI(:,:,I,J),
#endif 
     &       SNOWI(I,J),MSI1,MSI(I,J),HSI(:,I,J),SSI(:,I,J))

            FWSIM(I,J)=RSI(I,J)*(ACE1I+SNOWI(I,J)+MSI(I,J)-
     &           SUM(SSI(1:LMI,I,J)))
            END IF
          END DO
        END DO
C**** Set atmospheric arrays
        DO J=J_0,J_1
          DO I=I_0,I_1
            IF (FOCEAN(I,J).gt.0) THEN
C**** set total atmopsheric pressure anomaly in case needed by ocean
              APRESS(I,J) = 100.*(P(I,J)+PTOP-1013.25d0)+RSI(I,J)
     *             *(SNOWI(I,J)+ACE1I+MSI(I,J))*GRAV
              GTEMP(1:2,2,I,J)=((HSI(1:2,I,J)-SSI(1:2,I,J)*LHM)/
     *             (XSI(1:2)*(SNOWI(I,J)+ACE1I))+LHM)*BYSHI
#ifdef TRACERS_WATER
              GTRACER(:,2,I,J)=TRSI(:,1,I,J)/(MHS(1,I,J)+MHS(LMI+1,I,J)
     *             -SSI(1,I,J))
#endif
            END IF
          END DO
        END DO
      ELSE          ! fixed SST case, save implied heat convergence
        DO J=J_0,J_1
          DO I=I_0,I_1
            IF (FOCEAN(I,J).gt.0) THEN
             OA(I,J,13)=OA(I,J,13)+(RSI(I,J)*SUM(MHS(3+LMI:2+2*LMI,I,J))
     *             -RSISAVE(I,J)*SUM(HSI(1:LMI,I,J)))
C**** reset sea ice concentration
              RSI(I,J)=RSISAVE(I,J)
            END IF
          END DO
        END DO
      END IF
C****
      RETURN
      END SUBROUTINE ADVSI
