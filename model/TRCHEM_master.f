#include "rundeck_opts.h"
      SUBROUTINE masterchem
!@sum masterchem main chemistry routine
!@auth Drew Shindell (modelEifications by Greg Faluvegi)
!@ver  1.0 (based on masterchem000_M23p)   
!@calls fastj2_drv,checktracer,Crates,Oxfam,HOxfam,NOxfam,chemstep
C
C IF ALTERING THIS ROUTINE, PLEASE SEE THE WARNING ABOUT THE CHANGEL
C VARIABLE IN THE STRATOSPHERIC OVERWRITE SECTION.
c
C**** GLOBAL parameters and variables:
c
!!    use precision_mod, only : reduce_precision 
      USE PARAM, only : get_param, is_set_param
      USE SOMTQ_COM, only   : qmom
      USE DOMAIN_DECOMP_1D, only : PACK_DATA ! for DU_O3
      USE DOMAIN_DECOMP_ATM,only: GRID,GET,AM_I_ROOT,
     &                        GLOBALSUM,GLOBALMAX,
     &                        write_parallel,writet8_column,
     &                        writet_parallel
      USE MODEL_COM, only   : Q,JDAY,IM,JM,sig,ptop,psf,ls1,JYEAR,
     &                        COUPLED_CHEM,JHOUR, itime, itimeI, itime0
      USE CONSTANT, only    : radian,gasc,mair,mb2kg,pi,avog,rgas,pO2,
     &                        bygrav,lhe,undef,teeny
      USE DYNAMICS, only    : pedn,LTROPO
      USE FILEMANAGER, only : openunit,closeunit,nameunit
      USE RAD_COM, only     : COSZ1,alb,rcloudfj=>rcld,
     &                        rad_to_chem,chem_tracer_save,H2ObyCH4,
     &                        SRDN,rad_to_file,ghg_yr,clim_interact_chem
#if (defined SHINDELL_STRAT_EXTRA) && (defined ACCMIP_LIKE_DIAGS)
     &                        ,stratO3_tracer_save
#endif
      USE GEOM, only        : BYAXYP, AXYP, LAT2D_DG, IMAXJ, LAT2D
      USE FLUXES, only      : tr3Dsource
      USE TRACER_COM, only  : n_Ox,n_NOx,n_N2O5,n_HNO3,n_H2O2,n_CH3OOH,
     &                        n_HCHO,n_HO2NO2,n_CO,n_CH4,n_PAN,
     &                        n_Isoprene,n_AlkylNit,n_Alkenes,n_stratOx,
     &                        n_Terpenes,
     &                        n_Paraffin,ntm_chem,n_DMS,n_MSA,n_SO2,
     &                        tr_wd_type,nWater,trm,trmom,
     &                        n_isopp1g,n_isopp1a,n_isopp2g,n_isopp2a,
     &                        n_apinp1g,n_apinp1a,n_apinp2g,n_apinp2a,
#ifdef TRACERS_AEROSOLS_Koch
     &                        n_SO4,
#endif
#ifdef TRACERS_AMP
     &                        n_M_AKK_SU,n_M_ACC_SU,n_M_DD1_SU,
     &                        n_M_DS1_SU,n_M_DD2_SU,n_M_DS2_SU,
     &                        n_M_SSA_SU,n_M_OCC_SU,n_M_BC1_SU,
     &                        n_M_BC2_SU,n_M_BC3_SU,n_M_DBC_SU,
     &                        n_M_BOC_SU,n_M_BCS_SU,n_M_MXX_SU,
#endif
#ifdef TRACERS_TOMAS
     &                        n_ASO4,nbins,  
#endif
     &                        n_H2O2_s,oh_live,no3_live,
     &                        nChemistry,nOverwrite,rsulf1,rsulf2,
     &                        rsulf3,rsulf4,TR_MM,trname
     *                        ,mass2vol,vol2mass
#ifdef SHINDELL_STRAT_CHEM
     &                        ,n_HBr,n_HOCl,n_HCl,n_ClONO2,n_ClOx,
     &                        n_BrOx,n_BrONO2,n_CFC,n_N2O,n_HOBR
#endif
#ifdef TRACERS_HETCHEM
     &                        ,krate,n_N_d1,n_N_d2,n_N_d3
#endif
#ifdef TOMAS_HETCHEM
     &                        ,krate,n_NO3p,nOther
#endif
#ifdef INTERACTIVE_WETLANDS_CH4
      USE TRACER_SOURCES, only: avg_model,n__sw
#endif
      USE TRDIAG_COM, only    : taijs=>taijs_loc,taijls=>taijls_loc
     &     ,ijlt_NO3,jls_COp,jls_COd,jls_Oxp,jls_N2O5sulf,jls_O3vmr
     &     ,jls_Oxd,jls_OxpT,jls_OxdT,ijs_NO2_1030,ijs_NO2_1030c
     &     ,ijlt_COp,ijlt_COd,ijlt_Oxd,ijlt_Oxp,ijlt_pO1D,ijs_O3mass
     &     ,ijlt_pOH,ijlt_OxpHO2,ijlt_OxpCH3O2,ijlt_OxlHO2,ijlt_OxlALK
     &     ,ijlt_OxlOH,ijs_NO2_1330,ijs_NO2_1330c,ijlt_NO2vmr,ijlt_NOvmr
     &     ,ijlt_JO1D,ijlt_JNO2,ijlt_JH2O2,ijlt_O3ppbv,ijlt_O3cmatm
#ifdef SHINDELL_STRAT_CHEM
     &     ,jls_ClOcon,jls_H2Ocon
#endif
      USE TRCHEM_Shindell_COM
#ifdef TRACERS_AEROSOLS_SOA
      USE TRACERS_SOA, only: soa_aerosolphase,voc2nox,soa_apart,
     &                       whichsoa,apartmolar,LM_soa
#endif  /* TRACERS_AEROSOLS_SOA */
      use zonalmean_mod, only : zonalmean_ij2ij

#ifdef CACHED_SUBDD
      use subdd_mod, only : subdd_groups,subdd_type,subdd_ngroups
     &     ,inc_subdd,find_groups, LmaxSUBDD
#endif
      use photolysis, only: fastj2_drv,o3_fastj
     &                     ,sza,szamax,zj,jppj,jpnl,sf3_fact,sf2_fact

      IMPLICIT NONE

#ifdef CACHED_SUBDD
      integer :: igrp,ngroups,grpids(subdd_ngroups),k
      type(subdd_type), pointer :: subdd
      real*8, dimension(grid%i_strt_halo:grid%i_stop_halo,
     &                  grid%j_strt_halo:grid%j_stop_halo,
     &                  LM) :: mrno,mrno2
#endif
C**** Local parameters and variables and arguments:
!@param by35 1/35 used for spherical geometry constant
      REAL*8, PARAMETER  :: by35=1.d0/35.d0
      REAL*8, PARAMETER  :: bymair = 1.d0/mair
!@var FASTJ_PFACT temp factor for vertical pressure-weighting
!@var FACT1,2,3 temp variable for strat overwrite
!@var bydtsrc reciprocal of the timestep dtsrc
!@var local logical for error checking 
!@var byam75 the reciprocal air mass near 75 hPa level
!@var average tropospheric ch4 value near 569 hPa level
!@var PRES2 local nominal pressure for verticle interpolations
!@var thick thickness of each layer in various units
!@var bythick recipricol thickness of each layer (1/m) saved on
!@+ model layers. Should replace thick with 1/this too?
!@var photO2_glob for O2->Ox diagnostic
!@var ClOx_old total ClOx at start of chemical timestep
!@var ClTOT total chlorine in all forms (reactive and reservoir)
!@var BrTOT total bromine in all forms (reactive and reservoir)
!@var colmO2, colmO3 are overhead oxygen and ozone columns
!@var CH4FACT, r179 for setting CH4 ICs and strat distribution
!@var changeClONO2,changeClOx,changeHOCl,changeHCl nighttime changes
!@+ also changeBrOx,changeBrONO2,changeBrOx2,changeHBr
!@var changehetClONO2 nighttime het change in ClONO2 (on sulfate)
!@var chgHT3,chgHT4,chgHT5 reaction rates for het rxns on pscs
!@var rmrClOx,rmrBrOx dummy vars with mixing ratios of halogens
!@var rmv dummy variable for halogne removal in trop vs height
!@var changeL 2D array holds the local change due to chem or strat 
!@+   overwrite until adding to tr3Dsource (replaces 4D "change")
!@var PIfact strat-overwrite scaling
!@var pfactor to convert units on species chemical changes
!@var bypfactor to convert units on species chemical changes
!@var dNO3,gwprodHNO3,gprodHNO3,gwprodN2O5,changeAldehyde,
!@+   changeAlkenes,changeIsoprene,changeHCHO,changeAlkylNit,
!@+   changeTerpenes,changeisopp1g,changeisopp2g,changeapinp1g,
!@+   changeapinp2g,changeHNO3,changeNOx,changeN2O5,wprodHCHO
!@+   working variables to calculate nighttime chemistry changes
!@var rlossN,rprodN,ratioN variables for nitrogen conservation
!@var I,J,L,N,igas,inss,LL,Lqq,JJ,L2,n2 dummy loop variables
!@var avgTT_CH4 Itime avg CH4 # density at LTROPO between 20N and 20S
!@var avgTT_H2O Itime avg H2O # density at LTROPO between 20N and 20S
!@var countTT # of points between 20N and 20S on LTROPO plane
!@var aero yes(1) or no(0) tag of non-zero rkext from Crates
!@var maxl chosen tropopause 0=LTROPO(I,J), 1=LS1-1
!@var sumOx for summing regional Ox tracers
!@var bysumOx reciprocal of sum of regional Ox tracers
!@var maxPSC a limit placed on some PSC reactions to prevent sudden overflows
      REAL*8, DIMENSION(LM,NTM) :: changeL
      REAL*8, DIMENSION(NTM)    :: PIfact
      REAL*8, DIMENSION(LM)     :: PRES2,rh,bythick
      REAL*8 :: tempChangeNOx
      REAL*8 :: FACT1,FACT2,FACT3,FACT4,FACT5,FACT6,FACT7,fact_so4,
     &  FASTJ_PFACT,bydtsrc,byam75,byavog,CH4FACT,r179,rlossN,maxPSC,
     &  rprodN,ratioN,pfactor,bypfactor,gwprodHNO3,gprodHNO3,
     &  gwprodN2O5,wprod_sulf,wprodCO,dNO3,wprodHCHO,prod_sulf,
     &  RVELN2O5,changeAldehyde,changeAlkenes,changeAlkylNit,
     &  changeIsoprene,changeHCHO,changeHNO3,changeNOx,changeN2O5,
     &  rNO3plusNO2,rN2O5decomp,rHCHOplusNO3,rAldplusNO3,rIsopplusNO3,
     &  rClOplusNO2,rDMSplusNO3,wlossN2O5,wlossNOx,
     &  wprodNOx,rlossNO3,rprodNO3,wprodN2O5,pNO3temp,rBrOplusNO2,
     &  wlossClONO2,wprodClONO2,changeN,sphericalCorrectionReg1,
     &  sphericalCorrectionReg2,sphericalCorrectionReg3,
     &  sphericalCorrectionReg4,
#ifdef TRACERS_TERP
     &  changeTerpenes,rTerpplusNO3,
#endif  /* TRACERS_TERP */
#ifdef TRACERS_AEROSOLS_SOA
     & changeisopp1g,changeisopp2g,
#ifdef TRACERS_TERP
     & changeapinp1g,changeapinp2g,
#endif  /* TRACERS_TERP */
#endif  /* TRACERS_AEROSOLS_SOA */
     &  changeOx,fraQ,CH4_569,count_569,thick, changeCO,tempChangeOx
#ifdef TRACERS_HETCHEM
      REAL*8 :: changeN_d1,changeN_d2,changeN_d3
#endif
#ifdef TOMAS_HETCHEM
      REAL*8 :: changeNO3p
#endif
#ifdef INTERACTIVE_WETLANDS_CH4
      REAL*8 :: temp_SW
#endif
#ifdef SHINDELL_STRAT_CHEM
      REAL*8, DIMENSION(LM)     :: ClOx_old  
      REAL*8 :: BRTOT,CLTOT,colmO2,colmO3,changeClONO2,changeClOx,
     & changeHOCl,changeHCl,changehetClONO2,chgHT3,
     & chgHT4,chgHT5,rmrClOx,rmrBrOx,rmv,rmrOx,avgTT_H2O,avgTT_CH4,
     & countTT,bHNO3,mHNO3,HNO3_thresh,Ttemp,changeBrOx,changeBrONO2,
     & changeBrOx2,changeHBr
      INTEGER, DIMENSION(LM)    :: aero
#endif
      INTEGER                   :: igas,LL,I,J,L,N,inss,Lqq,L2,n2,
     &                          ierr,ierr_loc,Jqq,Iqq,maxl,iu,
     & itemp_iter,ih1330e,ih1030e,ih1030,ih1330,m,istep,index1,index2
      LOGICAL                   :: error, jay
      CHARACTER*4               :: ghg_name
      CHARACTER*80              :: ghg_file,title
      character(len=300)        :: out_line

      real*8 :: ghg_out(LM,GRID%I_STRT_HALO:GRID%I_STOP_HALO,
     &                     GRID%J_STRT_HALO:GRID%J_STOP_HALO)

      real*8 :: ss27x2,ss27x2_c,OHpptv,HO2pptv,ObyO3,NO2byNO,ClbyClO
      real*8, dimension(JM)         :: DU_O3_glob
#ifdef SHINDELL_STRAT_CHEM
c      real*8, dimension(JM,LM)      :: photO2_glob
#ifdef TRACERS_TERP
      integer, parameter :: iN2O5plusH2O=109,iNO3plusNO2=103,
     &                      iN2O5decomp=96,iClOplusNO2=107,
     &                      iClONO2plusH2O=110,iClONO2plusHCl=111,
     &                      iHOClplusHCl=112,iN2O5plusHCl=113,
     &                      iTerpenesO3=93,iTerpenesNO3=94
#else
      integer, parameter :: iN2O5plusH2O=106,iNO3plusNO2=100,
     &                      iN2O5decomp=93,iClOplusNO2=104,
     &                      iClONO2plusH2O=107,iClONO2plusHCl=108,
     &                      iHOClplusHCl=109,iN2O5plusHCl=110
#endif  /* TRACERS_TERP */
#else
#ifdef TRACERS_TERP
      integer, parameter :: iNO3plusNO2=56,iN2O5decomp=50,
     &                      iTerpenesO3=47,iTerpenesNO3=48
#else
      integer, parameter :: iNO3plusNO2=53,iN2O5decomp=47
#endif  /* TRACERS_TERP */
#endif
      real*8, dimension(GRID%I_STRT_HALO:GRID%I_STOP_HALO,
     &                  GRID%J_STRT_HALO:GRID%J_STOP_HALO) ::
     &     avgTT_CH4_part,avgTT_H2O_part,countTT_part,
     &     CH4_569_part,count_569_part

      real*8, dimension(GRID%I_STRT_HALO:GRID%I_STOP_HALO,
     &                  GRID%J_STRT_HALO:GRID%J_STOP_HALO) ::
     &     surfIsop,zonalIsop
      
      INTEGER :: J_0, J_1, J_0S, J_1S, J_0H, J_1H, I_0, I_1
      integer :: initial_GHG_setup
      LOGICAL :: HAVE_SOUTH_POLE, HAVE_NORTH_POLE     
      real*8 :: qsat ! this is a function in UTILDBL.f
#ifdef TRACERS_AEROSOLS_SOA
      real*8 voc2nox_denom
#endif  /* TRACERS_AEROSOLS_SOA */
#ifdef TRACERS_TOMAS
      integer :: nb
#endif

      CALL GET(grid, J_STRT    =J_0,  J_STOP    =J_1,
     &               I_STRT    =I_0,  I_STOP    =I_1,
     &               J_STRT_SKP=J_0S, J_STOP_SKP=J_1S,
     &               J_STRT_HALO=J_0H, J_STOP_HALO=J_1H,
     &               HAVE_SOUTH_POLE = have_south_pole,
     &               HAVE_NORTH_POLE = have_north_pole)
      
      byavog = 1.d0/avog

! calculate what longitudes to accumulate for 10:30am/1:30pm NO2 diags:
! Um... Does use of Jhour here assume starting the model at midnight?
      istep = NINT(real(IM)/24.) ! number of boxes per hour
      ! ih1030/1330 are westmost I index that hour (careful: int arith.)
      ih1030 = istep*(10-Jhour)+IM/2+NINT(real(istep)/2.)-(istep-1)/2
      ih1330 = istep*(13-Jhour)+IM/2+NINT(real(istep)/2.)-(istep-1)/2
      if(ih1030 < 0) ih1030 = IM+ih1030  
      if(ih1330 < 0) ih1330 = IM+ih1330  
      if(ih1030 > IM) ih1030 = ih1030-IM
      if(ih1330 > IM) ih1330 = ih1330-IM
      ih1330e=ih1330+istep-1
      ih1030e=ih1030+istep-1
      if(ih1030e < 0) ih1030e = IM+ih1030e  
      if(ih1330e < 0) ih1330e = IM+ih1330e  
      if(ih1030e > IM) ih1030e = ih1030e-IM
      if(ih1330e > IM) ih1330e = ih1330e-IM

! meanwhile, initialize the instantaneous SUBDD of NO2 column:
! I think it will be overwritten for all i,j, so this can be
! a temporary check:
      save_NO2column(I_0:I_1,J_0:J_1)=undef
 
      if (is_set_param('initial_GHG_setup')) then
        call get_param('initial_GHG_setup', initial_GHG_setup)
        if (initial_GHG_setup == 1 .and. itime == itimeI) then
C--------special section for ghg runs ---------
          write(out_line,*)'Warning: INITIAL_GHG_SETUP is on!'
          call write_parallel(trim(out_line))
          if(use_rad_ch4>0 .or. use_rad_n2o>0 .or. use_rad_cfc>0)then
            rad_to_file(1,:,I_0:I_1,J_0:J_1)=
     &           rad_to_chem(1,:,I_0:I_1,J_0:J_1)
            rad_to_file(2,:,I_0:I_1,J_0:J_1)=
     &           rad_to_chem(2,:,I_0:I_1,J_0:J_1)
            do j=J_0,J_1
              do i=I_0,I_1
                rad_to_file(3,:,i,j)=rad_to_chem(3,:,i,j)*2.69e20*
     &               byavog*axyp(i,j)*tr_mm(n_N2O) ! i.e. in trm units now!
                rad_to_file(4,:,i,j)=rad_to_chem(4,:,i,j)*2.69e20*
     &               byavog*axyp(i,j)*tr_mm(n_CH4) ! i.e. in trm units now!
                rad_to_file(5,:,i,j)=rad_to_chem(5,:,i,j)*2.69e20*
     &               byavog*axyp(i,j)*tr_mm(n_CFC)*fact_CFC ! i.e. in trm units now!
              enddo
            enddo 
            if(ghg_yr/=0)then; write(ghg_name,'(I4.4)')ghg_yr
            else; write(ghg_name,'(I4.4)')jyear; endif
            ghg_file='GHG_IC_'//ghg_name
            call openunit(ghg_file,iu,.true.,.false.)
            do m=1,5
             ghg_out(:,I_0:I_1,J_0:J_1)=rad_to_file(m,:,I_0:I_1,J_0:J_1)
             CALL WRITET8_COLUMN(grid,iu,NAMEUNIT(iu),GHG_OUT,ghg_file)
            enddo
            call closeunit(iu)          
            if(AM_I_ROOT( ))then
              write(6,*)'Kludge in masterchem to output inital'
              write(6,*)'conditions for ghgs.: ',trim(ghg_file)
              write(6,*)'First time step has used default values.'
              write(6,*)'If you wish to produce a correct'
              write(6,*)'first time step, then redo setup for this'
              write(6,*)' rundeck with initial_GHG_setup set to 0.'
              write(6,*)'Address questions to G. Faluvegi or T. Clune.'
              write(6,*)'Thanks.'
            endif
          end if
        end if
      end if

C Some INITIALIZATIONS :
      byavog  = 1.d0/avog
      bydtsrc = 1.d0/dtsrc
      BYFJM   = 1.d0/real(JM)
      PRES2(:)= SIG(:)*(PSF-PTOP)+PTOP
#ifdef SHINDELL_STRAT_CHEM
      if(H2ObyCH4 /= 0. .and. clim_interact_chem > 0)                
     &call stop_model('H2ObyCH4.ne.0 .and. clim_interact_chem > 0',13)
#endif
#ifdef INTERACTIVE_WETLANDS_CH4
C Not really related to chemistry, but convenient place to update
C running-averages for interactive wetlands CH4:
      do J=J_0,J_1; do I=I_0,IMAXJ(J)
        temp_SW=ALB(I,J,1)*(SRDN(I,J)+1.d-20)*COSZ1(I,J)
        call running_average(temp_SW,I,J,1.d0,n__sw)
      end do      ; end do
#endif

C Calculation of gas phase reaction rates for sulfur chemistry:
C Now called from tracer_3Dsource
c      CALL GET_SULF_GAS_RATES
      
#if (defined TRACERS_HETCHEM) || (defined TOMAS_HETCHEM)
c Calculation of removal rates on dust surfaces:
      CALL HETCDUST
#endif

#ifndef SHINDELL_STRAT_CHEM 
      if(which_trop==0 .or. (JPNL<LS1-1.and.which_trop==1))
     &call stop_model('Are you sure JPNL is right?',255)
#else
      ! Note to self: move all Itime==ItimeI things to TRCHEM_init.f
      if(Itime==ItimeI)then
        if(use_rad_n2o > 0)then
          write(out_line,*) 'Warning:use_rad_n2o overrides PIfact_N2O'
          call write_parallel(trim(out_line))
        endif
        if(use_rad_cfc > 0)then
          write(out_line,*) 'Warning:use_rad_cfc overrides PIfact_CFC'
          call write_parallel(trim(out_line))
        endif
      endif
#endif

c Set "chemical time step". Really this is a method of applying only
c a fraction of the chemistry change to the tracer mass for the first
c 30 hours.  That fraction is: dt2/dtscr.  E.g. in the first hour it
c is (dtsrc/24)/dtsrc = 1/24th of the chemistry change is applied.
c This is to work around initial instabilities.

      if(allowSomeChemReinit == 1)then
        if(Itime-ItimeI <= 3)then
          dt2=dtsrc/24.d0          ! e.g. 150.
        elseif(Itime-ItimeI > 3 .and. Itime-ItimeI <= 6)then
          dt2=dtsrc/12.d0          ! e.g. 300.
        elseif(Itime-ItimeI > 6 .and. Itime-ItimeI <= 11)then
          dt2=dtsrc/6.d0           ! e.g. 600.
        elseif(Itime-ItimeI > 11 .and. Itime-ItimeI <= 30)then
          dt2=dtsrc/2.4d0          ! e.g. 1500.
        elseif(Itime-ItimeI > 30)then
          dt2=dtsrc                ! e.g. 3600
        endif
      else 
        dt2=dtsrc
      endif
      maxPSC=0.2d0/dt2

c Calculate new photolysis rates every n_phot main timesteps:
      MODPHOT= 0 ! old days was: MOD(Itime-ItimeI,n_phot)

C CALCULATE TX, THE REAL TEMPERATURE:
C (note this section is already done in DIAG.f)
      IF(HAVE_SOUTH_POLE) THEN
        DO L=1,LM
          TX(1,1,L)=T(1,1,L)*PK(L,1,1)
          TX(I_0:I_1,1,L)=TX(1,1,L)
        END DO
      ENDIF  
      IF(HAVE_NORTH_POLE) THEN
        DO L=1,LM
          TX(1,JM,L)=T(1,JM,L)*PK(L,1,JM)
          TX(I_0:I_1,JM,L)=TX(1,JM,L)
        END DO
      ENDIF
      DO L=1,LM
        DO J=J_0,J_1
          TX(I_0:I_1,J,L)=T(I_0:I_1,J,L)*PK(L,I_0:I_1,J)
        END DO
      END DO

#ifdef SHINDELL_STRAT_CHEM
C info to set strat H2O based on tropical tropopause H2O and CH4:
      if(Itime == ItimeI .and. allowSomeChemReinit == 1 )then
        avgTT_H2O_part(I_0:I_1,J_0:J_1)=0.d0
        avgTT_CH4_part(I_0:I_1,J_0:J_1)=0.d0
        countTT_part(I_0:I_1,J_0:J_1)=0.d0
        do J=J_0,J_1
          do I=I_0,IMAXJ(J)
            if(LAT2D_DG(I,J) >= -20. .and. LAT2D_DG(I,J) <= 20.)then
              avgTT_H2O_part(I,J) = Q(I,J,LTROPO(I,J))*MWabyMWw
              if(use_rad_ch4 > 0) then
                avgTT_CH4_part(I,J) =
     &          rad_to_chem(4,LTROPO(I,J),I,J)
     &          *2.69d20*byavog*mair*BYAM(LTROPO(I,J),I,J)
              else
                avgTT_CH4_part(I,J) =
     &          trm(I,J,LTROPO(I,J),n_CH4)
     &          *mass2vol(n_CH4)*BYAXYP(I,J)*BYAM(LTROPO(I,J),I,J)
              endif
              countTT_part(I,J) = 1.d0
            end if
          end do
        end do
        CALL GLOBALSUM(grid, avgTT_CH4_part, avgTT_CH4, all=.true.)
        CALL GLOBALSUM(grid, avgTT_H2O_part, avgTT_H2O, all=.true.)
        CALL GLOBALSUM(grid, countTT_part,   countTT,   all=.true.)
        if(countTT <= 0.)call stop_model('countTT.le.0',255)
      end if
#endif

! Define acetone in terms of Isoprene:
!kt Terpenes should also be included here in the future
      do j=J_0,J_1
        do i=I_0,IMAXJ(j)
          surfIsop(i,j)=trm(i,j,1,n_Isoprene)*mass2vol(n_Isoprene)*
     &         byaxyp(i,j)*byAM(1,i,j)
        enddo
      enddo
      call zonalmean_ij2ij(surfIsop,zonalIsop)
      do j=J_0,J_1
        do i=I_0,IMAXJ(j)
          select case(which_trop)
          case(0); maxl=ltropo(I,J)
          case(1); maxl=ls1-1
          case default; call stop_model('which_trop problem',255)
          end select
          do L=1,maxl
            acetone(i,j,L)=max(0.d0, ! in molec/cm3
     &      (1.25d0*(
     &        zonalIsop(i,j)-trm(i,j,L,n_Isoprene)*mass2vol(n_Isoprene)*
     &        byaxyp(i,j)*byAM(L,i,j)))*PMID(L,i,j)/(TX(i,j,L)*cboltz))
          enddo
          do L=maxl+1,LM
            acetone(i,j,L)=0.d0
          enddo
        enddo
      enddo

      ierr_loc = 0

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      j_loop: DO J=J_0,J_1          ! ---> MAIN J LOOP BEGINS <---

      i_loop: DO I=I_0,IMAXJ(J)     ! ---> MAIN I LOOP BEGINS <---
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      if(checktracer_on==1) call checktracer(I,J)
      select case(which_trop)
      case(0); maxl=ltropo(I,J)
      case(1); maxl=ls1-1
      case default; call stop_model('which_trop problem 4',255)
      end select

      DO L=1,LM
c Initialize the 2D change variable:
       changeL(L,:)=0.d0  ! (LM,NTM)
c Save presure, temperature, thickness, rel. hum. in local arrays:
       pres(L)=PMID(L,I,J)
       ta(L)=TX(I,J,L)
       rh(L)=Q(i,j,l)/min(1.d0,QSAT(ta(L),lhe,pres(L)))
       bythick(L)=1.d0/
     & (rgas*bygrav*TX(i,j,L)*LOG(pedn(L,i,j)/pedn(L+1,i,j)))
c Calculate M and set fixed ratios for O2 & H2:
       y(nM,L)=pres(L)/(ta(L)*cboltz)
       y(nO2,L)=y(nM,L)*pO2
#ifdef SHINDELL_STRAT_CHEM
       if(pres2(l) > 20.d0)then
         y(nH2,L)=y(nM,L)*pfix_H2
       else
         ! Was: y(nH2,L)=y(nM,L)*pfix_H2*7.d1/(7.d1+L-maxl+1)
         ! Now: a drop of 0.3 molec/cm3 per 12 hPa decrease:
         y(nH2,L)=y(nM,L)*pfix_H2 + 2.5d-2*(pres2(L)-20.d0)
       endif
       CLTOT=0.d0
#else
       y(nH2,L)=y(nM,L)*pfix_H2
#endif

c Tracers (converted from mass to number density):
       do igas=1,ntm_chem
         y(igas,L)=trm(I,J,L,igas)*y(nM,L)*mass2vol(igas)*
     &   BYAXYP(I,J)*BYAM(L,I,J)
       enddo

! If we are fixing methane for chemistry purposes set it's y here:
! 0.55866d0 below is 1/1.79 (HALOE observations)
      if(fix_CH4_chemistry == 1) then
        if(lat2d_dg(i,j) < 0.) then         ! Southern Hemisphere
          y(n_CH4,1:LS1-1)=y(nM,1:LS1-1)*ch4_init_sh*1.d-6 !troposphere
          if(abs(lat2d_dg(i,j)) > 30.) then ! extratropics
            y(n_CH4,LS1:LM)=                ! stratosphere
     &      y(nM,LS1:LM)*ch4_init_sh*0.55866d0*1.d-6*CH4altX(LS1:LM)
          else                              ! tropics
            y(n_CH4,LS1:LM)=                ! stratosphere
     &      y(nM,LS1:LM)*ch4_init_sh*0.55866d0*1.d-6*CH4altT(LS1:LM)
          end if
        else                                ! Northern Hemisphere
          y(n_CH4,1:LS1-1)=y(nM,1:LS1-1)*ch4_init_sh*1.d-6 !troposphere
          if(abs(lat2d_dg(i,j)) > 30.) then ! extratropics
            y(n_CH4,LS1:LM)=                ! stratosphere
     &      y(nM,LS1:LM)*ch4_init_nh*0.55866d0*1.d-6*CH4altX(LS1:LM)
          else                              ! tropics
            y(n_CH4,LS1:LM)=                ! stratosphere
     &      y(nM,LS1:LM)*ch4_init_nh*0.55866d0*1.d-6*CH4altT(LS1:LM)
          end if
        end if
      end if

C Concentrations of DMS and SO2 for sulfur chemistry:
       if (coupled_chem == 1) then
         ydms(i,j,l)=trm(i,j,l,n_dms)*y(nM,L)*(28.0D0/62.0D0)*
     &   BYAXYP(I,J)*BYAM(L,I,J)
         yso2(i,j,l)=trm(i,j,l,n_so2)*y(nM,L)*(28.0D0/64.0D0)*
     &   BYAXYP(I,J)*BYAM(L,I,J)
       else
         ! Convert from pptv to molecule cm-3:
         ydms(i,j,l)=dms_offline(i,j,l)*1.0D-12*y(nM,L)
         yso2(i,j,l)=so2_offline(i,j,l)*1.0D-12*y(nM,L)
       endif

#ifndef SHINDELL_STRAT_CHEM
c If desired, fix the methane concentration used in chemistry:
C WHY IS THIS NECESSARY ANY MORE, NOW THAT GET_CH4_IC IS CALLED?
       if(fix_CH4_chemistry == 1) THEN
         if(LAT2D_DG(I,J) < 0.) THEN ! SH
           y(n_CH4,L)=y(nM,L)*ch4_init_sh*1.d-6
         else                        ! NH
           y(n_CH4,L)=y(nM,L)*ch4_init_nh*1.d-6
         endif
       end if
#endif
#ifdef SHINDELL_STRAT_CHEM
c Save initial ClOx amount for use in ClOxfam:
       ClOx_old(L)=trm(I,J,L,n_ClOx)*y(nM,L)*mass2vol(n_ClOx)*
     & BYAXYP(I,J)*BYAM(L,I,J)
#endif

c Limit N2O5 number density:
       if(y(n_N2O5,L) < 1.) y(n_N2O5,L)=1.d0
c Set H2O, based on Q:
       y(nH2O,L)=Q(I,J,L)*MWabyMWw*y(nM,L)
#ifdef SHINDELL_STRAT_CHEM
c Initialize stratospheric y(H2O) & GCM Q variable (!),
c based on tropical tropopause H2O and CH4:
       if(allowSomeChemReinit == 1) then
         if(Itime == ItimeI .and. L > LTROPO(I,J)) then
           y(nH2O,L) =  y(nM,L)*(avgTT_H2O/countTT +
     &     2.d0*(avgTT_CH4/countTT-y(n_CH4,L)/y(nM,L)))
           if(clim_interact_chem > 0)then 
             fraQ=(y(nH2O,L)/(y(nM,L)*MWabyMWw))/Q(I,J,L)
             Q(I,J,L)=y(nH2O,L)/(y(nM,L)*MWabyMWw)
             if(fraQ < 1.)qmom(:,i,j,l)=qmom(:,i,j,l)*fraQ
#ifdef TRACERS_WATER
C**** Add water to relevant tracers as well
              do n=1,ntm
                select case (tr_wd_type(n))
                case (nWater)       ! water: initialise tracers
                  trm(i,j,l,n) = trm(i,j,l,n)*fraQ
                  if(fraQ < 1.)trmom(:,i,j,l,n) = trmom(:,i,j,l,n)*fraQ
                end select
              end do
#endif
           end if
         end if
       end if
#endif

c Initialize various other species:
c - set [NO]=0 (where?) for first HOx calc, NO2 = NOx:
c - set reactive species for use in family chemistry & nighttime NO2:

       if(Itime == ItimeI .and. allowSomeChemReinit == 1)then 
         y(nAldehyde,L)=y(nM,L)*pfix_Aldehyde
       else
         y(nAldehyde,L)=yAldehyde(I,J,L)
       endif
       yNO3(I,J,L)=pNO3(I,J,L)*y(n_NOx,L)
       y(nNO2,L)     =y(n_NOx,L)*pNOx(I,J,L)
       y(nNO,L)      =y(n_NOx,L)-(y(nNO2,L)+yNO3(I,J,L))
       if(y(nNO,L).lt.1.d0)y(nNO,L)=1.d0
       y(nO3,L)      =pOx(I,J,L)*y(n_Ox,L)
       y(nCH3O2,L)   =yCH3O2(I,J,L)
       y(nC2O3,L)    =yC2O3(I,J,L)
       y(nXO2,L)     =yXO2(I,J,L)
       y(nXO2N,L)    =yXO2N(I,J,L)
       y(nRXPAR,L)   =yRXPAR(I,J,L)
       y(nROR,L)     =yROR(I,J,L)
#ifdef SHINDELL_STRAT_CHEM
       y(nCl2,L)     =yCl2(I,J,L)
       y(nCl2O2,L)   =yCl2O2(I,J,L)
       y(nOClO,L)    =y(n_ClOx,L)*pOClOx(I,J,L)
       y(nClO,L)     =y(n_ClOx,L)*pClOx(I,J,L)
       y(nCl,L)      =y(n_ClOx,L)*pClx(I,J,L)
       y(nBr,L)      =y(n_BrOx,L)*(1.d0-pBrOx(I,J,L))
       y(nBrO,L)     =y(n_BrOx,L)*pBrOx(I,J,L)
#endif
      END DO ! L

C For solar zenith angle, we use the arccosine of the COSZ1
C from the radiation code, which is the cosine of the solar zenith 
C angle averaged over the physics time step.
C If the solar zenith angle (sza) from the radiation code is > 90 deg,
C (and hence COSZ1 is set to 0), recalculate it with get_sza routine:
      IF(COSZ1(I,J) == 0.d0) THEN
        call get_sza(I,J,sza)
      ELSE
        sza = acos(COSZ1(I,J))*byradian
      END IF

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                 BEGIN PHOTOLYSIS                               C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      if(MODPHOT == 0)then ! i.e. not (necessarily) every time step

       ! additional SUNLIGHT criterion (see also fam chem criterion):
       if((ALB(I,J,1) /= 0.d0).AND.(sza < szamax))then
       
c Pass O3 array (in ppmv) to fastj. Above these levels fastj2 uses
C Nagatani climatological O3, read in by chem_init: 
        DO LL=1,LM
          O3_FASTJ(LL)=y(nO3,LL)/y(nM,LL)
        END DO

! calculate photolysis rates
        call fastj2_drv(I, J, ta, rh)
        call photo_acetone(I,J,sza*radian,Jacet) ! simpler calculation for acetone

C Define and alter resulting photolysis coefficients (zj --> ss):
#ifdef SHINDELL_STRAT_CHEM
        colmO2=5.6d20 
        colmO3=5.0d16 
#endif
        ! Using MAX() here because
        ! Letting these spherical corrections get too small (0?) causes NaNs
        ! e.g. in the SOA code. Also sza can be > 90, so COS could go negative:
        sphericalCorrectionReg1=
     &    (MAX(1.d-2,DCOS(sza*radian)))**reg1Power_SpherO2andN2Ocorr
        sphericalCorrectionReg2=
     &    (MAX(1.d-2,DCOS(sza*radian)))**reg2Power_SpherO2andN2Ocorr
        sphericalCorrectionReg3=
     &    (MAX(1.d-2,DCOS(sza*radian)))**reg3Power_SpherO2andN2Ocorr
        sphericalCorrectionReg4=
     &    (MAX(1.d-2,DCOS(sza*radian)))**reg4Power_SpherO2andN2Ocorr

        DO L=min(JPNL,LM),1,-1
          do inss=1,JPPJ
            ss(inss,L,I,J)=zj(L,inss)
#ifdef SHINDELL_STRAT_CHEM 
#ifndef SHINDELL_SKIP_WINDOW_TUNE
            !reduce rates for gases that photolyze in window region (~200nm):
            if(inss == 27 .or. inss == 28) then ! for O2 and N2O reactions:
              ! Apply spherical corrections:
              if(pres2(L)>reg1TopPres_SpherO2andN2Ocorr)then
                ss(inss,L,I,J)=ss(inss,L,I,J)*sphericalCorrectionReg1
              else if(pres2(L)>reg2TopPres_SpherO2andN2Ocorr .and.
     &                pres2(L).le.reg1TopPres_SpherO2andN2Ocorr) then
                ss(inss,L,I,J)=ss(inss,L,I,J)*sphericalCorrectionReg2
              else if(pres2(L)>reg3TopPres_SpherO2andN2Ocorr .and.
     &                pres2(L).le.reg2TopPres_SpherO2andN2Ocorr) then
                ss(inss,L,I,J)=ss(inss,L,I,J)*sphericalCorrectionReg3
              else
                ss(inss,L,I,J)=ss(inss,L,I,J)*sphericalCorrectionReg4
              end if
            end if
            ! Then apply linear corrections for same reactions:
            if(inss == 27) then
              ss(inss,L,I,J)=ss(inss,L,I,J)*windowO2corr
            end if
            if(inss == 28) then
              ss(inss,L,I,J)=ss(inss,L,I,J)*windowN2Ocorr
            end if
#endif /* not defined to skip */
#endif /* SHINDELL_STRAT_CHEM */
          enddo
          taijls(i,j,l,ijlt_JO1D)=taijls(i,j,l,ijlt_JO1D)+ss(2,l,i,j)
          taijls(i,j,l,ijlt_JNO2)=taijls(i,j,l,ijlt_JNO2)+ss(1,l,i,j)
          taijls(i,j,l,ijlt_JH2O2)=taijls(i,j,l,ijlt_JH2O2)+ss(4,l,i,j)
#ifdef SHINDELL_STRAT_CHEM
          thick=
     &    1.d-3*rgas*bygrav*TX(I,J,L)*LOG(PEDN(L,i,j)/PEDN(L+1,i,j))
          colmO2=colmO2+y(nO2,L)*thick*1.d5
          colmO3=colmO3+y(nO3,L)*thick*1.d5
! SF3 is photolysis of water in Schumann-Runge bands based on:
! Nicolet, Pl. Space Sci., p 871, 1983.
! SF3_fact is, if x[ ] = bin4_flux[ ]:
! {(x[present] - x[1988]) / (x[1991] - x[1988])} * 0.1E-6
! This gets ADDED to the 1.3E-6 factor in the SF3 calculation. Here,
! bin4_flux is a proxy for the flux from all 175-200nm bins. 
! (Drew says the ratio would be the same.)
          if(SF2_fact == 0.)call stop_model('SF2_fact=0 in master',255)
          if(pres2(L) <= 10.)then
            if((SF3_FACT+1.3d-6) < 0.)call stop_model
     &      ('(SF3_FACT+1.3d-6) < 0 in master',255)
            SF3(I,J,L)=6.d0*(SF3_FACT+1.3d-6)*EXP(-1.d-7*colmO2**.35)
     &      *by35*SQRT(1.224d3*COSZ1(I,J)**2.+1d0)
            SF3(I,J,L)=SF3(I,J,L)*5.d-2
          else
            SF3(I,J,L)=0.d0
          endif
! SF2 is photlysis of NO in bands (0-0) and (1-0) based on Nicolet,
! Pl. Space Sci., p 111, 1980. SF2_fact is a ratio representative of
! bands (0-0) and (1-0); =  bin5_flux[present] / bin5_flux[1988] :
          if(colmO2 > 2.d19)then
            SF2(I,J,L)=4.5d-6*EXP(-(1.d-8*colmO2**.38+5.d-19*colmO3))
          else
            SF2(I,J,L)=4.75d-6*EXP(-1.5d-20*colmO2)
          endif
          SF2(I,J,L)=SF2(I,J,L)*SF2_fact*
     &    by35*SQRT(1.224d3*COSZ1(I,J)**2.+1.d0)
#endif
        END DO

       endif ! (sunlight)
      endif  ! (modphot)                           
      
CCCCCCCCCCCCCCCCC END PHOTOLYSIS SECTION CCCCCCCCCCCCCCCCCCCCCCCCC

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                 BEGIN CHEMISTRY                                C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


#ifdef SHINDELL_STRAT_CHEM
      ! Check for PSC's existance. 
      pscX(:)=.false.
      do L=1,LM 
        if(pres2(L) <= 250.d0 .and. pres2(L) >= 30.d0)then    ! pres crit for now
          if(lat2D_dg(I,J)<=PSClatS.or.lat2D_dg(I,J)>=PSClatN)then! lat crit for now
            if(lat2d_dg(I,J)<=PSClatS)then
              Ttemp=ta(L)+Tpsc_offset_S
            else if(lat2d_dg(I,J)>=PSClatN)then
              Ttemp=ta(L)+Tpsc_offset_N
            endif
            if(Ttemp <= T_thresh)then                      ! cold enough forya?
              bHNO3=38.9855d0-11397.d0/Ttemp+0.009179d0*Ttemp ! H2O and HNO3
              mHNO3= -2.7836d0 - 8.8d-4*Ttemp                 ! criteria from
              HNO3_thresh=2.69d19/760.d0*10.d0**              ! Hanson+Mauersberger
     &        (mHNO3*log10(y(nH2O,L)*760.d0/2.69d19)+bHNO3)   ! 1988
              if(y(n_HNO3,L) >= HNO3_thresh) pscX(L)=.true.! <-- yes PSC
            endif ! temperature
          endif   ! lat
        endif     ! pressure
      enddo       ! L
#endif
     
c Calculate the chemical reaction rates:
      call Crates(I,J,rh,bythick
#ifdef SHINDELL_STRAT_CHEM
     &                ,aero
#endif
     &           )      

#ifdef TRACERS_AEROSOLS_SOA
      voc2nox(:)=0.d0
#endif  /* TRACERS_AEROSOLS_SOA */

      if((ALB(I,J,1) /= 0.d0).AND.(sza < szamax))then
CCCCCCCCCCCCCCCCCCCC   SUNLIGHT   CCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

c      When PSCs present, ensure heterogenous reactions do not destroy
c      more than exists to prevent non-conserving reaction overflows at
c      end of polar night:
       do L=1,LM
         if(pscX(L))then
           if(rr(iN2O5plusH2O,L)>=maxPSC) rr(iN2O5plusH2O,L)=maxPSC
           if(rr(iClONO2plusH2O,L)>=maxPSC) rr(iClONO2plusH2O,L)=maxPSC
           if(rr(iClONO2plusHCl,L)>=maxPSC) rr(iClONO2plusHCl,L)=maxPSC
           if(rr(iClONO2plusHCl,L)>=maxPSC*(y(n_HCl,L)/y(n_ClONO2,L)))
     &       rr(iClONO2plusHCl,L)=maxPSC*(y(n_HCl,L)/y(n_ClONO2,L))
           if(rr(iHOClplusHCl,L)>=maxPSC) rr(iHOClplusHCl,L)=maxPSC
           if(rr(iHOClplusHCl,L)>=maxPSC*(y(n_HCl,L)/y(n_HOCl,L)))
     &       rr(iHOClplusHCl,L)=maxPSC*(y(n_HCl,L)/y(n_HOCl,L))
           if(rr(iN2O5plusHCl,L)>=maxPSC) rr(iN2O5plusHCl,L)=maxPSC
           if(rr(iN2O5plusHCl,L)>=maxPSC*(y(n_HCl,L)/y(n_N2O5,L)))
     &       rr(iN2O5plusHCl,L)=maxPSC*(y(n_HCl,L)/y(n_N2O5,L))
         end if
       end do


CCCCCCCCCCCCCCCCC FAMILY PARTITIONING CCCCCCCCCCCCCCCCCCCCCCCCCC

#ifdef SHINDELL_STRAT_CHEM
       call Oxfam(LM,I,J)
       call HOxfam(LM,I,J)
       call NOxfam(LM,I,J)
       call BrOxfam(LM,I,J)
#else
       call Oxfam(maxl,I,J)
       call HOxfam(maxl,I,J)
       call NOxfam(maxl,I,J)
#endif
CCCCCCCCCCCCCCCCC NON-FAMILY CHEMISTRY CCCCCCCCCCCCCCCCCCCCCCCC
#ifdef TRACERS_AEROSOLS_SOA
! calculate voc2nox for SOA precursor chemistry
      do L=1,LM
        voc2nox_denom=(4.2d-12*exp(180.d0/ta(L))*y(nNO,L)+
     &                 rr(43,L)*y(nHO2,L)+
     &                 1.7d-14*exp(1300.d0/ta(L))*yXO2(I,J,L))
        if (voc2nox_denom==0.d0) then
          voc2nox(L)=0.d0
        else
          voc2nox(L)=4.2d-12*exp(180.d0/ta(L))*y(nNO,L)/
     &               voc2nox_denom
        endif
      enddo
#endif  /* TRACERS_AEROSOLS_SOA */

      call chemstep(I,J,changeL,ierr_loc)
      if(ierr_loc > 0) cycle i_loop

C Save 3D radical arrays to pass to aerosol code:
      if(coupled_chem == 1) then
        do l=1,LM
          oh_live(i,j,l)=y(nOH,L)
          no3_live(i,j,l)=yNO3(i,j,l)
        end do
      endif

#ifdef SHINDELL_STRAT_CHEM
      call ClOxfam(LM,I,J,ClOx_old) ! needed something from chemstep.
#endif

      call printDaytimeChemistryDiags()
      
      else

CCCCCCCCCCCCCCCCCCCC END SUNLIGHT CCCCCCCCCCCCCCCCCCCCCCCCCCCCC

CCCCCCCCCCCCCCCCCCCC   DARKNESS  CCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C*****************************************************************
C               ABOUT: N2O5 sink on sulfate aerosol 
C                REACTION PROBABLITY FORMULATION:
C
C To evaluate 'k' for N2O5 + H2O -> 2HNO3, assume k = GAMMA*SA*v/4
C (Dentener and Crutzen, 1993). GAMMA = RGAMMASULF defined below,
c SA = given, v = molecular velocity (cm/s) where v = SQRT (8.Kb.T / PI*M);
C Kb = 1.38062E-23; T = Temp (K); M = mass N2O5 (Kg)
C
C Off-line sulfate fields to run in uncoupled mode give SO4 in
C cm2/cm3 'surface area density'.
C
C On-line sulfate in coupled mode must be converted to aerosol
C surface (cm2 aerosol/cm3 air) via aerosol volume fraction
C (cm3 aerosol/cm3 air). Assume a monodispersed aerosol with diameter
C of 0.078 microns. Specific aerosol density = 1.1g/cm3 (1.7g/cm3 ?)
C See Dentener and Crutzen, 1993 for details. Mr[sulfate] = 96.0g;
C Mr[(NH4)HSO4] = 115.0gC
C*****************************************************************

CCCCCCCCCCCCCCCC NIGHTTIME CCCCCCCCCCCCCCCCCCCCCC
#ifdef SHINDELL_STRAT_CHEM
      LL=LM
#else
      LL=maxl
#endif
#ifdef TRACERS_AEROSOLS_SOA
      call soa_apart ! calculate current apartmolar factors
#ifdef SOA_DIAGS
     &              (I,J)
#endif  /* SOA_DIAGS */
#endif  /* TRACERS_AEROSOLS_SOA */
      DO L=1,LL  

! ----------------- RGAMMASULF ---------------------------------------
! Until more sophisticated method arrives, or when aerosol tracers are
! off, use method recommended by Faye, based on Kane et al., JPC, 2001
! and Hallquist et al., PCCP, 2003:
! For RH>50%, gamma=0.015. For RH <=50%:
! T<=290K, gamma=0.052 - 2.79d-4*RH [RH in percent]
! T>290K, gamma=above - log10(T-290)*0.05 [minimum gamma=0.001]
!
        if(rh(l)>0.5)then
          rgammasulf = 1.5d-2
        else
          rgammasulf = 5.2d-2 - 2.79d-4*100.d0*rh(l)
          if(ta(l)>290.) rgammasulf=
     &    max(1.d-3,rgammasulf-log10(ta(l)-290.d0)*5.d-2)
        endif
! --------------------------------------------------------------------
        if (coupled_chem == 1) then
          ! Convert SO4 from mass (kg) to aerosol surface per grid box:
          ! Here there is a factor of 1d-3  that converts kg/m3 to g/cm3
          ! and 1.76d5 is cm2/g from Dentener and Crutzen, 1993.
          ! So 1.d-3*1.76d5=1.76d2, and that value is for a relative
          ! humidity of 0.75 (1/0.75 = 1.33333d0 below). Reciprocal
          ! layer thickness below is in 1/m units:
#ifndef TRACERS_TOMAS
          sulfate(i,j,L)=
#ifdef TRACERS_AMP
     &      (trm(i,j,L,n_M_AKK_SU)+trm(i,j,L,n_M_ACC_SU)+
     &       trm(i,j,L,n_M_DD1_SU)+trm(i,j,L,n_M_DS1_SU)+
     &       trm(i,j,L,n_M_DD2_SU)+trm(i,j,L,n_M_DS2_SU)+
     &       trm(i,j,L,n_M_SSA_SU)+trm(i,j,L,n_M_OCC_SU)+
     &       trm(i,j,L,n_M_BC1_SU)+trm(i,j,L,n_M_BC2_SU)+
     &       trm(i,j,L,n_M_BC3_SU)+trm(i,j,L,n_M_DBC_SU)+
     &       trm(i,j,L,n_M_BOC_SU)+trm(i,j,L,n_M_BCS_SU)+
     &       trm(i,j,L,n_M_MXX_SU))
#else
     &      trm(i,j,L,n_SO4)
#endif
     &      *1.76d2*byaxyp(i,j)*bythick(L)
     &      *max(0.1d0,rh(L)*1.33333d0)
#endif
#ifdef TRACERS_TOMAS
           sulfate(i,j,l)=0.0
           do nb=1,nbins
              sulfate(i,j,L)=sulfate(i,j,l)+trm(i,j,L,n_ASO4(nb))*1.76d2
     &            *byaxyp(i,j)*bythick(L)*max(0.1d0,rh(L)*1.33333d0)
           enddo
#endif
        endif

        pfactor=axyp(I,J)*AM(L,I,J)/y(nM,L)
        bypfactor=1.D0/pfactor
        RVELN2O5=SQRT(TX(I,J,L)*RKBYPIM)*100.d0
C       Calculate sulfate sink, and cap it at 20% of N2O5:
c       in troposphere loss is rxn on sulfate, in strat rxn w PSC or sulfate
        wprod_sulf=
     &  dt2*sulfate(I,J,L)*y(n_N2O5,L)*RGAMMASULF*RVELN2O5*0.25d0
#ifdef SHINDELL_STRAT_CHEM
        if(pres2(L)>5.d0)then
c       if there is reaction on strat particulate (in Crates), use that
          if(rr(iN2O5plusH2O,L)>1.0d-25)wprod_sulf=
     &                       DT2*y(n_N2O5,L)*rr(iN2O5plusH2O,L)
        else
          wprod_sulf=0.d0
        endif
#endif
        ! we used to limit (w)prod_sulf here w.r.t. N2O5
#ifndef SHINDELL_STRAT_CHEM /* that is a not def */
        CALL INC_TAJLS(I,J,L,jls_N2O5sulf,-prod_sulf*vol2mass(n_N2O5))
#endif

C*****************************************************************
c g signifies gas phase
c while prod_sulf and wprod_sulf are sulfate rxn
c wprods are in molecules/cm3/s
c prods/mass2vol are in mass units to add to tracers
c
c NO3 amounts are a function of reaction 7 (NO2 + O3 -> NO3),
c 24, 25 (leave out 28, 0.9*32, iN2O5decomp (47), iNO3plusNO2 (53)
c outside NOx family) NO2, similarly leave out 29, 45, and
c iN2O5decomp (47). Keep NOx unchanged as this is only intrafamily.
C*****************************************************************

c       calculate NO3 vs NO2 (assume no NO at night)
        do itemp_iter=1,5
          rprodNO3=rr(7,L)*y(nNO2,L)*y(n_Ox,L)*pOx(I,J,L)
          rlossNO3=(2.d0*rr(25,L)*yNO3(I,J,L))*yNO3(I,J,L)-
     &    (rr(36,L)*y(n_Alkenes,L))*yNO3(I,J,L)
          if((rlossNO3+rprodNO3)==0.)
     &    call stop_model('(rlossNO3+rprodNO3)=0',255)
          pNO3temp=rprodNO3/(rlossNO3+rprodNO3)
          if(pNO3temp.gt.0.99d0)pNO3temp=0.95d0
          if(pNO3temp.lt.0.01d0)pNO3temp=0.05d0
          yNO3(I,J,L)=pNO3temp*y(n_NOx,L)
          y(nNO2,L)=y(n_NOx,L)-yNO3(I,J,L)
        end do

c       set reaction rates, then limit any uniformly across all
c       paths if lead to negative conc:
        rNO3plusNO2=rr(iNO3plusNO2,L)*y(nNO2,L)*yNO3(I,J,L)*dt2
        rN2O5decomp=rr(iN2O5decomp,L)*y(n_N2O5,L)*dt2
        chgHT5=rr(iN2O5plusHCl,L)*y(n_N2O5,L)*dt2
        rHCHOplusNO3=y(n_HCHO,L)*rr(28,L)*yNO3(I,J,L)*dt2
        rAldplusNO3=2.5d-15*yAldehyde(I,J,L)*yNO3(I,J,L)*dt2
        rIsopplusNO3=rr(32,L)*y(n_Isoprene,L)*yNO3(I,J,L)*dt2
#ifdef TRACERS_TERP
        rTerpplusNO3=rr(iTerpenesNO3,L)*y(n_Terpenes,L)
     &   *yNO3(I,J,L)*dt2
#endif  /* TRACERS_TERP */
        rClOplusNO2=y(nClO,L)*rr(iClOplusNO2,L)*y(nNO2,L)*dt2
        rDMSplusNO3=ydms(i,j,l)*rsulf3(i,j,l)*yNO3(I,J,L)*dt2
        rBrOplusNO2=rr(108,L)*y(n_NOx,L)*pNOx(I,J,L)
     &   *y(n_BrOx,L)*pBrOx(I,J,L)*dt2
        chgHT3=rr(iClONO2plusHCl,L)*y(n_ClONO2,L)*dt2
        changehetClONO2=-1.d0*(rr(iClONO2plusH2O,L)*y(n_ClONO2,L))*dt2

c       Examine prod and loss of N2O5
        wlossN2O5=rN2O5decomp+wprod_sulf+chgHT5
        wprodN2O5=rNO3plusNO2
        if(wlossN2O5>0.99d0*(y(n_N2O5,L)+wprodN2O5))then
          if(wlossN2O5==0.)call stop_model('wprodN2O5=0',255)
          ratioN=0.99d0*(y(n_N2O5,L)+wprodN2O5)/wlossN2O5
          rN2O5decomp=rN2O5decomp*ratioN
          wprod_sulf=wprod_sulf*ratioN
          chgHT5=chgHT5*ratioN
        endif

c       Examine prod and loss of NOx
        wlossNOx=rHCHOplusNO3+rAldplusNO3+2.d0*rNO3plusNO2
     &   +0.9d0*rIsopplusNO3+rClOplusNO2+rDMSplusNO3+rBrOplusNO2
#ifdef TRACERS_TERP
     &   +0.9d0*rTerpplusNO3
#endif  /* TRACERS_TERP */
        wprodNOx=2.d0*rN2O5decomp
        if(wlossNOx>0.99d0*(y(n_NOx,L)+wprodNOx))then
          if(wlossNOx==0.)call stop_model('wlossNOx=0',255)
          ratioN=0.99d0*(y(n_NOx,L)+wprodNOx)/wlossNOx
          rHCHOplusNO3=rHCHOplusNO3*ratioN
          rAldplusNO3=rAldplusNO3*ratioN
          rNO3plusNO2=rNO3plusNO2*ratioN
          rIsopplusNO3=rIsopplusNO3*ratioN
#ifdef TRACERS_TERP
          rTerpplusNO3=rTerpplusNO3*ratioN
#endif  /* TRACERS_TERP */
          rClOplusNO2=rClOplusNO2*ratioN
          rDMSplusNO3=rDMSplusNO3*ratioN
          rBrOplusNO2=rBrOplusNO2*ratioN
        endif

#ifdef SHINDELL_STRAT_CHEM
         if(rClOplusNO2 >= y(nClO,L))rClOplusNO2=0.8d0*y(nClO,L)
         changeClONO2=rClOplusNO2
         changeClOx=-1.d0*changeClONO2
#endif

c       Examine prod and loss of ClONO2
        wlossClONO2=changehetClONO2+chgHT3
        wprodClONO2=rClOplusNO2
        if(wlossClONO2>0.99d0*(y(n_ClONO2,L)+wprodClONO2))then
          if(wlossClONO2==0.)call stop_model('wprodClONO2=0',255)
          ratioN=0.99d0*(y(n_ClONO2,L)+wprodClONO2)/wlossClONO2
          changehetClONO2=changehetClONO2*ratioN
          chgHT3=chgHT3*ratioN
        endif

        prod_sulf=wprod_sulf*pfactor

C       LOWER LIMIT ON N2O5:
        if(y(n_N2O5,L) <= 1.d0) y(n_N2O5,L)=1.d0

C Calculate and limit gaseous changes to HNO3, HCHO, N2O5, Aldehyde,
C Alkenes, Isoprene, Terpenes (if used) and AlkylNit:

        gwprodHNO3=rHCHOplusNO3+rAldplusNO3
        if(gwprodHNO3 > y(n_HCHO,L))gwprodHNO3=y(n_HCHO,L)
        gprodHNO3=gwprodHNO3*pfactor

        changeAldehyde=(rr(36,L)*y(n_Alkenes,L)
     &                 +rr(32,L)*y(n_Isoprene,L)*0.12d0
#ifdef TRACERS_TERP
     &                 +rr(iTerpenesNO3,L)*y(n_Terpenes,L)*0.12d0
#endif  /* TRACERS_TERP */
     &                 -2.5d-15*yAldehyde(I,J,L))*
     &              yNO3(I,J,L)*dt2
        if(-changeAldehyde > 0.75d0*yAldehyde(I,J,L))changeAldehyde=
     &  -0.75d0*yAldehyde(I,J,L)

        changeAlkenes=(rr(32,L)*y(n_Isoprene,L)*0.45d0
#ifdef TRACERS_TERP
     &                +rr(iTerpenesNO3,L)*y(n_Terpenes,L)*0.45d0
#endif  /* TRACERS_TERP */
     &                -rr(36,L)*y(n_Alkenes,L))*
     &               yNO3(I,J,L)*dt2
     &               +(rr(31,L)*y(n_Isoprene,L)
#ifdef TRACERS_TERP
     &                +rr(iTerpenesO3,L)*y(n_Terpenes,L)
#endif  /* TRACERS_TERP */
     &                -rr(35,L)*y(n_Alkenes,L))*
     &               y(nO3,L)*dt2
        if(-changeAlkenes > 0.75d0*y(n_Alkenes,L))changeAlkenes=
     &  -0.75d0*y(n_Alkenes,L)

#ifdef TRACERS_AEROSOLS_SOA
! WARNING!!!
! No nighttime production from OH reactions, since no nighttime OH exist.
! This should be improved in the future.
      if (L<=LM_soa) then
        changeisopp1g=apartmolar(L,whichsoa(n_isopp1a))*
     &                (rr(31,L)*y(nO3,L))*y(n_Isoprene,L)*dt2
        if(-changeisopp1g > 0.75d0*y(n_isopp1g,L))changeisopp1g=
     &  -0.75d0*y(n_isopp1g,L)
        changeisopp2g=apartmolar(L,whichsoa(n_isopp2a))*
     &                (rr(31,L)*y(nO3,L))*y(n_Isoprene,L)*dt2
        if(-changeisopp2g > 0.75d0*y(n_isopp2g,L))changeisopp2g=
     &  -0.75d0*y(n_isopp2g,L)
#ifdef TRACERS_TERP
        changeapinp1g=apartmolar(L,whichsoa(n_apinp1a))*
     &                (rr(iTerpenesO3,L)*y(nO3,L))*y(n_Terpenes,L)*dt2
        if(-changeapinp1g > 0.75d0*y(n_apinp1g,L))changeapinp1g=
     &  -0.75d0*y(n_apinp1g,L)
        changeapinp2g=apartmolar(L,whichsoa(n_apinp2a))*
     &                (rr(iTerpenesO3,L)*y(nO3,L))*y(n_Terpenes,L)*dt2
        if(-changeapinp2g > 0.75d0*y(n_apinp2g,L))changeapinp2g=
     &  -0.75d0*y(n_apinp2g,L)
#endif  /* TRACERS_TERP */
      else
        changeisopp1g=0.d0
        changeisopp2g=0.d0
#ifdef TRACERS_TERP
        changeapinp1g=0.d0
        changeapinp2g=0.d0
#endif  /* TRACERS_TERP */
      endif
#endif  /* TRACERS_AEROSOLS_SOA */

        changeIsoprene=-(rr(32,L)*yNO3(I,J,L)
     &                 +rr(31,L)*y(nO3,L))*y(n_Isoprene,L)*dt2
        if(-changeIsoprene > 0.75d0*y(n_Isoprene,L))changeIsoprene=
     &  -0.75d0*y(n_Isoprene,L)

#ifdef TRACERS_TERP
        changeTerpenes=-(rr(iTerpenesNO3,L)*yNO3(I,J,L)
     &                 +rr(iTerpenesO3,L)*y(nO3,L))*y(n_Terpenes,L)*dt2
        if(-changeTerpenes > 0.75d0*y(n_Terpenes,L))changeTerpenes=
     &  -0.75d0*y(n_Terpenes,L)
#endif  /* TRACERS_TERP */

        changeHCHO=(rr(36,L)*y(n_Alkenes,L)
     &             +rr(32,L)*y(n_Isoprene,L)*0.03d0
#ifdef TRACERS_TERP
     &             +rr(iTerpenesNO3,L)*y(n_Terpenes,L)*0.03d0
#endif  /* TRACERS_TERP */
     &             )*yNO3(I,J,L)*dt2
     &             -gwprodHNO3+(rr(31,L)*y(n_Isoprene,L)*0.9d0
#ifdef TRACERS_TERP
     &                         +rr(iTerpenesO3,L)*y(n_Terpenes,L)*0.9d0
#endif  /* TRACERS_TERP */
     &             +rr(35,L)*y(n_Alkenes,L))*y(nO3,L)*0.64d0*dt2

        changeAlkylNit=rIsopplusNO3*0.9d0
#ifdef TRACERS_TERP
     &                +rTerpplusNO3*0.9d0
#endif  /* TRACERS_TERP */
        ! removed limit here on changeAlkylNit

c Convert some changes to molecules/cm3/s:
        changeHNO3=gwprodHNO3+2*wprod_sulf  !always positive

        wlossNOx=rHCHOplusNO3+rAldplusNO3+2.d0*rNO3plusNO2
     &   +0.9d0*rIsopplusNO3+rClOplusNO2+rDMSplusNO3
#ifdef TRACERS_TERP
     &   +0.9d0*rTerpplusNO3
#endif  /* TRACERS_TERP */
        wprodNOx=2.d0*rN2O5decomp

        changeNOx=wprodNOx-wlossNOx

        wlossN2O5=rN2O5decomp+wprod_sulf
        gwprodN2O5=rNO3plusNO2
        changeN2O5=gwprodN2O5-wlossN2O5

#ifdef SHINDELL_STRAT_CHEM

c       Nighttime changes in Bromine-containing species
        changeBrOx=-1.d0*rBrOplusNO2
        if(rBrOplusNO2>0.5d0*y(n_BrOx,L))
     &   changeBrOx=-0.5d0*y(n_BrOx,L)
        changeBrONO2=-changeBrOx
        changeNOx=changeNOx+changeBrOx
        if(-1.d0*changeNOx>y(n_NOx,L))changeNOx=-1.d0*y(n_NOx,L)

c        Br+H2O2 converts to HBr+HO2. HO2 assumed to revert to H2O2
        changeBrOx2=-rr(79,L)*y(n_H2O2,L)
     &    *y(n_BrOx,L)*(1.d0-pBrOx(I,J,l))*dt2
        if(-changeBrOx2>0.2d0*y(n_BrOx,L))changeBrOx2=-0.2d0*y(n_BrOx,L)
        changeBrOx=changeBrOx+changeBrOx2
        changeHBr=-changeBrOx2

c Heterogeneous reaction ClONO2+H2O on sulfate (and PSCs if present):
         if(rr(iClONO2plusH2O,L) > 2.d-35)then
           changeClONO2=changeClONO2+changehetClONO2
           changeHOCl=-changehetClONO2
           changeHNO3=changeHNO3-changehetClONO2
         else
           changeHOCl=0.d0
         endif
         
         changeHCl=0.d0

c Polar Stratospheric Clouds (PSC) Chemistry:
c 106 N2O5    +H2O     -->HNO3    +HNO3  (calculated above)
c 107 ClONO2  +H2O     -->HOCl    +HNO3  (calculated above)
c 108 ClONO2  +HCl     -->Cl      +HNO3  !really makes Cl2
c 109 HOCl    +HCl     -->Cl      +H2O   !really makes Cl2
c 110 N2O5    +HCl     -->Cl      +HNO3  !really makes ClNO2 (calc above)
         if(pscX(L)) then  ! PSCs exist
           if(chgHT3 >= 0.2d0*y(n_HCl,L))chgHT3=0.2d0*y(n_HCl,L)
           chgHT4=rr(iHOClplusHCl,L)*y(n_HOCl,L)*dt2
           if(chgHT4 >= 0.2d0*y(n_HOCl,L))chgHT4=0.2d0*y(n_HOCl,L)
           if(chgHT4 >= 0.2d0*y(n_HCl,L))chgHT4=0.2d0*y(n_HCl,L)
           if(chgHT5 >= 0.5d0*y(n_HCl,L))chgHT5=0.5d0*y(n_HCl,L)
           changeClONO2=changeClONO2-chgHT3
           changeHOCl=changeHOCl-chgHT4
           changeN2O5=changeN2O5-chgHT5
c          Note that really the following 3 produce Cl2, not ClOx, and Cl2
C          at night is stable and doesn't go back into ClONO2, so
C          should eventually keep track of Cl2/ClOx partitioning!
           changeHCl=changeHCl-chgHT3-chgHT4-chgHT5
           changeHNO3=changeHNO3+chgHT3+chgHT5
           changeClOx=changeClOx+chgHT3+chgHT4+chgHT5
           ! Here we USED TO remove some of the HNO3 formed heterogeneously,
           ! as it doesn't come back to the gas phase. No longer.
         endif
#endif

#if(defined TRACERS_HETCHEM) || (defined TOMAS_HETCHEM)
C       Include reactions on dust for HNO3:
        changeHNO3 = changeHNO3 - krate(i,j,l,1,1)*y(n_HNO3,l)*dt2

#ifdef TRACERS_HETCHEM
        changeN_d1 = krate(i,j,l,2,1) * y(n_HNO3,l) *dt2
        changeN_d2 = krate(i,j,l,3,1) * y(n_HNO3,l) *dt2
        changeN_d3 = krate(i,j,l,4,1) * y(n_HNO3,l) *dt2
#endif

#ifdef TOMAS_HETCHEM
        changeNO3p = krate(i,j,l,1,1) * y(n_HNO3,l) *dt2
#endif
#endif

C Apply Alkenes, AlkyNit, and Aldehyde changes here:
        y(n_Alkenes,L)  =y(n_Alkenes,L)  +changeAlkenes
        y(n_AlkylNit,L) =y(n_AlkylNit,L) +changeAlkylNit
        yAldehyde(I,J,L)=yAldehyde(I,J,L)+changeAldehyde

#ifdef TRACERS_AEROSOLS_SOA
        y(n_isopp1g,L)  =y(n_isopp1g,L)  +changeisopp1g
        y(n_isopp2g,L)  =y(n_isopp2g,L)  +changeisopp2g
#ifdef TRACERS_TERP
        y(n_apinp1g,L)  =y(n_apinp1g,L)  +changeapinp1g
        y(n_apinp2g,L)  =y(n_apinp2g,L)  +changeapinp2g
#endif  /* TRACERS_TERP */
#endif  /* TRACERS_AEROSOLS_SOA */

C Note: the lower limit of minKG placed on the resulting tracer mass
C from the following changes is to prevent negative tracer mass:

C -- HCHO --
c       Gas phase NO3 + HCHO -> HNO3 + CO yield of HCHO & CO:
        changeL(L,n_HCHO)=changeHCHO*pfactor*vol2mass(n_HCHO)
        if(-changeL(L,n_HCHO) > trm(I,J,L,n_HCHO))then
          changeL(L,n_HCHO)=-.95d0*trm(I,J,L,n_HCHO)
          changeHCHO=changeL(L,n_HCHO)*mass2vol(n_HCHO)*bypfactor
        endif
        IF((trm(i,j,l,n_HCHO)+changeL(l,n_HCHO)) < minKG) THEN
          changeL(l,n_HCHO) = minKG - trm(i,j,l,n_HCHO)
          changeHCHO=changeL(L,n_HCHO)*mass2vol(n_HCHO)*bypfactor
        ENDIF
        wprodHCHO=changeHCHO
C -- CO --
        changeL(L,n_CO)=rHCHOplusNO3*pfactor*vol2mass(n_CO)
        changeCO=changeL(L,n_CO)*mass2vol(n_CO)*bypfactor
        if((trm(i,j,l,n_CO)+changeL(l,n_CO)) < minKG)then
          changeL(l,n_CO) = minKG - trm(i,j,l,n_CO)
          changeCO=changeL(L,n_CO)*mass2vol(n_CO)*bypfactor
        endif
        wprodCO=rHCHOplusNO3   ! <---note
        if(changeL(L,n_CO) >= 0.) then  
          CALL INC_TAJLS(I,J,L,jls_COp,changeL(L,n_CO))
#ifdef HTAP_LIKE_DIAGS
          taijls(i,j,l,ijlt_COp)=taijls(i,j,l,ijlt_COp)+changeCO*cpd
#endif
        else
          CALL INC_TAJLS(I,J,L,jls_COd,changeL(L,n_CO))
#ifdef HTAP_LIKE_DIAGS
          taijls(i,j,l,ijlt_COd)=taijls(i,j,l,ijlt_COd)+changeCO*cpd
#endif
        endif       
C -- HNO3 --  (HNO3 from gas and het phase rxns )
        changeL(L,n_HNO3)=changeHNO3*pfactor*vol2mass(n_HNO3)
        IF((trm(i,j,l,n_HNO3)+changeL(l,n_HNO3)) < minKG) THEN
          changeL(l,n_HNO3) = minKG - trm(i,j,l,n_HNO3)
          changeHNO3=changeL(L,n_HNO3)*mass2vol(n_HNO3)*bypfactor
        END IF
#ifdef TRACERS_HETCHEM
#ifdef TRACERS_NITRATE
        changeL(L,n_N_d1)=changeN_d1*pfactor*vol2mass(n_N_d1)
!       if(i==36.and.j==28.and.l==1) then
!         write(out_line,*)'Mchange L 2 ', changeL(L,n_N_d1),changeN_d1
!         call write_parallel(trim(out_line),crit=.true.)
!       endif
        IF((trm(i,j,l,n_N_d1)+changeL(l,n_N_d1)) < minKG) THEN
          changeL(l,n_N_d1) = minKG - trm(i,j,l,n_N_d1)
          changeN_d1=changeL(L,n_N_d1)*mass2vol(n_N_d1)*bypfactor
        END IF
        changeL(L,n_N_d2)=changeN_d2*pfactor*vol2mass(n_N_d2)
        IF((trm(i,j,l,n_N_d2)+changeL(l,n_N_d2)) < minKG) THEN
          changeL(l,n_N_d2) = minKG - trm(i,j,l,n_N_d2)
          changeN_d2=changeL(L,n_N_d2)*mass2vol(n_N_d2)*bypfactor
        END IF
        changeL(L,n_N_d3)=changeN_d3*pfactor*vol2mass(n_N_d3)
        IF((trm(i,j,l,n_N_d3)+changeL(l,n_N_d3)) < minKG) THEN
          changeL(l,n_N_d3) = minKG - trm(i,j,l,n_N_d3)
          changeN_d3=changeL(L,n_N_d3)*mass2vol(n_N_d3)*bypfactor
        END IF
#endif  /* TRACERS_NITRATE */
#endif  /* TRACERS_HETCHEM */

#ifdef TOMAS_HETCHEM
        changeL(L,n_NO3p)=changeNO3p*pfactor*vol2mass(n_NO3p)
        IF((trm(i,j,l,n_NO3p)+changeL(l,n_NO3p)) < minKG) THEN
          changeL(l,n_NO3p) = minKG - trm(i,j,l,n_NO3p)
          changeNO3p=changeL(L,n_NO3p)*mass2vol(n_NO3p)*bypfactor
        END IF

#endif
C -- N2O5 --  (N2O5 from gas and het phase rxns)
        changeL(L,n_N2O5)=changeN2O5*pfactor*vol2mass(n_N2O5)
        IF((trm(i,j,l,n_N2O5)+changeL(l,n_N2O5)) < minKG) THEN
          changeL(l,n_N2O5) = minKG - trm(i,j,l,n_N2O5)
          changeN2O5=changeL(L,n_N2O5)*mass2vol(n_N2O5)*bypfactor
        END IF
c -- NOx --   (NOx from gas phase rxns)
        changeL(L,n_NOx)=changeNOx*pfactor*vol2mass(n_NOx)
        IF((trm(i,j,l,n_NOx)+changeL(l,n_NOx)) < minKG) THEN
          changeL(l,n_NOx) = minKG - trm(i,j,l,n_NOx)
          changeNOx=changeL(L,n_NOx)*mass2vol(n_NOx)*bypfactor
        END IF
C -- Alkenes --  (Alkenes from gas phase rxns)
        changeL(L,n_Alkenes)=
     &  changeAlkenes*pfactor*vol2mass(n_Alkenes)
        IF((trm(i,j,l,n_Alkenes)+changeL(l,n_Alkenes)) < minKG)THEN
          changeL(l,n_Alkenes) = minKG - trm(i,j,l,n_Alkenes)
          changeAlkenes=changeL(L,n_Alkenes)*mass2vol(n_Alkenes)
     &    *bypfactor
        END IF
#ifdef TRACERS_AEROSOLS_SOA
C -- isopp1g --  (isopp1g from gas phase rxns)
        changeL(L,n_isopp1g)=
     &  changeisopp1g*pfactor*vol2mass(n_isopp1g)
        IF((trm(i,j,l,n_isopp1g)+changeL(l,n_isopp1g)) < minKG)THEN
          changeL(l,n_isopp1g) = minKG - trm(i,j,l,n_isopp1g)
          changeisopp1g=changeL(L,n_isopp1g)*mass2vol(n_isopp1g)
     &    *bypfactor
        END IF
C -- isopp2g --  (isopp2g from gas phase rxns)
        changeL(L,n_isopp2g)=
     &  changeisopp2g*pfactor*vol2mass(n_isopp2g)
        IF((trm(i,j,l,n_isopp2g)+changeL(l,n_isopp2g)) < minKG)THEN
          changeL(l,n_isopp2g) = minKG - trm(i,j,l,n_isopp2g)
          changeisopp2g=changeL(L,n_isopp2g)*mass2vol(n_isopp2g)
     &    *bypfactor
        END IF
#ifdef TRACERS_TERP
C -- apinp1g --  (apinp1g from gas phase rxns)
        changeL(L,n_apinp1g)=
     &  changeapinp1g*pfactor*vol2mass(n_apinp1g)
        IF((trm(i,j,l,n_apinp1g)+changeL(l,n_apinp1g)) < minKG)THEN
          changeL(l,n_apinp1g) = minKG - trm(i,j,l,n_apinp1g)
          changeapinp1g=changeL(L,n_apinp1g)*mass2vol(n_apinp1g)
     &    *bypfactor
        END IF
C -- apinp2g --  (apinp2g from gas phase rxns)
        changeL(L,n_apinp2g)=
     &  changeapinp2g*pfactor*vol2mass(n_apinp2g)
        IF((trm(i,j,l,n_apinp2g)+changeL(l,n_apinp2g)) < minKG)THEN
          changeL(l,n_apinp2g) = minKG - trm(i,j,l,n_apinp2g)
          changeapinp2g=changeL(L,n_apinp2g)*mass2vol(n_apinp2g)
     &    *bypfactor
        END IF
#endif  /* TRACERS_TERP */
#endif  /* TRACERS_AEROSOLS_SOA */
c -- Isoprene -- (Isoprene from gas phase rxns)
        changeL(L,n_Isoprene)=
     &  changeIsoprene*pfactor*vol2mass(n_Isoprene)
        IF((trm(i,j,l,n_Isoprene)+changeL(l,n_Isoprene)) < minKG)
     &  THEN
          changeL(l,n_Isoprene) = minKG - trm(i,j,l,n_Isoprene)
          changeIsoprene=changeL(L,n_Isoprene)*mass2vol(n_Isoprene)
     &    *bypfactor
        END IF
#ifdef TRACERS_TERP
c -- Terpenes -- (Terpenes from gas phase rxns)
        changeL(L,n_Terpenes)=
     &  changeTerpenes*pfactor*vol2mass(n_Terpenes)
        IF((trm(i,j,l,n_Terpenes)+changeL(l,n_Terpenes)) < minKG)
     &  THEN
          changeL(l,n_Terpenes) = minKG - trm(i,j,l,n_Terpenes)
          changeTerpenes=changeL(L,n_Terpenes)*mass2vol(n_Terpenes)
     &    *bypfactor
        END IF
#endif  /* TRACERS_TERP */
c -- AlkylNit -- (AlkylNit from gas phase rxns)
        changeL(L,n_AlkylNit)=
     &  changeAlkylNit*pfactor*vol2mass(n_AlkylNit)
        IF((trm(i,j,l,n_AlkylNit)+changeL(l,n_AlkylNit)) < minKG)
     &  THEN
          changeL(l,n_AlkylNit) = minKG - trm(i,j,l,n_AlkylNit)
          changeAlkylNit=changeL(L,n_AlkylNit)*mass2vol(n_AlkylNit)
     &    *bypfactor
        END IF

C Save 3D radical arrays to pass to aerosol code:
C Make sure we get the nightime values; Set OH to zero for now:
        if(coupled_chem == 1) then
          oh_live(i,j,l)=0.d0
          no3_live(i,j,l)=yNO3(i,j,l)
        endif

#ifdef SHINDELL_STRAT_CHEM
c --  Ox --   ( Ox from gas phase rxns)
        changeOx=-1.d0*rr(7,L)*y(nNO2,L)*y(n_Ox,L)*pOx(I,J,L)*dt2
        changeL(L,n_Ox)=changeOx*pfactor*vol2mass(n_Ox)
        IF((trm(i,j,l,n_Ox)+changeL(l,n_Ox)) < minKG) THEN
          changeL(l,n_Ox) = minKG - trm(i,j,l,n_Ox)
          changeOx=changeL(L,n_Ox)*mass2vol(n_Ox)*bypfactor
        END IF
#if (defined SHINDELL_STRAT_EXTRA) && (defined ACCMIP_LIKE_DIAGS)
        if(trm(i,j,l,n_Ox)==0.)call stop_model('zero Ox denom',255)
        changeL(L,n_stratOx)=changeL(L,n_Ox)*
     &  trm(i,j,l,n_stratOx)/trm(i,j,l,n_Ox)
        if((trm(i,j,l,n_stratOx)+changeL(l,n_stratOx)) < minKG)
     &  changeL(l,n_stratOx) = minKG - trm(i,j,l,n_stratOx)
#endif
        ! then come diags:
        if(changeL(L,n_Ox) >= 0.) then
          CALL INC_TAJLS(I,J,L,jls_Oxp,changeL(L,n_Ox))
          if(L<=maxl)CALL INC_TAJLS(I,J,L,jls_OxpT,changeL(L,n_Ox))
#ifdef HTAP_LIKE_DIAGS
          taijls(i,j,l,ijlt_Oxp)=taijls(i,j,l,ijlt_Oxp)+changeOx*cpd
#endif
        else
          CALL INC_TAJLS(I,J,L,jls_Oxd,changeL(L,n_Ox))
          if(L<=maxl)CALL INC_TAJLS(I,J,L,jls_OxpT,changeL(L,n_Ox))
#ifdef HTAP_LIKE_DIAGS
          taijls(i,j,l,ijlt_Oxd)=taijls(i,j,l,ijlt_Oxd)+changeOx*cpd
#endif
        endif
c -- ClONO2 --   (ClONO2 from gas and het phase rxns)
        changeL(L,n_ClONO2)=changeClONO2*pfactor*
     &  vol2mass(n_ClONO2)
        IF((trm(i,j,l,n_ClONO2)+changeL(l,n_ClONO2)) < minKG) THEN
          changeL(l,n_ClONO2) = minKG - trm(i,j,l,n_ClONO2)
          changeClONO2=changeL(L,n_ClONO2)*mass2vol(n_ClONO2)*
     &    bypfactor
        END IF
c -- ClOx --   (ClOx from gas and het phase rxns)
        changeL(L,n_ClOx)=changeClOx*pfactor*vol2mass(n_ClOx)
        IF((trm(i,j,l,n_ClOx)+changeL(l,n_ClOx)) < minKG) THEN
          changeL(l,n_ClOx) = minKG - trm(i,j,l,n_ClOx)
          changeClOx=changeL(L,n_ClOx)*mass2vol(n_ClOx)*bypfactor
        END IF
        if(pscX(L))then
c -- HOCl --   (HOCl from het phase rxns)
          changeL(L,n_HOCl)=changeHOCl*pfactor*vol2mass(n_HOCl)
          IF((trm(i,j,l,n_HOCl)+changeL(l,n_HOCl)) < minKG) THEN
            changeL(l,n_HOCl) = minKG - trm(i,j,l,n_HOCl)
            changeHOCl=changeL(L,n_HOCl)*mass2vol(n_HOCl)*bypfactor
          END IF
c -- HCl --   (HCl from het phase rxns)
          changeL(L,n_HCl)=changeHCl*pfactor*vol2mass(n_HCl)
          IF((trm(i,j,l,n_HCl)+changeL(l,n_HCl)) < minKG) THEN
            changeL(l,n_HCl) = minKG - trm(i,j,l,n_HCl)
            changeHCl=changeL(L,n_HCl)*mass2vol(n_HCl)*bypfactor
          END IF
c -- HBr --   (HBr from gas phase rxns)
          changeL(L,n_HBr)=changeHBr*pfactor*vol2mass(n_HBr)
          IF((trm(i,j,l,n_HBr)+changeL(l,n_HBr)) < minKG) THEN
            changeL(l,n_HBr) = minKG - trm(i,j,l,n_HBr)
            changeHBr=changeL(L,n_HBr)*mass2vol(n_HBr)*bypfactor
          END IF
c -- BrOx --   (BrOx from gas phase rxns)
          changeL(L,n_BrOx)=changeBrOx*pfactor*vol2mass(n_BrOx)
          IF((trm(i,j,l,n_BrOx)+changeL(l,n_BrOx)) < minKG) THEN
            changeL(l,n_BrOx) = minKG - trm(i,j,l,n_BrOx)
            changeBrOx=changeL(L,n_BrOx)*mass2vol(n_BrOx)*bypfactor
          END IF
c -- BrONO2 --   (BrONO2 from gas phase rxns)
          changeL(L,n_BrONO2)=changeBrONO2*pfactor*vol2mass(n_BrONO2)
          IF((trm(i,j,l,n_BrONO2)+changeL(l,n_BrONO2)) < minKG) THEN
            changeL(l,n_BrONO2) = minKG - trm(i,j,l,n_BrONO2)
           changeBrONO2=changeL(L,n_BrONO2)*mass2vol(n_BrONO2)*bypfactor
          END IF
        endif  ! PSCs exist
#endif

        call printNightChemistryDiags()
        call checkNighttimeTolerances()

C       ACCUMULATE 3D NO3 diagnostic: 
        if (yNO3(I,J,L) > 0.d0 .and. yNO3(I,J,L) < 1.d20)
     &  taijls(i,j,l,ijlt_NO3)=taijls(i,j,l,ijlt_NO3)+yNO3(i,j,l)

#ifdef SHINDELL_STRAT_CHEM
        if (y(nClO,L) > 0.d0 .and. y(nClO,L) < 1.d20)
     &  CALL INC_TAJLS2(I,J,L,jls_ClOcon,y(nClO,L)/y(nM,L))
        if (y(nH2O,L) > 0.d0 .and. y(nH2O,L) < 1.d20)
     &  CALL INC_TAJLS2(I,J,L,jls_H2Ocon,y(nH2O,L)/y(nM,L))
#endif
     
       enddo  ! L loop

CCCCCCCCCCCCCCCC END NIGHTTIME CCCCCCCCCCCCCCCCCCCC

      endif
CCCCCCCCCCCCCCCCCCCC END DARKNESS CCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      save_NO2column(i,j)=0.d0 ! initialize sum outside L loop.

#ifdef SHINDELL_STRAT_CHEM
      LL=LM
#else
      LL=maxl
#endif
      DO L=1,LL  

C Lower limit on HO2NO2 : 
        if(trm(i,j,l,n_HO2NO2)+changeL(l,n_HO2NO2) < minKG)
     &  changeL(l,n_HO2NO2) = minKG - trm(i,j,l,n_HO2NO2)

#ifdef SHINDELL_STRAT_CHEM
c Tropospheric halogen sink Br & Cl :
        if(y(nH2O,L)/y(nM,L) > 10.d-6)then ! sink by wet removal in trop
          rmv=0.5d0
          changeL(L,n_ClOx)  =-trm(I,J,L,n_ClOx)  *rmv
          changeL(L,n_HCl)   =-trm(I,J,L,n_HCl)   *rmv
          changeL(L,n_HOCl)  =-trm(I,J,L,n_HOCl)  *rmv
          changeL(L,n_ClONO2)=-trm(I,J,L,n_ClONO2)*rmv
          changeL(L,n_BrOx)  =-trm(I,J,L,n_BrOx)  *rmv
          changeL(L,n_HBr)   =-trm(I,J,L,n_HBr)   *rmv
          changeL(L,n_HOBr)  =-trm(I,J,L,n_HOBr)  *rmv
          changeL(L,n_BrONO2)=-trm(I,J,L,n_BrONO2)*rmv
        else
c Set CLTOT based on CFCs (2.4 ppbv yield from complete oxidation of
c 1.8 ppbv CFC plus 0.8 ppbv background which is tied to methane) :
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          !WARNING: RESETTING SOME Y's HERE; SO DON'T USE THEM BELOW!     
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          y(n_ClOx,L)=(trm(I,J,L,n_ClOx)+changeL(L,n_ClOx))*y(nM,L)*
     &    mass2vol(n_ClOx)*BYAXYP(I,J)*BYAM(L,I,J)
          y(n_HCl,L)= (trm(I,J,L,n_HCl)+changeL(L,n_HCl))*y(nM,L)*
     &    mass2vol(n_HCl)*BYAXYP(I,J)*BYAM(L,I,J)
          y(n_HOCl,L)=(trm(I,J,L,n_HOCl)+changeL(L,n_HOCl))*y(nM,L)*
     &    mass2vol(n_HOCl)*BYAXYP(I,J)*BYAM(L,I,J)
          y(n_ClONO2,L)=(trm(I,J,L,n_ClONO2)+changeL(L,n_ClONO2))*
     &    y(nM,L)*mass2vol(n_ClONO2)*BYAXYP(I,J)*BYAM(L,I,J)
          CLTOT=((y(n_CFC,1)/y(nM,1)-y(n_CFC,L)/y(nM,L))*(3.0d0/1.8d0)*
     &    y(n_CFC,1)/(1.8d-9*y(nM,1)))
          CLTOT=CLTOT+0.8d-9*(y(n_CH4,1)/y(nM,1)-y(n_CH4,L)/y(nM,L))/
     &    (y(n_CH4,1)/y(nM,1))
          CLTOT=CLTOT*y(nM,L)/
     &    (y(n_ClOx,L)+y(n_HCl,L)+y(n_HOCl,L)+y(n_ClONO2,L))
          if(prnchg.and.J == jprn.and.I == iprn.and.L == lprn)then  
            write(out_line,66) CLTOT
            call write_parallel(trim(out_line),crit=jay)
 66         format ('CLTOT = ',F20.5)
          endif
          IF(CLTOT <= 0.999d0 .OR. CLTOT >= 1.001d0) THEN
            changeL(L,n_ClOx)=changeL(L,n_ClOx)*CLTOT+
     &      trm(I,J,L,n_ClOx)*(CLTOT-1.D0)
            changeL(L,n_HCl)=changeL(L,n_HCl)*CLTOT+
     &      trm(I,J,L,n_HCl)*(CLTOT-1.D0)
            changeL(L,n_HOCl)=changeL(L,n_HOCl)*CLTOT+
     &      trm(I,J,L,n_HOCl)*(CLTOT-1.D0)
c           Conserve N wrt ClONO2 once inital Cl changes past:
            if(Itime-ItimeI >= 6 .OR. allowSomeChemReinit .NE. 1)then ! note logic
              changeL(L,n_NOx)=changeL(L,n_NOx)-
     &        (trm(I,J,L,n_ClONO2)+changeL(L,n_ClONO2))*
     &        (CLTOT-1.D0)*tr_mm(n_NOx)/tr_mm(n_ClONO2)
              if(-changeL(L,n_NOx) > trm(I,J,L,n_NOx))changeL(L,n_NOx)=
     &        -0.8d0*trm(I,J,L,n_NOx)
            endif
            changeL(L,n_ClONO2)=changeL(L,n_ClONO2)*CLTOT+
     &      trm(I,J,L,n_ClONO2)*(CLTOT-1.D0)
          ENDIF

c Set Total Bromine based on CFCs (4.5 pptv yield
C from complete oxidation of 1.8 ppbv CFC plus 0.5 pptv background) :
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          !WARNING: RESETTING SOME Y's HERE; SO DON'T USE THEM BELOW!     
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          y(n_BrOx,L)=(trm(I,J,L,n_BrOx)+changeL(L,n_BrOx))*y(nM,L)*
     &    mass2vol(n_BrOx)*BYAXYP(I,J)*BYAM(L,I,J)
          y(n_HBr,L)= (trm(I,J,L,n_HBr)+changeL(L,n_HBr))*y(nM,L)*
     &    mass2vol(n_HBr)*BYAXYP(I,J)*BYAM(L,I,J)
          y(n_HOBr,L)=(trm(I,J,L,n_HOBr)+changeL(L,n_HOBr))*y(nM,L)*
     &    mass2vol(n_HOBr)*BYAXYP(I,J)*BYAM(L,I,J)
          y(n_BrONO2,L)=(trm(I,J,L,n_BrONO2)+changeL(L,n_BrONO2))*
     &    y(nM,L)*mass2vol(n_BrONO2)*BYAXYP(I,J)*BYAM(L,I,J)
     
          BRTOT=((y(n_CFC,1)/y(nM,1)-y(n_CFC,L)/y(nM,L))*(4.5d-3/1.8d0)
     &    *y(n_CFC,1)/(1.8d-9*y(nM,1)))
          BRTOT=BRTOT+0.5d-12*(y(n_CH4,1)/y(nM,1)-y(n_CH4,L)/y(nM,L))/
     &    (y(n_CH4,1)/y(nM,1))
          BRTOT=BRTOT*y(nM,L)/
     &    (y(n_BrOx,L)+y(n_HBr,L)+y(n_HOBr,L)+y(n_BrONO2,L))
          if(prnchg.and.J == jprn.and.I == iprn.and.L == lprn)then  
            write(out_line,67) BRTOT
            call write_parallel(trim(out_line),crit=jay)
 67         format ('BrTOT = ',F20.5)
          endif
          IF(BRTOT <= 0.999d0 .OR. BRTOT >= 1.001d0) THEN
            changeL(L,n_BrOx)=changeL(L,n_BrOx)*BRTOT+
     &      trm(I,J,L,n_BrOx)*(BRTOT-1.D0)
            changeL(L,n_HBr)=changeL(L,n_HBr)*BRTOT+
     &      trm(I,J,L,n_HBr)*(BRTOT-1.D0)
            changeL(L,n_HOBr)=changeL(L,n_HOBr)*BRTOT+
     &      trm(I,J,L,n_HOBr)*(BRTOT-1.D0)
c           Conserve N wrt BrONO2 once inital Br changes past:
            if(Itime-ItimeI >= 6 .OR. allowSomeChemReinit .NE. 1)then ! note logic
              changeL(L,n_NOx)=changeL(L,n_NOx)-
     &        (trm(I,J,L,n_BrONO2)+changeL(L,n_BrONO2))*
     &        (BRTOT-1.D0)*tr_mm(n_NOx)/tr_mm(n_BrONO2)
              if(-changeL(L,n_NOx) > trm(I,J,L,n_NOx))changeL(L,n_NOx)=
     &        -0.8d0*trm(I,J,L,n_NOx)
            endif
            changeL(L,n_BrONO2)=changeL(L,n_BrONO2)*BRTOT+
     &      trm(I,J,L,n_BrONO2)*(BRTOT-1.D0)
          ENDIF
        endif ! i.e. y(nH2O,L)/y(nM,L) <= 10.d-6 
#endif
#ifdef TRACERS_AEROSOLS_SOA
        pfactor=axyp(I,J)*AM(L,I,J)/y(nM,L)
        bypfactor=1.D0/pfactor
        call soa_aerosolphase(I,J,L,changeL,bypfactor)
#endif  /* TRACERS_AEROSOLS_SOA */

        tempChangeNOx= ! this needed for several diags below:
     &  changeL(L,n_NOx)*mass2vol(n_NOx)*y(nM,L)/(axyp(I,J)*AM(L,I,J))

        tempChangeOx=
     &  changeL(L,n_Ox)*mass2vol(n_Ox)*y(nM,L)/(axyp(I,J)*AM(L,I,J))

! Accumulate NO2 10:30am/1:30pm tropo column diags:
! -- moved from sunlight/darkness sections because needed changeNOx
! -- saved here in molecules/cm2
        if(L<=min(maxl,LTROPO(I,J)))then
          if((ALB(I,J,1) /= 0.d0).AND.(sza < szamax))then

            index1=0 ; index2=0

            if(ih1030 < ih1030e)then ! normal case
              if(i>=ih1030.and.i<=ih1030e)then
                index1=ijs_NO2_1030; index2=ijs_NO2_1030c
              end if 
            else                     ! crossing date line
              if(i<=ih1030.or.i>=ih1030e)then
                index1=ijs_NO2_1030; index2=ijs_NO2_1030c
              end if 
            end if
            if(ih1330 < ih1330e)then ! normal case
              if(i>=ih1330.and.i<=ih1330e)then
                index1=ijs_NO2_1330; index2=ijs_NO2_1330c
              end if
            else                     ! crossing date line
              if(i<=ih1330.or.i>=ih1330e)then
                index1=ijs_NO2_1330; index2=ijs_NO2_1330c
              end if
            end if

            if(index1/=0 .and. index2/=0)then
              thick= ! layer thickness in cm
     &        1.d2*rgas*bygrav*TX(I,J,L)*LOG(PEDN(L,i,j)/PEDN(L+1,i,j))
              taijs(i,j,index1)=taijs(i,j,index1)+thick*
     &        pNOx(i,j,L)*(y(n_NOx,L)+tempChangeNOx)
              if(L==1)taijs(i,j,index2)=taijs(i,j,index2)+1.d0
            end if

          end if ! sunlight criteria
        end if ! troposphere criterion

! Also save instantaneous NO2 tropospheric column for SUBDDiag:
! Conversion is only from molecules/cm3 to molecules/cm2:
! save_NO2column is initialized to 0 outside this L loop.
! [note: we should consolodate all these "thick/byThick" guys.]
        if(L<=min(maxl,LTROPO(I,J)))then
          thick= ! layer thickness in cm
     &    1.d2*rgas*bygrav*TX(I,J,L)*LOG(PEDN(L,i,j)/PEDN(L+1,i,j))
          save_NO2column(i,j) = save_NO2column(i,j)+
     &    thick*pNOx(i,j,L)*(y(n_NOx,L)+tempChangeNOx)
        end if

#ifdef ACCMIP_LIKE_DIAGS
! accumulate some 3D diagnostics in moles/m3/s units:
        ! chemical_production_of_O1D_from_ozone:
        taijls(i,j,l,ijlt_pO1D)=taijls(i,j,l,ijlt_pO1D)+
     &  ss(2,l,i,j)*y(nO3,l)*cpd

        ! chemical_production_of_OH_from_O1D_plus_H2O:
        taijls(i,j,l,ijlt_pOH)=taijls(i,j,l,ijlt_pOH)+
     &  2.d0*rr(10,l)*y(nH2O,l)*y(nO1D,l)*cpd

        ! chemical_production_rate_of_ozone_by_HO2_plus_NO:
        taijls(i,j,l,ijlt_OxpHO2)=taijls(i,j,l,ijlt_OxpHO2)+
     &  rr(6,l)*y(nHO2,l)*y(nNO,l)*cpd
   
        ! chemical_production_rate_of_ozone_by_CH3O2_plus_NO:
        taijls(i,j,l,ijlt_OxpCH3O2)=taijls(i,j,l,ijlt_OxpCH3O2)+
     &  rr(20,l)*y(nCH3O2,l)*y(nNO,l)*cpd
    
        ! chemical_destruction_rate_of_ozone_by_OH:
        taijls(i,j,l,ijlt_OxlOH)=taijls(i,j,l,ijlt_OxlOH)+
     &  rr(2,l)*y(nOH,l)*y(nO3,l)*cpd ! (positive)

        !chemical_destruction_rate_of_ozone_by_HO2:
        taijls(i,j,l,ijlt_OxlHO2)=taijls(i,j,l,ijlt_OxlHO2)+
     &  rr(4,l)*y(nOH,l)*y(nO3,l)*cpd ! (positive)

        !chemical_destruction_rate_of_ozone_by_Alkenes:
        taijls(i,j,l,ijlt_OxlALK)=taijls(i,j,l,ijlt_OxlALK)+
     &  rr(35,l)*y(n_Alkenes,l)*y(nO3,l)*cpd ! (positive)

        !Save 3D NO2 separately from NOx (pppv here):
        ! need to add NOx change to match the NOx tracer diag:
        taijls(i,j,l,ijlt_NO2vmr)=taijls(i,j,l,ijlt_NO2vmr)+
     &  pNOx(i,j,l)*(y(n_NOx,l)+tempChangeNOx)/y(nM,l)

        !Save 3D NO separately from NOx (pppv here):
        ! need to add NOx change to match the NOx tracer diag:
        taijls(i,j,l,ijlt_NOvmr)=taijls(i,j,l,ijlt_NOvmr)+
     &  (1.d0-pNOx(i,j,l))*(y(n_NOx,l)+tempChangeNOx)/y(nM,l)
#endif

        ! Below there is a 3D O3 diagnostic in cm-atm units for more
        ! direct NINT input. Here try to save it in vmr(ppbv) for humans,
        ! (both JL and IJL).
        taijls(i,j,L,ijlt_O3ppbv)=taijls(i,j,L,ijlt_O3ppbv)+
     &  1.e9*pOx(i,j,L)*(y(n_Ox,L)+tempChangeOx)/y(nM,L)
        CALL INC_TAJLS2  ! (V/V air)
     &  (I,J,L,jls_O3vmr,pOx(i,j,L)*(y(n_Ox,L)+tempChangeOx)/y(nM,L))

#ifdef SHINDELL_STRAT_CHEM
      call printSS27x2Etc()
#endif
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Save chemistry changes for applying in apply_tracer_3Dsource.  C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

        DO N=1,NTM_CHEM
          tr3Dsource(i,j,l,nChemistry,n) = changeL(l,n) * bydtsrc
        END DO
#if (defined SHINDELL_STRAT_EXTRA) && (defined ACCMIP_LIKE_DIAGS)
        tr3Dsource(i,j,l,nChemistry,n_stratOx)=
     &  changeL(l,n_stratOx)*bydtsrc
#endif

        ! save NO2 volume mixing ratio for sub-daily diagnosic:
        mNO2(i,j,L)=pNOx(i,j,L)*y(n_NOx,L)/y(nM,L)

#ifdef CACHED_SUBDD
        ! above line with mNO2 should really be changed to include 
        ! +tempChangeNOx, but I am leaving it for backwards compatability. 
        ! However, for newer 3D SUBDDiags, better make it match NOx tracer:
        mrno2(i,j,L)=pNOx(i,j,L)*(y(n_NOx,L)+tempChangeNOx)/y(nM,L)
        mrno(i,j,L)=(1d0-pNOx(i,j,L))*(y(n_NOx,L)+tempChangeNOx)/y(nM,L)
#endif
     
#ifdef TRACERS_HETCHEM
#ifdef TRACERS_NITRATE
        tr3Dsource(i,j,l,nChemistry,n_N_d1) = changeL(l,n_N_d1) *bydtsrc
        tr3Dsource(i,j,l,nChemistry,n_N_d2) = changeL(l,n_N_d2) *bydtsrc
        tr3Dsource(i,j,l,nChemistry,n_N_d3) = changeL(l,n_N_d3) *bydtsrc
#endif  /* TRACERS_NITRATE */
#endif  /* TRACERS_HETCHEM */

#ifdef TOMAS_HETCHEM
! Save nitrate formed from HETCHEM to NO3p for now.
        tr3Dsource(i,j,l,nOther,n_NO3p) = changeL(l,n_NO3p) *bydtsrc
#endif

      END DO ! end current altitude loop

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      END DO i_loop ! ---> MAIN I LOOP ENDS <---

      END DO j_loop ! ---> MAIN J LOOP ENDS <---
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCC!$OMP END PARALLEL DO  

      ! check if there was that error in certain section of chemstep
      ! anywhere in the world; if so, stop the model (all processors):
      
      ! Currently the section where ierr could become non-zero is 
      ! commented. If it remains so, remove it from calls/sums/checks
      ! in this program:
      call globalmax(grid,ierr_loc,ierr)
      if(ierr > 0) then ! all processors call stop_model
        if(am_i_root()) write(6,*) 'chemstep Oxcorr fault'  
        call stop_model('chemstep Oxcorr fault',255)
      endif

#ifdef CACHED_SUBDD
      call find_groups('taijlh',grpids,ngroups)
      do igrp=1,ngroups
        subdd => subdd_groups(grpids(igrp))
        do k=1,subdd%ndiags
          select case (subdd%name(k))
          case ('MRNO2')
            call inc_subdd(subdd,k,mrno2)
          case ('MRNO')
            call inc_subdd(subdd,k,mrno)
          end select
        enddo ! k
      enddo ! igroup

      call find_groups('taijph',grpids,ngroups)
      do igrp=1,ngroups
        subdd => subdd_groups(grpids(igrp))
        do k=1,subdd%ndiags
          select case (subdd%name(k))
          case ('MRNO2cp')
            call inc_subdd(subdd,k,mrno2)
          case ('MRNOcp')
            call inc_subdd(subdd,k,mrno)
          end select
        enddo ! k
      enddo ! igroup
#endif

CCCCCCCCCCCCCCCCCC END CHEMISTRY SECTION CCCCCCCCCCCCCCCCCCCCCCCCC

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                 BEGIN OVERWRITING                              C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C If fix_CH4_chemistry is turned on, reset the CH4 tracer everywhere
C to initial conditions and set the chemistry change to zero...
      if(fix_CH4_chemistry == 1)then
        tr3Dsource(:,:,:,nChemistry,n_CH4) = 0.d0 
        call get_CH4_IC(1)
      end if 

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Special cases of overwriting, when doing stratospheric chemistry C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC     
#ifdef SHINDELL_STRAT_CHEM      
C N2O, CFC, and optional CH4 L=1 overwriting: with all these "fact"s
C this looks complicated, but basically, you are either converting
C from mixing ratio to KG (normal case) or from cm*atm to KG  
C (interactive radiation case - for more on that conversion, see
C the notes on O3MULT in the TRCHEM_Shindell_COM program):
      PIfact(:)=1.d0     
      if(PI_run == 1)then
        PIfact(n_NOx)=PIratio_N
        if(use_rad_n2o == 0) PIfact(n_N2O)=PIratio_N2O
        if(use_rad_cfc == 0) PIfact(n_CFC)=PIratio_CFC
      endif
      fact2=n2o_pppv  ! default N2O mixing ratio overwrite
      fact3=cfc_pppv  ! default CFC mixing ratio overwrite
      fact7=fact_cfc
      if(use_rad_cfc == 0)fact7=1.d0
      do j=J_0,J_1
       do i=I_0,IMAXJ(j)
        fact6=2.69d20*axyp(i,j)*byavog
        fact1=bymair*am(1,i,j)*axyp(i,j)
        fact5=fact6 
        fact4=fact6
        if(use_rad_n2o == 0)fact4=fact1 
        if(use_rad_cfc == 0)fact5=fact1
        if(use_rad_n2o > 0)fact2=rad_to_chem(3,1,i,j)
        if(use_rad_cfc > 0)fact3=rad_to_chem(5,1,i,j)
        tr3Dsource(i,j,1,nOverwrite,n_N2O)=(fact2*fact4*
     &  tr_mm(n_N2O)*PIfact(n_N2O) - (trm(i,j,1,n_N2O)+ 
     &  tr3Dsource(i,j,1,nChemistry,n_N2O)*dtsrc))*bydtsrc
        tr3Dsource(i,j,1,nOverwrite,n_CFC)=(fact3*fact5*fact7*
     &  tr_mm(n_CFC)*PIfact(n_CFC) - (trm(i,j,1,n_CFC)+
     &  tr3Dsource(i,j,1,nChemistry,n_CFC)*dtsrc))*bydtsrc
        if(use_rad_ch4 > 0)then
          tr3Dsource(i,j,1,nOverwrite,n_CH4)=(rad_to_chem(4,1,i,j)*
     &    fact6*tr_mm(n_CH4)-(trm(i,j,1,n_CH4)+
     &    tr3Dsource(i,j,1,nChemistry,n_CH4)*dtsrc))*bydtsrc
        endif
       end do
      end do

! Ox, stratOx, NOx, BrOx and ClOx, have overwriting where P<PltOx:
      !(Interpolate BrOx & ClOx altitude-dependence to model resolution)
      CALL LOGPINT(LCOalt,PCOalt,BrOxaltIN,LM,PRES2,BrOxalt,.true.)
      CALL LOGPINT(LCOalt,PCOalt,ClOxaltIN,LM,PRES2,ClOxalt,.true.)
      do L=LS1,LM
       if(pres2(L) < pltOx)then
        do j=J_0,J_1         
          do i=I_0,IMAXJ(j)
            ! -- Ox --
            tr3Dsource(i,j,L,nOverwrite,n_Ox)=(rad_to_chem(1,L,i,j)*
     &      axyp(i,j)*O3MULT - (trm(i,j,L,n_Ox)+
     &      tr3Dsource(i,j,L,nChemistry,n_Ox)*dtsrc))*bydtsrc
#if (defined SHINDELL_STRAT_EXTRA) && (defined ACCMIP_LIKE_DIAGS)
            ! -- stratOx --
            tr3Dsource(i,j,L,nOverwrite,n_stratOx)=
     &      (rad_to_chem(1,L,i,j)*axyp(i,j)*O3MULT - (
     &      trm(i,j,L,n_stratOx)+tr3Dsource(i,j,L,nChemistry,n_stratOx)
     &      *dtsrc))*bydtsrc
#endif
            ! -- ClOx --
            tr3Dsource(i,j,L,nOverwrite,n_ClOx)=(1.d-11*ClOxalt(l)
     &      *vol2mass(n_CLOx)*am(L,i,j)*axyp(i,j) - (
     &      trm(i,j,L,n_ClOx)+tr3Dsource(i,j,L,nChemistry,n_ClOx)*dtsrc
     &      ))*bydtsrc    
            ! -- BrOx --
            tr3Dsource(i,j,L,nOverwrite,n_BrOx)=(1.d-11*BrOxalt(l)
     &      *vol2mass(n_BrOx)*am(L,i,j)*axyp(i,j) - (
     &      trm(i,j,L,n_BrOx)+tr3Dsource(i,j,L,nChemistry,n_BrOx)*dtsrc
     &      ))*bydtsrc
            ! -- NOx --
            tr3Dsource(i,j,L,nOverwrite,n_NOx)=(75.d-11 !75=1*300*2.5*.1
     &      *am(L,i,j)*axyp(i,j)*PIfact(n_NOx)-(trm(i,j,L,n_NOx)+ 
     &      tr3Dsource(i,j,L,nChemistry,n_NOx)*dtsrc))*bydtsrc
          end do ! I 
        end do   ! J
       end if    ! pressure
      end do     ! L

#else

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C No stratospheric chemistry; overwrite all tracers in stratosphere C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C W A R N I N G : This section should never be used if there is
C chemistry done in the stratosphere, because the stratospheric
C changes below assume that the changeL variable is zero for L>maxl
C at this point in the code. Not the case if chemistry was done.
C To put it another way, the overwritings below are explicitly 
C functions of tracer mass UNCHANGED by chemistry !
#ifdef SHINDELL_STRAT_EXTRA
      call stop_model("SHINDELL_STRAT_EXTRA w/o strat chem",255)
#endif

C determine scaling factors, if any:
      PIfact(:)=1.d0
      if(PI_run == 1) then
        do N=1,NTM
          select case(trname(n))
          case('NOx','HNO3','N2O5','HO2NO2')
            PIfact(n)=PIratio_N
          case('CO')
            PIfact(n)=PIratio_CO_S
          case('PAN','Isoprene','AlkylNit','Alkenes','Paraffin'
#ifdef TRACERS_TERP
     &        ,'Terpenes'
#endif  /* TRACERS_TERP */
     &        )
            PIfact(n)=PIratio_other
          end select
        end do
      endif

C Calculate an average tropical CH4 value near 569 hPa::
      CH4_569_part(I_0:I_1,J_0:J_1)=0.d0   
      count_569_part(I_0:I_1,J_0:J_1)=0.d0
      DO J=J_0,J_1
        do I=I_0,IMAXJ(J)
        if(LAT2D_DG(I,J) >= -30. .and. LAT2D_DG(I,J) <= 30.)then
            count_569_part(I,J)=1.d0
            CH4_569_part(I,J)=1.d6*byaxyp(i,j)*
     &      (F569M*trm(I,J,L569M,n_CH4)*byam(L569M,I,J)+
     &      F569P*trm(I,J,L569P,n_CH4)*byam(L569P,I,J))
        end if
        end do
      END DO
      CALL GLOBALSUM(grid,  CH4_569_part,  CH4_569, all=.true.)
      CALL GLOBALSUM(grid,count_569_part,count_569, all=.true.)
      if(count_569 <= 0.)call stop_model('count_569.le.0',255)
      CH4_569 = CH4_569 / count_569

      r179=1.d0/1.79d0 ! 1.79 is observed trop. CH4

      do j=J_0,J_1
       do i=I_0,IMAXJ(J)
         select case(which_trop)
         case(0); maxl=ltropo(I,J)
         case(1); maxl=ls1-1
         case default; call stop_model('which_trop problem 5',255)
         end select
         changeL(:,:)=0.d0 ! initilize the change (LM,NTM)
         
         do L=maxl+1,LM    ! => BEGIN LOOP OVER STRATOSPHERE <=
c         Update stratospheric ozone to amount set in radiation:
          changeL(L,n_Ox)=rad_to_chem(1,L,I,J)*AXYP(I,J)*O3MULT
     &                                         - trm(I,J,L,n_Ox)
          byam75=F75P*byam(L75P,I,J)+F75M*byam(L75M,I,J)
          FACT1=2.0d-9*AXYP(I,J)*am(L,I,J)*byam75
C         We think we have too little stratospheric NOx, so, to
C         increase the flux into the troposphere, increasing
C         previous stratospheric value by 70% here: GSF/DTS 9.15.03: 
          changeL(L,n_NOx)=trm(I,J,L,n_Ox)*2.3d-4*1.7d0*PIfact(n_NOx)
     &                                         - trm(I,J,L,n_NOx)
          changeL(L,n_N2O5)=  FACT1*PIfact(n_N2O5)- trm(I,J,L,n_N2O5)
          changeL(L,n_HNO3)=trm(I,J,L,n_Ox)*1.d-3*PIfact(n_HNO3)
     &                                         - trm(I,J,L,n_HNO3)
          changeL(L,n_H2O2)=  FACT1            - trm(I,J,L,n_H2O2)
          changeL(L,n_CH3OOH)=FACT1            - trm(I,J,L,n_CH3OOH)
          changeL(L,n_HCHO)=  FACT1            - trm(I,J,L,n_HCHO)
          ! here 70. = 1.4E-7/2.0E-9
          changeL(L,n_HO2NO2)=FACT1*70.d0*PIfact(n_HO2NO2)
     &                                         - trm(I,J,L,n_HO2NO2)
          ! Don't know if we still want this 40% hardcode in here:
          changeL(L,n_CO)=COIC(I,J,L)*0.4d0*PIfact(n_CO)-trm(I,J,L,n_CO)
          changeL(L,n_PAN)= FACT1*1.d-4*PIfact(n_PAN)
     &                           - trm(I,J,L,n_PAN)
          changeL(L,n_Isoprene)= FACT1*1.d-4*PIfact(n_Isoprene)
     &                           - trm(I,J,L,n_Isoprene)
          changeL(L,n_AlkylNit)= FACT1*1.d-4*PIfact(n_AlkylNit)
     &                           - trm(I,J,L,n_AlkylNit)
          changeL(L,n_Alkenes) = FACT1*1.d-4*PIfact(n_Alkenes)
     &                           - trm(I,J,L,n_Alkenes)
          changeL(L,n_Paraffin)= FACT1*1.d-4*PIfact(n_Paraffin)
     &                           - trm(I,J,L,n_Paraffin)
#ifdef TRACERS_TERP
!kt Terpenes start from zero. Some other tracers might need that too.
          changeL(L,n_Terpenes)= 0.d0*FACT1*1.d-4*PIfact(n_Terpenes)
     &                           - trm(I,J,L,n_Terpenes)
#endif  /* TRACERS_TERP */

c Overwrite stratospheric ch4 based on HALOE obs for tropics and
C extratropics and scale by the ratio of near-569hPa mixing ratios
C to 1.79:
          if (fix_CH4_chemistry == 0) then ! -------------------------
          CH4FACT=CH4_569*r179
          IF(ABS(LAT2D_DG(I,J)) > 30.) THEN ! extratropics
            DO L2=L,maxl+1,-1
              IF(CH4altX(L2) /= 0.d0) THEN
                CH4FACT=CH4FACT*CH4altX(L2)
                EXIT
              END IF
            END DO
          ELSE                              ! tropics
            DO L2=L,maxl+1,-1
              IF(CH4altT(L2) /= 0.d0) THEN
                CH4FACT=CH4FACT*CH4altT(L2)
                EXIT
              END IF
            END DO
          END IF
          select case(PI_run)
          case(1) ! yes, use scalings
            if(lat2d(i,j).lt.0.)then; PIfact(n_CH4)=ch4_init_sh*1.d-6
            else                    ; PIfact(n_CH4)=ch4_init_nh*1.d-6
            end if
            changeL(L,n_CH4)=am(l,i,j)*axyp(i,j)*vol2mass(n_CH4)
     &      *PIfact(n_CH4)  - trm(I,J,L,n_CH4)
          case default
            ! also ensure that strat overwrite is only a sink:
            changeL(L,n_CH4)=-MAX(0d0,trm(I,J,L,n_CH4)-
     &      (AM(L,I,J)*AXYP(I,J)*CH4FACT*1.d-6))
          end select
          else   ! ------ fixed CH4 set in get_CH4_IC(1) -------------
            changeL(L,n_CH4)=0.d0
          end if ! ---------------------------------------------------

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Save overwrite changes for applying in apply_tracer_3Dsource.  C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

          DO N=1,NTM_CHEM
            tr3Dsource(i,j,l,nOverwrite,n) = changeL(l,n)*bydtsrc
          END DO

        end do ! => END LOOP OVER STRATOSPHERE <=
       end do  ! i
      end do   ! j
#endif

CCCCCCCCCCCCCCCCCC END OVERWRITE SECTION CCCCCCCCCCCCCCCCCCCCCC

c Save new tracer O3 and CH4 fields for use in radiation or elsewhere:
c (radiation code wants atm-cm units):
      do j=J_0,J_1
        if(prnchg)DU_O3(J)=0.d0 ! Drew's diagnostic...
        do i=I_0,imaxj(j) 
#ifdef SHINDELL_STRAT_CHEM
         maxl = LM
#else
         select case(which_trop)
         case(0); maxl=ltropo(I,J)
         case(1); maxl=ls1-1
         case default; call stop_model('maxl error',255)
         end select
#endif
         do L=1,maxl
           ! Pass O3 to the rad code:
           chem_tracer_save(1,L,i,j)=pOx(i,j,L)*(trm(i,j,L,n_Ox) +
     &     (tr3Dsource(i,j,L,nChemistry,n_Ox) +
     &     tr3Dsource(i,j,L,nOverwrite,n_Ox))*dtsrc)
     &     *byaxyp(i,j)*byO3MULT
           chem_tracer_save(2,L,i,j)=(trm(i,j,L,n_CH4) +
     &     (tr3Dsource(i,j,L,nChemistry,n_CH4) +
     &     tr3Dsource(i,j,L,nOverwrite,n_CH4))*dtsrc)
     &     *byaxyp(i,j)*avog/(tr_mm(n_CH4)*2.69e20)
           if(prnchg)DU_O3(J)=DU_O3(J)+chem_tracer_save(1,L,i,j)
           ! Above 3D O3 diagnostic in ppbv units is saved (for humans to see).
           ! Here do it in atm-cm units for direct NINT input for rad code.
           taijls(i,j,L,ijlt_O3cmatm)=taijls(i,j,L,ijlt_O3cmatm)+
     &     chem_tracer_save(1,L,i,j)
         end do

         if(maxl < LM) then
           do l=maxl+1,LM
             chem_tracer_save(1,l,i,j)=rad_to_chem(1,l,i,j)
             chem_tracer_save(2,l,i,j)=rad_to_chem(4,l,i,j)
             if(prnchg)DU_O3(J)=DU_O3(J)+chem_tracer_save(1,l,i,j)
           end do
         end if

#if (defined SHINDELL_STRAT_EXTRA) && (defined ACCMIP_LIKE_DIAGS)
         strato3_tracer_save(1:maxl,i,j)=(trm(i,j,1:maxl,n_stratOx) +
     &   (tr3Dsource(i,j,1:maxl,nChemistry,n_stratOx) +
     &   tr3Dsource(i,j,1:maxl,nOverwrite,n_stratOx))*dtsrc)
     &   *byaxyp(i,j)*byO3MULT
         if(maxl < LM)strato3_tracer_save(maxl+1:LM,i,j)=
     &   rad_to_chem(1,maxl+1:LM,i,j)
#endif
         ! accumulate diag for the column sum of O3 mass hopefully similarly to
         ! how taijn Ox_Total_mass is done:
         taijs(i,j,ijs_O3mass)=taijs(i,j,ijs_O3mass)+sum(
     &   pOx(i,j,:)*(trm(i,j,:,n_Ox)+(tr3Dsource(i,j,:,nChemistry,n_Ox)+
     &   tr3Dsource(i,j,:,nOverwrite,n_Ox))*dtsrc))*byaxyp(i,j)

        end do ! i
        if(prnchg)DU_O3(J)=1.d3*DU_O3(J)/IMAXJ(J)
      end do   ! j
      
      if(prnchg)then
       call PACK_DATA( grid, DU_O3, DU_O3_glob )
       IF(AM_I_ROOT()) THEN
         write(6,*) 'Ozone column fm -90 to +90'
         write(6,'(46(f4.0,1x))') (DU_O3_glob(J),J=1,JM)
       END IF
      endif

      RETURN

      contains

      subroutine printSS27x2Etc()
      if(prnchg .and. J == jprn .and. I == iprn) then
        jay = (J >= J_0 .and. J <= J_1)
        write(out_line,*)'O3pO2 means O3prof from O2 Herz & SRB:'
        call write_parallel(trim(out_line),crit=jay)
        write(out_line,*)
     &  'L, O3pO2, O3pO2*C, OHpptv, HO2pptv, O/O3, NO2/NO, Cl/ClO:'
        call write_parallel(trim(out_line),crit=jay)
        do L=LS1,LM
          if((ALB(I,J,1) /= 0.d0).AND.(sza < szamax))then
!old:       ss27x2=2.d0*ss(27,L,i,j)*y(nO2,L)
            ss27x2=2.d0*ss(27,L,i,j)*y(nO2,L) * (rr(98,L)*y(nO2,L))/
     &           (rr(98,L)*y(nO2,L)+rr(88,L)*y(nO3,L))
          else
            ss27x2=0.d0
          end if
          ss27x2_c=ss27x2*DCOS(SZA*radian) ! prob. no longer wanted
          OHpptv=1.d12*y(nOH,L)/y(nM,L)
          HO2pptv=1.d12*y(nHO2,L)/y(nM,L)
          ObyO3=y(nO,L)/y(nO3,L)
          NO2byNO=y(nNO2,L)/y(nNO,L)
          ClbyClO=y(nCl,l)/y(nClO,L)
          write(out_line,'(I3,7(1X,E20.5))')
     &    L,ss27x2,ss27x2_c,OHpptv,HO2pptv,ObyO3,NO2byNO,ClbyClO
          call write_parallel(trim(out_line),crit=jay)
        end do
      end if
      end subroutine printSS27x2Etc


      subroutine printDaytimeChemistryDiags()
CCCCCCCCCCCCC PRINT SOME CHEMISTRY DIAGNOSTICS CCCCCCCCCCCCCCCC
      if(prnchg .and. J == jprn .and. I == iprn) then
       jay = (J >= J_0 .and. J <= J_1) 
       l=lprn
       write(out_line,*) ' '
       call write_parallel(trim(out_line),crit=jay)
       write(out_line,*) 'Family ratios at I,J,L: ',i,j,l
       call write_parallel(trim(out_line),crit=jay)
       write(out_line,*) 'OH/HO2 = ',y(nOH,l)/y(nHO2,l)
       call write_parallel(trim(out_line),crit=jay)
       write(out_line,*) 'O/O3 = ',y(nO,l)/y(nO3,l)
       call write_parallel(trim(out_line),crit=jay)
       write(out_line,*) 'O1D/O3 = ',y(nO1D,l)/y(nO3,l),
     &  '  J(O1D) = ',ss(2,l,I,J)
       call write_parallel(trim(out_line),crit=jay)
       write(out_line,*) 'NO/NO2 = ',y(nNO,l)/y(nNO2,l),
     &  '   J(NO2) = ',ss(1,l,I,J)
       call write_parallel(trim(out_line),crit=jay)
       write(out_line,*) 'conc OH = ',y(nOH,l)
       call write_parallel(trim(out_line),crit=jay)
#ifdef SHINDELL_STRAT_CHEM
       write(out_line,*) 'Cl,ClO,Cl2O2,OClO,Cl2 = ',y(nCl,l),
     &  y(nClO,l),y(nCl2O2,l),y(nOClO,l),y(nCl2,l)
       call write_parallel(trim(out_line),crit=jay)
       write(out_line,*) 'Br,BrO = ',y(nBr,l),y(nBrO,l)
       call write_parallel(trim(out_line),crit=jay)
       write(out_line,*) 'pCl,pClO,pOClO,pBrO = ',pClx(I,J,l),
     &  pClOx(I,J,l),pOClOx(I,J,l),pBrOx(I,J,l)
       call write_parallel(trim(out_line),crit=jay)
#endif
       write(out_line,*)
     & 'sun, SALBFJ,sza,I,J,Itime= ',ALB(I,J,1),sza,I,J,Itime
       call write_parallel(trim(out_line),crit=jay)
!c     do inss=1,JPPJ
!r      write(out_line,195) ' J',inss,ay(ks(inss)),' = ',
!a   &  (ss(inss,Lqq,I,J),Lqq=1,LS1-1)
!s      call write_parallel(trim(out_line),crit=jay)
!h     enddo
!e     write(out_line,196) ' RCloud',(RCLOUDFJ(Lqq,I,J),Lqq=1,LS1-1)
!s     call write_parallel(trim(out_line),crit=jay)
!?     write(out_line,196) ' Ozone ',(y(nO3,Lqq),Lqq=1,LS1-1)
!?     call write_parallel(trim(out_line),crit=jay)
!?     write(out_line,*) ' '
!?     call write_parallel(trim(out_line),crit=jay)
      endif
 195  format(a2,i2,1x,a8,a3,11(1x,e9.2))
 196  format(a7,9x,11(1x,e9.2))
CCCCCCCCCCCCCCCCCCCC END CHEM DIAG SECT CCCCCCCCCCCCCCCCCCCCCCC
      end subroutine printDaytimeChemistryDiags

      subroutine printNightChemistryDiags()
CCCCCCCCCCCCC PRINT SOME CHEMISTRY DIAGNOSTICS CCCCCCCCCCCCCCCC
        if(prnchg.and.J == jprn.and.I == iprn.and.L == lprn)then
          jay = (J >= J_0 .and. J <= J_1)
          write(out_line,*)
     &    'dark, SALBFJ,sza,I,J,L,Itime= ',ALB(I,J,1),sza,I,J,L,Itime
          call write_parallel(trim(out_line),crit=jay)
#ifdef SHINDELL_STRAT_CHEM
          if(pscX(L))then
            write(out_line,*) 'There are PSCs, T =',ta(L)
            call write_parallel(trim(out_line),crit=jay)
          else
            write(out_line,*) 'There are no PSCs, T =',ta(L)
            call write_parallel(trim(out_line),crit=jay)
          endif
#endif
          write(out_line,198) ay(n_NOx),': ',
     &    changeNOx,' molecules produced; ',
     &    100.d0*(changeNOx)/y(n_NOx,L),' percent of'
     &    ,y(n_NOx,L),'(',1.d9*y(n_NOx,L)/y(nM,L),' ppbv)'
          call write_parallel(trim(out_line),crit=jay)
          write(out_line,198) ay(n_HNO3),': ',
     &    changeHNO3,' molecules produced; ',
     &    100.d0*(changeHNO3)/y(n_HNO3,L),' percent of'
     &    ,y(n_HNO3,L),'(',1.d9*y(n_HNO3,L)/y(nM,L),' ppbv)'
          call write_parallel(trim(out_line),crit=jay)
#if (defined TRACERS_HETCHEM) || (defined TOMAS_HETCHEM)
          write(out_line,198) ay(n_HNO3),': ',
     &    (-krate(i,j,l,1,1)*y(n_HNO3,l)*dt2),' molecules dest dust ',
     &    (100.d0*(-krate(i,j,l,1,1)*y(n_HNO3,l)*dt2))/y(n_HNO3,L),
     &    ' percent of'
     &    ,y(n_HNO3,L),'(',1.d9*y(n_HNO3,L)/y(nM,L),' ppbv)'
          call write_parallel(trim(out_line),crit=jay)
#endif
          write(out_line,198) ay(n_N2O5),': ',
     &    changeN2O5,' net molec produced; ',
     &    100.d0*(changeN2O5)/y(n_N2O5,L),' percent of'
     &    ,y(n_N2O5,L),'(',1.d9*y(n_N2O5,L)/y(nM,L),' ppbv)'
          call write_parallel(trim(out_line),crit=jay)
          write(out_line,198) ay(n_N2O5),': ',
     &    gwprodN2O5,' molec prod fm gas;  ',
     &    100.d0*(gwprodN2O5)/y(n_N2O5,L),' percent of'
     &    ,y(n_N2O5,L),'(',1.d9*y(n_N2O5,L)/y(nM,L),' ppbv)'
          call write_parallel(trim(out_line),crit=jay)
          write(out_line,198) ay(n_N2O5),': ',
     &    -wprod_sulf,' molec prod fm sulf; ',
     &    -100.d0*(wprod_sulf)/y(n_N2O5,L),' percent of'
     &    ,y(n_N2O5,L),'(',1.d9*y(n_N2O5,L)/y(nM,L),' ppbv)'
          call write_parallel(trim(out_line),crit=jay)
          write(out_line,198) ay(n_HCHO),': ',
     &    wprodHCHO,' molecules produced; ',
     &    100.d0*(wprodHCHO)/y(n_HCHO,L),' percent of'
     &    ,y(n_HCHO,L),'(',1.d9*y(n_HCHO,L)/y(nM,L),' ppbv)'
          call write_parallel(trim(out_line),crit=jay)
          write(out_line,198) 'Aldehyde',': ',
     &    changeAldehyde,' molecules produced; ',
     &    100.d0*(changeAldehyde)/yAldehyde(I,J,L),' percent of'
     &    ,yAldehyde(I,J,L),'(',1.d9*yAldehyde(I,J,L)/y(nM,L),' ppbv)'
          call write_parallel(trim(out_line),crit=jay)
          write(out_line,198) 'Alkenes ',': ',
     &    changeAlkenes,' molecules produced; ',
     &    100.d0*(changeAlkenes)/y(n_Alkenes,L),' percent of'
     &    ,y(n_Alkenes,L),'(',1.d9*y(n_Alkenes,L)/y(nM,L),' ppbv)'
          call write_parallel(trim(out_line),crit=jay)
#ifdef TRACERS_AEROSOLS_SOA
          write(out_line,198) 'isopp1g ',': ',
     &    changeisopp1g,' molecules produced; ',
     &    100.d0*(changeisopp1g)/y(n_isopp1g,L),' percent of'
     &    ,y(n_isopp1g,L),'(',1.d9*y(n_isopp1g,L)/y(nM,L),' ppbv)'
          call write_parallel(trim(out_line),crit=jay)
          write(out_line,198) 'isopp2g ',': ',
     &    changeisopp2g,' molecules produced; ',
     &    100.d0*(changeisopp2g)/y(n_isopp2g,L),' percent of'
     &    ,y(n_isopp2g,L),'(',1.d9*y(n_isopp2g,L)/y(nM,L),' ppbv)'
          call write_parallel(trim(out_line),crit=jay)
#ifdef TRACERS_TERP
          write(out_line,198) 'apinp1g ',': ',
     &    changeapinp1g,' molecules produced; ',
     &    100.d0*(changeapinp1g)/y(n_apinp1g,L),' percent of'
     &    ,y(n_apinp1g,L),'(',1.d9*y(n_apinp1g,L)/y(nM,L),' ppbv)'
          call write_parallel(trim(out_line),crit=jay)
          write(out_line,198) 'apinp2g ',': ',
     &    changeapinp2g,' molecules produced; ',
     &    100.d0*(changeapinp2g)/y(n_apinp2g,L),' percent of'
     &    ,y(n_apinp2g,L),'(',1.d9*y(n_apinp2g,L)/y(nM,L),' ppbv)'
          call write_parallel(trim(out_line),crit=jay)
#endif  /* TRACERS_TERP */
#endif  /* TRACERS_AEROSOLS_SOA */
          write(out_line,198) 'Isoprene',': ',
     &    changeIsoprene,' molecules produced; ',
     &    100.d0*(changeIsoprene)/y(n_Isoprene,L),' percent of'
     &    ,y(n_Isoprene,L),'(',1.d9*y(n_Isoprene,L)/y(nM,L),' ppbv)'
          call write_parallel(trim(out_line),crit=jay)
#ifdef TRACERS_TERP
          write(out_line,198) 'Terpenes',': ',
     &    changeTerpenes,' molecules produced; ',
     &    100.d0*(changeTerpenes)/y(n_Terpenes,L),' percent of'
     &    ,y(n_Terpenes,L),'(',1.d9*y(n_Terpenes,L)/y(nM,L),' ppbv)'
          call write_parallel(trim(out_line),crit=jay)
#endif  /* TRACERS_TERP */
          write(out_line,198) 'AlkylNit',': ',
     &    changeAlkylNit,' molecules produced; ',
     &    100.d0*(changeAlkylNit)/y(n_AlkylNit,L),' percent of'
     &    ,y(n_AlkylNit,L),'(',1.d9*y(n_AlkylNit,L)/y(nM,L),' ppbv)'
          call write_parallel(trim(out_line),crit=jay)
#ifdef SHINDELL_STRAT_CHEM
          write(out_line,198) ay(n_ClONO2),': ',
     &    changeClONO2,' molecules produced; ',
     &    100.d0*(changeClONO2)/y(n_ClONO2,L),' percent of'
     &    ,y(n_ClONO2,L),'(',1.d9*y(n_ClONO2,L)/y(nM,L),' ppbv)'
          call write_parallel(trim(out_line),crit=jay)
          write(out_line,198) ay(n_ClOx),': ',
     &    changeClOx,' molecules produced; ',
     &    100.d0*(changeClOx)/y(n_ClOx,L),' percent of'
     &    ,y(n_ClOx,L),'(',1.d9*y(n_ClOx,L)/y(nM,L),' ppbv)'
          call write_parallel(trim(out_line),crit=jay)
          write(out_line,198) ay(n_HOCl),': ',
     &    changeHOCl,' molecules produced; ',
     &    100.d0*(changeHOCl)/y(n_HOCl,L),' percent of'
     &    ,y(n_HOCl,L),'(',1.d9*y(n_HOCl,L)/y(nM,L),' ppbv)'
          call write_parallel(trim(out_line),crit=jay)
          write(out_line,198) ay(n_HCl),': ',
     &    changeHCl,' molecules produced; ',
     &    100.d0*(changeHCl)/y(n_HCl,L),' percent of'
     &    ,y(n_HCl,L),'(',1.d9*y(n_HCl,L)/y(nM,L),' ppbv)'
          call write_parallel(trim(out_line),crit=jay)

          write(out_line,198) ay(n_BrONO2),': ',
     &    changeBrONO2,' molecules produced; ',
     &    100.d0*(changeBrONO2)/y(n_BrONO2,L),' percent of'
     &    ,y(n_BrONO2,L),'(',1.d9*y(n_BrONO2,L)/y(nM,L),' ppbv)'
          call write_parallel(trim(out_line),crit=jay)
          write(out_line,198) ay(n_BrOx),': ',
     &    changeBrOx,' molecules produced; ',
     &    100.d0*(changeBrOx)/y(n_BrOx,L),' percent of'
     &    ,y(n_BrOx,L),'(',1.d9*y(n_BrOx,L)/y(nM,L),' ppbv)'
          call write_parallel(trim(out_line),crit=jay)
          write(out_line,198) ay(n_HBr),': ',
     &    changeHBr,' molecules produced; ',
     &    100.d0*(changeHBr)/y(n_HBr,L),' percent of'
     &    ,y(n_HBr,L),'(',1.d9*y(n_HBr,L)/y(nM,L),' ppbv)'
          call write_parallel(trim(out_line),crit=jay)

C         write(out_line,199) 'NO2, NO3  = ',y(nNO2,L),yNO3(I,J,L)
C         call write_parallel(trim(out_line),crit=jay)
C         write(out_line,198) ay(n_Ox),': ',
C    &    changeOx,' molecules produced; ',
C    &    100.d0*(changeOx)/y(n_Ox,L),' percent of'
C    &    ,y(n_Ox,L),'(',1.d9*y(n_Ox,L)/y(nM,L),' ppbv)'
C         call write_parallel(trim(out_line),crit=jay)
#endif
          write(out_line,199) 'NO2, NO3  = ',y(nNO2,L),yNO3(I,J,L)
          call write_parallel(trim(out_line),crit=jay)
c     diagnose N conservation
          changeN=changeNOx+2.d0*changeN2O5+changeHNO3+
     &    changeClONO2+changeBrONO2+changeAlkylNit
          write(out_line,197) '===> N conservation: changeN= ',changeN
          call write_parallel(trim(out_line),crit=jay)
        endif
 197    format(a30,e13.3)
 198    format(1x,a8,a2,e13.3,a21,f10.0,a11,2x,e13.3,3x,a1,f12.5,a6)
 199    format(1x,a20,2(2x,e13.3))
CCCCCCCCCCCCCCCCCCCC END CHEM DIAG SECT CCCCCCCCCCCCCCCCCCCCCCC
        end subroutine printNightChemistryDiags


        subroutine checkNighttimeTolerances()
C Make sure nighttime chemistry changes are not too big:
        error=.false.
        if(changeNOx < -1.d15.OR.changeNOx > 1.d15) then
          write(6,*) 'Big chg@ Itime,I,J,L,NOx ',Itime,I,J,L,changeNOx
          write(6,*) 'rlossN,rprodN,ratioN =',rlossN,rprodN,ratioN
          error=.true.
        end if
        if(changeHNO3 < -1.d15.OR.changeHNO3 > 1.d15) then
          write(6,*) 'Big chg@ Itime,I,J,L,HNO3',Itime,I,J,L,changeHNO3
          error=.true.
        end if
        if(changeN2O5 < -1.d15.OR.changeN2O5 > 1.d15) then
          write(6,*) 'Big chg@ Itime,I,J,L,N2O5',Itime,I,J,L,changeN2O5
          error=.true.
        end if
        if(wprodHCHO < -1.d15.OR.wprodHCHO > 1.d15) then
          write(6,*)'Big chg@ Itime,I,J,L,HCHO',Itime,I,J,L,wprodHCHO
          error=.true.
        endif
        if(error)call stop_model('nighttime chem: big changes',255)
        end subroutine checkNighttimeTolerances

      END SUBROUTINE masterchem



      subroutine photo_acetone(I,J,sza,Jacet)
!@sum calculate photolysis rate for acetone geometrically
!@+ taken from the UK Harwell Model
!@auth Drew Shindell (modelEifications by Greg Faluvegi)
!@ver  1.0 (based on dsf_master_chem_GCM.S)
!
      use model_com, only: jday,jhour,LS1,LM
      use dynamics, only: pmid
      use geom, only:  lat2d ! lat is in radians
      use constant, only: twopi,pi,radian,teeny

!@var dec declination angle of the earth
!@var lha local hour angle
!@var CC COS(lat)*COS(dec)
!@var SS SIN(lat)*SIN(dec)

      real*8, parameter :: C1=9.269d-7,C2=1.563d0,C3=0.301d0
      real*8 :: dec,lha,CC,SS,sec_func,Jacet0
      real*8, intent(in):: sza ! passed in radians
      real*8, dimension(LM), intent(out) :: Jacet
      integer, intent(in) :: i,j
      integer :: L

      dec=radian*23.455d0*COS( ((jday-173)*twopi)/365.d0 )
      lha=twopi*real(jhour)/24.d0

      CC=COS(lat2d(I,J))*COS(dec)
      SS=SIN(lat2d(I,J))*SIN(dec)

      sec_func=1.d0/max(teeny,COS(lha)*CC+SS)
 
      Jacet0=max(0.d0,C1*(COS(sza)*C2)*EXP(-1.*C3*sec_func))
      Jacet(:)=0.d0
      do L=1,LS1-1
        Jacet(L)=3.d0*Jacet0/LOG(pmid(L,i,j))
      enddo
      
      return
      end subroutine photo_acetone



      SUBROUTINE Crates(I,J,rh,bythick
#ifdef SHINDELL_STRAT_CHEM
     &   ,aero
#endif
     & )
!@sum Crates calculate chemical reaction rates for each altitude,
!@+   using JPL 00.  Includes special calculations for pressure
!@+   dependent reactions. Specifically:
!@+   #13 CO+OH->HO2+CO2, #15 HO2+HO2->H2O2+O2, #16 OH+HNO3->H2O+NO3,
!@+   and reactions #29, and #42.
!@auth Drew Shindell (modelEifications by Greg Faluvegi)
!@ver  1.0 (based on masterchem000_M23p)

C**** GLOBAL parameters and variables:
      USE MODEL_COM, only: LM,JM,LS1,ptop,psf,sig,Itime,ItimeI,
     & coupled_chem
      USE RAD_COM, only  : rad_to_chem
      USE CONSTANT, only : PI
      USE DYNAMICS, only : LTROPO, am
      USE TRACER_COM, only : trm, n_N2O5, vol2mass
#ifdef TRACERS_AEROSOLS_Koch
     &                       ,n_SO4
#endif
#ifdef TRACERS_AMP
     &                       ,n_M_AKK_SU,n_M_ACC_SU,n_M_DD1_SU
     &                       ,n_M_DS1_SU,n_M_DD2_SU,n_M_DS2_SU
     &                       ,n_M_SSA_SU,n_M_OCC_SU,n_M_BC1_SU
     &                       ,n_M_BC2_SU,n_M_BC3_SU,n_M_DBC_SU
     &                       ,n_M_BOC_SU,n_M_BCS_SU,n_M_MXX_SU
#endif
#ifdef TRACERS_TOMAS
     &                       ,n_ASO4,nbins
#endif
      USE TRDIAG_COM, only : jls_N2O5sulf
      USE TRCHEM_Shindell_COM, only: nr2,nr3,nmm,nhet,ta,ea,rr,pe,
     & cboltz,r1,sb,nst,y,nM,nH2O,ro,sn,which_trop,sulfate,RKBYPIM,dt2,
     & RGAMMASULF
#ifdef SHINDELL_STRAT_CHEM
     &        ,pscX
#endif

#ifdef TRACERS_AEROSOLS_SOA
      USE TRACER_COM, only: n_isopp1a,n_isopp2a
#ifdef TRACERS_TERP
     &                     ,n_apinp1a,n_apinp2a
#endif  /* TRACERS_TERP */
      USE TRACERS_SOA, only: KpCALC,kpart,kpart_ref,kpart_temp_ref,
     &                       whichsoa,dH_isoprene,dH_apinene
#endif  /* TRACERS_AEROSOLS_SOA */
      USE GEOM, only : lat2d_dg,byaxyp,axyp
      IMPLICIT NONE

C**** Local parameters and variables and arguments:
!@var I,J passed horizontal position indicies
!@var dd,pp,fw,rkp,rk2,rk3M,nb,rrrr,temp dummy "working" variables
!@var L,jj dummy loop variables
!@var byta reciprocal of the local temperature
!@var rkext aerosol extinction from SAGE obs
!@var pscEx NAT PSC surface conc per unit volume (cm^2/cm^3)
!@var Ltop is number of levels with chemistry
!@var beta branching ratio for (HO2+NO) reactions
!@var pcon variable for some pressure conversions
      REAL*8:: byta,dd,pp,fw,rkp,rk2,rk3M,rrrr,temp,beta,pcon,waterPPMV
      real*8 :: associationReaction, activationReaction,pfactor,
     & bypfactor,k0T,k0TM,kinfT,kinfTbyM,RVELN2O5,wprod_sulf,prod_sulf
      INTEGER             :: L,jj,nb,Ltop
      INTEGER, INTENT(IN) :: I,J
      REAL*8, DIMENSION(LM):: rh,bythick
#ifdef SHINDELL_STRAT_CHEM
!@var PRES local nominal pressure
!@var LAXt,LAXb lowest and highest levels to have nonzero 
!@+   RAD-code aerosol extinction 
!@var aero array =1 for nonzero rkext, otherwise 0.
      REAL*8, DIMENSION(LM) :: PSCEX,PRES,rkext
      INTEGER               :: LAXt,LAXb
      INTEGER, INTENT(OUT), dimension(LM) :: aero
#endif

#ifdef SHINDELL_STRAT_CHEM
      aero(:)=0
      Ltop=LM
      PRES(:)=SIG(:)*(PSF-PTOP)+PTOP
      rkext(:)=0.d0 ! initialize over L
      if(rad_to_chem(2,1,I,J) /= 0.)call stop_model('kext prob 0',255)
      do L=2,Ltop
        if(rad_to_chem(2,L,I,J) /= 0..and.rad_to_chem(2,L-1,I,J) == 0.)
     &  LAXb=L
        if(rad_to_chem(2,L,I,J) == 0..and.rad_to_chem(2,L-1,I,J) /= 0.)
     &  LAXt=L-1
      end do
#else
      select case(which_trop)
      case(0); Ltop=ltropo(I,J)
      case(1); Ltop=ls1-1
      case default; call stop_model('which_trop problem 6',255)
      end select
#endif
      do L=1,Ltop            !  --> BEGIN ALTITUDE LOOP <--
        byta=1.d0/ta(L)
        pcon=y(nM,L)*ta(L)*cboltz/1013.d0
        do jj=1,nr2             ! bimolecular rates start
          IF(ea(jj) /= 0.d0) THEN
            rr(jj,L)=pe(jj)*exp(-ea(jj)*byta)
          ELSE
            rr(jj,L)=pe(jj)
          END IF
c         for #12, k based on three-parameters from JPL2011
          if(jj == 12) rr(jj,L)=rr(jj,L)*(ta(L)**0.667)
c         for #13, k= based on termolecular reaction from JPL2011
c         (see paged 185-188 and note D1)
          if(jj == 13) then
            k0TM=y(nM,L)*pe(jj)*((300.d0*byta)**1.4)
            kinfT=1.1d-12*(300.d0*byta)**(-1.3)
            dd=k0TM/kinfT
            pp=0.6d0**(1.d0/(1.d0+(log10(dd))**2.))
            associationReaction=(k0TM/(1.d0+dd))*pp
            k0T=1.5d-13*((300.d0*byta)**(-0.6))
            kinfTbyM=(2.1d9*((300.d0*byta)**(-6.1)))/y(nM,L)
            dd=k0T/kinfTbyM
            pp=0.6d0**(1.d0/(1.d0+(log10(dd))**2.))
            activationReaction=(k0T/(1.d0+dd))*pp
            rr(jj,L)=associationReaction+activationReaction
          end if
c         for reaction #15, k=(kc+kp)fw, kc=rr
          if(jj == 15)then
            rkp=2.1d-33*y(nM,L)*exp(920.d0*byta) 
            fw=(1.d0+1.4d-21*y(nH2O,L)*exp(2200.d0*byta))
            rr(jj,L)=(rr(jj,L)+rkp)*fw
          endif
c         for #16, k=[pe*exp(-e(jj)/ta(l))]+k3[M]/(1+k3[M]/k2)
          if(jj == 16)then
            rk3M=y(nM,l)*6.5d-34*exp(1335.d0*byta)
            rk2=2.7d-17*exp(2199.d0*byta) 
            rr(jj,L)=rr(jj,L)+rk3M/(1.d0+(rk3M/rk2))
          endif
          if(jj == 29)rr(jj,L)=rr(jj,L)/y(nM,L)!PAN+M really PAN
          if(jj == 42)rr(jj,L)=rr(jj,L)/y(nM,L)!ROR+M really ROR
!         for #6 & 91 (HO2+NO) calculate branching ratio here          
!         Butkovskaya et al J.Phys.Chem 2007         
          waterPPMV=1.d6*y(nH2O,L)/y(nM,L)
          if(ta(L)<298.d0 .and. waterPPMV > 100.)then
            beta=(530.d0*byta + 6.4d-4*pcon*760.d0 - 1.73d0)*1.d-2
          else
            beta=0.d0
          endif
#ifdef SHINDELL_STRAT_CHEM
          if(jj == 91)rr(jj,L)=rr(jj,L)*beta
          if(jj ==  6)rr(jj,L)=rr(jj,L)*(1.d0-beta)
#else
          if(jj == 45)rr(jj,L)=rr(jj,L)*beta
          if(jj ==  6)rr(jj,L)=rr(jj,L)*(1.d0-beta)
#endif
        end do                ! bimolecular rates end

        ! here we used to scale rr(86 & 87) N2O+O(1D)-->NO+NO and N2O+O(1D)-->N2+O2

        do jj=1,nr3           ! trimolecular rates start
          rr(nr2+jj,L)=y(nM,L)*ro(jj)*(300.d0*byta)**sn(jj)
          if(r1(jj) .ne. 0.d0)then 
            dd=rr(nr2+jj,L)/(r1(jj)*(300.d0*byta)**sb(jj))
            pp=0.6d0**(1.d0/(1.d0+(log10(dd))**2.))
            rr(nr2+jj,L)=(rr(nr2+jj,L)/(1.d0+dd))*pp
          endif
        end do                ! trimolecular rates end

        nb=nr2-nmm
        if(nmm >= 1) then
          do jj=1,nmm         ! monomolecular rates start
           ! 0.5 for precision,correct following line:
           rrrr=exp(0.5d0*ea(jj+nb)*byta)
           rr(jj+nb,L)=rr(nst(jj),L)/(rrrr*pe(jj+nb)*rrrr*y(nM,l))     
          end do              ! monomolecular rates end
        end if

#ifdef SHINDELL_STRAT_CHEM
c Calculate rates for heterogeneous reactions (Divided by solid
C in Chem1). sticking coefficients from JPL '02:
c       1=N2O5 + H2O --> 2HNO3          gamma=0.2, 0.0004 (PSC)
c       2=ClONO2 + H2O --> HOCl + HNO3  gamma=0.8d-2 (aero), 0.004 (PSC)
c       3=ClONO2 + HCl --> Cl2 + HNO3   gamma=0.2
c       4=HOCl + HCl --> Cl2 + H2O      gamma=0.1
c       5=N2O5 + HCl --> ClNO2 + HNO3   gamma=0.003
C
C Aerosols (14-33 km) & PSCs 14-22 km.
C
c Aerosol profiles and latitudinal distribution of extinction 
c coefficients(in km**-1) are from SAGE II data on GISS web site:

        if(pres(l) >= 245.d0 .or. pres(l) <= 5.d0)then 
          do jj=nr2+nr3+2,nr2+nr3+nhet
            rr(jj,L)=1.0d-35
          enddo 
          ! Add rxn of N2O5 on sulfate analogous to what is done in darkness:
          if(rh(L)>0.5)then
            rgammasulf = 1.5d-2
          else
            rgammasulf = 5.2d-2 - 2.79d-4*100.d0*rh(L)
            if(ta(L)>290.) rgammasulf=
     &      max(1.d-3,rgammasulf-log10(ta(L)-290.d0)*5.d-2)
          endif
          if (coupled_chem == 1) then
            ! Convert SO4 from mass (kg) to aerosol surface per grid box:
            ! Here there is a factor of 1d-3  that converts kg/m3 to g/cm3
            ! and 1.76d5 is cm2/g from Dentener and Crutzen, 1993.
            ! So 1.d-3*1.76d5=1.76d2, and that value is for a relative
            ! humidity of 0.75 (1/0.75 = 1.33333d0 below). Reciprocal
            ! layer thickness below is in 1/m units:
#ifndef TRACERS_TOMAS
          sulfate(i,j,L)=
#ifdef TRACERS_AMP
     &      (trm(i,j,L,n_M_AKK_SU)+trm(i,j,L,n_M_ACC_SU)+
     &       trm(i,j,L,n_M_DD1_SU)+trm(i,j,L,n_M_DS1_SU)+
     &       trm(i,j,L,n_M_DD2_SU)+trm(i,j,L,n_M_DS2_SU)+
     &       trm(i,j,L,n_M_SSA_SU)+trm(i,j,L,n_M_OCC_SU)+
     &       trm(i,j,L,n_M_BC1_SU)+trm(i,j,L,n_M_BC2_SU)+
     &       trm(i,j,L,n_M_BC3_SU)+trm(i,j,L,n_M_DBC_SU)+
     &       trm(i,j,L,n_M_BOC_SU)+trm(i,j,L,n_M_BCS_SU)+
     &       trm(i,j,L,n_M_MXX_SU))
#else
     &      trm(i,j,L,n_SO4)
#endif
     &      *1.76d2*byaxyp(i,j)*bythick(L)
     &      *max(0.1d0,rh(L)*1.33333d0)
#endif
#ifdef TRACERS_TOMAS
           sulfate(i,j,l)=0.0
           do nb=1,nbins
              sulfate(i,j,L)=sulfate(i,j,l)+trm(i,j,L,n_ASO4(nb))*1.76d2
     &            *byaxyp(i,j)*bythick(L)*max(0.1d0,rh(L)*1.33333d0)
           enddo
#endif
          endif
          pfactor=axyp(I,J)*AM(L,I,J)/y(nM,L)
          bypfactor=1.D0/pfactor
          RVELN2O5=SQRT(ta(L)*RKBYPIM)*100.d0
C         Calculate sulfate sink, and cap it at 20% of N2O5:
c         in troposphere loss is rxn on sulfate, in strat rxn w PSC or sulfate
          wprod_sulf=
     &    dt2*sulfate(I,J,L)*y(n_N2O5,L)*RGAMMASULF*RVELN2O5*0.25d0

          if(wprod_sulf > 0.2d0*y(n_N2O5,L))wprod_sulf=0.2d0*y(n_N2O5,L)
          prod_sulf=wprod_sulf*pfactor
          CALL INC_TAJLS(I,J,L,jls_N2O5sulf,-prod_sulf*vol2mass(n_N2O5))
          rr(nr2+nr3+1,L)=wprod_sulf/(dt2*y(n_N2O5,L))
        else  
          if((pres(l) < 245.d0.and.pres(l) > 150.d0) .or. 
     &    LAXb < 1.or.LAXb > ltop.or.LAXt < 1.or.LAXt > ltop)then 
            rkext(l)=0.d0
          else
            if(pres(l) <= 150..and.pres(l) > 31.60)then
              if(l < LAXb) then
                rkext(l)=5.d-2*rad_to_chem(2,LAXb,i,j)
              else if(l > LAXt) then
                rkext(l)=0.33d0*rkext(l-1)
              else
                rkext(l)=5.d-2*rad_to_chem(2,l,i,j)
              endif
            endif
            if(pres(l) <= 31.6d0.and.pres(l) >= 17.8d0)then
              if(l < LAXb) then
                call stop_model('kext problem 1',255)
              else if(l > LAXt) then
                rkext(l)=2.0d0*rkext(l-1)
              else
                rkext(l)=5.d-2*rad_to_chem(2,l,i,j)
              endif
            endif
            if(pres(l) <= 17.8d0.and.pres(l) >= 10.0d0)then
              if(l < LAXb) then
                call stop_model('kext problem 2',255)
              else if(l > LAXt) then
                rkext(l)=16.d0*8.33333d-2*rkext(l-1)
              else
                rkext(l)=5.d-2*rad_to_chem(2,l,i,j)
              endif
            endif
            if(pres(l) <= 10.0d0.and.pres(l) >= 4.6d0)then
              if(l < LAXb) then
                call stop_model('kext problem 3',255)
              else if(l > LAXt) then
                rkext(l)=0.4d0*6.6667d-1*rkext(l-1)
              else
                rkext(l)=0.5d-2*rad_to_chem(2,l,i,j)
              endif
            endif
          endif

          ! here for certain pressures in the tropics we used to scale rkext
           
          if(rkext(l) /= 0.)aero(l) = 1

          ! NAT PSC surface conc per unit volume (cm^2/cm^3)
          if(pscX(l))then
            pscEx(l)=2.d-6
          else
            pscEx(l)=0.d0
          endif

c         Reaction 1 on sulfate and PSCs:      
          temp=sqrt(8.d0*1.38d-16*ta(l)*6.02d23/(PI*108.d0))
          rr(nr2+nr3+1,L)=0.5d0*rkext(l)*1.d-5*temp*0.2d0
          if(pres(l) > 31.6d0) rr(nr2+nr3+1,L)=
     &    rr(nr2+nr3+1,L)+0.25d0*pscEx(l)*temp*0.0004d0

c         Reaction 2 on sulfate and PSCs:      
          temp=sqrt(8.d0*1.38d-16*ta(l)*6.02d23/(PI*97.d0))
          rr(nr2+nr3+2,L)=0.5d0*rkext(l)*1.d-5*temp*0.8d-2
          if(pres(l) > 31.6d0) rr(nr2+nr3+2,L)=
     &    rr(nr2+nr3+2,L)+0.25d0*pscEx(l)*temp*4.d-3 

          if(pres(l) > 31.6d0) then
            rr(nr2+nr3+3,L)=0.25d0*pscEx(l)*temp*0.2d0
            rr(nr2+nr3+4,L)=
     &      sqrt(8.d0*1.38d-16*ta(l)*6.02d23/(PI*52.d0))
            rr(nr2+nr3+4,L)=0.25d0*pscEx(l)*rr(nr2+nr3+4,L)*0.1d0
            rr(nr2+nr3+5,L)=
     &      sqrt(8.d0*1.38d-16*ta(l)*6.02d23/(PI*108.d0))
            rr(nr2+nr3+5,L)=0.25d0*pscEx(l)*rr(nr2+nr3+5,L)*0.003d0
          endif
        endif  

        if(pres(l) < 245.d0 .and. pres(l) > 5.d0)then
          wprod_sulf=dt2*y(n_N2O5,L)*rr(nr2+nr3+1,L)
          prod_sulf=wprod_sulf*pfactor
          CALL INC_TAJLS(I,J,L,jls_N2O5sulf,-prod_sulf*vol2mass(n_N2O5))
        end if
#endif
      end do                  !  --> END ALTITUDE LOOP <--

#ifdef TRACERS_AEROSOLS_SOA
      do L=1,LM ! this should be up to LM, no matter if strat chem is on or off
        kpart(L,whichsoa(n_isopp1a))=
     &       KpCALC(dH_isoprene,kpart_ref(whichsoa(n_isopp1a)),ta(L),
     &              kpart_temp_ref(whichsoa(n_isopp1a)))
        kpart(L,whichsoa(n_isopp2a))=
     &       KpCALC(dH_isoprene,kpart_ref(whichsoa(n_isopp2a)),ta(L),
     &              kpart_temp_ref(whichsoa(n_isopp2a)))
#ifdef TRACERS_TERP
        kpart(L,whichsoa(n_apinp1a))=
     &       KpCALC(dH_apinene,kpart_ref(whichsoa(n_apinp1a)),ta(L),
     &              kpart_temp_ref(whichsoa(n_apinp1a)))
        kpart(L,whichsoa(n_apinp2a))=
     &       KpCALC(dH_apinene,kpart_ref(whichsoa(n_apinp2a)),ta(L),
     &              kpart_temp_ref(whichsoa(n_apinp2a)))
#endif  /* TRACERS_TERP */
      enddo
#endif  /* TRACERS_AEROSOLS_SOA */
 
      RETURN
      END SUBROUTINE Crates



#ifdef TRACERS_SPECIAL_Shindell
      SUBROUTINE checktracer(I,J)
!@sum checktracer for various debugging of tracer chemistry
!@auth Drew Shindell (modelEifications by Greg Faluvegi)
!@ver  1.0 (based on masterchem000_M23p)

C**** GLOBAL parameters and variables:
      USE MODEL_COM, only  : Itime, LM, LS1
      USE DYNAMICS, only   : LTROPO
      USE TRACER_COM, only : ntm, trname, n_Ox, ntm_chem
      USE TRCHEM_Shindell_COM, only: y, nM, which_trop

      IMPLICIT NONE

C**** Local parameters and variables and arguments:
!@var I,J passed horizontal position indicies
!@var L,igas dummy loop variables
!@var tlimit if tracer goes above this limit, model stops
!@var checkOx logical: should I check for large tropospheric Ox?
!@var checkmax logical: should I check for large tracers throughout?
!@var checkNeg logical: should I check for negative tracers?
!@var checkNaN logical: should I check for unreal tracers?
!@var maxL LTROPO(I,J) or LS1-1, depending upon which_trop variable

      INTEGER                  :: L, igas, maxL
      INTEGER, INTENT(IN)      :: I,J
      REAL*8, DIMENSION(ntm_chem) :: tlimit
    
      LOGICAL :: checkOx, checkmax, checkNeg, checkNan
      DATA checkNeg /.true./
      DATA checkNan /.true./
      DATA checkOx  /.true./
      DATA checkmax /.true./
      integer,dimension(4) ::  error
      integer :: is_error
      character*80, dimension(4) :: message

      tlimit(:)=1.d-5
      message(1)='Ox too big.'
      message(2)='A tracer is too big.'
      message(3)='A tracer is negative.'
      message(4)='A tracer is NaN.'
      error(:)=0

      IF(i == 1.and.j == 1)
     & WRITE(6,*) 'WARNING: checktracer call is active.'
      select case(which_trop)
      case(0); maxl=ltropo(I,J)
      case(1); maxl=ls1-1
      case default; call stop_model('which_trop problem 7',255)
      end select

C check if ozone gets really big:
       IF(checkOx) THEN
       do L=1,maxL
         if(y(n_Ox,L)/y(nM,L) > 1.d-5) then
           write(6,*)'Ox @ I,J,L,Ox,Itime:',I,J,L,y(n_Ox,L),Itime
           error(1)=1
         end if
       end do
#ifdef SHINDELL_STRAT_CHEM
       do L=maxL+1,LM
         if(y(n_Ox,L)/y(nM,L) > 1.5d-5) then
           write(6,*)'Ox @ I,J,L,Ox,Itime:',I,J,L,y(n_Ox,L),Itime
           error(1)=1
         end if
       end do
#endif
       END IF
       
c general check on maximum of tracers:
!!! note, this is only as useful as the limits you set are appropriate!
      IF(checkmax) THEN
      do L=1,LM
       do igas=1,ntm_chem
        if(y(igas,L)/y(nM,L) > tlimit(igas)) then
          write(6,*) trname(igas),'@ I,J,L,Y :',
     &    I,J,L,y(igas,L)/y(nM,L)
          error(2)=1
        end if
       end do
      end do
      END IF

c check for negative tracers:
      IF(checkNeg) THEN
      do L=1,LM
       do igas=1,ntm_chem
        if(y(igas,L) < 0.d0) THEN
          write(6,*)trname(igas),
     &    'negative @ tau,I,J,L,y:',Itime,I,J,L,y(igas,L)
          error(3)=1
        end if
       enddo
      end do
      END IF

c check for unreal (not-a-number) tracers (maybe SGI only?):
      IF(checkNaN) THEN
      do L=1,LM
       do igas=1,ntm_chem
        if(.NOT.(y(igas,L) > 0.d0.OR.y(igas,L) <= 0.d0)) THEN
         write(6,*)trname(igas),
     &   'is not a number @ tau,I,J,L,y:',Itime,I,J,L,y(igas,L)
         error(4)=1
        end if
       enddo
      end do
      END IF

      is_error=0
      do L=1,4
        if (error(L)>0)then
          is_error=1
          write(6,*)trim(message(L))
        endif
      enddo
      if(is_error==1)call stop_model('error in checktracer',255)
   
      RETURN

      END SUBROUTINE checktracer
#endif
      
      
