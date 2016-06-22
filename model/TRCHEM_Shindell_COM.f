#include "rundeck_opts.h"

      MODULE TRCHEM_Shindell_COM
!@sum  TRCHEM_Shindell_COM declares variables for tracer chemistry
!@+    and sources.
!@auth Drew Shindell (modelEifications by Greg Faluvegi)
!@ver  1.0 (based on various chemistry modules of B436Tds3YM23 model)
c
#ifdef TRACERS_ON
      USE MODEL_COM, only  : im,jm,lm,psf,ptop,sig,sige,dsig,bydsig,
     &                       dtsrc,Itime,ItimeI,T
      USE CONSTANT, only   : pi, mair, mwat, radian,avog,pO2
      USE DYNAMICS, only   : am, byam, PMID, PK
      USE TRACER_COM, only : ntm, trm, TR_MM, ntm_soa, ntm_terp

      IMPLICIT NONE
      SAVE

C**************  P  A  R  A  M  E  T  E  R  S  *******************
!@param p_1 number of reactants per reaction
!@param p_2 number of rxns in assembled lists (check with print rxn list)
!@param p_3 number of rxns in assembled lists (check with print rxn list)
!@param p_4 number of rxns in assembled lists (check with print rxn list)
!@param p_5 number of levels from top down with SRB flux
!@param n_rx maximum number of chemical reactions
!@param n_bi maximum number of bimolecular reactions
!@param n_tri maximum number of trimolecular reactions
!@param n_nst maximum number of monomolecular decompositions
!@param n_fam maximum number of chemical families
!@param JPPJ_Shindell number of photolysis reactions in the Shindell chemistry
!@param luselb Use reflective photolysis boundary treatment
!@param zlbatm Optical depth above which to set lower boundary
!@param CMEQ1 ?
!@param nc total number of molecules included (incl. O2 and N2)
!@param ny number of chemically calculated gases (no O2 or N2)
!@param numfam number of chemical families
!@param n_phot how often to do photolysis (in increments of DTsrc)
!@param O3MULT =2.14d-2 This is the conversion from (atm*cm) units
!@+     (i.e. 1000 Dobson Units) to KG/m2. It is: 
!@+     1.E4*2.69E19*48./6.02E26 where 1.E4 is cm2/m2, 2.69E19 is 
!@+     molecules/cm3 at 1 atm pressure, 48. is molecular wt of O3,
!@+     and 6.02E26 is Avogadro's number in molecules/Kmol.
!@param cpd conversion from molecules/cm3 to mole/m3
!@param BYO3MULT = 1/O3MULT
!@param pfix_H2 fixed ratio of H2/M
!@param pfix_Aldehyde fixed ratio of Aldehyde/M for initial conditions
!@param checktracer_on integer to turn on the checktracer call
!@param MWabyMWw ratio of molecular weights of air/water
!@param RKBYPIM=8.*RBOLTZ/pi/MASSN2O55=8.*1.38062D-23/3.14159/1.793D-25
!@param cboltz Boltzman's Constant = 1.3806d-19
!@param byradian 1/radian = conversion from radians to degrees
!@param LCOalt number of levels in the several tracer IC arrays
!@param LCH4alt number of levels in the CH4altIN array
!@param PCOalt pressures at LCOalt levels
!@param PCH4alt pressures at LCH4alt levels
!@param T_thresh threshold temperature used in master chem
!@param n2o_pppv default N2O L=1 overwriting in pppv
!@param cfc_pppv default CFC L=1 overwriting in pppv
!@param cfc_rad95 the average L=1 radiation code CFC11+CFC12 value
!@+     for 1995 (pppv). 
!@param fact_cfc ratio of our default CFC L=1 overwriting to the 
!@+     radiation's 1995 L=1 CFC11+CFC12 value.
!@param PSClatS SH latitude limit for PSCs
!@param PSClatN NH latitude limit for PSCs
!@param minKG minimum kg for trm before we set to this post-change
      INTEGER, PARAMETER ::
     & LCOalt =   23,
     & LCH4alt=    6,
     & p_1   =     2 
      INTEGER, PARAMETER ::
     & p_2   =   209,
     & p_3   =   500,
     & p_4   =   209,
#ifdef TRACERS_TERP
     & n_rx  =   113,
     & n_bi  =    97,
#else
     & n_rx  =   110,
     & n_bi  =    94,
#endif  /* TRACERS_TERP */
     & n_tri =    11,
     & n_nst =     3,
     & nc     =   53+ntm_terp+ntm_soa,     !formerly in param sub
     & ny     =   51+ntm_terp+ntm_soa,     !formerly in param sub  
     & numfam =    4,     !formerly in param sub  
     & nC2O3=     26+ntm_terp+ntm_soa,
     & nXO2=      27+ntm_terp+ntm_soa,
     & nXO2N=     28+ntm_terp+ntm_soa,
     & nRXPAR=    29+ntm_terp+ntm_soa,
     & nROR=      30+ntm_terp+ntm_soa,
     & nAldehyde= 31+ntm_terp+ntm_soa,
     & nH2O=      32+ntm_terp+ntm_soa,
     & nCH3O2=    33+ntm_terp+ntm_soa,
     & nH2=       34+ntm_terp+ntm_soa,
     & nOH=       35+ntm_terp+ntm_soa,
     & nHO2=      36+ntm_terp+ntm_soa,
     & nO3=       37+ntm_terp+ntm_soa,
     & nO=        38+ntm_terp+ntm_soa,
     & nO1D=      39+ntm_terp+ntm_soa,
     & nNO=       40+ntm_terp+ntm_soa,
     & nNO2=      41+ntm_terp+ntm_soa,
     & nNO3=      42+ntm_terp+ntm_soa,
     & nHONO=     43+ntm_terp+ntm_soa,
     & nCl2O2=    44+ntm_terp+ntm_soa,
     & nClO=      45+ntm_terp+ntm_soa,
     & nOClO=     46+ntm_terp+ntm_soa,
     & nCl2=      47+ntm_terp+ntm_soa,
     & nCl=       48+ntm_terp+ntm_soa,
     & nBrCl=     49+ntm_terp+ntm_soa,
     & nBrO=      50+ntm_terp+ntm_soa,
     & nBr=       51+ntm_terp+ntm_soa,
     & nO2=       52+ntm_terp+ntm_soa,
     & nM=        53+ntm_terp+ntm_soa,     !you must always put nM last (highest number)
     & JPPJ_Shindell = 28,
     & n_fam =     5
      INTEGER, PARAMETER ::
     & p_5   =    14,
C ----------------------------------------------     
c     & n_Ox=        1,    ! note, these
c     & n_NOx=       2,    ! first 15 species are
c     & n_N2O5=      3,    ! tracers, and therefore
c     & n_HNO3=      4,    ! these parameters are
c     & n_H2O2=      5,    ! to be defined in 
c     & n_CH3OOH=    6,    ! TRACER_COM.f.
c     & n_HCHO=      7,    ! Note the UNDERSCORE!
c     & n_HO2NO2=    8,    !  T
c     & n_CO=        9,    !  R
c     & n_CH4=      10,    !  A
c     & n_PAN=      11,    !  C
c     & n_Isoprene= 12,    !  E
c     & n_AlkylNit= 13,    !  R
c     & n_Alkenes=  14,    !  S
c     & n_Paraffin= 15,    !
c     & n_Terpenes= 16,    ! ---------------
C ----------------------------------------------   
     & n_phot=     2  
      INTEGER, PARAMETER, DIMENSION(12) :: MDOFM =
     & (/31,59,90,120,151,181,212,243,273,304,334,365/)
     
      REAL*8, PARAMETER ::  O3MULT       = 2.14d-2,
     &                      BYO3MULT     = 1./O3MULT,
     &                      T_thresh     = 200.d0,
     &                      pfix_H2      = 560.d-9,
     &                      pfix_Aldehyde= 2.d-9,
     &                      MWabyMWw     = mair/mwat,
     &                      RKBYPIM      = 1.961d2,
     &                      cboltz       = 1.3806d-19,
     &                      zlbatm       = 4.d0,
     &                      CMEQ1        = 0.25d0,
     &                      byradian     = 1.d0/radian,
     &                      cpd          = 1.d6/avog,
     &                      minKG        = 0.d0
#ifdef SHINDELL_STRAT_CHEM
     &                     ,cfc_pppv     = 1722.d-12
     &                     ,n2o_pppv     = 316.3d-9
     &                     ,cfc_rad95    = 794.d-12 
     &                     ,fact_cfc     = cfc_pppv/cfc_rad95
#endif
C Please note: since PCOalt is essentially the nominal 
C pressures for the 23-level GCM, I'm going to use it
C to define BrOx,ClOx,ClONOs,HCL,COIC,OxIC,CFCIC,N2OICX,CH4ICX too:
      REAL*8, PARAMETER, DIMENSION(LCOalt) :: PCOalt = (/
     & 0.9720D+03,0.9445D+03,0.9065D+03,
     & 0.8515D+03,0.7645D+03,0.6400D+03,0.4975D+03,0.3695D+03,
     & 0.2795D+03,0.2185D+03,0.1710D+03,0.1335D+03,0.1016D+03,
     & 0.7120D+02,0.4390D+02,0.2470D+02,0.1390D+02,0.7315D+01,
     & 0.3045D+01,0.9605D+00,0.3030D+00,0.8810D-01,0.1663D-01/)
#ifdef SHINDELL_STRAT_CHEM
      REAL*8, PARAMETER, DIMENSION(LCOalt) ::  
     &     BrOxaltIN = (/1.d-2,1.d-2,1.d-2,1.d-2,1.d-2,1.d-2,1.d-2,
     &     1.d-2,1.d-2,1.d-2,1.d-2,0.12d0,0.12d0,0.12d0,0.12d0,0.06d0,
     &     0.06d0,0.06d0,0.06d0,0.06d0,0.06d0,0.06d0,0.06d0/)
     &     ,ClOxaltIN = (/1.d0,1.d0,1.d0,1.d0,1.d0,1.d0,1.d0,1.d0,1.d0,
     &     1.d0,1.d0,8.d0,8.d0,8.d0,8.d0,8.d1,8.d1,8.d1,8.d1,8.d0,8.d0,
     &     8.d0,8.d0/)
     &     ,ClONO2altIN = (/1.d0,1.d0,1.d0,1.d0,1.d0,1.d0,1.d0,1.d0,
     &     1.d0,1.d0,1.d0,1.d0,1.d0,5.d1,5.d1,5.d1,5.d1,5.d1,5.d1,5.d1,
     &     5.d1,5.d1,5.d1/)
     &     ,HClaltIN = (/1.d0,1.d0,1.d0,1.d0,1.d0,1.d0,1.d0,1.d0,1.d0,
     &     1.d0,1.d0,2.5d1,4.0d1,9.0d1,1.7d2,1.9d2,2.5d2,2.5d2,2.5d2,
     &     2.5d2,2.5d2,2.5d2,2.5d2/)
#endif    
      REAL*8, PARAMETER, DIMENSION(LCH4alt) :: PCH4alt = 
     &                     (/569d0, 150d0, 100d0, 32d0, 3.2d0, 0.23d0/)
      REAL*8, PARAMETER, DIMENSION(LCH4alt) ::   
     &   CH4altINT =(/1.79d0, 1.75d0, 1.620d0,1.460d0,0.812d0,0.230d0/),
     &   CH4altINX =(/1.79d0, 1.75d0, 1.440d0,1.130d0,0.473d0,0.202d0/)    
     
!@dbparam Tpsc_offset_N NH offset for the above T_thresh
!@dbparam Tpsc_offset_S SH offset for the above T_thresh
!@dbparam reg1Power_SpherO2andN2Ocorr first from surface region power of
!@+ cos(sza)^x of spherical correction to ss(27) and ss(28)
!@dbparam reg2Power_SpherO2andN2Ocorr second from surface region power of
!@+ cos(sza)^x of spherical correction to ss(27) and ss(28)
!@dbparam reg3Power_SpherO2andN2Ocorr third from surface region power of
!@+ cos(sza)^x of spherical correction to ss(27) and ss(28)
!@dbparam reg4Power_SpherO2andN2Ocorr fourth and last from surface region power of
!@+ cos(sza)^x of spherical correction to ss(27) and ss(28)
!@dbparam reg1TopPres_SpherO2andN2Ocorr pressure at top of first from surface region
!@+ for spherical correction to ss(27) and ss(28) (hPa)
!@dbparam reg2TopPres_SpherO2andN2Ocorr pressure at top of second from surface region
!@+ for spherical correction to ss(27) and ss(28) (hPa)
!@dbparam reg3TopPres_SpherO2andN2Ocorr pressure at top of third from surface region
!@+ for spherical correction to ss(27) and ss(28) (hPa)
! (fourth = top region needs no upper pressure)
!@dbparam windowO2corr linear correction to ss(27) O2 in window region (in addition to spherical)
!@dbparam windowN2Ocorr linear correction to ss(28) N2O in window region (in addition to spherical)
!@dbparam ch4_init_sh,ch4_init_nh initial methane conc. (ppmv) 
!@+       defaults are for 1990
!@dbparam allowSomeChemReinit (1=YES) to allow some chemistry variables
!@+       to cold-start even if tracers don't. Warning: this includes
!@+       model strat Q( ) spec. hum. reinitialization, and default =1!
!@dbparam fix_CH4_chemistry (1=YES 0=NO) whether or not to used a fixed
!@+       value for methane in the chemistry code. USE -1 for initial
!@+       conditions from file CH4_IC (but L>LS1-1 only now!)
!@+       but use_rad_CH4=1 overrides this.
!@dbparam scale_ch4_IC_file multiplicative factor of CH4 IC if 
!@+       fix_CH4_chemistry=-1 (but only above LS1-1 !)
!@dbparam use_rad_ch4 =1 replaces CH4 surface sources with L=1
!+        overwriting with radiation code values.
!@dbparam use_rad_n2o =1 as ch4 case above
!@dbparam use_rad_cfc =1 as ch4 case above
!@dbparam Lmax_rad_O3 model levels to use tracer Ox in rad code (if on)
!@dbparam Lmax_rad_CH4 model levels to use tracer CH4 in rad code(if on)
!@dbparam which_trop 1=ls1-1 is tropopause, 0=LTROPO(I,J) is tropopause
!@dbparam PI_run used to turn on (1) and off (0) use of PI_ratio*
!@dbparam PIratio_N to scale NOx, HNO3, N2O5, HO2NO2
!@+       initial conditions and stratospheric overwriting.
!@dbparam PIratio_CO_T to scale tropospheric CO IC and overwrite
!@dbparam PIratio_CO_S to scale stratospheric CO IC and overwrite
!@dbparam PIratio_other to scale PAN,Isoprene,AlkyNit,Alkenes,Paraffin
!@+       ,Terpenes
!@+       initial conditions and stratospheric overwriting.
!@dbparam PIratio_N2O preindustrial ratio for N2O ICs and L=1 overwrite
!@dbparam PIratio_CFC preindustrial ratio for CFC ICs and L=1 overwrite
!@+       with model time (JYEAR, JMON, JDAY) 
!@dbparam PltOx for pres<PltOx Ox, NOx, ClOx, and BrOx get overwritten

      INTEGER ::        fix_CH4_chemistry = 0
     &                 ,which_trop        = 0
     &                 ,PI_run            = 0
     &                 ,use_rad_ch4       = 0
     &                 ,use_rad_n2o       = 0
     &                 ,use_rad_cfc       = 0
     &                 ,Lmax_rad_O3       = LM
     &                 ,Lmax_rad_CH4      = LM
     &                 ,checktracer_on    = 0
     &                 ,allowSomeChemReinit = 1
      REAL*8 ::             ch4_init_sh   = 1.750d0,
     &                      ch4_init_nh   = 1.855d0,
     &                      scale_ch4_IC_file= 1.d0, 
     &                      PIratio_N     = 0.667d0,
     &                      PIratio_CO_T  = 0.667d0,
     &                      PIratio_CO_S  = 0.500d0,
     &                      PIratio_other = 0.500d0
#ifdef SHINDELL_STRAT_CHEM
     &                     ,PIratio_N2O   = 0.896d0
     &                     ,PIratio_CFC   = 0.000d0
     &                     ,PltOx         = 0.000d0
     &                     ,Tpsc_offset_N = -10.d0
     &                     ,Tpsc_offset_S = -10.d0
     &                     ,reg1Power_SpherO2andN2Ocorr = 2.0d0
     &                     ,reg2Power_SpherO2andN2Ocorr = 2.0d0
     &                     ,reg3Power_SpherO2andN2Ocorr = 1.0d0
     &                     ,reg4Power_SpherO2andN2Ocorr = 0.5d0
     &                     ,reg1TopPres_SpherO2andN2Ocorr = 50.d0
     &                     ,reg2TopPres_SpherO2andN2Ocorr = 10.d0
     &                     ,reg3TopPres_SpherO2andN2Ocorr = 5.d0
     &                     ,windowN2Ocorr = 0.8d0
     &                     ,windowO2corr  = 0.8d0
     &                     ,PSClatS       = -50.d0
     &                     ,PSClatN       =  50.d0
#endif

      LOGICAL, PARAMETER :: luselb            = .false.

C**************  V  A  R  I  A  B  L  E  S *******************  
!@var nn reactant's number in mol list, first index reactant 1 or 2,
!@+      second - reaction number
!@var nnr reaction product's number in mol list, indicies like nn
!@var nps reaction numbers by molecule, photolytic production
!@var nds reaction numbers by molecule, photolytic destruction
!@var npnr reaction numbers by molecule, photolytic production
!@var ndnr reaction numbers by molecule, photolytic destruction
!@var kps reaction numbers by molecule, chemical production
!@var kds reaction numbers by molecule, chemical destruction
!@var kpnr reaction numbers by molecule, chemical production
!@var kdnr reaction numbers by molecule, chemical destruction
!@var fam ___?
!@var nst reverse reaction number for dissociation reactions
!@var lprn,jprn,iprn l, j, and i point for chemistry debugging
!@var ay name of gas being considered
!@var y concentration of gas, 1st index=gas number, 2nd=verticle level
!@var rr rate constant of chemical reaction, first index - reaction
!@+   number, 2nd is verticle level
!@var ss photodissociation coefficient, indicies; rxn #,L,I,J
!@var pe rate constant for bimolecular chemical reaction
!@var ea activation energy constant for bimolecular chemical reactions
!@var ro,r1,sn,sb rate parameters for trimolecular reactions
!@var conc concentration of optically important gases (O2 & O3), first
!@+   vertivle level, second=gas number (1=O2,2=O3)
!@var TXL temperature profile
!@var prnrts logical: print rate of each chemical reaction?
!@var prnchg logical: print chemical changes?
!@var prnls logical: print reaction lists by species?
!@var yNO3,pHOx,pNOx,pOx,yCH3O2,yC2O3,yROR,yXO2,yAldehyde,yXO2N,yRXPAR?
!@var mNO2 3D vol mixing ratio of NO2 saved for subdaily diagnostics
!@var yCl2,yCl2O2 3D arrays to remember some non-tracer species...
!@var NCFASTJ number of levels in the fastj atmosphere
!@var MIEDX Type of aerosol scattering, currently 6 set up:
!@+   1=Rayly 2=iso 3=iso-equiv 4=bkgrd-sulf,5=volc-sulf,6=liq water
!@var PFASTJ pressure sent to FASTJ
!@var   Rayleigh parameters (effective cross-section) (cm2)
!@var odtmp Optical depth (temporary array)
!@var XLTAU    TTAU along the slant path
!@var XL      Slant path between points
!@var nfam number of beginning molecule of each chemical family
!@var SALBFJ surface albedo parameter from radiation to fastj
!@var OxICIN Ox initial conditions (unit=PPPM,LCOalt levels)
!@var OxICINL column version of OxICIN
!@var COICIN CO initial conditions (unit=PPPM,LCOalt levels)
!@var COICINL column version of OxICIN
!@var N2OICIN N2O initial conditions (unit=PPPM,LCOalt levels)
!@var N2OICINL column version of N2OICIN
!@var CH4ICIN CH4 initial conditions (unit=PPPM,LCOalt levels)
!@var CH4ICINL column version of CH4ICIN
!@var CFCICIN CFC initial conditions (unit=PPPM,LCOalt levels)
!@var CFCICINL column version of CFCICIN
!@var BrOxaltIN altitude dependence BrOx (unitless,LCOalt levels)
!@var ClOxaltIN altitude dependence ClOx (unitless,LCOalt levels)
!@var ClONO2altIN altitude dependence ClONO2 (unitless,LCOalt levels)
!@var HClaltIN altitude dependence HCl (unitless,LCOalt levels)
!@var CH4altINT tropical strat adjustments to CH4 (LCH4alt levels)
!@var CH4altINX xtra-tropical strat adjustments to CH4 LCH4alt levels)
!@var OxIC Ox initial conditions (unit=KG,LM levels)
!@var OxICL column version of OxIC
!@var COIC CO initial conditions (unit=KG,LM levels)
!@var COICL column version of COIC
!@var N2OICX N2O initial conditions (unit=KG,LM levels) X=not Jean's
!@var N2OICL column version of N2OICX
!@var CH4ICX CH4 initial conditions (unit=KG,LM levels) X=not Jean's
!@var CH4ICL column version of CH4ICX
!@var CFCIC CFC initial conditions (unit=KG,LM levels)
!@var CFCICL column version of CFCIC
!@var BrOxalt altitude dependence BrOx (unitless,LM levels)
!@var ClOxalt altitude dependence ClOx (unitless,LM levels)
!@var ClONO2alt altitude dependence ClONO2 (unitless,LM levels)
!@var HClalt altitude dependence HCl (unitless,LM levels)
!@var CH4altT tropical strat adjustments to CH4 (unitless, LM levels)
!@var CH4altX xtra-tropical strat adjustments to CH4 (LM levels)
!@var BYFJM = 1/JM
!@var MODPHOT if MODPHOT=0 do photolysis, else skip it
!@var TX temperature variable for master chem
!@var ta, pres local arrays to hold temperature,pressure
!@var FASTJLAT,FASTJLON latitude & LONGITUDE (degrees) for use in fastj
!@var sulfate N2O5 sulfate sink (formerly SRC(I,J,L,20) variable)   
!@var dms_offline DMS concentration for HOx sink reactions
!@var so2_offline SO2 concentration for HOx conversion reactions
!@var prod_sulfate  N2O5 change by sulfate reactions in mass units
!@var wprod_sulf N2O5 change by sulfate reactions in molecules/cm3/s
!@var DT2 variable chemical time step, set in masterchem
!@var nr total number of        reactions read in from gs_jpl00_trop_15
!@var nr3 #of trimolecular      reactions read in from gs_jpl00_trop_15
!@var nr2 #of mono+bi-molecular reactions read in from gs_jpl00_trop_15
!@var nmm #of monomolecular     reactions read in from gs_jpl00_trop_15
!@var nhet #of heterogenous     reactions read in from gs_jpl00_trop_15
!@var ratioNs,ratioN2,rNO2frac,rNOfrac,rNOdenom variables for nitrogen
!@+   conservation (strat)
!@var chemrate,photrate ?   
!@var MDOFM cumulative days at end of each month
!@var L75P first model level above nominal 75 hPa
!@var L75M first model level below nominal 75 hPa
!@var F75P interpolation coeff. of higher altitude value (units ln(P))
!@var F75M interpolation coeff. of lower altitude value (units ln(P))
!@var L569P first model level above nominal 569 hPa
!@var L569M first model level below nominal 569 hPa
!@var F569P interpolation coeff. of higher altitude value (units ln(P))
!@var F569M interpolation coeff. of lower altitude value (units ln(P))
!@var DU_O3 total column ozone in latitude band
!@var SF3 is H2O photolysis in Schumann-Runge Bands
!@var SF2 is NO photolysis in Schumann-Runge Bands
!@var Jacet photolysis rate for acetone (not done through fastj)
!@var acetone 3D acetone mixing ratio (static for now)
!@var pscX column logical for the existance of polar strat clouds(PSCs)
!@var sOx_acc accumulated SURFACE ozone (Ox) (special for SUBDD)
!@var sNOx_acc accumulated SURFACE NOx (special for SUBDD)
!@var sCO_acc accumulated SURFACE CO (special for SUBDD)
!@var l1Ox_acc accumulated L=1 ozone (Ox) (special for SUBDD)
!@var l1NO2_acc accumulated L=1 NO2 (special for SUBDD)
!@var save_NO2column instantaneous NO2 column (for SUBDD exporting)
!@var RGAMMASULF N2O5-->HNO3 conversion on aerosols?
      INTEGER :: nr,nr2,nr3,nmm,nhet,MODPHOT,L75P,L75M,L569P,L569M,
     &lprn,jprn,iprn,MIEDX,NCFASTJ
#ifdef SHINDELL_STRAT_CHEM
      INTEGER, DIMENSION(n_fam)        :: nfam = 
     &     (/37+ntm_terp+ntm_soa,40+ntm_terp+ntm_soa,
     &       44+ntm_terp+ntm_soa,50+ntm_terp+ntm_soa,0/)
#else
      INTEGER, DIMENSION(n_fam)        :: nfam = 
     &     (/27+ntm_terp+ntm_soa,30+ntm_terp+ntm_soa,0,0/)
#endif
      INTEGER, DIMENSION(p_1,p_2)      :: nn, nnr
      INTEGER, DIMENSION(p_3)          :: nps, nds, npnr, ndnr
      INTEGER, DIMENSION(p_4)          :: kps, kds, kpnr, kdnr
      INTEGER, DIMENSION(n_nst)        :: nst

C**************  Latitude-Dependant (allocatable) *******************
      REAL*8, ALLOCATABLE, DIMENSION(:)       :: DU_O3
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:)   :: acetone
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:,:) :: ss
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:)   :: yNO3,pHOx,pNOx,pOx,
     & yCH3O2,yC2O3,yROR,yXO2,yAldehyde,yXO2N,yRXPAR,TX,sulfate,OxIC,
     & CH4ICX,dms_offline,so2_offline,yso2,ydms,mNO2,COIC,pNO3
#ifdef SHINDELL_STRAT_CHEM
     & ,pClOx,pClx,pOClOx,pBrOx,yCl2,yCl2O2,N2OICX,CFCIC,SF3,SF2
#endif
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:):: COICIN,OxICIN,CH4ICIN
#ifdef SHINDELL_STRAT_CHEM
     &                                       ,N2OICIN,CFCICIN
#endif
      REAL*8, ALLOCATABLE, DIMENSION(:,:):: sOx_acc,sNOx_acc,sCO_acc,
     & l1Ox_acc,l1NO2_acc,save_NO2column

C**************  Not Latitude-Dependant ****************************      
      REAL*8 :: XLTAU,BYFJM,
     & FASTJLAT,FASTJLON,DT2,F75P,F75M,F569P,F569M,RGAMMASULF
#ifdef SHINDELL_STRAT_CHEM
     & ,ratioNs,ratioN2,rNO2frac,rNOfrac,rNOdenom
#endif
      REAL*8, DIMENSION(nc,LM)         :: y
      REAL*8, DIMENSION(n_rx,LM)       :: rr
      REAL*8, DIMENSION(n_bi)          :: pe, ea
      REAL*8, DIMENSION(n_tri)         :: ro, r1, sn, sb
      REAL*8, DIMENSION(LM)            :: odtmp,ta,pres,Jacet
      REAL*8, DIMENSION(p_2,LM)        :: chemrate, photrate
      REAL*8, DIMENSION(ny,LM)         :: dest, prod
      REAL*8, DIMENSION(LCOalt)        :: COICINL,OxICINL,CH4ICINL
#ifdef SHINDELL_STRAT_CHEM
     &                                   ,N2OICINL,CFCICINL
#endif
      REAL*8, DIMENSION(LM)  :: CH4altT,CH4altX,COICL,OxICL,CH4ICL
#ifdef SHINDELL_STRAT_CHEM
     &                        ,BrOxalt,ClOxalt,ClONO2alt,HClalt
     &                        ,N2OICL,CFCICL,OxlossbyH
#endif

      LOGICAL                      :: fam,prnrts,prnchg,prnls      
#ifdef SHINDELL_STRAT_CHEM
      LOGICAL, DIMENSION(LM)       :: pscX
#endif

      CHARACTER*8, DIMENSION(nc)   :: ay
      
#endif
      END MODULE TRCHEM_Shindell_COM
      
      
      
      subroutine alloc_trchem_shindell_com(grid)
!@SUM  To allocate arrays whose sizes now need to be determined
!@+    at run-time
!@auth G.Faluvegi
!@ver  1.0
      use domain_decomp_atm, only : dist_grid, get
      use model_com, only     : im,lm
      use TRCHEM_Shindell_COM, only: DU_O3,ss,yNO3,sOx_acc,l1Ox_acc,
     & pHOx,pNOx,pOx,yCH3O2,yC2O3,yROR,yXO2,yAldehyde,yXO2N,yRXPAR,
     & TX,sulfate,COIC,OxIC,CH4ICX,dms_offline,so2_offline,yso2,ydms,
     & COICIN,OxICIN,CH4ICIN,JPPJ_Shindell,LCOalt,acetone,mNO2,
     & l1NO2_acc,sNOx_acc,sCO_acc,save_NO2column,pNO3
#ifdef SHINDELL_STRAT_CHEM
     & ,pClOx,pClx,pOClOx,pBrOx,yCl2,yCl2O2,N2OICX,CFCIC,SF3,SF2,
     & N2OICIN,CFCICIN
#endif

      IMPLICIT NONE

      type (dist_grid), intent(in) :: grid
      integer :: ier, J_1H, J_0H, I_1H, I_0H
      logical :: init = .false.

      if(init)return
      init=.true.
    
      call get( grid , J_STRT_HALO=J_0H, J_STOP_HALO=J_1H )
      I_0H=GRID%I_STRT_HALO
      I_1H=GRID%I_STOP_HALO
 
      allocate(          ss(JPPJ_Shindell,LM,I_0H:I_1H,J_0H:J_1H) )
      allocate(     acetone(I_0H:I_1H,J_0H:J_1H,LM)      )
      allocate(        yNO3(I_0H:I_1H,J_0H:J_1H,LM)      )
      allocate(        mNO2(I_0H:I_1H,J_0H:J_1H,LM)      )
      allocate(        pHOx(I_0H:I_1H,J_0H:J_1H,LM)      )
      allocate(        pNOx(I_0H:I_1H,J_0H:J_1H,LM)      )
      allocate(        pNO3(I_0H:I_1H,J_0H:J_1H,LM)      )
      allocate(         pOx(I_0H:I_1H,J_0H:J_1H,LM)      )
      allocate(      yCH3O2(I_0H:I_1H,J_0H:J_1H,LM)      )
      allocate(       yC2O3(I_0H:I_1H,J_0H:J_1H,LM)      )
      allocate(        yROR(I_0H:I_1H,J_0H:J_1H,LM)      )
      allocate(        yXO2(I_0H:I_1H,J_0H:J_1H,LM)      )
      allocate(   yAldehyde(I_0H:I_1H,J_0H:J_1H,LM)      )
      allocate(       yXO2N(I_0H:I_1H,J_0H:J_1H,LM)      )
      allocate(      yRXPAR(I_0H:I_1H,J_0H:J_1H,LM)      )
      allocate(          TX(I_0H:I_1H,J_0H:J_1H,LM)      )
      allocate(     sulfate(I_0H:I_1H,J_0H:J_1H,LM)      )
      allocate(        OxIC(I_0H:I_1H,J_0H:J_1H,LM)      )
      allocate(        COIC(I_0H:I_1H,J_0H:J_1H,LM)      )
      allocate(      CH4ICX(I_0H:I_1H,J_0H:J_1H,LM)      )
      allocate( dms_offline(I_0H:I_1H,J_0H:J_1H,LM)      )
      allocate( so2_offline(I_0H:I_1H,J_0H:J_1H,LM)      )
      allocate(        yso2(I_0H:I_1H,J_0H:J_1H,LM)      )
      allocate(        ydms(I_0H:I_1H,J_0H:J_1H,LM)      )
      allocate(      OxICIN(I_0H:I_1H,J_0H:J_1H,LCOalt)  )
      allocate(      COICIN(I_0H:I_1H,J_0H:J_1H,LCOalt)  )
      allocate(     CH4ICIN(I_0H:I_1H,J_0H:J_1H,LCOalt)  )
      allocate(     sOx_acc(I_0H:I_1H,J_0H:J_1H)         )
      allocate(     sNOx_acc(I_0H:I_1H,J_0H:J_1H)        )
      allocate(     sCO_acc(I_0H:I_1H,J_0H:J_1H)         )
      allocate(    l1Ox_acc(I_0H:I_1H,J_0H:J_1H)         )
      allocate(    l1NO2_acc(I_0H:I_1H,J_0H:J_1H)        )
      allocate(save_NO2column(I_0H:I_1H,J_0H:J_1H)       )

      sOx_acc=0.; sNOx_acc=0.; sCO_acc=0.; l1Ox_acc=0. ; l1NO2_acc=0.

      allocate(       DU_O3(J_0H:J_1H)                   )
#ifdef SHINDELL_STRAT_CHEM
      allocate(       pClOx(I_0H:I_1H,J_0H:J_1H,LM)      )
      allocate(        pClx(I_0H:I_1H,J_0H:J_1H,LM)      )
      allocate(      pOClOx(I_0H:I_1H,J_0H:J_1H,LM)      ) 
      allocate(       pBrOx(I_0H:I_1H,J_0H:J_1H,LM)      )
      allocate(        yCl2(I_0H:I_1H,J_0H:J_1H,LM)      ) 
      allocate(      yCl2O2(I_0H:I_1H,J_0H:J_1H,LM)      )
      allocate(      N2OICX(I_0H:I_1H,J_0H:J_1H,LM)      )
      allocate(       CFCIC(I_0H:I_1H,J_0H:J_1H,LM)      )
      allocate(         SF3(I_0H:I_1H,J_0H:J_1H,LM)      )
      allocate(         SF2(I_0H:I_1H,J_0H:J_1H,LM)      )
      allocate(     N2OICIN(I_0H:I_1H,J_0H:J_1H,LCOalt)  )
      allocate(     CFCICIN(I_0H:I_1H,J_0H:J_1H,LCOalt)  )
#endif
      
      return
      end subroutine alloc_trchem_shindell_com
      
