E000TRACER.R GISS Model E                                gas 06/00


Preprocessor Options
#define TRACERS_ON                  ! include tracers code
#define TRACERS_WATER               ! tracers can interact with water
#define TRACERS_DRYDEP              ! include tracer dry deposition
#define TRACERS_SPECIAL_Shindell    ! includes drew's chemical tracers
#define TRACERS_DUST
#define TRACERS_NITRATE
#define TRACERS_HETCHEM
#define TRACERS_AEROSOLS_Koch
#define EDGAR_HYDE_SOURCES       ! use EDGAR-HYDE tracers sources instead
!  OFF #define SHINDELL_STRAT_CHEM      ! turns on stratospheric chemistry
!  OFF #define regional_Ox_tracers      ! turns on regional Ox tracers   
!  OFF #define INTERACTIVE_WETLANDS_CH4 ! turns on interactive CH4 wetland source
End Preprocessor Options

Object modules: (in order of decreasing priority)
RES_M23                             ! horiz/vert resolution
MODEL_COM GEOM_B IORSF              ! model variables and geometry
MODELE                              ! Main and model overhead
PARAM PARSER                        ! parameter database
DOMAIN_DECOMP ALLOC_DRV             ! domain decomposition, allocate global distributed arrays
ATMDYN_COM ATMDYN MOMEN2ND          ! atmospheric dynamics
STRATDYN STRAT_DIAG                 ! strospheric dynamics (incl. gw drag)
QUS_COM QUSDEF QUS_DRV              ! advection of tracers
TQUS_DRV                            ! advection of Q
! ---TRACER SPECIFIC CODES----------
TRACER_COM TRACERS_DRV              ! common and driver for tracers
TRACERS_SPECIAL_Shindell            ! routines specific to drew's 15-tracers
TRCHEM_Shindell_COM                 ! Drew Shindell's tracers common
TRCHEM_calc                         ! chemical reaction calculations
TRCHEM_init                         ! chemistry initialization, I/O
TRCHEM_family                       ! tracer family chemistry
TRCHEM_fastj                        ! tracer chem photlysis code/rad transf
TRCHEM_master                       ! trop chem "driver"/strat prescrioption
TRACERS_AEROSOLS_Koch_e4            ! Koch aerosols
TRDUST TRDUST_COM TRACERS_DUST      ! Dust aerosols
TRACER_NITRATE                      ! Nitrate aerosol chemistry
TRACER_HETCHEM                      ! Heterogeneous aerosol chemistry 
! ----------------------------------
TRDRYDEP                            ! tracer dry deposition from Harvard CTM
TRACERS                             ! generic tracer code
TRDIAG_COM TRACER_PRT               ! tracer diagnostic printout
CLOUDS2 CLOUDS2_DRV CLOUDS_COM      ! clouds modules
SURFACE FLUXES                      ! surface calculation and fluxes
GHY_COM GHY_DRV GHY                 ! land surface and soils
VEG_DRV VEG_COM VEGETATION          ! vegetation
PBL_COM PBL_DRV PBL                 ! atmospheric pbl
! pick exactly one of the next 2 choices: ATURB or DRYCNV
ATURB                               ! turbulence in whole atmosphere
! DRYCNV                            ! drycnv
LAKES_COM LAKES                     ! lake modules
SEAICE SEAICE_DRV                   ! seaice modules
LANDICE LANDICE_DRV                 ! land ice modules
ICEDYN_DUM ! or: ICEDYN ICEDYN_DRV  ! dynamic ice modules
OCEAN OCNML                         ! ocean modules
SNOW_DRV SNOW                       ! snow model
RAD_COM RAD_DRV RADIATION           ! radiation modules
DIAG_COM DIAG DEFACC DIAG_PRT       ! diagnostics
CONST FFT72 UTILDBL SYSTEM          ! utilities
POUT_netcdf                         ! post-processing output

Data input files:
!-----------------------------------------------
!
!    GCM INPUT
!
AIC=AIC.RES_M23.D771201
GIC=GIC.E046D3M20A.1DEC1955
OCNML=Z1O.B4X5.cor ! needed for post-processing only
OSST=OST4X5.B.1975-84avg.Hadl1.1
SICE=SICE4X5.B.1975-84avg.Hadl1.1
CDN=CD4X500S
VEG=V72X46.1.cor2_no_crops ! VEG=V72X46.1.cor2
CROPS=CROPS_72X46N.cor4
SOIL=S4X50093 TOPO=Z72X46N.cor4_nocasp ! bdy.cond
REG=REG4X5           ! special regions-diag
RVR=RD4X525.RVR.2    ! river direction file
ZVAR=ZVAR4X5         ! topographic variation for gwdrag
RADN1=sgpgxg.table8    ! rad.tables
RADN2=radfil33k                   !     8/2003 version
RADN3=miescatpar.abcdv2
MSU_wts=MSU.RSS.weights.data
TAero_PRE=sep2003_PRE_Koch_kg_m2_ChinSEA_Liao_1850 ! pre-industr trop. aerosols
TAero_SUI=sep2003_SUI_Koch_kg_m2_72x46x9_1875-1990 ! industrial sulfates
TAero_OCI=sep2003_OCI_Koch_kg_m2_72x46x9_1875-1990 ! industrial organic carbons
TAero_BCI=sep2003_BCI_Koch_kg_m2_72x46x9_1875-1990 ! industrial black carbons
RH_QG_Mie=oct2003.relhum.nr.Q633G633.table
RADN6=dust8.tau9x8x13
RADN7=STRATAER.VOL.1850-1999.Apr02
RADN8=cloud.epsilon4.72x46
RADN9=solar.lean02.ann.uvflux    ! need KSOLAR=2
RADNE=topcld.trscat8
! new ozone files (minimum 1, maximum 9 files)
O3file_01=mar2004_o3_shindelltrop_72x46x49x12_1850
O3file_02=mar2004_o3_shindelltrop_72x46x49x12_1890
O3file_03=mar2004_o3_shindelltrop_72x46x49x12_1910
O3file_04=mar2004_o3_shindelltrop_72x46x49x12_1930
O3file_05=mar2004_o3_shindelltrop_72x46x49x12_1950
O3file_06=mar2004_o3_shindelltrop_72x46x49x12_1960
O3file_07=mar2004_o3_shindelltrop_72x46x49x12_1970
O3file_08=mar2004_o3_shindelltrop_72x46x49x12_1980
O3file_09=mar2004_o3_shindelltrop_72x46x49x12_1990
O3trend=mar2004_o3timetrend_46x49x2412_1850_2050
GHG=GHG.Mar2004.txt
dH2O=dH2O_by_CH4_monthly
BC_dep=BC.Dry+Wet.depositions.ann
TOP_INDEX=top_index_72x46.ij.ext
!-----------------------------------------------
!
!  TROPOSPHERIC CHEMISTRY INPUT	
!
! choose these for trop-only chem model
!-----------------------------------------------
MOLEC=chem_files/su_dust_moleculesE
JPLRX=chem_files/gs_jpl00_trop_15_fix
JPLPH=chem_files/ds_photlist_trop_15
RATJ=chem_files/ratj.giss_15
SPECFJ=chem_files/jv_spec00_15.dat
!-----------------------------------------------
ATMFJ=chem_files/jv_atms.dat
DRYCOEFF=chem_files/drydep.coef
VEGTYPE=chem_files/vegtype.global
OLSON=chem_files/drydep.table
LAI01=chem_files/lai01.global
LAI02=chem_files/lai02.global
LAI03=chem_files/lai03.global
LAI04=chem_files/lai04.global
LAI05=chem_files/lai05.global
LAI06=chem_files/lai06.global
LAI07=chem_files/lai07.global
LAI08=chem_files/lai08.global
LAI09=chem_files/lai09.global
LAI10=chem_files/lai10.global
LAI11=chem_files/lai11.global
LAI12=chem_files/lai12.global
Ox_IC=gsin/Ox_init_cond_M23_4x5 !see README in /usr/people/cmrun/gsin
! next one needed only if correct_strat_Ox=.true.
Ox_corr=gsin/corrOx_modelE_v4
!
!----------NORMAL-CASE--------------------------
!-----------------------------------------------
! keep next two on for edgar-hyde and normal case:
Isoprene_VEGETATION=gsin/isoprene_vegetation_4x5 !monthly  KG/hr/grid
SULFATE_SA=NOy_sinks/sulfate_fakeM23_M_SA
DMS_FIELD=dms_conc
SO2_FIELD=so2_conc
! SULFATE_SA=NOy_sinks/sulfate_mass_EA07TnbM23
!---------EDGAR-HYDE-case------------------------
CO_FC_EH=gsin/EDGAR_HYDE/CO_FC1990.dat
CO_IN_EH=gsin/EDGAR_HYDE/CO_IN1990.dat
CO_AW_EH=gsin/EDGAR_HYDE/CO_AW1990_mon.dat
CO_SB_EH=gsin/EDGAR_HYDE/CO_SB1990_mon.dat! these for edgar-hyde
CO_BF_EH=gsin/EDGAR_HYDE/CO_BC1990_mon.dat
CO_DF_EH=gsin/EDGAR_HYDE/CO_DF1990_mon.dat
Alkenes_VEGETATION=gsin/alkenes_vegetation_4x5   !monthly  KG/hr/grid
Alkenes_FC_EH=gsin/EDGAR_HYDE/ALKFC1990.dat
Alkenes_IN_EH=gsin/EDGAR_HYDE/ALKIN1990.dat
Alkenes_FP_EH=gsin/EDGAR_HYDE/ALKFP1990.dat  !these for edgar-hyde
Alkenes_AW_EH=gsin/EDGAR_HYDE/ALKAW1990_mon.dat
Alkenes_SB_EH=gsin/EDGAR_HYDE/ALKSB1990_mon.dat
Alkenes_BF_EH=gsin/EDGAR_HYDE/ALKBC1990_mon.dat
Alkenes_DF_EH=gsin/EDGAR_HYDE/ALKDF1990_mon.dat
Paraffin_VEGETATION=gsin/paraffin_vegetation_4x5 !monthly  KG/hr/grid
Paraffin_FC_EH=gsin/EDGAR_HYDE/PARFC1990.dat
Paraffin_IN_EH=gsin/EDGAR_HYDE/PARIN1990.dat
Paraffin_FP_EH=gsin/EDGAR_HYDE/PARFP1990.dat  !these for edgar-hyde
Paraffin_AW_EH=gsin/EDGAR_HYDE/PARAW1990_mon.dat 
Paraffin_SB_EH=gsin/EDGAR_HYDE/PARSB1990_mon.dat
Paraffin_BF_EH=gsin/EDGAR_HYDE/PARBC1990_mon.dat
Paraffin_DF_EH=gsin/EDGAR_HYDE/PARDF1990_mon.dat
NOx_FC_EH=gsin/EDGAR_HYDE/NOXFC1990.dat 
NOx_IN_EH=gsin/EDGAR_HYDE/NOXIN1990.dat
NOx_AL_EH=gsin/EDGAR_HYDE/NOXAL1990_mon.dat
NOx_AW_EH=gsin/EDGAR_HYDE/NOXAW1990_mon.dat
NOx_BF_EH=gsin/EDGAR_HYDE/NOXBC1990_mon.dat
NOx_DF_EH=gsin/EDGAR_HYDE/NOXDF1990_mon.dat
NOx_SB_EH=gsin/EDGAR_HYDE/NOXSB1990_mon.dat
NOx_AIRCRAFT=NOy_sources/aircraft_4x5_1990 !!!!remember this one
CH4_ANIMALS_EH=gsin/EDGAR_HYDE/CH4AN1990.dat
CH4_F_FUEL_C_EH=gsin/EDGAR_HYDE/CH4FC1990.dat
CH4_LANDFILL_EH=gsin/EDGAR_HYDE/CH4LF1990.dat
CH4_F_FUEL_P_EH=gsin/EDGAR_HYDE/CH4FP1990.dat
CH4_SOIL_ABS=methane/gcm_data/CH4SOILABS_4X5   !Annual 1.194
CH4_TERMITES=methane/gcm_data/CH4TRMITE_4X5    !Annual 0.999
CH4_WETL=methane/gcm_data/CH4WETL+TUNDRA_4X5  !Monthly 0.9818; zonal also
CH4_AW_EH=gsin/EDGAR_HYDE/CH4AW1990_mon.dat
CH4_DF_EH=gsin/EDGAR_HYDE/CH4DF1990_mon.dat
CH4_SB_EH=gsin/EDGAR_HYDE/CH4SB1990_mon.dat
CH4_BF_EH=gsin/EDGAR_HYDE/CH4BC1990_mon.dat
CH4_AL_EH=gsin/EDGAR_HYDE/CH4AL1990_mon.dat
!-----------------------------------------------
!
! AEROSOL MODELE INPUT
!
DMS_SEA=DMS.dat
AER_CHEM=Sulf_chem_drewE_20
AER_OH_STRAT=Strat_OH_drewE_20
SO2_IND=SO2_ind2000.AEROCOM_DEC03
SO2_INDh=SO2.1875-1990.Lefohn  !vanAardenne
SO2_BIOMASS=bioburn
SO2_VOLCANO=SO2_volc_conti2000.AEROCOM_FEB12
AIRCRAFT=MM_fuel_2015_subsonic_M20
BC_FOSSIL_FUEL=BC_ff2000.AEROCOM_DEC03
BC_BIOMASS=BC_Biomass_1997-2001
OC_BIOMASS=OM_Biomass_1997-2001
BC_INDh=BC_OC.1875-1990.Jan05
OC_BIOFUEL=OC_bf2000.AEROCOM_DEC03
OC_FOSSIL_FUEL=OC_ff2000.AEROCOM_DEC03
TERPENE=terp
NH3SOURCE=GISS_NH3.4x5
O3_FIELD=Ox_3D_field_bell
!-----------------------------------------------
!
! DUST MODEL INPUT 
!
VTRSH=vtr-mod-o0.mean-pb
FRCLAY=claygcm-f
FRSILT=siltgcm-f
DRYHR=text5hr-f
GIN=Ginoux_dstsrc
LKTAB=table_emission
ERS=ERS1_1993_MONTHLY
LKTAB1=table_wspdf
OFFDUST=dust.1296.c4s4.12months-bg   ! Monthly offline dust
dust_bin1=DUST_bin1_2000_new.nc      ! AEROCOM prescribed dust emission; 4 bins
dust_bin2=DUST_bin2_2000_new.nc
dust_bin3=DUST_bin3_2000_new.nc
dust_bin4=DUST_bin4_2000_new.nc

!----------------------------------------------
!      nudging input 
!u0.nc=uwnd.1995.GISS4x5.nc
!v0.nc=vwnd.1995.GISS4x5.nc
!u1.nc=uwnd.2000.GISS4x5_MANIP.nc
!v1.nc=vwnd.2000.GISS4x5_MANIP.nc
!u2.nc=uwnd.2000.GISS4x5_MANIP.nc
!v2.nc=vwnd.2000.GISS4x5_MANIP.nc

Label and Namelist:
E000TRACER (PD control, corrected Alkenes, Paraffin doubling)
R=00BG/B
DTFIX=300
&&PARAMETERS

ocn_cycl=1      ! =0 if ocean varies from year to year

X_SDRAG=.00025,.000025  ! used for lin. sdrag above P_SDRAG mb
C_SDRAG=0.     ! no constant sdrag
P_SDRAG=.1     ! lin. sdrag above .1mb (top 2 layers) except near poles
PP_SDRAG=.1    ! lin. sdrag above 1.mb near poles (top 4 layers)
ANG_SDRAG=1    ! if =1: sdrag conserves ang mom.
PBREAK = 200.  ! The level for GW breaking above.
DEFTHRESH=0.000035 !the default is 15d-6
PCONPEN=500.   ! penetrating convection defn for GWDRAG
CMC = 0.0000003

KOCEAN=0
U00ice  = 0.59   ! tune this first to get reas.alb/cldcvr (range: .4-.6), then
! u00wtrx = 1.4  ! 1880 conditions 
u00wtrx = 1.404    ! 1980 conditions (NEEDS TO BE VERIFIED!)
cond_scheme=2  ! more elaborate conduction scheme (GHY, Nancy Kiang)

H2ObyCH4=1.  ! activates strat.H2O generated by CH4
KSIALB=0        ! 6-band albedo (Hansen) (=1 A.Lacis orig. 6-band alb)
KSOLAR=2

! for setting fixed methane value for chemistry:
fix_CH4_chemistry=0
pfix_CH4_S=1.700d-6
pfix_CH4_N=1.850d-6

! parameters that control the Shapiro filter
DT_XUfilter=450. ! Shapiro filter on U in E-W direction; usually same as DT (below)
DT_XVfilter=450. ! Shapiro filter on V in E-W direction; usually same as DT (below)
DT_YVfilter=0.   ! Shapiro filter on V in N-S direction
DT_YUfilter=0.   ! Shapiro filter on U in N-S direction

LMCM=16              ! max level of moist convection
XCDNST=300.,10000.   ! strat. gw drag parameters
DTsrc = 1800.        ! half-hour physics time step (default: DTsrc=3600.)
DT=450.,             ! dynamics time step
NIsurf=1,            ! number of surface time steps

NSUBDD=0        ! saving sub-daily diags
Kvflxo=0        ! saving VFLXO (daily)
Ndisk=240       ! use =240 on halem 24 on ra/daley
! pick sub-daily frequency diags, next line ::
!SUBDD='UALL VALL TALL RALL SAT PS STX STY QSEN QLAT Z1000 SWD SWU LWD'
SUBDD=' '
!SUBDD1='QS SNOWD ICEF PREC '
SUBDD1=' '
LmaxSUBDD=18    ! only need to save up to level LmaxSUBDD
KCOPY=2         ! saving acc + rsf
isccp_diags=0   ! use =0 to save cpu time
!-----------------
COUPLED_CHEM=1     ! to couple chemistry and aerosols
rad_interact_tr=0 ! 1=use calculated Ox in radiation, 0=use climatology
                   ! (either case does the rad-forcing calculation)
imAER=0            ! MODEL AS IT IS AEROSOL SOURCES = 0 , AEROCOM SOURCES = 1
imDust=0           ! =1 AEROCOM dust emissions (default=0)
diag_rad=0         ! switches on (=1) comprehensive rad diags (default 0:off)
adiurn_dust = 0    ! no daily dust output
rad_forc_lev=1     ! use LTROPO(I,J) level for rad forcing diags.
NIPRNT=1           ! number of initial print-outs
prather_limits=1   ! to avoid some negative tracers in sub-gridscale
!-----------------
! parameters that control the atmospheric/boundary conditions
! if set to 0, the current (day/) year is used: transient run
crops_yr=-1 ! if -1, crops in VEG-file is used
s0_yr=1979 ! 1990
s0_day=182
ghg_yr=1979 ! 1990
ghg_day=182
volc_yr=1979 ! 1990
volc_day=182
aero_yr=1979 ! 1990
o3_yr=1995 ! 1990 !1979
&&END_PARAMETERS

 &INPUTZ
   QCHECK=.false.
   kdiag = 9,9,9,9,9,9,9,0,9,9,9,9,9
   YEARI=1981,MONTHI=1,DATEI=29,HOURI=0,
   YEARE=1987,MONTHE=3,DATEE=1,HOURE=1,
   ISTART=2,IRANDI=0, YEARE=1981,MONTHE=1,HOURE=1,DATEE=29,IWRITE=1,JWRITE=1,
 &END
