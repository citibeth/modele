E_AR5_V2_CAMP_F79.R GISS Model E  1850 ocn/atm     Larissa    07/15/2014

E_AR5_V2_CAMP_F79: E_AR5_V2_CADI + MATRIX aerosol microphysics;
               U00a=0.61, U00b=1.0, WMUI_multiplier=1.

E4F40: modelE as frozen (or not yet) in July 2009
modelE4 2x2.5 hor. grid with 40 lyrs, top at .1 mb (+ 3 rad.lyrs)
atmospheric composition from year 1850
ocean data: prescribed, 1876-1885 climatology
uses turbulence scheme (no dry conv), grav.wave drag
time steps: dynamics 3.75 min leap frog; physics 30 min.; radiation 2.5 hrs
filters: U,V in E-W and N-S direction (after every physics time step)
         U,V in E-W direction near poles (after every dynamics time step)
         sea level pressure (after every physics time step)

Preprocessor Options
#define NEW_IO
#define RAD_O3_GCM_HRES
#define RAD_O3_DECADAL_INPUT
#define TRAC_ADV_CPU
#define USE_ENT                  ! include dynamic vegetation model
#define TRACERS_ON               ! include tracers code
#define TRACERS_WATER            ! wet deposition and water tracer
!  OFF  #define TRACERS_DUST             ! include dust tracers
!  OFF  #define TRACERS_DUST_Silt4       ! include 4th silt size class of dust
#define TRACERS_DRYDEP           ! default dry deposition
#define TRDIAG_WETDEPO           ! additional wet deposition diags for tracers
#define NO_HDIURN                ! exclude hdiurn diagnostics
#define TRACERS_SPECIAL_Shindell    ! includes drew's chemical tracers
#define SHINDELL_STRAT_CHEM         ! turns on stratospheric chemistry
!  OFF #define AUXILIARY_OX_RADF ! radf diags for climatology or tracer Ozone
#define TRACERS_TERP                ! include terpenes in gas-phase chemistry
#define BIOGENIC_EMISSIONS       ! turns on interactive isoprene emissions
!  OFF #define TRACERS_AEROSOLS_Koch    ! Dorothy Koch's tracers (aerosols, etc)
!  OFF #define TRACERS_AEROSOLS_SOA     ! Secondary Organic Aerosols
!  OFF #define SOA_DIAGS                ! Additional diagnostics for SOA
!  OFF #define TRACERS_NITRATE
!  OFF #define TRACERS_HETCHEM
#define BC_ALB                      ! optional tracer BC affects snow albedo
#define TRACERS_AMP
#define TRACERS_AMP_M1
#define CLD_AER_CDNC                ! aerosol-cloud interactions
#define BLK_2MOM                    ! aerosol-cloud interactions
!  OFF #define WATER_MISC_GRND_CH4_SRC ! adds lake, ocean, misc. ground sources for CH4
!  OFF #define CALCULATE_FLAMMABILITY  ! activated code to determine flammability of surface veg
!  OFF #define DYNAMIC_BIOMASS_BURNING  ! alter biomas burning my flammability
!  OFF #define CALCULATE_LIGHTNING ! turn on Colin Price lightning when TRACERS_SPECIAL_Shindell off
!  OFF #define SHINDELL_STRAT_EXTRA     ! non-chemistry stratospheric tracers
!  OFF #define INTERACTIVE_WETLANDS_CH4 ! turns on interactive CH4 wetland source
!  OFF #define NUDGE_ON                 ! nudge the meteorology
!  OFF #define GFED_3D_BIOMASS          ! turns on IIASA AR4 GFED biomass burning
!  OFF #define HTAP_LIKE_DIAGS    ! adds many diags, changes OH diag, adds Air tracer
!  OFF #define ACCMIP_LIKE_DIAGS  ! adds many diags as defined by ACCMIP project
End Preprocessor Options

Object modules:
     ! resolution-specific source codes
RES_F79                             ! horiz/vert resolution, 2x2.5, top at 0.1mb, 40 layers
DIAG_RES_F                          ! diagnostics
FFT144                              ! Fast Fourier Transform

    ! lat-lon grid specific source codes
GEOM_B                              ! model geometry
DIAG_ZONAL GCDIAGb                  ! grid-dependent code for lat-circle diags
DIAG_PRT POUT_netcdf                ! diagn/post-processing output
IO_DRV                              ! new i/o

     ! GISS dynamics with gravity wave drag
ATMDYN MOMEN2ND                     ! atmospheric dynamics
QUS_DRV                             ! advection of T

QUS3D                               ! advection of Q and tracers
TRDUST_COM TRDUST TRDUST_DRV        ! dust tracer specific code
! Codes common to atmospheric tracer sets
TRACER_COM                          ! configurable tracer code
TRACERS_DRV |-O0|                   ! O0 speeds compilation and no difference for this file
TRACERS                             ! generic tracer code
TRDRYDEP                            ! dry deposition of tracers
TRDIAG_COM TRACER_PRT               ! tracer diagnostic printout
! ---TRACER SPECIFIC CODES----------
TRACERS_SPECIAL_Shindell            ! routines specific to drew's 15-tracers
TRCHEM_Shindell_COM                 ! Drew Shindell's tracers common
TRCHEM_calc                         ! chemical reaction calculations
TRCHEM_init                         ! chemistry initialization, I/O
TRCHEM_family                       ! tracer family chemistry
TRCHEM_fastj2                       ! used for trop+strat chem version
TRCHEM_master                       ! trop chem "driver"/strat prescrioption
BIOGENIC_EMISSIONS                  ! old N.Unger interactive isoprene
!Aerosol Micro Physics
TRAMP_drv        |-extend_source  |  
TRAMP_actv       |-extend_source  |  
TRAMP_diam       |-extend_source  | 
TRAMP_nomicrophysics |-extend_source  |  
TRAMP_subs       |-extend_source  |  
TRAMP_coag       |-extend_source  |  
TRAMP_depv       |-extend_source  | 
TRAMP_param_GISS |-extend_source  |    
TRAMP_config  
TRAMP_dicrete    |-extend_source  |    
TRAMP_init       |-extend_source  |  
TRAMP_quad       |-extend_source  |  
TRAMP_matrix     |-extend_source  |        
TRAMP_setup      |-extend_source  |  
TRAMP_npf        |-extend_source  |  
TRAMP_rad        |-extend_source  |
! When using ISORROPIA Thermodynamics
!TRAMP_thermo_isorr2 |-extend_source  | 
!TRAMP_isocom2              
!TRAMP_isofwd2          
!TRAMP_isorev2
!isrpia.inc
! When using EQSAM Thermodynamics
TRAMP_thermo_eqsam |-extend_source  |  
TRAMP_eqsam_v03d
TRACERS_AEROSOLS_Koch_e4
TRDIAG

MODEL_COM                           ! model variables 
MODELE                              ! Main and model overhead
ALLOC_DRV                           ! domain decomposition, allocate global distributed arrays
ATMDYN_COM                          ! atmospheric dynamics
ATM_UTILS                           ! utilities for some atmospheric quantities
QUS_COM QUSDEF                      ! T/Q moments, 1D QUS
CLOUDS2 CLOUDS2_DRV CLOUDS_COM      ! clouds modules
SURFACE FLUXES                      ! surface calculation and fluxes
GHY_COM GHY_DRV    ! + giss_LSM     ! land surface and soils + snow model
VEG_DRV                             ! vegetation
! VEG_COM VEGETATION                ! old vegetation
ENT_DRV  ENT_COM   ! + Ent          ! new vegetation
PBL_COM PBL_DRV PBL                 ! atmospheric pbl
ATURB                               ! turbulence in whole atmosphere
LAKES_COM LAKES                     ! lake modules
SEAICE SEAICE_DRV                   ! seaice modules
LANDICE LANDICE_DRV                 ! land ice modules
ICEDYN_DRV ICEDYN                   ! ice dynamics modules
RAD_COM RAD_DRV RADIATION           ! radiation modules
  RAD_native_O3_AR5_V2
RAD_UTILS ALBEDO READ_AERO          ! radiation, albedo, prescribed aerosols
DIAG_COM DIAG DEFACC                ! diagnostics
CLD_AEROSOLS_Menon_MBLK_MAT_E29q BLK_DRV ! aerosol-cloud interactions
lightning                                ! Colin Price lightning model

OCEAN OCNML                         ! ocean modules

Components:
shared ESMF_Interface solvers giss_LSM 
dd2d
Ent

Component Options:
OPTS_Ent = ONLINE=YES PS_MODEL=FBB PFT_MODEL=ENT    
OPTS_dd2d = NC_IO=PNETCDF

Data input files:
    ! start from the restart file of an earlier run ...                 ISTART=8
! AIC=1....rsfE... ! initial conditions, no GIC needed, use
!! AIC=1JAN1961.rsfE4F40.MXL65m  ! end of run with KOCEAN=0

    ! start from observed conditions AIC(,OIC), model ground data GIC   ISTART=2
AIC=AIC.RES_F79.D771201_Z2HX2fromZ1QX1N_B2.00   ! observed init cond (atm. only)
GIC=GIC.144X90.DEC01.1.ext.nc   ! initial ground conditions
OSST=OST_144x90.1876-1885avg.HadISST1.1         ! prescr. climatological ocean
SICE=SICE_144x90.1876-1885avg.HadISST1.1        ! prescr. climatological sea ice
!! q-flux ocean: use the next line instead,       set KOCEAN=1
!! OHT=OTSPEC.E4F40.MXL65m.1956-1960            ! ocean horizontal heat transports
!! OCNML is not used if KOCEAN=0, but needed in and to prepare for q-flux model
OCNML=Z1O.B144x90                               ! mixed layer depth
TOPO=Z2HX2fromZ1QX1N
VEG_DENSE=gsin/veg_dense_2x2.5 ! vegetation density for flammability calculations
RVR=RD_Fb.RVR.bin              ! river direction file

CDN=CD144X90.ext
VEG=V144x90_EntMM16_lc_max_trimmed_scaled_nocrops_4.ij
CROPS=CROPS_and_pastures_Pongratz_to_Hurtt_144X90N_nocasp
SOIL=S144X900098M.ext
TOP_INDEX=top_index_144x90_a.ij.ext
Z4var=Z4var144x89  
! probably need these (should convert to 144x90)
soil_textures=soil_textures_top30cm_2x2.5
SOILCARB_global=soilcarb_top30cm_nmaps_2x2.5bin.dat
GLMELT=GLMELT_144X90_gas.OCN
    ! resolution independent files
RADN1=sgpgxg.table8               ! rad.tables and history files
RADN2=LWTables33k.1a              ! rad.tables and history files
RADN4=LWCorrTables33k              ! rad.tables and history files
RADN5=H2Ocont_MT_CKD  ! Mlawer/Tobin_Clough/Kneizys/Davies H2O continuum table
! other available H2O continuum tables:
!    RADN5=H2Ocont_Ma_2004
!    RADN5=H2Ocont_Roberts
!    RADN5=H2Ocont_MT_CKD  ! Mlawer/Tobin_Clough/Kneizys/Davies
RADN3=miescatpar.abcdv2
! updated aerosols need MADAER=3
RH_QG_Mie=oct2003.relhum.nr.Q633G633.table
RADN6=dust_mass_CakmurMillerJGR06_72x46x20x7x12
RADN7=STRATAER.VOL.1850-2012.May13_hdr
RADN8=cloud.epsilon4.72x46
RADN9=solar.lean2015.ann1610-2014_hdr ! need KSOLAR=2
RADNE=topcld.trscat8
ISCCP=ISCCP.tautables
! ozone files (minimum 1, maximum 9 files + 1 trend file)
GHG=GHG_RCP45.txt 
dH2O=dH2O_by_CH4_monthly
BC_dep=BC.Dry+Wet.depositions.ann
! updated aerosols need MADAER=3
TAero_SUL=SUL_Koch2008_kg_m2_72x46x20_1890-2000h
TAero_SSA=SSA_Koch2008_kg_m2_72x46x20h
TAero_NIT=NIT_Bauer2008_kg_m2_72x46x20_1890-2000h
TAero_OCA=OCA_Koch2008_kg_m2_72x46x20_1890-2000h
TAero_BCA=BCA_Koch2008_kg_m2_72x46x20_1890-2000h
TAero_BCB=BCB_Koch2008_kg_m2_72x46x20_1890-2000h
! ozone files (minimum 1, maximum 17 files )
! this set for #defined RAD_O3_GCM_HRES
O3file_01=Ox/jan2012_o3_shindell_144x90x49x12_AVG1850
O3file_02=Ox/jan2012_o3_shindell_144x90x49x12_AVG1860
O3file_03=Ox/jan2012_o3_shindell_144x90x49x12_AVG1870
O3file_04=Ox/jan2012_o3_shindell_144x90x49x12_AVG1880
O3file_05=Ox/jan2012_o3_shindell_144x90x49x12_AVG1890
O3file_06=Ox/jan2012_o3_shindell_144x90x49x12_AVG1900
O3file_07=Ox/jan2012_o3_shindell_144x90x49x12_AVG1910
O3file_08=Ox/jan2012_o3_shindell_144x90x49x12_AVG1920
O3file_09=Ox/jan2012_o3_shindell_144x90x49x12_AVG1930
O3file_10=Ox/jan2012_o3_shindell_144x90x49x12_AVG1940
O3file_11=Ox/jan2012_o3_shindell_144x90x49x12_AVG1950
O3file_12=Ox/jan2012_o3_shindell_144x90x49x12_AVG1960
O3file_13=Ox/jan2012_o3_shindell_144x90x49x12_AVG1970
O3file_14=Ox/jan2012_o3_shindell_144x90x49x12_AVG1980
O3file_15=Ox/jan2012_o3_shindell_144x90x49x12_AVG1990
O3file_16=Ox/jan2012_o3_shindell_144x90x49x12_AVG2000
O3file_17=Ox/jan2012_o3_shindell_144x90x49x12_AVG2010
Ox_ref=jan2010_o3_shindell_144x90x49x12_April1850 ! for radiative forcing reference

!-----------------------------------------------
!  resolution-independent chemistry input files:
!-----------------------------------------------
MOLEC=chem_files/ds4_moleculesE_terp_soa
JPLRX=chem_files/JPL2011_JUN13_fastterp_y2
JPLPH=chem_files/ds4_photlist_T25
RATJ=chem_files/ratj.giss_25
SPECFJ=chem_files/jv_spec_AV_X68d.dat ! for define AR5_FASTJ_XSECS case, use jv_spec_AV.dat
ATMFJ=chem_files/jv_atms.dat
! fltran file used if rad_FL.ne.0:
FLTRAN=solar.lean2015.ann1610-2014_hdr_fastj2 ! KSOLAR=2
!-----------------------------------------------
MOLEC=chem_files/ds4_moleculesE_terp!_soa
!-----------------------------------------------
!  3D chemistry input files:
!-----------------------------------------------
N2O_IC=gsin/N2O_IC_M23_4x5_6.17_conc_2x2.5_conc
CFC_IC=gsin/CFC_IC_M23_4x5_6.17_conc_2x2.5_conc
CH4_IC=gsin/CH4_IC_M23_4x5_6.17_conc_2x2.5_conc
Ox_IC=gsin/Ox_init_cond_M23_4x5_conc_2x2.5_conc
CO_IC=gsin/CO_init_cond_M23_conc_2x2.5_conc
SULFATE_SA=temp_2x2.5/sulfate_pi_fakeM23_M_SA_2x2.5gf ! really 4x5 and 9-layer
DMS_FIELD=temp_2x2.5/dms_conc_2x2.5gf ! really 4x5
SO2_FIELD=temp_2x2.5/so2_conc_2x2.5gf ! really 4x5

! files for dust tracers
ERS=ERS1_1993_MONTHLY.144x90.threshold-13 ! ERS data
DSRC=Ginoux_source_v2009_VegMask_144x90   ! preferred dust sources
! alternative preferred dust source files:
! DSRC=Ginoux2001_source_VegMask_144x90
! DSRC=Ginoux_source_v2009_NoVegMask_144x90
! DSRC=GriniZender_DustSources_144x90
! DSRC=Tegen_DustSources_144x90
LKTAB=log_dust_emission_60ms-1 ! look up table for emission calculations
LKTAB1=table_wspdf             ! look up table for wind speed probabilities
! AEROCOM year 2000 dust emissions (uncomment for AEROCOM simulations)
! also set imDust=1 further down
!dust_bin1=DUST_bin1_2000_2x2.5.nc
!dust_bin2=DUST_bin2_2000_2x2.5.nc
!dust_bin3=DUST_bin3_2000_2x2.5.nc
!dust_bin4=DUST_bin4_2000_2x2.5.nc
!------- Needed for dry deposition ---------
VEGTYPE=chem_files/vegtype.global_2x2.5gf ! really 4x5
OLSON=chem_files/drydep.table
DRYCOEFF=chem_files/drydep.coef
LAI01=chem_files/lai01.global_2x2.5gf ! really 4x5
LAI02=chem_files/lai02.global_2x2.5gf ! really 4x5
LAI03=chem_files/lai03.global_2x2.5gf ! really 4x5
LAI04=chem_files/lai04.global_2x2.5gf ! really 4x5
LAI05=chem_files/lai05.global_2x2.5gf ! really 4x5
LAI06=chem_files/lai06.global_2x2.5gf ! really 4x5
LAI07=chem_files/lai07.global_2x2.5gf ! really 4x5
LAI08=chem_files/lai08.global_2x2.5gf ! really 4x5
LAI09=chem_files/lai09.global_2x2.5gf ! really 4x5
LAI10=chem_files/lai10.global_2x2.5gf ! really 4x5
LAI11=chem_files/lai11.global_2x2.5gf ! really 4x5
LAI12=chem_files/lai12.global_2x2.5gf ! really 4x5

!---------- mostly transient, mostly AR5 gas tracer emissions ------------------
CO_01=AR5_emis/F/T/CO_ind_AR5_1850-2000_2x2.5_h
CO_02=AR5_emis/F/T/CO_tra_AR5_1850-2000_2x2.5_h
CO_03=AR5_emis/F/T/CO_wst_AR5_1850-2000_2x2.5_h
CO_04=AR5_emis/F/T/CO_awb_AR5_1850-2000_2x2.5_h
CO_05=AR5_emis/F/T/CO_dom_AR5_1850-2000_2x2.5_h
CO_06=AR5_emis/F/T/m_CO_shp_AR5_1850-2000_2x2.5_h
CO_07=AR5_emis/F/T/CO_slv_AR5_1990-2000_2x2.5_h
CO_08=AR5_emis/F/T/CO_ene_AR5_1850-2000_2x2.5_h
CO_09=AR5_emis/F/T/CO_agr_AR5_1990-2000_2x2.5_h
CO_10=AR5_emis/F/T/CO_forestfire_AR5_1900-2000_2x2.5_h
CO_11=AR5_emis/F/T/CO_grassfire_AR5_1900-2000_2x2.5_h
NOx_AIRC=AR5_emis/F/T/NOx_air_AR5_1910-2000_2x2.5
NOx_01=AR5_emis/F/NAT/NOx_Soil_GEIA_2x2.5_HALF_h ! half because we have ag source
NOx_02=AR5_emis/F/T/NOx_awb_AR5_1850-2000_2x2.5_h
NOx_03=AR5_emis/F/T/NOx_dom_AR5_1850-2000_2x2.5_h
NOx_04=AR5_emis/F/T/NOx_ene_AR5_1850-2000_2x2.5_h
NOx_05=AR5_emis/F/T/NOx_ind_AR5_1850-2000_2x2.5_h
NOx_06=AR5_emis/F/T/m_NOx_shp_AR5_1850-2000_2x2.5_h
NOx_07=AR5_emis/F/T/NOx_tra_AR5_1850-2000_2x2.5_h
NOx_08=AR5_emis/F/T/NOx_wst_AR5_1850-2000_2x2.5_h
NOx_09=AR5_emis/F/T/NOx_agr_AR5_1850-2000_2x2.5_h
NOx_10=AR5_emis/F/T/NOx_forestfire_AR5_1900-2000_2x2.5_h
NOx_11=AR5_emis/F/T/NOx_grassfire_AR5_1900-2000_2x2.5_h
! Note that the Isoprene emis file is ignored when BIOGENIC_EMISSIONS
! directive is on. But I am commenting it anyway.
! if BIOGENIC_EMISSIONS or PS_BVOC are defined, and only one
! Terpenes file is available, the Isoprene file is needed too
! Isoprene_01=ORCHIDEE_Isoprene_1990_2x2.5_h
Terpenes_01=ORCHIDEE_Terpenes_1990_2x2.5_h
Terpenes_02=ORCHIDEE_ORVOC_1990_2x2.5_h
! ========= please remember that Alkenes =================
! ========= and Paraffin emissions files =================
! ========= must now be in Kmole units,  =================
! ========= not Kg units ...             =================
Alkenes_01=AR5_emis/F/NAT/Alkenes_vegetation_GEIA_2x2.5_h_1
Alkenes_02=AR5_emis/F/T/m_Alkenes_shp_AR5_1850-2000_2x2.5_h
Alkenes_03=AR5_emis/F/T/Alkenes_wst_AR5_1850-2000_2x2.5_h
Alkenes_04=AR5_emis/F/T/Alkenes_dom_AR5_1850-2000_2x2.5_h
Alkenes_05=AR5_emis/F/T/Alkenes_ind_AR5_1850-2000_2x2.5_h
Alkenes_06=AR5_emis/F/T/Alkenes_tra_AR5_1850-2000_2x2.5_h
Alkenes_07=AR5_emis/F/T/Alkenes_ene_AR5_1850-2000_2x2.5_h
Alkenes_08=AR5_emis/F/T/Alkenes_awb_AR5_1890-2000_2x2.5_h
Alkenes_09=AR5_emis/F/T/Alkenes_agr_AR5_1990-2000_2x2.5_h
Alkenes_10=AR5_emis/F/T/Alkenes_forestfire_AR5_1900-2000_2x2.5_h
Alkenes_11=AR5_emis/F/T/Alkenes_grassfire_AR5_1900-2000_2x2.5_h
Paraffin_01=AR5_emis/F/NAT/Paraffin_vegetation_GEIA_2x2.5_h_1
Paraffin_02=AR5_emis/F/T/m_Paraffin_shp_AR5_1850-2000_2x2.5_h
Paraffin_03=AR5_emis/F/T/Paraffin_wst_AR5_1850-2000_2x2.5_h
Paraffin_04=AR5_emis/F/T/Paraffin_dom_AR5_1850-2000_2x2.5_h
Paraffin_05=AR5_emis/F/T/Paraffin_ind_AR5_1850-2000_2x2.5_h
Paraffin_06=AR5_emis/F/T/Paraffin_tra_AR5_1850-2000_2x2.5_h
Paraffin_07=AR5_emis/F/T/Paraffin_slv_AR5_1860-2000_2x2.5_h
Paraffin_08=AR5_emis/F/T/Paraffin_ene_AR5_1850-2000_2x2.5_h
Paraffin_09=AR5_emis/F/T/Paraffin_awb_AR5_1890-2000_2x2.5_h
Paraffin_10=AR5_emis/F/T/Paraffin_agr_AR5_1990-2000_2x2.5_h
Paraffin_11=AR5_emis/F/T/Paraffin_forestfire_AR5_1900-2000_2x2.5_h
Paraffin_12=AR5_emis/F/T/Paraffin_grassfire_AR5_1900-2000_2x2.5_h
! OFF CH4_01=AR5_emis/F/NAT/CH4gsfMGOLjal_blank_2x2.5_h
! OFF CH4_02=AR5_emis/F/T/CH4_agr_AR5_1850-2000_2x2.5_h
! OFF CH4_03=AR5_emis/F/T/CH4_awb_AR5_1850-2000_2x2.5_h
! OFF CH4_04=AR5_emis/F/T/CH4_dom_AR5_1850-2000_2x2.5_h  
! OFF CH4_05=AR5_emis/F/T/CH4_ene_AR5_1850-2000_2x2.5_h
! OFF CH4_06=AR5_emis/F/T/CH4_ind_AR5_1850-2000_2x2.5_h
! OFF CH4_07=AR5_emis/F/T/CH4_wst_AR5_1850-2000_2x2.5_h
! OFF CH4_08=AR5_emis/F/T/m_CH4_shp_AR5_1850-2000_2x2.5_h
! OFF CH4_09=AR5_emis/F/T/CH4_tra_AR5_1850-2000_2x2.5_h
! OFF CH4_10=AR5_emis/F/NAT/CH4SOILABS_2x2.5_h
! OFF CH4_11=AR5_emis/F/NAT/CH4TRMITE_2x2.5_h
! OFF CH4_12=AR5_emis/F/NAT/CH4WETL+TUNDRA_2x2.5_h
! OFF CH4_13=AR5_emis/F/T/CH4_forestfire_AR5_1900-2000_2x2.5_h
! OFF CH4_14=AR5_emis/F/T/CH4_grassfire_AR5_1900-2000_2x2.5_h
codirect_01=AR5_emis/F/HTAP_codirect_emissions_2x2.5_h
!------------ end of chem emissions files ---------------

!-------Aerosol inputs--------------
!----oxidants needed if not coupled to chem ----
!
AER_CHEM=OXID_E__1TgfF40_2x2.5
OFFLINE_HNO3.nc=HNO3_dummy_2000_GISS2x2.nc
!OFFLINE_HNO3.nc=HNO3_dummy_1850_GISS2x2.nc
!-----------------------------------------------
!       AMP Radiation Input
!
AMP_MIE_TABLES=AMP_MIE_Q_G_A_S.nc
AMP_CORESHELL_TABLES=AMP_CORESHELL_TABLES.nc
!------nitrate inputs--------------
NH3_01=NH3hCON_OCEANflux_Jan10_2x2.5_h
NH3_02=NH3_agr_AR5_1850-2000_2x2.5_h
NH3_03=NH3_awb_AR5_1850-2000_2x2.5_h
NH3_04=NH3_dom_AR5_1850-2000_2x2.5_h
NH3_05=NH3_ind_AR5_1850-2000_2x2.5_h
NH3_06=NH3_ene_AR5_1850-2000_2x2.5_h
NH3_07=NH3_tra_AR5_1850-2000_2x2.5_h
NH3_08=NH3_forestfire_AR5_1900-2000_2x2.5_h
NH3_09=NH3_grassfire_AR5_1900-2000_2x2.5_h
! ------- aerosol -----------
PSREF=ANN1960.E70F40pi.prsurf  ! time avg. surf. pres. on model grid
SO2_VOLCANO=SO2_volc_2000_2x2.5_pres.nc
! -------- aerosol --------------
DMS_SEA=DMS_Kettle_Andreae_2x2.5 ! needed when OFFLINE_DMS_SS=0
!DMS_FLUX=DMS.2x2.5_AEROCOM.nc ! needed when OFFLINE_DMS_SS=1
SALT1=SALT_bin1_2000_2x2.5.nc
SALT2=SALT_bin1_4_2000_2x2.5.nc  ! valis for coarse size range 1 - 4 um (r80=1.7)
!SALT2=SALT_bin2_2000_2x2.5.nc ! valid for coarse size range 1 - 10 ym (r80=5)
dust_bin1=DUST_bin1_2000_2x2.5.nc
dust_bin2=DUST_bin2_2000_2x2.5.nc
dust_bin3=DUST_bin3_2000_2x2.5.nc
dust_bin4=DUST_bin4_2000_2x2.5.nc
!--------AR5 inputs --------------------
M_BC1_BC_AIRC=AR5_emis/F/T/BCIA_air_AR5_1910-2000_2x2.5
M_BC1_BC_01=BCII_awb_AR5_1850-2000_2x2.5_h
M_BC1_BC_02=BCII_dom_AR5_1850-2000_2x2.5_h
M_BC1_BC_03=BCII_ene_AR5_1850-2000_2x2.5_h
M_BC1_BC_04=BCII_ind_AR5_1850-2000_2x2.5_h
M_BC1_BC_05=BCII_tra_AR5_1850-2000_2x2.5_h
M_BC1_BC_06=BCII_wst_AR5_1850-2000_2x2.5_h
M_BC1_BC_07=BCII_shp_AR5_1850-2000_2x2.5_h
M_OCC_OC_01=OCII_awb_AR5_1850-2000_2x2.5_h
M_OCC_OC_02=OCII_dom_AR5_1850-2000_2x2.5_h
M_OCC_OC_03=OCII_ene_AR5_1850-2000_2x2.5_h
M_OCC_OC_04=OCII_ind_AR5_1850-2000_2x2.5_h
M_OCC_OC_05=OCII_tra_AR5_1850-2000_2x2.5_h
M_OCC_OC_06=OCII_wst_AR5_1850-2000_2x2.5_h
M_OCC_OC_07=OCII_shp_AR5_1850-2000_2x2.5_h
!M_BOC_BC_01=BCB_forestfire_AR5_1900-2000_2x2.5_h
!M_BOC_BC_02=BCB_grassfire_AR5_1900-2000_2x2.5_h
!M_BOC_OC_01=OCB_forestfire_AR5_1900-2000_2x2.5_h
!M_BOC_OC_02=OCB_grassfire_AR5_1900-2000_2x2.5_h

M_BC1_BC_08=BCB_forestfire_AR5_1900-2000_2x2.5_h
M_BC1_BC_09=BCB_grassfire_AR5_1900-2000_2x2.5_h
M_OCC_OC_08=OCB_forestfire_AR5_1900-2000_2x2.5_h
M_OCC_OC_09=OCB_grassfire_AR5_1900-2000_2x2.5_h
 
SO2_01=SO2_shp_AR5_1850-2000_2x2.5_h
SO2_02=SO2_awb_AR5_1850-2000_2x2.5_h
SO2_03=SO2_dom_AR5_1850-2000_2x2.5_h
SO2_04=SO2_ind_AR5_1850-2000_2x2.5_h
SO2_05=SO2_tra_AR5_1850-2000_2x2.5_h
SO2_06=SO2_wst_AR5_1850-2000_2x2.5_h
SO2_07=SO2_ene_AR5_1850-2000_2x2.5_h
SO2_08=SO2_forestfire_AR5_1900-2000_2x2.5_h
SO2_09=SO2_grassfire_AR5_1900-2000_2x2.5_h

MSU_wts=MSU.RSS.weights.data      ! MSU-diag
REG=REG2X2.5                      ! special regions-diag

Label and Namelist:  (next 2 lines)
E_AR5_V2_CAMP_F79 (E_AR5_V2_CADI + MATRIX aerosol microphysics)


&&PARAMETERS
! parameters set for choice of ocean model:
KOCEAN=0        ! ocean is prescribed
!! KOCEAN=1        ! ocean is computed
Kvflxo=0        ! usually set to 1 only during a prescr.ocn run by editing "I"
!  Kvflxo=1     ! saves VFLXO files to prepare for q-flux runs (mkOTSPEC)
ocn_cycl=1      ! =0 if ocean varies from year to year

variable_lk=1   ! variable lakes

USE_UNR_DRAG=1      !if 1 => SDRAG is turned off and unresolved drag is applied.
                    !if 0 => SDRAG is intact and alternative gwd is not employed.

! cond_scheme=2   ! newer conductance scheme (N. Kiang) ! not used with Ent

FS8OPX=1.,1.,1.,1.,1.5,1.5,1.,1.
FT8OPX=1.,1.,1.,1.,1.,1.,1.,1.

! Increasing U00a decreases the high cloud cover; increasing U00b decreases net rad at TOA
U00a=0.61   ! above 850mb w/o MC region;  tune this first to get 30-35% high clouds
U00b=1.00   ! below 850mb and MC regions; tune this last  to get rad.balance
WMUI_multiplier = 1.

PTLISO=15.       ! press(mb) above which rad. assumes isothermal layers
H2ObyCH4=0.      ! activates strat.H2O generated by CH4
KSOLAR=2         ! 2: use long annual mean file ; 1: use short monthly file

initial_GHG_setup = 1 ! Set to 0 after initial setup.

! parameters that control the atmospheric/boundary conditions
! if set to 0, the current (day/) year is used: transient run
master_yr=1850
!crops_yr=1850  ! if -1, crops in VEG-file is used
!s0_yr=1850
!s0_day=182
!ghg_yr=1850
!ghg_day=182
volc_yr=-1
!volc_day=182
!aero_yr=1850
od_cdncx=0.        ! don't include 1st indirect effect
cc_cdncx=0.        ! don't include 2nd indirect effect (used 0.0036)
!albsn_yr=1850
dalbsnX=.024
!o3_yr=-1850
!aer_int_yr=1850    !select desired aerosol emissions year or 0 to use JYEAR
! atmCO2=368.6          !uatm for year 2000 - enable for CO2 tracer runs

!variable_orb_par=0
!orb_par_year_bp=100  !  BP i.e. 1950-orb_par_year_bp AD = 1850 AD
!!!!!!!!!!!!!!!!!!!!!!!
! Please note that making o3_yr non-zero tells the model
! to override the transient chemistry tracer emissions'
! use of model year and use abs(o3_yr) instead!
!!!!!!!!!!!!!!!!!!!!!!!
madaer=3         ! 3: updated aerosols          ; 1: default sulfates/aerosols
!--------- general aerosol parameters-------------
rad_forc_lev=1     ! 0: for TOA, 1: for tropopause for rad forcing diags.
rad_interact_aer=1 ! 1: couples aerosols to radiation, 0: use climatology
prather_limits=1   ! 1: to avoid some negative tracers in sub-gridscale
diag_rad=1         ! 1: additional radiation diagnostics
diag_wetdep=1      ! 1: additional wet deposition diagnostics
OCB_om2oc=1.4      ! biomass burning organic matter to organic carbon ratio (default is 1.4)
!BBinc=1.4          ! enhancement factor for carbonaceous aerosols (1.4 for AR5 emissions, 1.0 elsewhere)
!--- number of biomass burning sources (per tracer)
NH3_nBBsources=2
SO2_nBBsources=2
M_BC1_BC_nBBsources=2
M_OCC_OC_nBBsources=2
OFFLINE_DMS_SS=0   ! 1= read offline DMS and dust from aerocom file
!------------------  AMP parameters
AMP_DIAG_FC=  1    ! 2=nmode radiation calls  ||  1=one radiation call
AMP_RAD_KEY = 1    ! 1=Volume Mixing || 2=Core-Shell || 3=Maxwell Garnett
!--------- dust aerosol parameters----------------
imDust=0                          ! 0: PDF emission scheme, 1: AEROCOM
adiurn_dust=0                     ! 1: daily dust diagnostics for selected grid boxes
fracClayPDFscheme=0.04D0          ! clay emission parameter from calibration
fracSiltPDFscheme=0.08D0          ! silt emission parameter from calibration
!to_conc_soildust=1   ! three-dimensional dust output as concentration [kg/m^3]
!-------------------------------------------------
!-----------------------------------------------
!  Start tracer code parameters:
!-----------------------------------------------
!--- define emission sectors above files belong to ---
! example: CH4_13_sect='WET'

!      (careful; they're allowed to overlap):
!       ---------define-REGIONS------------
!        global S.Asia E.Asia Europe N.Amer
REG_S=    -90.,    5.,   15.,   25.,   15.
REG_N=     90.,   35.,   50.,   65.,   55.
REG_W=   -180.,   50.,   95.,  -10., -125.
REG_E=    180.,   95.,  160.,   50.,  -60.
!       ---define-regions-names/order------
REGIONS_ARE='global S_Asia E_Asia Europe N_America'
!-fit-here--|                                                              |---
!       ---define-factors-by-sector--------
!        global S.Asia E.Asia Europe N.Amer
SECT_01= 1.000, 1.000, 1.000, 1.000, 1.000 ! WET (for example)
!       ---define-sectors-names/order------
SECTORS_ARE='WET'
!-fit-here--|                                                              |---
!-----
aircraft_Tyr1=1910 ! regardless of the type of run, if you have non-transient
aircraft_Tyr2=2000 ! aircraft emission files, set these two equal or omit them.

! Colin Price lightning model needs resolution-dependant tuning:
tune_lt_land=2.4920d0 ! =2.2d0*2.17d0, *0.5*1.2*0.5 for 2x2.5, * 2 for bugfix, *0.87
tune_lt_sea= 5.5221d0 ! =3.9d0*2.17d0, *0.25*1.5 for 2x2.5, * 2 for bugfix, *0.87

! -----------------------------------
! Pressure above which Ox, NOx, BrOx, and ClOx will be
! overwritten with climatology based on NINT ozone input.
PltOx=0.0
! NH and SH polar stratospheric cloud formation temperature offsets:
Tpsc_offset_N=-10.d0 
Tpsc_offset_S=-10.d0
! -----------------------------------
! 40-layer model with VMP clouds chemistry tuning parameters: 
! (If you do *not* have VMP clouds comment these out and uncomment the next section):
reg1Power_SpherO2andN2Ocorr=2.0
reg1TopPres_SpherO2andN2Ocorr=50.
reg2Power_SpherO2andN2Ocorr=2.0
reg2TopPres_SpherO2andN2Ocorr=10.
reg3Power_SpherO2andN2Ocorr=1.0
reg3TopPres_SpherO2andN2Ocorr=5.
reg4Power_SpherO2andN2Ocorr=0.5
windowN2Ocorr=0.8
windowO2corr=0.8
! -----------------------------------
! Older, 40-layer model without VMP clouds chemistry tuning parameters:
! reg1Power_SpherO2andN2Ocorr=0.5
! reg1TopPres_SpherO2andN2Ocorr=50.
! reg2Power_SpherO2andN2Ocorr=0.5
! reg2TopPres_SpherO2andN2Ocorr=10.
! reg3Power_SpherO2andN2Ocorr=0.5
! reg3TopPres_SpherO2andN2Ocorr=5.
! reg4Power_SpherO2andN2Ocorr=0.5
! windowN2Ocorr=0.6
! windowO2corr=0.6
! -----------------------------------

COUPLED_CHEM=1     ! to couple chemistry and aerosols
use_sol_Ox_cycle=0 ! (=1) apply ozone changes in radiation, based on solar cycle
clim_interact_chem=1 ! 1=use calculated Ox/CH4 in radiation, 0=use climatology
                   ! If = 0, consider turning on AUXILIARY_OX_RADF CPP directive.
                   ! Note: 0 also turns off chemistry(H2O)-->Q(humidity) feedback
                   ! if you want humidity feedback on but radiation feedback off
                   ! you could do: clim_interact_chem=1, Lmax_rad_{O3,CH4}=0...
! Lmax_rad_O3=0    ! Ox levels used in rad code default is LM
! Lmax_rad_CH4=0   ! CH4 levels used in rad code default is LM
use_rad_n2o=1      ! use the radiation code's N2O
use_rad_cfc=1      ! use rad code cfc11+cfc12, adjusted
use_rad_ch4=1      ! use rad code CH4, shut off sfc sources
rad_FL=1           ! use rad code insolation getting fastj2 photon flux
which_trop=0       ! choose tropopause for chemistry purposes:
                   ! 0=LTROPO(I,J), 1=LS1-1
fix_CH4_chemistry=0    ! for setting fixed methane value for chemistry:
ch4_init_sh=0.791      ! init cond/fixed conditions SH CH4 ppmv
ch4_init_nh=0.791      ! init cond/fixed conditions NH CH4 ppmv
scale_ch4_IC_file=1.d0 ! multiplicative factor on CH4 IC file (fix_CH4_chemistry=-1)

! For altering tracer initial conditions and overwriting by a factor:
! set PI_run=1 and change the corresponding factors below. [For altering
! emissions, use the sectors above in the rundeck.
!PI_run        = 1       ! =1 turns on below factors:
PIratio_N     = 0.667d0 ! {NOx, HNO3, N2O5, HO2NO2}
PIratio_CO_T  = 0.333d0 ! CO in troposphere
PIratio_CO_S  = 0.500d0 ! CO in stratosphere
PIratio_other = 0.500d0 ! {PAN,Isoprene,AlkyNit,Alkenes,Paraffin}
PIratio_N2O   = 1.000d0 ! {N2O ICs, L=1 overwrit}, set to 1 for use_rad_n2o=1
PIratio_CFC   = 1.000d0 ! {CFC ICs, L=1 overwrit}, set to 1 for use_rad_cfc=1
!--- number of biomass burning sources (per tracer)
Alkenes_nBBsources=2
CO_nBBsources=2
NOx_nBBsources=2
Paraffin_nBBsources=2
! OFF FOR NOW: CH4_nBBsources=2

DTsrc=1800.      ! cannot be changed after a run has been started
DT=180.
! parameters that control the Shapiro filter
DT_XUfilter=180. ! Shapiro filter on U in E-W direction; usually same as DT
DT_XVfilter=180. ! Shapiro filter on V in E-W direction; usually same as DT
DT_YVfilter=0.   ! Shapiro filter on V in N-S direction
DT_YUfilter=0.   ! Shapiro filter on U in N-S direction

NIsurf=2         ! (surf.interaction NIsurf times per physics time step)
NRAD=5           ! radiation (every NRAD'th physics time step)
! parameters that affect at most diagn. output:  standard if DTsrc=1800. (sec)
aer_rad_forc=0   ! if set =1, radiation is called numerous times - slow !!
cloud_rad_forc=1 ! calls radiation twice; use =0 to save cpu time
SUBDD=' '        ! no sub-daily frequency diags
NSUBDD=0         ! saving sub-daily diags every NSUBDD-th physics time step (1/2 hr)
KCOPY=2          ! saving acc + rsf
isccp_diags=1    ! use =0 to save cpu time, but you lose some key diagnostics
nda5d=13         ! use =1 to get more accurate energy cons. diag (increases CPU time)
nda5s=13         ! use =1 to get more accurate energy cons. diag (increases CPU time)
ndaa=13
nda5k=13
nda4=48          ! to get daily energy history use nda4=24*3600/DTsrc

Nssw=2           ! until diurnal diags are fixed, Nssw has to be even
Ndisk=960
&&END_PARAMETERS

 &INPUTZ
 YEARI=1949,MONTHI=12,DATEI=1,HOURI=0, ! pick IYEAR1=YEARI (default) or < YEARI
 YEARE=1961,MONTHE=1,DATEE=1,HOURE=0,     KDIAG=12*0,9,
 ISTART=2,IRANDI=0, YEARE=1949,MONTHE=12,DATEE=1,HOURE=1,
 &END
