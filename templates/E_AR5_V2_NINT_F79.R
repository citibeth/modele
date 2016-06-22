E_AR5_V2_NINT_F79.R     GISS Model E                Tiehan      05/13/2015 

!E_AR5_V2_NINT.R GISS Model E  1850 ocn/atm          Larissa     07/15/2014

!E_AR5_V2_NINT is for AR5 Prime production runs; 
!                 U00a=0.63, U00b=1.0, WMUI_multiplier=1.; 
!                 grav.wave adjustment: CMTN=0.1, CDEF=1.6 

!! delete lines starting with '!!' unless E4F40 prepares a q-flux ocean run
!! E4qsF40.R GISS Model E  1850 atm, ocn: q-flux 65m             rar 07/15/2009

!! E4qsF40 = E4F40 with 65m q-flux ocean
E4F40 = modelE as frozen in April 2010:
modelE1 (3.0) 2x2.5 hor. grid with 40 lyrs, top at .1 mb (+ 3 rad.lyrs)
atmospheric composition from year 1850
ocean data: prescribed, 1876-1885 climatology
uses turbulence scheme (no dry conv), grav.wave drag
time steps: dynamics 3.75 min leap frog; physics 30 min.; radiation 2.5 hrs
filters: U,V in E-W and N-S direction (after every physics time step)
         U,V in E-W direction near poles (after every dynamics time step)
         sea level pressure (after every physics time step)

Preprocessor Options
!#define TRACERS_ON                  ! include tracers code
#define USE_ENT
#define NEW_IO
#define RAD_O3_GCM_HRES
#define RAD_O3_DECADAL_INPUT
End Preprocessor Options

Object modules:
     ! resolution-specific source codes
RES_F79                        ! horiz/vert resolution, 2x2.5, top at 0.1mb, 79 layers
DIAG_RES_F                          ! diagnostics
FFT144                              ! Fast Fourier Transform

    ! lat-lon grid specific source codes
GEOM_B                              ! model geometry
DIAG_ZONAL GCDIAGb                  ! grid-dependent code for lat-circle diags
DIAG_PRT POUT                       ! diagn/post-processing output
IO_DRV                              ! old i/o

     ! GISS dynamics with gravity wave drag
ATMDYN MOMEN2ND                     ! atmospheric dynamics
QUS_DRV QUS3D                       ! advection of Q/tracers

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

RVR=RD_Fb.RVR.bin          ! river direction file

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
RADN6=DUST_Tcadi2012_Bauer_kg_m2_144x90x40
RADN7=STRATAER.VOL.1850-2012.May13_hdr
RADN8=cloud.epsilon4.72x46
RADN9=solar.lean2015.ann1610-2014_hdr ! need KSOLAR=2
RADNE=topcld.trscat8
ISCCP=ISCCP.tautables
! ozone files (minimum 1, maximum 9 files + 1 trend file)
GHG=GHG_RCP45.txt
dH2O=dH2O_by_CH4_monthly
BC_dep=BC_dep_Tcadi2012_Bauer_kg_m2_s_144x90_1850-2100
! updated aerosols need MADAER=3
TAero_SUL=SUL_Tcadi2012_Bauer_kg_m2_144x90x40_1850-2100
TAero_SSA=SSA_Tcadi2012_Bauer_kg_m2_144x90x40
TAero_NIT=NIT_Bauer2008_kg_m2_72x46x40_1850-2100h
TAero_OCA=OCA_Tcadi2012_Bauer_kg_m2_144x90x40_1850-2100
TAero_BCA=BCA_Tcadi2012_Bauer_kg_m2_144x90x40_1850-2100
TAero_BCB=BCB_Tcadi2012_Bauer_kg_m2_144x90x40_1850-2100
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

MSU_wts=MSU.RSS.weights.data      ! MSU-diag
REG=REG2X2.5                      ! special regions-diag

Label and Namelist:  (next 2 lines)
E_AR5_V2_NINT_F79 (E_AR5_NINT + U00a=0.63; U00b=1.0; CMTN=0.1; CDEF=1.6) 

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

xCDpbl=1.
cond_scheme=2   ! newer conductance scheme (N. Kiang) ! not used with Ent

FS8OPX=1.,1.,1.,1.,1.5,1.5,1.,1.
FT8OPX=1.,1.,1.,1.,1.,1.,1.,1.

! Increasing U00a decreases the high cloud cover; increasing U00b decreases net rad at TOA
U00a=0.63      ! above 850mb w/o MC region; tune this first to get 30-35% high clouds
U00b=1.00      ! below 850mb and MC regions; then tune this to get rad.balance

WMUI_multiplier = 1.

PTLISO=15.       ! press(mb) above which rad. assumes isothermal layers
H2ObyCH4=1.      ! activates strat.H2O generated by CH4
KSIALB=0         ! 6-band albedo (Hansen) (=1 A.Lacis orig. 6-band alb)
KSOLAR=2         ! 2: use long annual mean file ; 1: use short monthly file

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
madaer=3         ! 3: updated aerosols          ; 1: default sulfates/aerosols

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
Ndisk=480
&&END_PARAMETERS

 &INPUTZ
 YEARI=1949,MONTHI=12,DATEI=1,HOURI=0, ! pick IYEAR1=YEARI (default) or < YEARI
 YEARE=1971,MONTHE=1,DATEE=1,HOURE=0,     KDIAG=12*0,9,
 ISTART=2,IRANDI=0, YEARE=1949,MONTHE=12,DATEE=1,HOURE=1,
!! suggested settings for E4qsF40:
!! YEARI=1901,MONTHI=1,DATEI=1,HOURI=0,
!! YEARE=1931,MONTHE=1,DATEE=1,HOURE=0,     KDIAG=12*0,9,
!! ISTART=8,IRANDI=0, YEARE=1901,MONTHE=1,DATEE=1,HOURE=1,
 &END
