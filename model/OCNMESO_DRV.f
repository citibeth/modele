#include "rundeck_opts.h"

      module ocnmeso_com
      implicit none

!@dbparam use_tdmix whether to use thickness-diffusion mesoscale code
!@dbparam use_gmscz whether to use exponential vertical structure in mesoscale K
      integer :: use_tdmix=0,use_gmscz=0

!@dbparam zsmult multiplier for exp decay scale when use_gmscz>0
!@dbparam kvismult multiplier for Visbeck K when using simple_mesodiff
!@dbparam enhance_shallow_kmeso whether to enhance shallow-ocean diffusivity
      real*8 :: zsmult=1d0
      real*8 :: kvismult=1d0
      integer :: enhance_shallow_kmeso=0

!@dbparam kbg (m2/s) minimum mesoscale diffusivity
      real*8 :: kbg=100d0

!@var g3d potential enthalpy (J/kg)
!@var s3d salinity (kg/kg)
!@var v3d specific volume (ref to mid point pressure)
!@var r3d density (ref to mid point pressure)
!@var p3d mid point pressure
!@var rhox along-layer x-gradient of potential density (ref. to local pres)
!@var rhoy along-layer y-gradient of potential density (ref. to local pres)
!@var rhomz minus the z-gradient of potential density (ref. to local pres)
!@var byrhoz 1/rhomz
!@var dze3d approx distance between layer midpoints
!@var bydze3d 1/dze3d
      real*8, allocatable, dimension(:,:,:) ::
     &     g3d,s3d,p3d,r3d,v3d,
     &     rhox,rhoy,rhomz,byrhoz,
     &     dze3d,bydze3d

      end module ocnmeso_com

      subroutine alloc_ocnmeso_com
      use ocean, only : im,lmo
      use ocnmeso_com
      use param, only : sync_param
      use oceanr_dim, only : ogrid
      use domain_decomp_1d, only : get
      implicit none
      integer :: j_0h,j_1h

      call get(ogrid, j_strt_halo=j_0h, j_stop_halo=j_1h)

      allocate( g3d (lmo,im,j_0h:j_1h) )
      allocate( s3d (lmo,im,j_0h:j_1h) )
      allocate( p3d (lmo,im,j_0h:j_1h) )
      allocate( r3d (lmo,im,j_0h:j_1h) )
      allocate( v3d (lmo,im,j_0h:j_1h) )
      allocate( rhox(lmo,im,j_0h:j_1h) )
      allocate( rhoy(lmo,im,j_0h:j_1h) )

      allocate( rhomz  (im,j_0h:j_1h,lmo) )
      allocate( byrhoz (im,j_0h:j_1h,lmo) )
      allocate( dze3d  (im,j_0h:j_1h,lmo) )
      allocate( bydze3d(im,j_0h:j_1h,lmo) )

      call sync_param('ocean_use_tdmix',use_tdmix)
      if(use_tdmix==1) then
        call alloc_tdmix
      endif
      call sync_param('ocean_use_gmscz',use_gmscz)

      if(use_gmscz<0 .or. use_gmscz>2) then
        call stop_model('bad value of ocean_use_gmscz',255)
      endif

      if(use_gmscz>0) then
        call sync_param('ocean_zsmult',zsmult)
      endif

#ifdef SIMPLE_MESODIFF
      call sync_param('ocean_kvismult',kvismult)
      call sync_param('ocean_enhance_shallow_kmeso',
     &     enhance_shallow_kmeso)
#endif

      call sync_param('ocean_kmeso_bg',kbg)

      end subroutine alloc_ocnmeso_com

      subroutine ocnmeso_drv
      use ocean, only : im,jm,lmo,lmm,mo,g0m,s0m,dts,dxypo
      use ocean, only : use_qus,
     *     gxmo,gymo,gzmo, gxxmo,gyymo,gzzmo, gxymo,gyzmo,gzxmo,
     *     sxmo,symo,szmo, sxxmo,syymo,szzmo, sxymo,syzmo,szxmo
      use ocean, only :
     &     nbyzm,nbyzu,nbyzv, i1yzm,i2yzm, i1yzu,i2yzu, i1yzv,i2yzv
      use ocean_dyn, only : mmi
      use domain_decomp_1d, only : get, halo_update,
     &     south,north, hassouthpole, hasnorthpole
      use oceanr_dim, only : grid=>ogrid
      use odiag, only : oijl=>oijl_loc,oij=>oij_loc,
     &    ijl_ggmfl,ijl_sgmfl
#ifdef TDMIX_AUX_DIAGS
     &   ,ijl_gsymmf,ijl_ssymmf
#endif
     &   ,ijl_mfub,ijl_mfvb,ijl_mfwb
      use odiag, only : ij_gmsc,ij_gmscz
#ifdef TRACERS_OCEAN
      use ocn_tracer_com, only : t_qlimit,ntm
      use ocean, only : trmo,
     &     txmo,tymo,tzmo,txxmo,tyymo,tzzmo,txymo,tyzmo,tzxmo
      use odiag, only: toijl=>toijl_loc,
     *               toijl_conc,toijl_tflx,toijl_gmfl
#endif
      use ocnmeso_com, only : kbg,use_tdmix,use_gmscz
     &     ,enhance_shallow_kmeso
      implicit none
      integer i,j,l,n
!@var gmscz (m) eddy activity depth scale
      real*8, dimension(im,grid%j_strt_halo:grid%j_stop_halo) ::
     &     k2d,gmscz
!@var k3d[xy] diffusivity at x- and y- edges (m2/s)
!@var mokg rescaled version of mo to have kg units like other tracers
      real*8, dimension(:,:,:), allocatable :: k3dx,k3dy,mokg
!@var fl3d 3D fluxes (kg/s) for diagnostic accumulations (see notes in tdmix)
      real*8, dimension(:,:,:,:), allocatable :: fl3d
#ifdef TDMIX_AUX_DIAGS
!@var fl3ds the symmetric part of fl3d
     &     ,fl3ds
#endif
      integer :: ind1,ind2
#ifdef TRACERS_OCEAN
      integer :: itrac
#endif
      real*8, dimension(:,:,:), allocatable ::
     &     Xxxmo,Xyymo,Xzzmo, Xxymo,Xyzmo,Xzxmo

c**** Extract domain decomposition info

      integer :: j_0, j_1, j_0h,j_1h, j_0s,j_1s
      call get(grid, j_strt = j_0, j_stop = j_1,
     &     j_strt_skp = j_0s, j_stop_skp = j_1s,
     &     j_strt_halo = j_0h, j_stop_halo = j_1h)

C**** Calculate horizontal and vertical density gradients.
      call densgrad

C**** Calculate surface mesoscale diffusivity
#ifdef SIMPLE_MESODIFF
      call simple_mesodiff(k2d)
#else
      call orig_mesodiff(k2d)
#endif

C**** Apply GM + Redi tracer fluxes

      if(use_tdmix==1) then

        allocate(mokg(im,grid%j_strt_halo:grid%j_stop_halo,lmo))
        allocate(k3dx(lmo,im,grid%j_strt_halo:grid%j_stop_halo))
        allocate(k3dy(lmo,im,grid%j_strt_halo:grid%j_stop_halo))
        allocate(fl3d(im,grid%j_strt_halo:grid%j_stop_halo,lmo,3))
#ifdef TDMIX_AUX_DIAGS
        allocate(fl3ds(im,grid%j_strt_halo:grid%j_stop_halo,lmo,3))
#endif

        call halo_update(grid,mo)

        mokg = mmi

!        if(use_kmeso2) then
!          call get_kmeso2(kbg,k2d)
!        endif

        if(use_gmscz>0) then
          call get_gmscz(gmscz)
        else
          gmscz = 100000.
        endif

        if(enhance_shallow_kmeso==1) then
          call shallow_enhance_kmeso(k2d,gmscz)
        endif

        call make_k3d(k2d,gmscz,k3dx,k3dy)

        call tdmix_prep(dts,k3dx,k3dy)

        ! Water mass is transported as a "tracer" with unit concentration.
        ! Bolus velocity diagnostics are inferred from this tracer.
        ! See notes in tdmix_mod regarding post-hoc partitioning of the
        ! vertical remapping flux into resolved and bolus-induced components.
        call tdmix(mokg,.false.,fl3d
#ifdef TDMIX_AUX_DIAGS
     &       ,fl3ds  ! will be zero for water mass
#endif
     &       )
        ind1 = ijl_mfub; ind2 = ind1 + 2
        oijl(:,:,:,ind1:ind2) = oijl(:,:,:,ind1:ind2) + fl3d
        do l=1,lmo
        do j=j_0,j_1
        do n=1,nbyzm(j,l)
        do i=i1yzm(n,j,l),i2yzm(n,j,l)
          mo(i,j,l) = mokg(i,j,l)/dxypo(j)
        enddo
        enddo
        enddo
        enddo

        if(use_qus.ne.1) then
          allocate(Xxxmo(im,j_0h:j_1h,lmo)); Xxxmo = 0.
          allocate(Xyymo(im,j_0h:j_1h,lmo)); Xyymo = 0.
          allocate(Xzzmo(im,j_0h:j_1h,lmo)); Xzzmo = 0.
          allocate(Xxymo(im,j_0h:j_1h,lmo)); Xxymo = 0.
          allocate(Xyzmo(im,j_0h:j_1h,lmo)); Xyzmo = 0.
          allocate(Xzxmo(im,j_0h:j_1h,lmo)); Xzxmo = 0.
        endif

        ! Heat transport
        if(use_qus==1) then
          call relax_qusmoms(mokg,g0m,
     &         gxmo,gymo,gzmo, gxxmo,gyymo,gzzmo, gxymo,gyzmo,gzxmo)
        else
          call relax_qusmoms(mokg,g0m,
     &         gxmo,gymo,gzmo, Xxxmo,Xyymo,Xzzmo, Xxymo,Xyzmo,Xzxmo)
        endif
        call tdmix(g0m,.false.,fl3d
#ifdef TDMIX_AUX_DIAGS
     &       ,fl3ds
#endif
     &       )
        ind1 = ijl_ggmfl; ind2 = ind1 + 2
        oijl(:,:,:,ind1:ind2) = oijl(:,:,:,ind1:ind2) + fl3d
#ifdef TDMIX_AUX_DIAGS
        ind1 = ijl_gsymmf; ind2 = ind1 + 2
        oijl(:,:,:,ind1:ind2) = oijl(:,:,:,ind1:ind2) + fl3ds
#endif

        ! Salt transport
        if(use_qus==1) then
          call relax_qusmoms(mokg,s0m,
     &         sxmo,symo,szmo, sxxmo,syymo,szzmo, sxymo,syzmo,szxmo)
        else
          call relax_qusmoms(mokg,s0m,
     &         sxmo,symo,szmo, Xxxmo,Xyymo,Xzzmo, Xxymo,Xyzmo,Xzxmo)
        endif
        call tdmix(s0m,.true. ,fl3d
#ifdef TDMIX_AUX_DIAGS
     &       ,fl3ds
#endif
     &       )
        ind1 = ijl_sgmfl; ind2 = ind1 + 2
        oijl(:,:,:,ind1:ind2) = oijl(:,:,:,ind1:ind2) + fl3d
#ifdef TDMIX_AUX_DIAGS
        ind1 = ijl_ssymmf; ind2 = ind1 + 2
        oijl(:,:,:,ind1:ind2) = oijl(:,:,:,ind1:ind2) + fl3ds
#endif

#ifdef TRACERS_OCEAN
        ! Tracer transport
        ind1 = toijl_gmfl; ind2 = ind1 + 2
        do n=1,ntm
          if(use_qus==1) then
            call relax_qusmoms(mokg,trmo(1,j_0h,1,n),
     &       txmo (1,j_0h,1,n),tymo (1,j_0h,1,n),tzmo (1,j_0h,1,n),
     &       txxmo(1,j_0h,1,n),tyymo(1,j_0h,1,n),tzzmo(1,j_0h,1,n),
     &       txymo(1,j_0h,1,n),tyzmo(1,j_0h,1,n),tzxmo(1,j_0h,1,n)
     &       )
          else
            call relax_qusmoms(mokg,trmo(1,j_0h,1,n),
     &           txmo(1,j_0h,1,n),tymo(1,j_0h,1,n),tzmo(1,j_0h,1,n),
     &           Xxxmo,Xyymo,Xzzmo, Xxymo,Xyzmo,Xzxmo)
          endif
          call tdmix(trmo(1,j_0h,1,n),t_qlimit(n),fl3d
#ifdef TDMIX_AUX_DIAGS
     &         ,fl3ds
#endif
     &         )
          toijl(:,:,:,ind1:ind2,n) = toijl(:,:,:,ind1:ind2,n) + fl3d
#ifdef TDMIX_AUX_DIAGS
          call stop_model(
     &     'ocnmeso_drv: add tracer acc space for tdmix_aux_diags',255)
#endif
        enddo
#endif

      else
        gmscz = 0.
        call gmkdif(k2d)
        call gmfexp(g0m,gxmo,gymo,gzmo,.false.,oijl(1,j_0h,1,ijl_ggmfl))
        call gmfexp(s0m,sxmo,symo,szmo,.true. ,oijl(1,j_0h,1,ijl_sgmfl))
#ifdef TRACERS_OCEAN
        do n = 1,ntm
          call gmfexp(trmo(1,j_0h,1,n),
     &         txmo(1,j_0h,1,n),tymo(1,j_0h,1,n),tzmo(1,j_0h,1,n),
     &         t_qlimit(n),toijl(1,j_0h,1,toijl_gmfl,n))
        end do
#endif

      endif ! use_tdmix or not

      ! set diagnostics
      do j=j_0,j_1
      do i=1,im
        if(lmm(i,j).eq.0) cycle
        oij(i,j,ij_gmsc) = oij(i,j,ij_gmsc) + k2d(i,j)
        oij(i,j,ij_gmscz) = oij(i,j,ij_gmscz) + gmscz(i,j)
      enddo
      enddo

      end subroutine ocnmeso_drv

      SUBROUTINE DENSGRAD
!@sum  DENSGRAD calculates all horizontal and vertical density gradients
      use ocean_dyn, only : bydh
      use ocnmeso_com, only : g3d,s3d,p3d,rhox,rhoy,rhomz,byrhoz,
     &     rho=>r3d,vbar=>v3d,bydzv=>bydze3d,dzv=>dze3d
      use constant, only: grav
      !use ocean, only : bydyv,bydxp
      use ocean, only : dyvo,dxpo
      use ocean, only : G0M,GZM=>GZMO, S0M,SZM=>SZMO, OPRESS, FOCEAN, MO
      use ocean, only : im,jm,lmo,lmm,lmu,lmv,dxypo
      use ocean_dyn, only  : dh
      use oceanr_dim, only : grid=>ogrid
      use domain_decomp_1d, only : get,halo_update_column
      use ocean, only : nbyzm,i1yzm,i2yzm,lmm
      use ocean, only : nbyzu,i1yzu,i2yzu,lmu
      use ocean, only : nbyzv,i1yzv,i2yzv,lmv
      implicit none

c
      REAL*8  BYRHO,DZVLM1,ARHO,ARHOX,ARHOY,ARHOZ,DH0,R1,R2,P12
      INTEGER IM1,IMAX
      Real*8,External   :: VOLGSP

      integer :: i,j,k,l,n,il,ir,jl,jr

      INTEGER :: J_0, J_1, J_0S, J_1S
      LOGICAL :: HAVE_NORTH_POLE,HAVE_SOUTH_POLE

      Real*8, Parameter :: z12eH=.28867513d0  !  z12eH = 1/SQRT(12)
      Real*8, dimension(lmo) :: gup,gdn,sup,sdn,pm
!@var dVBARdZ specific volume vertical difference (ref to lower point pressure)
      Real*8 :: vup,vdn,vupu,vdnu,bym,dvbardz

      real*8, dimension(0:lmo) :: pe

c**** Extract domain decomposition info
      CALL GET(grid, J_STRT = J_0, J_STOP = J_1,
     &               J_STRT_SKP  = J_0S,   J_STOP_SKP  = J_1S,
     &               HAVE_NORTH_POLE = HAVE_NORTH_POLE,
     &               HAVE_SOUTH_POLE = HAVE_SOUTH_POLE)


      RHOMZ = -0.0; BYRHOZ = -0.0;
      RHO = -0.0; BYDH = -0.0;
      DZV = -0.0; BYDZV = -0.0;

      do j=j_0,j_1
      do n=1,nbyzm(j,1)
      do i=i1yzm(n,j,1),i2yzm(n,j,1)
        pe(0) = opress(i,j)
        do l=1,lmm(i,j)
          bym = 1d0/(mo(i,j,l)*dxypo(j))
          pe(l) = pe(l-1) + mo(i,j,l)*grav
          g3d(l,i,j) = g0m(i,j,l)*bym
          s3d(l,i,j) = s0m(i,j,l)*bym
          p3d(l,i,j) = .5*(pe(l)+pe(l-1))

            GUP(L)=(G0M(I,J,L)-2*z12eH*GZM(I,J,L))*BYM
            GDN(L)=(G0M(I,J,L)+2*z12eH*GZM(I,J,L))*BYM
            SUP(L)=(S0M(I,J,L)-2*z12eH*SZM(I,J,L))*BYM
            SDN(L)=(S0M(I,J,L)+2*z12eH*SZM(I,J,L))*BYM
            sup(l) = max(0d0,sup(l))
            sdn(l) = max(0d0,sdn(l))
C**** Calculate potential specific volume (ref to mid-point pr)
            PM(L) = P3D(L,I,J)
            VUP = VOLGSP (GUP(L),SUP(L),PM(L))
            VDN = VOLGSP (GDN(L),SDN(L),PM(L))
            VBAR(L,I,J) = (VUP + VDN)*.5
C**** Vertical gradient calculated using lower box mid-point pr
            IF (L.gt.1) then 
              VUPU = VOLGSP (GUP(L-1),SUP(L-1),PM(L))
              VDNU = VOLGSP (GDN(L-1),SDN(L-1),PM(L))
              dVBARdZ = .5* (VUP + VDN - VUPU - VDNU)
              DZVLM1     = 0.5* (DH(I,J,L) + DH(I,J,L-1))
              DZV(I,J,L-1) = DZVLM1
              BYDZV(I,J,L-1) = 1d0/DZVLM1
              RHOMZ(I,J,L-1)=MAX(0d0,
     *             -dVBARdZ*BYDZV(I,J,L-1)/VBAR(L-1,I,J)**2)
C**** minus vertical gradient
              IF(RHOMZ(I,J,L-1).ne.0.)
     *             BYRHOZ(I,J,L-1)=1./RHOMZ(I,J,L-1)
            end if
C**** RHO(I,J,L)  Density=1/specific volume
            RHO(L,I,J)  = 1d0/VBAR(L,I,J)
            BYDH(I,J,L) = 1d0/DH(I,J,L)

        enddo
      enddo
      enddo
      enddo

C**** Copy to all longitudes at poles
      If(have_north_pole) Then
        Do L=1,LMM(1,JM)
          RHO(L,2:IM,JM) = RHO(L,1,JM)
          VBAR(L,2:IM,JM) = VBAR(L,1,JM)
          G3D(L,2:IM,JM) = G3D(L,1,JM)
          S3D(L,2:IM,JM) = S3D(L,1,JM)
          P3D(L,2:IM,JM) = P3D(L,1,JM)
          DZV(2:IM,JM,L) = DZV(1,JM,L)
          BYDZV(2:IM,JM,L) = BYDZV(1,JM,L)
          BYDH(2:IM,JM,L) = BYDH(1,JM,L)
          RHOMZ(2:IM,JM,L) = RHOMZ(1,JM,L)
          BYRHOZ(2:IM,JM,L) = BYRHOZ(1,JM,L)
        EndDo
      EndIf
      If(have_south_pole) Then
        Do L=1,LMM(1,1)
          RHO(L,2:IM,1) = RHO(L,1,1)
          VBAR(L,2:IM,1) = VBAR(L,1,1)
          G3D(L,2:IM,1) = G3D(L,1,1)
          S3D(L,2:IM,1) = S3D(L,1,1)
          P3D(L,2:IM,1) = P3D(L,1,1)
          DZV(2:IM,1,L) = DZV(1,1,L)
          BYDZV(2:IM,1,L) = BYDZV(1,1,L)
          BYDH(2:IM,1,L) = BYDH(1,1,L)
          RHOMZ(2:IM,1,L) = RHOMZ(1,1,L)
          BYRHOZ(2:IM,1,L) = BYRHOZ(1,1,L)
        EndDo
      EndIf

      CALL HALO_UPDATE_COLUMN(grid,RHO)
      CALL HALO_UPDATE_COLUMN(grid,G3D)
      CALL HALO_UPDATE_COLUMN(grid,S3D)
      CALL HALO_UPDATE_COLUMN(grid,P3D)

C**** Calculate density gradients
      rhox = 0.
      rhoy = 0.
      do j=j_0s,j_1s
      do n=1,nbyzu(j,1)
      do i=i1yzu(n,j,1),i2yzu(n,j,1)
        il = i
        if(i.eq.im) then
          ir = 1
        else
          ir = i+1
        endif
        do l=1,lmu(i,j)
          p12 = .5d0*(p3d(l,il,j)+p3d(l,ir,j)) 
          rhox(l,i,j) = (
     &         1d0/volgsp(g3d(l,ir,j),s3d(l,ir,j),p12)
     &        -1d0/volgsp(g3d(l,il,j),s3d(l,il,j),p12)
c     &         )*bydxp(j)!/dxpo(j)
     &         )*(1d0/dxpo(j))
        enddo
      enddo
      enddo
      enddo
      do j=max(2,j_0-1),j_1s
      do n=1,nbyzv(j,1)
      do i=i1yzv(n,j,1),i2yzv(n,j,1)
        jl = j
        jr = j+1
        do l=1,lmv(i,j)
          p12 = .5d0*(p3d(l,i,jl)+p3d(l,i,jr)) 
          rhoy(l,i,j) = (
     &         1d0/volgsp(g3d(l,i,jr),s3d(l,i,jr),p12)
     &        -1d0/volgsp(g3d(l,i,jl),s3d(l,i,jl),p12)
c     &         )*bydyv(j)!/dyvo(j)
     &         )*(1d0/dyvo(j))
        enddo
      enddo
      enddo
      enddo

      END SUBROUTINE DENSGRAD

      subroutine simple_mesodiff(ainv)
      use domain_decomp_1d, only : get
      use oceanr_dim, only : grid=>ogrid
      use constant, only : grav
      use ocnmeso_com, only : rho=>r3d,rhomz,rhox,rhoy
      use ocnmeso_com, only : kvismult,kbg
      use ocean, only : im,jm,lmo,lmm,lmu,lmv,sinpo,ze,dzo
      use ocean, only : sinic,cosic
      use ocean, only : nbyzv,i1yzv,i2yzv
      implicit none
! A prescription of mesoscale K a la "Visbeck" that is being used for
! a GISS-MIT MIP (excepting the latitude dependence).
! It is like the default method in that
! K = factor*L*L*[upper-ocean average of 1/t_eady=N*isopycnal_slope]
! but with the following differences:
! 1. The length scale L is a constant rather than the deformation
!    radius R_d=NH/f (or beta-plane R_d=sqrt(NH/beta) at low latitudes)
! 2. Off the equator, the explicit latitude dependence is 1/sin(lat)
!    rather than the 1/sin**2(lat) arising from use of L=R_d.  Near the
!    equator, the 1/sin(lat) dependence is capped, while the use of
!    beta-plane R_d in the default method makes K independent of N.
! 3. As N goes to zero, K goes to zero via imposition of a maximum
!    isopycnal slope rather than via R_d.  As N goes to infinity
!    but 1/t_eady goes to zero, K goes to zero, while the N**2 in
!    R_d in the default method permits large K (off the equator)
! 4. N*isopycnal_slope is computed locally and averaged over the upper
!    ocean rather than computed using averages of density gradients.
!    Where the ocean is shallower than an imposed limit, the
!    upper-ocean average of 1/t_eady is computed as per the
!    comments above for the mindepth parameter.
! 5. Minimum and maximum K values are imposed.

      real*8, dimension(im,grid%j_strt_halo:grid%j_stop_halo) :: ainv

! 
!@var snavg vertical average of the product of isoneutral slope times
!@+         N (sqrt of B-V freq).  snavg = 1/t_eady
      real*8 :: snavg,snsum,delz,ksurf,arhoh,arho,arhox,arhoy,arhoz
      real*8, parameter ::
!@param maxslope maximum slope used in t_eady calculation
     &         maxslope=4d-4
!@param mindepth (m) minimum averaging depth in t_eady calculation.  If
!@+              the local depth is less than mindepth, the calculation
!@+              is performed as if the depth were mindepth and isoneutral
!@+              slopes were zero between the local depth and mindepth.
     &        ,mindepth=400d0
!@param kfac the value of AMU used for SIMPLE_MESODIFF
     &        ,kfac0=12d-3  !*2d0
!@param lscale (m) fixed length scale for calculating mesoscale diffusivity
     &        ,lscale=2.5d5
!@param maxk (m2/s) upper bound for mesoscale diffusivity
     &        ,maxk=6000d0
!@param minsinlat minimum for 1/sin(lat) scaling of mesoscale diffusivity
     &        ,minsinlat=.1d0 ! roughly sin(6 degrees)
!
!@var mink (m2/s) lower bound for mesoscale diffusivity
      real*8 :: mink=100d0

      integer :: i,j,l,n,im1,lav
      integer :: j_0, j_1, j_0s, j_1s
      logical :: have_north_pole
      integer :: l1k
      real*8 :: kfac

      mink = kbg

      kfac = kfac0*kvismult

C**** Calculate level at 1km depth
      do l1k=1,lmo-1
        if(ze(l1k+1) .ge. 1d3) exit
      enddo

      call get(grid, j_strt = j_0, j_stop = j_1,
     &               j_strt_skp  = j_0s,   j_stop_skp  = j_1s,
     &               have_north_pole = have_north_pole)

      do j=j_0s,j_1s
      do i=1,im
        if(lmm(i,j).eq.0) cycle
        im1 = i-1
        if(i.eq.1) im1=im
        arho = 0.
        lav = min(l1k,lmm(i,j))
        snsum = 0.
        do l=1,lav
          arho  = arho  + rho(l,i,j)
          arhoz = .5d0*(rhomz(i,j,max(1,l-1))
     &                + rhomz(i,j,min(l,lmm(i,j)-1)) )
          if(arhoz.le.0.) cycle
          arhox = .5d0*(rhox(l,im1,j  )+rhox(l,i,j))
          arhoy = .5d0*(rhoy(l,i  ,j-1)+rhoy(l,i,j))
          arhoh = sqrt(arhox*arhox+arhoy*arhoy)
          snsum = snsum + dzo(l)*sqrt(arhoz)*min(maxslope, arhoh/arhoz)
        enddo
        arho  = arho / real(lav,kind=8)
        delz = max(mindepth, .5*(ze(lav)+ze(lav-1)-ze(1)))
        snavg = snsum * sqrt(grav/arho) / delz
        ksurf = (kfac * (lscale**2)) *snavg
        ksurf = ksurf/max(abs(sinpo(j)),minsinlat)
        ksurf = max(mink, min(maxk, ksurf))
        ainv(i,j) = ksurf
      enddo
      enddo

      if(have_north_pole) then
        arho = 0.
        lav = min(l1k,lmm(1,jm))
        snsum = 0.
        do l=1,lav
          i = 1
          j = jm
          arho  = arho  + rho(l,i,j)
          arhoz = .5d0*(rhomz(i,j,max(1,l-1))
     &                + rhomz(i,j,min(l,lmm(i,j)-1)) )
          if(arhoz.le.0.) cycle
          arhox = 0.
          arhoy = 0.
          j = jm-1
          do n=1,nbyzv(j,l)
          do i=i1yzv(n,j,l),i2yzv(n,j,l)
            arhox = arhox - sinic(i)*rhoy(l,i,j)
            arhoy = arhoy + cosic(i)*rhoy(l,i,j)
          enddo
          enddo
          arhox = arhox*2/im
          arhoy = arhoy*2/im
          arhoh = sqrt(arhox*arhox+arhoy*arhoy)
          snsum = snsum + dzo(l)*sqrt(arhoz)*min(maxslope, arhoh/arhoz)
        enddo
        i = 1
        j = jm
        arho  = arho / real(lav,kind=8)
        delz = max(mindepth, .5*(ze(lav)+ze(lav-1)-ze(1)))
        snavg = snsum * sqrt(grav/arho) / delz
        ksurf = (kfac * (lscale**2)) *snavg
        ksurf = ksurf/max(abs(sinpo(j)),minsinlat)
        ksurf = max(mink, min(maxk, ksurf))
        ainv(i,j) = ksurf
        ainv(2:im,j) = ainv(1,j)
      endif

      end subroutine simple_mesodiff

      subroutine get_gmscz(gmscz)
!@sum get_gmscz obtains a characteristic depth scale over which
!@+   surface-connected baroclinic eddies are active.  This version
!@+   calculates it as:
!@+        z-integral ( |grad_h(rho)| * z )
!@+        --------------------------------
!@+        z-integral ( |grad_h(rho)|     )
!@+   where grad_h is the horizontal gradient operator.
      use ocean, only : im,jm,lmo
      use ocean, only : dzo,ze
      use ocean, only : nbyzm,i1yzm,i2yzm,lmm
      use ocnmeso_com, only : zsmult
      use ocnmeso_com, only : rhox,rhoy
      use oceanr_dim, only : grid=>ogrid
      use domain_decomp_1d, only : get
      implicit none
!@var gmscz (m) eddy activity depth scale
      real*8, dimension(im,grid%j_strt_halo:grid%j_stop_halo) ::
     &     gmscz
c
      integer :: i,j,l,n,lmaxscz,il,ir,jl,jr
      integer :: j_0s,j_1s
      logical :: have_north_pole
      real*8 :: rxysum,arhox,arhoy,arhohdz


      call get(grid, j_strt_skp=j_0s, j_stop_skp=j_1s,
     &               have_north_pole=have_north_pole)

      do lmaxscz=1,lmo
        if(ze(lmaxscz+1).gt.3000.) exit
      enddo
      do j=j_0s,j_1s
      do n=1,nbyzm(j,1)
      do i=i1yzm(n,j,1),i2yzm(n,j,1)
        if(i.eq.1) then
          il = im
        else
          il = i-1
        endif
        ir = i
        jl = j-1
        jr = j
        gmscz(i,j) = 0.
        rxysum = 0.
        do l=1,min(lmm(i,j),lmaxscz)
          arhox = .5d0*(rhox(l,il,j)+rhox(l,ir,j))
          arhoy = .5d0*(rhoy(l,i,jl)+rhoy(l,i,jr))
          arhohdz = sqrt(arhox*arhox+arhoy*arhoy)*dzo(l)
          gmscz(i,j) = gmscz(i,j) + arhohdz*.5d0*(ze(l-1)+ze(l))
          rxysum = rxysum + arhohdz
        enddo
        gmscz(i,j) = zsmult * gmscz(i,j)/(rxysum+1d-30)
      enddo
      enddo
      enddo

      if(have_north_pole) gmscz(:,jm) = gmscz(:,jm-1)

      end subroutine get_gmscz

      subroutine shallow_enhance_kmeso(k2d,gmscz)
!@sum shallow_enhance_kmeso increase diffusivity in shallow waters
!@+   to help disperse river input and prevent too-low salinities
      use ocean, only : im,ze
      use ocean, only : nbyzm,i1yzm,i2yzm,lmm
      use oceanr_dim, only : grid=>ogrid
      use domain_decomp_1d, only : get
      implicit none
      real*8, dimension(im,grid%j_strt_halo:grid%j_stop_halo) ::
     &     k2d,gmscz
c
      integer :: i,j,l,n

      integer :: j_0,j_1

      call get(grid, j_strt=j_0, j_stop=j_1)

      do j=j_0,j_1
      do n=1,nbyzm(j,1)
      do i=i1yzm(n,j,1),i2yzm(n,j,1)
        if(ze(lmm(i,j)).lt.500.) then ! hard-coded definition of shallow
          k2d(i,j) = max(k2d(i,j),1200.) ! 1200 m2/s
          gmscz(i,j) = 1d4 ! remove vertical dependence
        endif
      enddo
      enddo
      enddo

      end subroutine shallow_enhance_kmeso

      subroutine make_k3d(k2d,gmscz,k3dx,k3dy)
!@sum make_k3d construct 3D diffusivity as the product of a
!@+   column-characteristic diffusivity k2d and a vertical shape
!@+   function having a characteristic depth scale gmscz.
!@+   Both k2d and gmscz vary horizontally.
      use ocean, only : im,jm,lmo
      use ocean, only : nbyzu,i1yzu,i2yzu,lmu
      use ocean, only : nbyzv,i1yzv,i2yzv,lmv
      use oceanr_dim, only : grid=>ogrid
      use domain_decomp_1d, only : get,halo_update
      implicit none
      real*8, dimension(im,grid%j_strt_halo:grid%j_stop_halo) ::
     &     k2d,gmscz
      real*8, dimension(lmo,im,grid%j_strt_halo:grid%j_stop_halo) ::
     &     k3dx,k3dy
c
      integer :: i,j,l,n,ll,lr,il,ir,jl,jr

      integer :: j_0,j_1,j_0s,j_1s

      real*8 :: ksurf,gmscze

      call get(grid, j_strt=j_0, j_stop=j_1,
     &               j_strt_skp=j_0s, j_stop_skp=j_1s)


      call halo_update(grid,k2d)
      call halo_update(grid,gmscz)

      do j=j_0s,j_1s
      do n=1,nbyzu(j,1)
      do i=i1yzu(n,j,1),i2yzu(n,j,1)
        il = i
        if(i.eq.im) then
          ir = 1
        else
          ir = i+1
        endif
        ksurf = .5d0*(k2d(il,j)+k2d(ir,j))
        gmscze = .5d0*(gmscz(il,j)+gmscz(ir,j))
        call do_1d(ksurf,gmscze,k3dx(:,i,j))
      enddo
      enddo
      enddo

      do j=max(2,j_0-1),j_1s
      do n=1,nbyzv(j,1)
      do i=i1yzv(n,j,1),i2yzv(n,j,1)
        ksurf = .5d0*(k2d(i,j)+k2d(i,j+1))
        gmscze = .5d0*(gmscz(i,j)+gmscz(i,j+1))
        call do_1d(ksurf,gmscze,k3dy(:,i,j))
      enddo
      enddo
      enddo

      contains

      subroutine do_1d(ks,zs,k1d)
      use ocnmeso_com, only : kbg,use_gmscz
      use ocean, only : ze
      real*8 :: ks,zs,k1d(lmo)
      integer :: l
      real*8 :: zbyzs,expfac
      do l=1,lmo
        zbyzs = .5d0*(ze(l-1)+ze(l))/zs
        if(use_gmscz==0) then
          expfac = 1d0
        elseif(use_gmscz==1) then
          expfac = exp(-zbyzs)
        else !if(use_gmscz==2) then
          expfac = exp(-(zbyzs-1d0)**2)
        endif
        k1d(l) = (ks-kbg)*expfac+kbg
      enddo

      end subroutine do_1d

      end subroutine make_k3d

      SUBROUTINE ORIG_MESODIFF(AINV)
!@sum  DENSGRAD calculates all horizontal and vertical density gradients
!@auth Gavin Schmidt/Dan Collins
!@ver  1.0
      use constant, only : omega,grav,radius
      use domain_decomp_1d, only : get
      use oceanr_dim, only : grid=>ogrid
      use ocean, only : im,jm,lmo,lmm,lmu,lmv,ze,dzo,sinpo,cospo,dypo
      use ocnmeso_com, only : rho=>r3d,g3d,s3d,p3d,rhox,rhoy,vbar=>v3d
      implicit none
      real*8, dimension(im,grid%j_strt_halo:grid%j_stop_halo) :: ainv
c
      REAL*8  BYRHO,CORI,BETA,ARHO,ARHOX,ARHOY,ARHOZ,AN,RD
     *     ,BYTEADY,DZSUMX,DZSUMY,R1,R2,P12,RHOZ1K,VUP,VDN
      REAL*8 :: HUP
      INTEGER I,J,L,IM1,LAV,LAVM

      INTEGER :: J_0, J_1, J_0S, J_1S
      LOGICAL :: HAVE_NORTH_POLE
      Real*8, External :: VOLGSP

!@var LUP level corresponding to 1km depth
      INTEGER :: LUP

!@var AMU = Visbeck scheme scaling parameter (1)
      REAL*8, PARAMETER :: AMU = 0.13d0

C**** Calculate level at 1km depth
        LUP=0
   10   LUP=LUP + 1
        IF (ZE(LUP+1).lt.1d3) GOTO 10


      HUP = ZE(LUP)

c**** Extract domain decomposition info
      CALL GET(grid, J_STRT = J_0, J_STOP = J_1,
     &               J_STRT_SKP  = J_0S,   J_STOP_SKP  = J_1S,
     &               HAVE_NORTH_POLE = HAVE_NORTH_POLE)

      AINV = 0.

C**** Calculate VMHS diffusion = amu* min(NH/f,equ.rad)^2 /Teady
!$OMP PARALLEL DO  PRIVATE(J,CORI,BETA,IM1,I,ARHO,ARHOX,ARHOY,
!$OMP&  DZSUMX,DZSUMY,LAV,L,ARHOZ,AN,RD,BYTEADY)
      !DO J=2,JM-1
      DO J=J_0S,J_1S
        CORI = ABS(2d0*OMEGA*SINPO(J))
        BETA = ABS(2d0*OMEGA*COSPO(J)/RADIUS)
        IM1=IM
        DO I=1,IM
          IF (LMM(I,J).gt.0) THEN

C**** Calculate average density + gradients over [1,LUP]
            ARHO  = 0.
            ARHOX = 0.
            ARHOY = 0.
            DZSUMX = 0.
            DZSUMY = 0.
            LAV = MIN(LUP,LMM(I,J))
            DO L=1,LAV
              ARHO  = ARHO  + RHO(L,I,J)
              IF(LMU(IM1,J).ge.L) THEN
                ARHOX = ARHOX + RHOX(L,IM1,J)*DZO(L)
                DZSUMX = DZSUMX + DZO(L)
              END IF
              IF(LMU(I  ,J).ge.L) THEN
                ARHOX = ARHOX + RHOX(L,I  ,J)*DZO(L)
                DZSUMX = DZSUMX + DZO(L)
              END IF
              IF(LMV(I,J-1).ge.L) THEN
                ARHOY = ARHOY + RHOY(L,I,J-1)*DZO(L)
                DZSUMY = DZSUMY + DZO(L)
              END IF
              IF(LMV(I,J  ).ge.L) THEN
                ARHOY = ARHOY + RHOY(L,I,J  )*DZO(L)
                DZSUMY = DZSUMY + DZO(L)
              END IF
            END DO
            ARHO  = ARHO / REAL(LAV,KIND=8)
            IF (DZSUMX.gt.0.) ARHOX = ARHOX / DZSUMX
            IF (DZSUMY.gt.0.) ARHOY = ARHOY / DZSUMY
            IF (LAV.gt.1) THEN
              ARHOZ = 2.*(RHO(LAV,I,J)-RHO(1,I,J))/
     *             (ZE(LAV)+ZE(LAV-1)-ZE(1))
            ELSE
              ARHOZ = 0.
            END IF
C**** avoid occasional inversions. IF ARHOZ<=0 then GM is pure vertical
C**** so keep at zero, and let KPP do the work.
            IF (ARHOZ.gt.0) THEN
              AN = SQRT(GRAV * ARHOZ / ARHO)
              RD = AN * HUP / CORI
              IF (RD.gt.ABS(J-.5*(JM+1))*DYPO(J)) RD=SQRT(AN*HUP/BETA)
              BYTEADY = GRAV * SQRT(ARHOX*ARHOX + ARHOY*ARHOY) / (AN
     *             *ARHO)
              AINV(I,J) = AMU * RD**2 * BYTEADY ! was = AIN
            END IF
          END IF
          IM1=I
        END DO
      END DO
!$OMP  END PARALLEL DO

C**** North pole
      if ( HAVE_NORTH_POLE ) then
        IF (LMM(1,JM).gt.0) THEN
          I = 1
          J = JM
C**** Vertical potential density gradient in top 1km
            LAV = MIN(LUP,LMM(I,J))
            LAVM = MAX(LAV/2,1)   ! mid depth
            VUP = VOLGSP (G3D(1,I,J),S3D(1,I,J),P3D(LAVM,I,J))
            VDN = VOLGSP (G3D(LAV,I,J),S3D(LAV,I,J),P3D(LAVM,I,J))
            RHOZ1K = (VUP - VDN)/VBAR(LAVM,I,J)**2

C**** Calculate average density + gradients over [1,LUP]
          ARHO  = 0. ; ARHOY = 0. ;  DZSUMY = 0.
          LAV = MIN(LUP,LMM(1,JM))
          DO L=1,LAV
            ARHO  = ARHO  + RHO(L,1,JM)
            DO I=1,IM
              IF(LMV(I,JM-1).ge.L) THEN
! take abs to get a non-directional scale
                ARHOY = ARHOY + ABS(RHOY(L,I,JM-1))*DZO(L)
                DZSUMY = DZSUMY + DZO(L)
              END IF
            END DO
          END DO
          ARHO  = ARHO / REAL(LAV,KIND=8)
          IF (DZSUMY.gt.0.) ARHOY = ARHOY / DZSUMY
          IF (LAV.gt.1) THEN
            ARHOZ=2*RHOZ1K/(ZE(LAV)+ZE(LAV-1)-ZE(1))
          ELSE
            ARHOZ = 0.
          END IF
C**** avoid occasional inversions. IF ARHOZ<=0 then GM is pure vertical
C**** so keep at zero, and let KPP do the work.
          IF (ARHOZ.gt.0) THEN
            AN = SQRT(GRAV * ARHOZ / ARHO)
            CORI = ABS(2d0*OMEGA*SINPO(JM))
            RD = AN * HUP / CORI
            BYTEADY = GRAV * ARHOY / (AN*ARHO)
            AINV(1,JM) = AMU * RD**2 * BYTEADY ! was = AIN
          END IF
        END IF
        AINV(2:IM,JM)=AINV(1,JM)
      endif

C****
c      USE FILEMANAGER
c      INTEGER iu_ODIFF
c      INTEGER, SAVE :: IFIRST = 1
c      CHARACTER TITLE*80
c      REAL*8, DIMENSION(IM,JM) ::  AINV_glob
c      IF (IFIRST.eq.1) THEN  !output GM diffusion coefficient
c        CALL PACK_DATA(grid,    AINV  ,    AINV_glob)
c        if( AM_I_ROOT() ) then
c          call openunit('ODIFF',iu_ODIFF,.true.,.false.)
c          TITLE = "Visbeck scaling for GM coefficient m^2/s"
c          WRITE(iu_ODIFF) TITLE,((REAL(AINV_glob(I,J),KIND=4),i=1,im),
c     *                            j=1,jm)
c          call closeunit(iu_ODIFF)
c        endif
c        IFIRST = 0
c      END IF

      RETURN
C****
      END SUBROUTINE ORIG_MESODIFF
