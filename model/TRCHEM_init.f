#include "rundeck_opts.h"
      SUBROUTINE cheminit
!@sum cheminit initialize model chemistry
!@auth Drew Shindell (modelEifications by Greg Faluvegi)
!@ver  1.0 (based on cheminit0C5_M23p & ds4p_chem_init_M23.f)
!@calls jplrts,fastj2_init,reactn

C**** GLOBAL parameters and variables:
      USE FILEMANAGER, only: openunit,closeunit
      USE MODEL_COM, only: Itime, ItimeI, IM
      USE DOMAIN_DECOMP_ATM, only : GET,grid
      USE TRACER_COM, only: oh_live,no3_live
      USE TRCHEM_Shindell_COM, only:ny,numfam,nn,nps,nds,
     &    ndnr,kps,kds,kpnr,kdnr,nnr,nr,npnr,nr2,nr3,nmm,nhet,
     &    prnls,prnrts,prnchg,lprn,jprn,iprn,ay,pHOx,pOx,pNOx,
     &    yCH3O2,yC2O3,yROR,yXO2,yAldehyde,yNO3,yRXPAR,yXO2N,acetone,
     &    allowSomeChemReinit,pNO3
     &    ,pCLOx,pCLx,pOClOx,pBrOx,yCl2,yCl2O2

      IMPLICIT NONE

C**** Local parameters and variables and arguments:
!@var iu_data temporary unit number
!@var i,l loop dummy
      integer :: iu_data,i,l,j
      integer :: J_0,J_1,J_0S,J_1S,J_1H,J_0H,I_0,I_1
         
      CALL GET(grid, J_STRT    =J_0,  J_STOP    =J_1,
     &               I_STRT    =I_0,  I_STOP    =I_1,
     &               J_STRT_SKP=J_0S, J_STOP_SKP=J_1S,
     &               J_STRT_HALO=J_0H, J_STOP_HALO=J_1H)

C Read chem diagnostics parameters and molecule names
C from MOLEC file:
      call openunit('MOLEC',iu_data,.false.,.true.)
      read(iu_data,100)prnls,prnrts,prnchg,lprn,jprn,iprn
      read(iu_data,110)ay
      call closeunit(iu_data)

C Read JPL chemical reactions/rates from unit JPLRX:
      call jplrts

! initialize fastj
      call fastj2_init

c Set up arrays of reaction numbers involving each molecule:
      call reactn

C Initialize a few (IM,JM,LM) arrays, first hour only:
      IF(Itime == ItimeI .and. allowSomeChemReinit == 1) THEN
        ! allowSomeChemReinit condition b/c these are in RSF files:
        pHOx(I_0:I_1,J_0:J_1,:)     =1.d0
        pOx(I_0:I_1,J_0:J_1,:)      =1.d0
        pNOx(I_0:I_1,J_0:J_1,:)     =1.d0
        pNO3(I_0:I_1,J_0:J_1,:)     =0.d0
        yCH3O2(I_0:I_1,J_0:J_1,:)   =1.d0 ! 1.d5 ??
        yC2O3(I_0:I_1,J_0:J_1,:)    =0.d0
        yROR(I_0:I_1,J_0:J_1,:)     =0.d0
        yXO2(I_0:I_1,J_0:J_1,:)     =0.d0
        yAldehyde(I_0:I_1,J_0:J_1,:)=0.d0
        yNO3(I_0:I_1,J_0:J_1,:)     =0.d0
        yXO2N(I_0:I_1,J_0:J_1,:)    =0.d0
        yRXPAR(I_0:I_1,J_0:J_1,:)   =0.d0
        oh_live(I_0:I_1,J_0:J_1,:)  =0.d0
        no3_live(I_0:I_1,J_0:J_1,:) =0.d0
        acetone(I_0:I_1,J_0:J_1,:)  =0.d0
        pClOx(I_0:I_1,J_0:J_1,:)    =1.d0
        pClx(I_0:I_1,J_0:J_1,:)     =0.d0
        pOClOx(I_0:I_1,J_0:J_1,:)   =0.d0
        pBrOx(I_0:I_1,J_0:J_1,:)    =1.d0
        yCl2(I_0:I_1,J_0:J_1,:)     =0.d0
        yCl2O2(I_0:I_1,J_0:J_1,:)   =0.d0
      END IF

 100  format(/3(50x,l1/),3(50x,i8/))
#ifdef TRACERS_AEROSOLS_SOA
#ifdef TRACERS_TERP
 110  format(6(///10(a8)),(///2(a8)))
#else
 110  format(6(///10(a8)),(///1(a8)))
#endif
#else
#ifdef TRACERS_TERP
 110  format(5(///10(a8)),(///4(a8)))
#else
 110  format(5(///10(a8)),(///3(a8)))
#endif
#endif  /* TRACERS_AEROSOLS_SOA */


      return
      END SUBROUTINE cheminit



      SUBROUTINE jplrts
!@sum jplrts read/set up chemical reaction rates from JPL
!@auth Drew Shindell (modelEifications by Greg Faluvegi)
!@ver  1.0 (based on cheminit0C5_M23p & ds4p_chem_init_M23)
!@calls lstnumc

C**** GLOBAL parameters and variables:
      USE DOMAIN_DECOMP_ATM, only: write_parallel
      USE FILEMANAGER, only: openunit,closeunit
      USE TRCHEM_Shindell_COM, only: nr,nr2,nr3,nmm,nhet,pe,ea,nst,ro,
     &                               r1,sn,sb,nn,nnr

      IMPLICIT NONE

C**** Local parameters and variables and arguments:
C
!@var ate temporary reactants names array
!@var i,ii,j dummy loop variable
!@var iu_data temporary unit number
      CHARACTER*8, DIMENSION(4) :: ate
      character(len=300) :: out_line
      INTEGER :: i,ii,j,iu_data

C Read in the number of each type of reaction:
      call openunit('JPLRX',iu_data,.false.,.true.)
      read(iu_data,124)nr,nr2,nr3,nmm,nhet
      write(out_line,*)' '
      call write_parallel(trim(out_line))
      write(out_line,*) 'Chemical reactions used in the model: '
      call write_parallel(trim(out_line))

      do i=1,nr               ! >>> begin loop over total reactions <<<
        if(i <= nr-nhet) then !non-hetero
          if(i <= nr2) then   !mono or bi
            if(i > nr2-nmm) then ! read monomolecular reactions
              if(i == nr2-nmm+1) read(iu_data,22)ate
              read(iu_data,16)ate,pe(i),ea(i),nst(i-nr2+nmm)
              write(out_line,30) i,ate(1),' + ',ate(2),
     &        ' --> ',ate(3),' + ',ate(4)
              call write_parallel(trim(out_line))
            else                  ! read bimolecular reactions
   5          read(iu_data,16)ate,pe(i),ea(i)
              write(out_line,30) i,ate(1),' + ',ate(2),
     &        ' --> ',ate(3),' + ',ate(4)
              call write_parallel(trim(out_line))
            end if
          else                    ! read trimolecular reactions
 20         if(i == nr2+1) read(iu_data,22)ate
            ii=i-nr2
            read(iu_data,21)ate,ro(ii),sn(ii),r1(ii),sb(ii)
            write(out_line,30) i,ate(1),' + ',ate(2),
     *      ' --> ',ate(3),' + ',ate(4)
            call write_parallel(trim(out_line))
          end if
        else                     ! read heterogeneous reactions
          if(i == nr-(nhet-1)) read(iu_data,22)ate
          read(iu_data,31)ate
          write(out_line,30) i,ate(1),' + ',ate(2),
     *    ' --> ',ate(3),' + ',ate(4)
          call write_parallel(trim(out_line))
        end if ! (i <= nr-nhet)
c
        do j=1,2
          call lstnum(ate(j),nn(j,i))
          call lstnum(ate(j+2),nnr(j,i))
        end do
      end do                ! >>> end loop over total reactions <<<

 124  format(///5(/43x,i3)///)
  27  format(/(30x,i2))
  21  format(4x,a8,1x,a8,3x,a8,1x,a8,e8.2,f5.2,e9.2,f4.1)
  22  format(/10x,4a8/)
  25  format(//32x,2f7.1,i6)
  16  format(4x,a8,1x,a8,3x,a8,1x,a8,e8.2,f8.0,i4)
  31  format(4x,a8,1x,a8,3x,a8,1x,a8)
  30  format(1x,i3,2x,a8,a3,a8,a5,a8,a3,a8)
      call closeunit(iu_data)
      return
      end SUBROUTINE jplrts



      SUBROUTINE lstnum(at,ks)
!@sum lstnum find molecule number in param list of molecules
!@auth Drew Shindell (modelEifications by Greg Faluvegi)
!@ver  1.0 (based on cheminit0C5_M23p & ds4p_chem_init_M23)

C**** GLOBAL parameters and variables:
      USE TRCHEM_Shindell_COM, only: nc,ay

      IMPLICIT NONE

C**** Local parameters and variables and arguments:
!@var at local copy of species name
!@var ks local variable to be passed back to jplrts nnr or nn array.
!@var j dummy loop variable
      INTEGER                  :: j
      INTEGER,     INTENT(OUT) :: ks
      CHARACTER*8, INTENT(IN)  :: at
      
      j=1
      do while(j <= nc)
        if(at == ay(j))then
          ks = j
          return
        else
          j = j + 1
          cycle
        endif
      enddo 
      ks = nc + 1
      if (at /= 'N2' .and. at /= 'H')
     &  call stop_model('ERROR: Tracer '//trim(at)//
     &    ' does not exist in the MOLEC file',255)

      return
      end SUBROUTINE lstnum



      subroutine fastj2_init
!@sum fastj2_init initialize fastj2 based on the currently active
!@+ chemistry scheme. It is the driver between the photolysis and
!@+ rest of chemistry code, so it can't be (in its current setup)
!@+ in the photolysis module
!@auth Kostas Tsigaridis

      use TRCHEM_Shindell_COM, only: iprn,jprn,prnrts,JPPJ_shindell
     &                              ,p_1
      use photolysis, only: phtlst,inphot
     &                     ,j_iprn,j_jprn,j_prnrts,jpnl,jppj,jlabel
     &                     ,jind,ks,kss,jfacta,zj
      implicit none

      j_iprn=iprn
      j_jprn=jprn
      j_prnrts=prnrts
      jppj=jppj_shindell

      allocate(jlabel(jppj))
      allocate(jind(jppj))
      allocate(ks(jppj))
      allocate(kss(p_1,jppj))
      allocate(jfacta(jppj))
      allocate(zj(jpnl,jppj))

C Read photolysis parameters and reactions from unit JPLPH:
      call phtlst

c fastj initialization routine:
      call inphot

      end subroutine fastj2_init


      SUBROUTINE reactn
!@sum reactn read chemical and photochemical reaction lists
!@auth Drew Shindell (modelEifications by Greg Faluvegi)
!@ver  1.0 (based on cheminit0C5_M23p & ds4p_chem_init_M23)
!@calls guide,printls

C**** GLOBAL parameters and variables:
      USE TRCHEM_Shindell_COM, only: nps,nds,kps,kds,nn,nnr,nr,
     &                      npnr,ndnr,kpnr,kdnr,prnls
      use photolysis, only: jppj,ks,kss

      IMPLICIT NONE

c Chemical reaction lists:
      call guide(npnr,ndnr,kpnr,kdnr,nn,nnr,2,nr)
c Photolysis reaction lists:
      call guide( nps, nds, kps, kds,ks,kss,1,JPPJ)
C Print out some diagnostics:
      if(prnls) call printls
      
      return
      end SUBROUTINE reactn



      SUBROUTINE guide(npr,ndr,kpr,kdr,xx,nnn,ns,nre)
!@sum guide read chemical and photochemical reaction lists
!@auth Drew Shindell (modelEifications by Greg Faluvegi)
!@ver  1.0 (based on cheminit0C5_M23p & ds4p_chem_init_M23)
!@calls calcls

C**** GLOBAL parameters and variables:
      USE TRCHEM_Shindell_COM, only: p_1,p_2,p_3,p_4

      IMPLICIT NONE

C**** Local parameters and variables and arguments:
!@var xx   = either nn  or ks   from reactn sub
!@var nnn  = either nnr or kss  from reactn sub
!@var kpr  = either kps or kpnr from reactn sub
!@var kdr  = either kds or kdnr from reactn sub
!@var npr  = either nps or npnr from reactn sub
!@var ndr  = either nds or ndnr from reactn sub
!@var ns   = either 1   or    2 from reactn sub
!@var nre number of reactions
      INTEGER,  DIMENSION(p_4)     :: kpr, kdr
      INTEGER,  DIMENSION(p_3)     :: npr, ndr
      INTEGER, DIMENSION(p_1,p_2)  :: xx, nnn
      INTEGER                      :: ns, nre

c Chemical and photolytic destruction:
      call calcls(xx,ns,nnn,2,ndr,kdr,nre)
c Chemical and photolytic production:
      call calcls(nnn,2,xx,ns,npr,kpr,nre)
      
      return
      end SUBROUTINE guide



      SUBROUTINE calcls(nn,ns,nnn,nns,ndr,kdr,nre)
!@sum calcls Set up reaction lists for calculated gases (1 to ny)
!@auth Drew Shindell (modelEifications by Greg Faluvegi)
!@ver  1.0 (based on cheminit0C5_M23p & ds4p_chem_init_M23)

C**** GLOBAL parameters and variables:
      USE DOMAIN_DECOMP_ATM, only: write_parallel
      USE TRCHEM_Shindell_COM, only: ny, numfam, p_2, p_3, p_4, nfam,
     &                               prnls

      IMPLICIT NONE

C**** Local parameters and variables and arguments:
!@var kdr  = either kdr or kpr from guide sub
!@var ndr  = either ndr or npr from guide sub
!@var nre  number of reactions
!@var nns number of partic_ on opposite side of reaction
!@var ns   = either ns or   2 from guide sub
!@var nn   = either xx or nnn from guide sub
!@var nnn  = either xx or nnn from guide sub
!@var ii,k,j,i,ij,i2,newfam,ifam dummy variables
      INTEGER, DIMENSION(p_4)    :: kdr
      INTEGER, DIMENSION(p_3)    :: ndr
      INTEGER :: nre, nns, ns, k, j, i, ij, i2, newfam, ifam, ii
      INTEGER, DIMENSION(ns,p_2) :: nn 
      INTEGER, DIMENSION(nns,p_2):: nnn
      character(len=300) :: out_line

      k=1
      nfam(numfam+1)=ny+1
      do j=1,numfam      !families, list only interfamily reactions
        newfam=0
        kdr(j)=k
        i_loop: do i=1,nre    ! 1 to # chem or phot reactions
          ij_loop: do ij=1,ns !ns # partic (prod & chem dest=2,phot dest=1)
            ! check if molecule # nn() is element of family j:
            if(nn(ij,i) >= nfam(j).and.nn(ij,i) < nfam(j+1))then
              ! check if reaction is intrafamily:
              do i2=1,nns  ! nns # partic on opposite side of reac.
                if(nnn(i2,i) >= nfam(j).and.nnn(i2,i) < nfam(j+1))
     &          cycle i_loop
              enddo
              ! don't write same reaction twice:
              if(k /= 1)then
                if(ndr(k-1) == i.and.newfam /= 0) cycle ij_loop
              endif
              ndr(k)=i
              k=k+1
              newfam=1
            endif
          enddo ij_loop
        enddo i_loop
      enddo

      do j=numfam+1,nfam(1)-1     ! individual non-family molecules
        kdr(j)=k
        do i=1,nre                ! 1 to # chem or phot reactions
          do ij=1,ns  !ns # partic (prod & chem dest=2,phot dest=1)
            if(nn(ij,i) /= j) cycle ! nn is mol # of participant
            ndr(k)=i
            k=k+1
          enddo
        enddo
      enddo

      do 100 j=nfam(1),ny !indiv family mols.,list only intrafamily
        do ii=1,numfam-1
          if(j < nfam(ii+1))then
            ifam=ii
            goto 110
          endif
        enddo
        ifam=numfam
 110    kdr(j)=k
        do 100 i=1,nre          ! 1 to # chem or phot reactions
          do 100 ij=1,ns  !ns # partic (prod & chem dest=2,phot dest=1)
            if(nn(ij,i) /= j)goto100       ! nn is mol # of participant
c           check that reaction is intrafamily
            do i2=1,nns  ! nns # participants on opposite side of reac.
              if(nnn(i2,i) >= nfam(ifam).and.nnn(i2,i) < nfam(ifam+1))
     &        then
                 ndr(k)=i
                 k=k+1
                 goto100
              endif
            enddo
 100  continue
      kdr(ny+1)=k

      if(prnls)then
        write(*,*) 'nn array size :',k
        write(out_line,*) 'nn array size :',k
        call write_parallel(trim(out_line))
      endif
      
      return
      end SUBROUTINE calcls



      SUBROUTINE printls
!@sum printls print out some chemistry diagnostics (reaction lists)
!@auth Drew Shindell (modelEifications by Greg Faluvegi)
!@ver  1.0 (based on cheminit0C5_M23p & ds4p_chem_init_M23)

C**** GLOBAL parameters and variables:
      USE DOMAIN_DECOMP_ATM, only: write_parallel
      USE TRCHEM_Shindell_COM, only: kpnr,npnr,kdnr,ndnr,kps,nps,
     &                         ny,nn,nnr,ay,kds,nds,nc
      use photolysis, only: ks,kss

      IMPLICIT NONE

C**** Local parameters and variables and arguments:
!@var ireac,igas,ichange,ii dummy variables
      INTEGER :: ireac,igas,ichange,ii
      character(len=300) :: out_line

c Print reaction lists:
      write(out_line,*) ' '
      call write_parallel(trim(out_line))
      write(out_line,*)
     &'______________ CHEMICAL PRODUCTION _______________'
      call write_parallel(trim(out_line))
      ireac=0
      do igas=1,ny
        write(out_line,*) ' '
        call write_parallel(trim(out_line))
        write(out_line,10) ay(igas)
        call write_parallel(trim(out_line))
        ichange=kpnr(igas+1)-kpnr(igas)
        if(ichange >= 1) then
          do ii=1,ichange
            ireac=ireac+1
            if (nnr(2,npnr(ireac)) > nc) then
              write(out_line,20)
     &        ' Reaction # ',npnr(ireac),' produces ',
     &        ay(nnr(1,npnr(ireac))),' and  ','X'
              call write_parallel(trim(out_line))
            else
              write(out_line,20)
     &        ' Reaction # ',npnr(ireac),' produces ',
     &        ay(nnr(1,npnr(ireac))),' and  ',ay(nnr(2,npnr(ireac)))
              call write_parallel(trim(out_line))
            end if
          enddo
        end if
      end do
      write(out_line,*) ' '
      call write_parallel(trim(out_line))
      write(out_line,*)
     &'______________ CHEMICAL DESTRUCTION _______________'
      call write_parallel(trim(out_line))
      ireac=0
      do igas=1,ny
        write(out_line,*) ' '
        call write_parallel(trim(out_line))
        write(out_line,10) ay(igas)
        call write_parallel(trim(out_line))
        ichange=kdnr(igas+1)-kdnr(igas)
        if(ichange >= 1) then
          do ii=1,ichange
            ireac=ireac+1
            write(out_line,20)
     &      ' Reaction # ',ndnr(ireac),' destroys ',
     *      ay(nn(1,ndnr(ireac))),' and  ',ay(nn(2,ndnr(ireac)))
            call write_parallel(trim(out_line))
          enddo
        end if
      end do
      write(out_line,*)
      call write_parallel(trim(out_line))
      write(out_line,*)
     &'______________ PHOTOLYTIC PRODUCTION _______________'
      call write_parallel(trim(out_line))
      ireac=0
      do igas=1,ny
        write(out_line,*) ' '
        call write_parallel(trim(out_line))
        write(out_line,10) ay(igas)
        call write_parallel(trim(out_line))
        ichange=kps(igas+1)-kps(igas)
        if(ichange >= 1) then
          do ii=1,ichange
            ireac=ireac+1
            write(out_line,20) ' Reaction # ',nps(ireac),' produces ',
     *      ay(kss(1,nps(ireac))),' and  ', ay(kss(2,nps(ireac)))
            call write_parallel(trim(out_line))
          enddo
        end if
      end do
      write(out_line,*) ' '
      call write_parallel(trim(out_line))
      write(out_line,*)
     & '______________ PHOTOLYTIC DESTRUCTION _______________'
      call write_parallel(trim(out_line))
      ireac=0
      do igas=1,ny
        write(out_line,*) ' '
        call write_parallel(trim(out_line))
        write(out_line,10) ay(igas)
        call write_parallel(trim(out_line))
        ichange=kds(igas+1)-kds(igas)
        if(ichange >= 1) then
          do ii=1,ichange
            ireac=ireac+1
            write(out_line,30) ' Reaction # ',nds(ireac),' destroys ',
     *      ay(ks(nds(ireac)))
            call write_parallel(trim(out_line))
          enddo
        end if
      end do
  10  format(1x,a8)
  20  format(a12,i4,a10,a8,a6,a8)
  30  format(a12,i4,a10,a8) 
       
      return
      end SUBROUTINE printls
