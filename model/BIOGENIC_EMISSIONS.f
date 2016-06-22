        MODULE BIOGENIC_EMIS
        USE MODEL_COM, only: IM,JM
        PARAMETER (NPOLYNB=20, NTYPE=16,NVEGTYPE=74)
        REAL*8 ISOPCOEFF(NPOLYNB),BASEISOP(IM,JM,nvegtype)

       END MODULE BIOGENIC_EMIS

       SUBROUTINE ISOPRENE_EMISSION(i,j,itype,emisop)

c          USE MODEL_COM, only: IM,JM
          USE BIOGENIC_EMIS
          USE TRACER_SOURCES, only: Isoprene_src
        use veg_com, only:  alaif, vdata
                !alaif = array(8,im,jm) of lai by veg functional type
                !vdata = array(im,jm,11) of cover fraction by func type
                !  In vdata, the 3rd subscript values 2..9 correspond to
                !  the 8 vegetation types.
                !The 8 veg functional types are:
                !tundr  grass  shrub  trees  decid evrgr  rainf crops
        use ghycom, only: canopy_temp_ij  
                !canopy_temp_ij = array(im,jm) of canopy temperatures (C)
                !Same temperature for all veg types over grid cell.
          USE RADNCB,     only : COSZ1,cfrac
          USE GEOM, only : BYDXYP
          USE SOCPBL, only : TGV
          USE PBL_DRV, only : psurf
          USE CONSTANT, only : rgas
          USE tracers_DRYDEP, only : xylai,IJREG,IREG,IJLAND,IJUSE

          INTEGER inveg
          INTEGER, INTENT(IN) :: I,J,ITYPE  
          REAL*8 TLAI,EMBIO,CLIGHT,TMMP,XLTMMP,BIOFIT,TCORR
          REAL*8 TMMPK(i,j),sfdens      
      real*8 , intent(out) :: emisop

      SELECT CASE(ITYPE)

      CASE(1:3) ! OCEAN, OCEAN ICE, LANDICE

      emisop=0.0

      CASE(4)     ! LAND

c Temperature of canpoy 

       sfdens=100.*psurf/(rgas*tgv) ! estimated surface density
       TMMPK(i,j) = canopy_temp_ij(i,j)+273.15

       EMISOP=0.
       TLAI=0.

      DO INVEG=1,IJREG(I,J)
            TLAI=TLAI+XYLAI(I,J,inveg)*BASEISOP(i,j,INVEG)

      END DO
C Light correction

C Only calculate for grid cell with sunlight and isoprene-emitting vegetation

      IF ((COSZ1(I,J).GT.0.0).AND.(TLAI.GT.0.)) THEN

c Note vdata array has 1 -11 types and alaif 1 8

         EMBIO=0.
         DO INVEG=1,ijreg(i,j)
            IF (xylai(I,J,inveg)*BASEISOP(i,j,INVEG).GT.0.) THEN
               CLIGHT=BIOFIT(ISOPCOEFF,xylai(I,J,inveg),
     &                             COSZ1(I,J),CFRAC(I,J))
               EMBIO=EMBIO+BASEISOP(i,j,INVEG)*
     &                          CLIGHT*ijuse(I,J,INVEG)/1000


            ENDIF
         END DO

C Temperature correction

         IF (TMMPK(i,j).GT.273.) THEN
            EMISOP=TCORR(TMMPK(i,j))*EMBIO

         ELSE
            EMISOP=0.

         ENDIF

      ENDIF  

C EMISOP = kg C emitted from grid cell per second
C 2D interative isoprene source 
C Convert to units kg C/m2/s


        EMISOP=EMISOP*BYDXYP(J)

      END SELECT

       RETURN                                                          
       END                              


      SUBROUTINE RDISOPCF                                              
C********************************************************              
C     Read polynomial coefficients from 'isoprene.coef'                
C********************************************************              

      USE BIOGENIC_EMIS
      USE FILEMANAGER, only: openunit,closeunit

      CHARACTER*80 DUM

C--  polynomial coefficients for biogenic isoprene emissions  

      call openunit('ISOPCF',iu_data,.false.,.true.)


       READ(iu_data,'(A80)') DUM

       READ(iu_data,'(8(1PE10.2))') (ISOPCOEFF(I),I=1,NPOLYNB)
       call closeunit(iu_data)

       RETURN                                                            
       END                       
C                                                                       
        SUBROUTINE RDISOBASE                                               
C**************************************************                     
C Read baseline emissions factors from 'isopemis.table'                 
C Units are atoms C cm^-2 leaf s^-1                                     
C Construct the base emission for each grid box                         
C Output is baseisop in kg C cm^-2 * surface area of cell (cm^2)        
C emitted in 1 hour time step                                           
C***************************************************                    
c                                                                       
      USE MODEL_COM, only: IM,JM
      USE BIOGENIC_EMIS
      USE tracers_DRYDEP, only : IJREG,IJLAND
      USE CONSTANT, only   : avog
      USE FILEMANAGER, only: openunit,closeunit
      use veg_com, only:   vdata
                !vdata = array(im,jm,11) of cover fraction by func type
                !  In vdata, the 3rd subscript values 2..9 correspond to
                !  the 8 vegetation types.
                !The 8 veg functional types are:
                !tundr  grass  shrub  trees  decid evrgr  rainf crops
      USE GEOM, only : DXYP
      INTEGER I,J,K                                              
      REAL*8 CONVERT(nvegtype),STEPH, FACTOR                       
      CHARACTER*80 DUM                                                  
c                                                                       

      call openunit('ISOPBASE',iu_data,.false.,.true.)


       READ(iu_data,'(A80)') DUM

       DO I=1,NVEGTYPE
       READ(iu_data,*) j, CONVERT(I)
       ENDDO

       call closeunit(iu_data)

C Compute the baseline ISOPRENE emissions, which depend on veg type   
                                                                       
C Emissions time step is 1 hour (STEPH)                                 

c         STEPH = 3600                                         

C Set up BASEISOP -- baseline ISOPRENE emissions                        
C Now hardwire molecular weight for Carbon = 0.012 kg/mol               
C ISOPRENE is traced in terms of equivalent C atoms                     
C                                                                       

            DO J=1,JM                                                         
c                                                                       
c         FACTOR = 12d-3*STEPH*DXYP(J)*1.D+04/avog                        

         FACTOR = 12d-3*DXYP(J)*1.D+04/avog      
c                                                                       
c 1D+4 because the convert data from file is / 1000        
c                                                                       
            DO I=1,IM                                                      
            DO K=1,ijreg(i,j)

         BASEISOP(i,j,K) = CONVERT(ijland(i,j,k)+1)*FACTOR        

            ENDDO                                                       
         ENDDO                                                          
      ENDDO                                                             
c                                                                       
      RETURN                                                            
      END                                                               
c                     

C TCORR function
************************************************************
C Temperature correction for isoprene emissions
C Sensitivity to local air temperature
***********************************************************
	
       REAL*8 FUNCTION TCORR(TEMP)

      IMPLICIT NONE

c temperature correction for isoprene emissions, Guenther et al.(92)

       REAL*8 R,CT1,CT2,T1,T3,TEMP
       DATA R,CT1,CT2,T1,T3 /8.314,95000.,230000.,303.,314./
       TCORR =
     *    EXP( CT1/(R*T1*TEMP)*(TEMP-T1) ) /
     *    (1 + EXP( CT2/(R*T1*TEMP)*(TEMP-T3) ))
        RETURN
        END





