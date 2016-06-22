      SUBROUTINE GRDTRNS (MODEV,
     -                   QI, FI, DFI, IDIN, JDIN, KDIN, PKI, LRFI,
     -                   PO, FO,      IDOU, JDOU, KDOU, QKO, LRFO, XYIJ)
C
C
C
      DIMENSION  QI(IDIN,JDIN,KDIN),  FI(IDIN,JDIN,KDIN),
     -                               DFI(IDIN,JDIN,KDIN),
     -           PKI(KDIN)
      DIMENSION  PO(IDOU,JDOU,KDOU),  FO(IDOU,JDOU,KDOU),
     -           QKO(KDOU)
      DIMENSION  XYIJ(1)
      DIMENSION  LRFI(6), LRFO(6)
      DIMENSION  FAB(2), QAB(2), DFAB(2)
      EQUIVALENCE  (FAB(1),FIJA),  (FAB(2),FIJB),
     -             (QAB(1),QIJA),  (QAB(2),QIJB),
     -             (DFAB(1),DFA),  (DFAB(2),DFB)
      DATA  EPSC  /1.0E-8/
      COMMON  /HMTCM1/ XIJ,NXI,DX,IPTH,
     -                 YIJ,MYJ,DY,JPTH,
     -                 IPL,IPU
C
C
C        STATEMENT FUNCTIONS
C
      IFLOR(X) = INT(X) - INT(.5-SIGN(.5,X))
      HMTINT (PX,DX,DFX1,FX1,FX2,DFX2) = FX1 + PX*(PX*(2.*PX-3.)*
     -             (FX1-FX2) + (PX-1.)*DX*(PX*(DFX1+DFX2) - DFX1)   )
C
C
C
C  =====================================================================
C  =====================================================================
C                   ==                        ==
C                   ==                        ==
C                   ==   SUBROUTINE GRDTRNS   ==
C                   ==                        ==
C                   ==                        ==
C  =====================================================================
C  =====================================================================
C
C
C  *** PURPOSE ***
C
C        GIVEN VOLUMES OF A FUNCTION  F  AND NEW VERTICAL COORDINATE  P
C        ON AN OLD COORDINATE GRID  (XI,YJ,QK),  TO TRANSFORM TO VOLUMES
C        OF  F  AND THE OLD VERTICAL COORDINATE  Q  ON A NEW COORDINATE
C        GRID  (XII,YJJ,PKK).
C
C
C  *** METHOD ***
C
C        TRANSFORMATION TO NEW HORIZONTAL GRID POINTS ON OLD VERTICAL
C        COORDINATE SURFACES IS VIA PIECEWISE BICUBIC QUASI-HERMITE
C        INTERPOLATION -- QUASI IN THE SENSE THAT THE DERIVATIVES D*/DX
C        AND D*/DY ARE OBTAINED BY SECOND ORDER CENTERED FINITE DIFFER-
C        ENCES.  VERTICAL TRANSFORM IS THEN PERFORMED TO NEW VERTICAL
C        COORDINATE SURFACES AT THE NEW HORIZONTAL GRID POINT USING
C        EITHER LINEAR OR HERMITE CUBIC INTERPOLATION.
C
C
C  *** ARGUMENT LIST ***
C
C        ON INPUT
C
C          QI - A THREE DIMENSIONAL VOLUME OF NEW VERTICAL COORDINATE
C                  DATA ON THE OLD COORDINATE GRID.
C          FI - A THREE DIMENSIONAL VOLUME OF FUNCTION VALUES ON THE OLD
C                  COORDINATE GRID.
C          DFI - A VOLUME OF FUNCTION VALUE VERTICAL DERIVATIVES (I.E., OF
C                  THE FI FUNCTION), OPTIONAL, NEED ONLY BE SUPPLIED IF
C                  HERMITE INTERPOLATION IS TO BE USED IN THE VERTICAL.
C          IDIN - FIRST DIMENSION OF QI, FI, DFI IN THE CALLING PROGRAM.
C          JDIN - SECOND DIMENSION OF QI, FI, DFI IN CALLING PROGRAM.
C          KDIN - NUMBER OF VERTICAL LEVELS IN THE VOLUMES QI, FI, DFI
C                  FOR THIS CALL.  MUST BE AT LEAST 2.
C          PKI - VECTOR OF OLD VERTICAL COORDINATE VALUES ASSOCIATED
C                  WITH VERTICAL LEVELS  1,...,KDIN IN INPUT VOLUMES.
C                  MUST BE MONOTONE (INCREASING OR DECREASING).
C          MODEV - VERTICAL INTERPOLATION MODE FLAG,
C             = -1,  HERMITE, AND DFI CONTAINS DF/DQ, THE DERIVATIVE
C                  WITH RESPECT TO THE NEW VERTICAL COORDINATE.
C             =  0,  LINEAR.
C             = +1,  HERMITE, AND DFI CONTAINS DF/DP, THE DERIVATIVE
C                  WITH RESPECT TO THE OLD VERTICAL COORDINATE.
C          LRFI - WE DEFINE THE  'TOTAL OLD COORDINATE SPACE'  TO BE THE
C                  FULL THREE DIMENSIONAL GRID IN OLD COORDINATE SPACE
C                  FOR WHICH THE USER HAS NEW COORD AND FUNCTION DATA.
C                  ON THIS CALL TO THIS ROUTINE, THE DATA IN THE INPUT
C                  ARRAYS MAY ENTIRELY SPAN THE HORIZONTAL DIMENSIONS OF
C                  THE TOTAL OLD COORD SPACE, OR ALTERNATELY, THE DATA
C                  MAY ONLY SPAN A SUBSPACE IN ONE OR BOTH OF THE HORI-
C                  ZONTAL DIMENSIONS.  IN THE LATTER CASE, IT IS ASSUMED
C                  THAT THIS CALL IS ONE OF A SEQUENCE OF CALLS WHICH
C                  WILL CYCLE THROUGH THE TOTAL INPUT SPACE AND THUS
C                  SPAN THE DEFICIENT DIMENSION(S).  IN EITHER CASE, THE
C                  ACTUAL DATA FOR THIS CALL MAY NOT FILL THE USER
C                  DECLARED DIMENSIONS OF QI, FI, DFI -- IDIN AND JDIN.
C             LRFI(1),LRFI(2) - THE DATA FOR THIS CALL, FOR LEVEL *, IS IN
C                       ARRAY(1,1,*)  TO  ARRAY(LRFI(1),LRFI(2),*)  .
C                  IF .LE.0, IDIN AND/OR JDIN TAKEN AS DEFAULTS.  OTHER-
C                  WISE MUST BE .GE.4, OR AN ERROR IS DIAGNOSED.
C             LRFI(3),LRFI(4) - THE HORIZONTAL FIRST AND SECOND DIMENSIONS
C                  OF THE USERS TOTAL OLD COORDINATE SPACE.  IF .LE.0,
C                  ASSUMED TO BE SAME AS LRFI(1),LRFI(2)  (OR THEIR
C                  DEFAULTS IDIN,JDIN ).
C             LRFI(5),LRFI(6) - HORIZONTAL GRID POINT  (1,1,*)  OF THE
C                  INPUT VOLUMES CORRESPONDS TO HORIZONTAL GRID POINT
C                  (LRFI(5),LRFI(6),*) IN THE USERS TOTAL OLD COORD
C                  SPACE.  IF .LE.0, DEFAULT VALUE(S) 1 ARE USED.
C             ** NOTE **  THE INFORMATION IN LRFI ALLOWS US TO COMPLETELY
C             LOCATE THE HORIZONTAL GRID OF THE ACTUAL DATA FOR THIS CALL
C             AS A SUBGRID OF THE TOTAL HORIZONTAL GRID FROM WHICH THE
C             USER IS WORKING.  FROM THIS, WE CAN TELL WHEN AN APPARENT
C             BOUNDARY GRID SQUARE IS ACTUALLY A BOUNDARY SQUARE IN THE
C             TOTAL OLD COORD SPACE, AND INTERPOLATE WITH A BOUNDARY
C             FORMULA, OR NOT INTERPOLATE, AS APPROPRIATE.
C          IDOU - FIRST DIMENSION OF ARRAYS PO, FO IN CALLING PROGRAM.
C          JDOU - SECOND DIMENSION OF ARRAYS PO, FO IN CALLING PROGRAM.
C          KDOU - NUMBER OF VERTICAL LEVELS OF NEW COORDINATE SPACE TO
C                  CONSTRUCT ON THIS CALL.
C          QKO - VECTOR OF NEW VERTICAL COORDINATE VALUES OF THE KDOU
C                  OUTPUT LEVELS DESIRED ON THIS CALL.  MUST BE MONOTONE.
C          LRFO - SIMILAR TO LRFI, EXCEPT APPLIES TO THAT HORIZONTAL
C                  PORTION OF THE TOTAL NEW COORDINATE SPACE (OUTPUT
C                  GRID) WHICH IS TO BE CONSTRUCTED ON THIS CALL.
C             LRFO(1),LRFO(2) - ON THIS CALL, CONSTRUCT A LRFO(1) BY
C                  LRFO(2) SUBGRID OF TOTAL NEW COORD SPACE AND STORE IN
C                  ARRAY(1,1,*)  TO  ARRAY(LRFO(1),LRFO(2),*)  OF OUTPUT
C                  ARRAYS.  DEFAULT VALUES OF IDOU AND/OR JDOU ARE USED
C                  IF .LE.0.
C             LRFO(3),LRFO(4) - HORIZONTAL FIRST AND SECOND DIMENSIONS OF
C                  USERS TOTAL NEW COORD SPACE.  IF .LE.0, DEFAULT VALUES
C                  LRFO(1) AND/OR LRFO(2) (OR THEIR DEFAULTS) ARE USED.
C             LRFO(5),LRFO(6) - HORIZONTAL ORIGIN OF THIS-CALL-SUBGRID
C                  OF TOTAL NEW COORD SPACE IN HORIZONTAL GRID OF NEW
C                  COORD SPACE.  DEFAULT VALUES 1 AND/OR 1 IF .LE.0 .
C          XYIJ - ARRAY CONTAINING OLD HORIZONTAL GRID POINT COORDS OF
C                  NEW HORIZONTAL GRID POINTS.  THE HORIZONTAL FIRST AND
C                  SECOND DIMENSIONS SHOULD ALWAYS EXACTLY MATCH THOSE OF
C                  THE TOTAL NEW COORD SPACE (LRFO(3) AND LRFO(4) OR
C                  THEIR DEFAULTS), OR ELSE THE ARRAY SHOULD BE LOADED IN
C                  SUCH A WAY TO MIMIC THIS DIMENSIONALITY.  THUS FOR
C                  IO=1,..,LRFO(3)  AND  JO=1,..,LRFO(4)  WE SHOULD HAVE
C                       XYIJ(IO,JO,1) = I OLD GRID COORD OF POINT
C                                       (IO,JO,*) IN TOTAL NEW GRID
C                       XYIJ(IO,JO,2) = J OLD GRID COORD OF POINT
C                                       (IO,JO,*) IN TOTAL NEW GRID.
C                  ON THIS CALL WE WILL BE EXAMINING
C                       XYIJ( LRFO(5), LRFO(6), *)
C                         TO
C                       XYIJ( LRFO(5)+LRFO(1)-1, LRFO(6)+LRFO(2)-1, *) .
C             **NOTE**  THE USER MAY WRITE A  'REAL FUNCTION XYIJ(IDX)',
C                  DECLARE IT EXTERNAL IN THE CALLING PROGRAM, AND DELETE
C                  THE DIMENSION STATEMENT ABOVE.  THE SINGLE ARGUMENT
C                  IDX WILL BE THE LINEAR SUBSCRIPT EQUIVALENT OF THE
C                  THREE DIMENSIONAL SUBSCRIPT  (IO,JO,1)  OR  (IO,JO,2).
C
C        ON OUTPUT
C
C          PO - THREE DIMENSIONAL VOLUME OF OLD VERTICAL COORDINATE DATA
C                  P TRANSFORMED TO NEW COORDINATE SPACE.
C          FO - THREE DIMENSIONAL VOLUME OF FUNCTION VALUES ON THE NEW
C                  COORDINATE GRID.  IF ANY GRID POINT OF THE PORTION OF
C                  THE TOTAL NEW COORDINATE GRID TO BE CONSTRUCTED THIS
C                  CALL IS EXTERIOR TO THE PORTION OF THE OLD COORD GRID
C                  SUPPLIED FOR THIS CALL, THEN THE CORRESPONDING ENTRIES
C                  OF THE OUTPUT ARRAYS PO AND FO WILL BE LEFT UNDEFINED.
C
C
C  =====================================================================
C  =====================================================================
C
C
C
C                      **************************
C                      *                        *
C                      ****  INITIALIZATION  ****
C                      *                        *
C                      **************************
C
C        CHECK INPUT DIMENSIONS AND POINTERS, SET DEFAULTS
C
      KDI = KDIN
      IDI = LRFI(1)
      IF (IDI.LE.0)  IDI = IDIN
      JDI = LRFI(2)
      IF (JDI.LE.0)  JDI = JDIN
      IIT = LRFI(3)
      IF (IIT.LE.0)  IIT = IDI
      JIT = LRFI(4)
      IF (JIT.LE.0)  JIT = JDI
      IIL = MAX0 (1,LRFI(5))
      JIL = MAX0 (1,LRFI(6))
      IIU = IIL + IDI - 1
      JIU = JIL + JDI - 1
      IF (  MAX0(4,MIN0(IDI,IDIN))  .NE. IDI )  GO TO 9001
      IF (  MAX0(4,MIN0(JDI,JDIN))  .NE. JDI )  GO TO 9002
C
C        CHECK OUTPUT DIMENSIONS AND POINTERS, SET DEFAULTS
C
      KDO = KDOU
      IDO = LRFO(1)
      IF (IDO.LE.0)  IDO = IDOU
      JDO = LRFO(2)
      IF (JDO.LE.0)  JDO = JDOU
      IOT = LRFO(3)
      IF (IOT.LE.0)  IOT = IDO
      JOT = LRFO(4)
      IF (JOT.LE.0)  JOT = JDO
      IOL = MAX0 (1,LRFO(5))
      JOL = MAX0 (1,LRFO(6))
      IOU = IOL + IDO - 1
      JOU = JOL + JDO - 1
      IF (IDO.GT.IDOU  .OR.  JDO.GT.JDOU)  GO TO 9011
      IF (IOU.GT.IOT   .OR.   JOU.GT.JOT)  GO TO 9012
C
C        SET BOUNDARY INTERPOLATION AND VERTICAL MODE FLAGS
C
      IBDL = MIN0 (IIL,2)
      JBDL = MIN0 (JIL,2)
      IBDU = MIN0 (IIT-IIU+1, 2)
      JBDU = MIN0 (JIT-JIU+1, 2)
      JPHMT = MODEV
C
C        PRESET OTHER CONSTANTS
C
      IDQ = IDIN
      IJOT = IOT*JOT
      IDM2 = IDI - 2
      JDM2 = JDI - 2
      IDM1 = IDI - 1
      JDM1 = JDI - 1
      XNRM = FLOAT(IIL) - 1.0
      YNRM = FLOAT(JIL) - 1.0
      XMX = FLOAT(IDI)
      YMX = FLOAT(JDI)
C
C
C                     *****************************
C                     *                           *
C                     ***  GRID TRANSFORMATION  ***
C                     *                           *
C                     *****************************
C
C
C  ===== OUTPUT GRID HORIZONTAL LOOP =====
C
      IOS = 0
      DO 7001 IXO = IOL,IOU
      IOS = IOS + 1
      JOS = 0
      DO 7000 JYO = JOL,JOU
      JOS = JOS + 1
      IDX = IXO + (JYO-1)*IOT
      JDY = IDX + IJOT
      XI = XYIJ(IDX) - XNRM
      YJ = XYIJ(JDY) - YNRM
      IF (XI.GT.XMX  .OR.  YJ.GT.YMX)  GO TO 7000
      NXI = MIN0 (IDM1,IFLOR(XI))
      MYJ = MIN0 (JDM1,IFLOR(YJ))
C
C                  COMPUTE PROCESSING PATH SELECTORS
C
      IPTH = MAX0 (IBDU*(NXI-IDM2), IBDL*MIN0(0,NXI-2) )
      JPTH = MAX0 (JBDU*(MYJ-JDM2), JBDL*MIN0(0,MYJ-2) )
      IF (  MAX0(IABS(IPTH),IABS(JPTH))  .GE. 2 )  GO TO 7000
      DX = XI - FLOAT (NXI)
      DY = YJ - FLOAT(MYJ)
      IPL = MAX0 (NXI-1, 1)
      IPU = MIN0 (NXI+2, IDI)
C
C  ===== INPUT VERTICAL LEVELS LOOP =====
C
         LAST = -1
         LVI = 0
 1900    LVI = LVI + 1
         IF (LVI.GT.KDI) GO TO 7000
         QIJB = QIJA
         FIJB = FIJA
         DFB = DFA
         QIJA = QSIHMT (QI(1,1,LVI),IDQ)
         IF (LVI.EQ.1)  GO TO 1900
C
C  ===== OUTPUT VERTICAL LEVELS LOOP =====
C
            LVO = 0
 2000       LVO = LVO + 1
            IF (LVO.GT.KDO) GO TO 1900
            IF ( (QKO(LVO)-QIJA) * (QKO(LVO)-QIJB) .GT. 0.0)  GO TO 2000
            DQ = QIJA - QIJB
            IF (DQ.EQ.0.) THEN
C           WRITE (*,'(3I4'' QIJA=QIJB=''1PE11.3)') IOS,JOS,LV0,QIJA
            PQ = .5
            ELSE
            PQ = (QKO(LVO)-QIJB)/DQ
            ENDIF
            OMPQ = 1.0 - PQ
            LVBRN = MAX0 (-1, (LAST+1) - LVI )
            IF (PQ*OMPQ .GT. EPSC)  GO TO 2015
            IDXC = IFIX(LVI - .5 + PQ)
            PO(IOS,JOS,LVO) = PKI(IDXC)
            LF1 = LVI - IDXC + 1
            LF2 = LF1 + LVBRN + 1
            IF (LF2.GE.3) GO TO 2004
            FAB(LF1) = QSIHMT (FI(1,1,IDXC),IDQ)
 2004       FO(IOS,JOS,LVO) = FAB(LF1)
            GO TO (2000,2006,2000,2000), LF2
 2006       LAST = IDXC
            IF (JPHMT) 2008,2000,2008
 2008       DFAB(LF1) = QSIHMT(DFI(1,1,IDXC),IDQ)
            GO TO 2000
 2015       POK = PQ*PKI(LVI) + OMPQ*PKI(LVI-1)
            LAST = LVI
            PO(IOS,JOS,LVO) = POK
            IF (LVBRN) 2025,2030,2035
 2025       FIJB = QSIHMT (FI(1,1,LVI-1),IDQ)
 2030       FIJA = QSIHMT (FI(1,1,LVI),IDQ)
 2035       IF (JPHMT) 2050,2040,2060
 2040       FO(IOS,JOS,LVO) = PQ*FIJA + OMPQ*FIJB
            GO TO 2000
 2050       DC = DQ
            GO TO 2070
 2060       DC = PKI(LVI) - PKI(LVI-1)
 2070       PC = PQ
            IF (LVBRN) 2075,2080,2085
 2075       DFB = QSIHMT (DFI(1,1,LVI-1),IDQ)
 2080       DFA = QSIHMT (DFI(1,1,LVI),IDQ)
 2085       FO(IOS,JOS,LVO) = HMTINT (PC,DC,DFB,FIJB,FIJA,DFA)
            GO TO 2000
 7000 CONTINUE
 7001 CONTINUE
C
C        SUCCESSFUL COMPLETION, NORMAL RETURN
C
 7007 RETURN
C
C        ERROR CONDITIONS
C
 9001 STOP 99001
 9002 STOP 99002
 9011 STOP 99011
 9012 STOP 99012
C
      END
C
C
      REAL FUNCTION QSIHMT(FCN,IDF)
      DIMENSION  FCN(IDF,1)
      REAL  RW(4)
      COMMON  /HMTCM1/ XIJ,NXI,DX,IPTH,
     -                 YIJ,MYJ,DY,JPTH,
     -                 IPL,IPU
      HMTINT (DL,F1,F2,DF1,DF2) =
     -        F1 + DL*( DL*(2.*DL-3.)*(F1-F2) +
     -               (DL-1.)*(DL*(DF1+DF2)-DF1)  )
C
      IW = 0
      DO 1050 IR = IPL,IPU
      IW = IW + 1
      FCNL = FCN(IR,MYJ)
      FCNU = FCN(IR,MYJ+1)
      IF (JPTH) 1020,1015,1010
 1010 DRVU = .5*FCN(IR,MYJ-1) - 2.*FCN(IR,MYJ) + 1.5*FCN(IR,MYJ+1)
      DRVL = .5*(FCN(IR,MYJ+1) - FCN(IR,MYJ-1))
      GO TO 1025
 1015 DRVL = .5*(FCN(IR,MYJ+1) - FCN(IR,MYJ-1))
      DRVU = .5*(FCN(IR,MYJ+2) - FCN(IR,MYJ))
      GO TO 1025
 1020 DRVL = -1.5*FCN(IR,MYJ) + 2.*FCN(IR,MYJ+1) - .5*FCN(IR,MYJ+2)
      DRVU = .5*(FCN(IR,MYJ+2) - FCN(IR,MYJ))
 1025 CONTINUE
      RW(IW) = HMTINT(DY,FCNL,FCNU,DRVL,DRVU)
 1050 CONTINUE
C
      JYW = 2 + MIN0(IPTH,0)
      FCNL = RW(JYW)
      FCNU = RW(JYW+1)
      IF (IPTH) 2020,2015,2010
 2010 DRVU = .5*RW(JYW-1) - 2.*RW(JYW) + 1.5*RW(JYW+1)
      DRVL = .5*(RW(JYW+1) - RW(JYW-1))
      GO TO 2025
 2015 DRVL = .5*(RW(JYW+1) - RW(JYW-1))
      DRVU = .5*(RW(JYW+2) - RW(JYW))
      GO TO 2025
 2020 DRVL = -1.5*RW(JYW) + 2.*RW(JYW+1) - .5*RW(JYW+2)
      DRVU = .5*(RW(JYW+2) - RW(JYW))
 2025 CONTINUE
C
      QSIHMT = HMTINT(DX,FCNL,FCNU,DRVL,DRVU)
C
      RETURN
      END
