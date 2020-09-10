C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
C     Automatically generated code
C
C     Tianfeng Lu
C     University of Connecticut
C     191 Auditorium Road U-3139
C     Storrs, CT 06269, USA
C     Email: tlu@engr.uconn.edu
C
C     September 18, 2014
C
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
C      SUBROUTINE GETRATES  (P, T, Y, DIFF, DT, ICKWRK, RCKWRK, WDOT)
C      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C      DIMENSION Y(*),DIFF(*),WDOT(*),ICKWRK(*),RCKWRK(*)
C      DIMENSION RF(268),RB(268),RKLOW(26),C(35)
C      DIMENSION XQ(18)
C
C      CALL YTCP(P, T, Y, C)
C      CALL RATT(T, RF, RB, RKLOW)
C      CALL RATX(T, C, RF, RB, RKLOW)
C      CALL QSSA(RF, RB, XQ)
C      CALL STIF(RF, RB, DIFF, DT, C)
C      CALL RDOT(RF, RB, WDOT)
C
C      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE CKWC  (T, C, WDOT)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
      DIMENSION C(*),WDOT(*)
      DIMENSION RF(268),RB(268),RKLOW(26), DIFF(35)
      DIMENSION XQ(18)
C
      CALL RATT(T, RF, RB, RKLOW)
      CALL RATX(T, C, RF, RB, RKLOW)
      CALL QSSA(RF, RB, XQ)
C      DIFF(1:35) = 0.d0
C      CALL STIF(RF, RB, DIFF, DT, C)
      CALL RDOT(RF, RB, WDOT)
C
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
C      SUBROUTINE CKWYR  (RHO, T, Y, ICKWRK, RCKWRK, WDOT)
C      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C      DIMENSION Y(*),WDOT(*),ICKWRK(*),RCKWRK(*)
C      DIMENSION RF(268),RB(268),RKLOW(26),C(35)
C      DIMENSION XQ(18)
C
C      CALL YTCR(RHO, T, Y, C)
C      CALL RATT(T, RF, RB, RKLOW)
C      CALL RATX(T, C, RF, RB, RKLOW)
C      CALL QSSA(RF, RB, XQ)
C      CALL RDOT(RF, RB, WDOT)
C
C      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE YTCP (P, T, Y, C)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
      DIMENSION Y(*), C(*)
C
      C(1) = Y(1)/1.70341023683548D2
      C(2) = Y(2)/1.007969975471497D0
      C(3) = Y(3)/1.599940013885498D1
      C(4) = Y(4)/1.700737011432648D1
      C(5) = Y(5)/3.300677025318146D1
      C(6) = Y(6)/2.015939950942993D0
      C(7) = Y(7)/1.801534008979797D1
      C(8) = Y(8)/3.401474022865295D1
      C(9) = Y(9)/3.199880027770996D1
      C(10) = Y(10)/1.503506028652191D1
      C(11) = Y(11)/1.604303026199341D1
      C(12) = Y(12)/3.00264904499054D1
      C(13) = Y(13)/2.80105504989624D1
      C(14) = Y(14)/4.400995063781738D1
      C(15) = Y(15)/2.603824067115784D1
      C(16) = Y(16)/2.805418062210083D1
      C(17) = Y(17)/3.007012057304382D1
      C(18) = Y(18)/4.304561078548431D1
      C(19) = Y(19)/4.107330095767975D1
      C(20) = Y(20)/4.208127093315125D1
      C(21) = Y(21)/5.606473112106323D1
      C(22) = Y(22)/5.510039126873016D1
      C(23) = Y(23)/5.610836124420166D1
      C(24) = Y(24)/6.912748157978058D1
      C(25) = Y(25)/7.013545155525208D1
      C(26) = Y(26)/8.416254186630249D1
      C(27) = Y(27)/9.818963217735291D1
      C(28) = Y(28)/1.122167224884033D2
      C(29) = Y(29)/1.262438127994537D2
      C(30) = Y(30)/1.272517827749252D2
      C(31) = Y(31)/1.402709031105042D2
      C(32) = Y(32)/1.68325083732605D2
      C(33) = Y(33)/2.013318539857864D2
      C(34) = Y(34)/2.163232841491699D2
      C(35) = Y(35)/2.801339912414551D1
C
      SUM = 0D0
      DO K = 1, 35
         SUM = SUM + C(K)
      ENDDO
      SUM = P/(SUM*T*8.314510D7)
C
      DO K = 1, 35
         C(K) = C(K)*SUM
      ENDDO
C
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE YTCR (RHO, T, Y, C)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
      DIMENSION Y(*), C(*)
      DATA SMALL/1D-50/
C
C     NC12H26
      C(1) = Y(1)/1.70341023683548D2
C     H
      C(2) = Y(2)/1.007969975471497D0
C     O
      C(3) = Y(3)/1.599940013885498D1
C     OH
      C(4) = Y(4)/1.700737011432648D1
C     HO2
      C(5) = Y(5)/3.300677025318146D1
C     H2
      C(6) = Y(6)/2.015939950942993D0
C     H2O
      C(7) = Y(7)/1.801534008979797D1
C     H2O2
      C(8) = Y(8)/3.401474022865295D1
C     O2
      C(9) = Y(9)/3.199880027770996D1
C     CH3
      C(10) = Y(10)/1.503506028652191D1
C     CH4
      C(11) = Y(11)/1.604303026199341D1
C     CH2O
      C(12) = Y(12)/3.00264904499054D1
C     CO
      C(13) = Y(13)/2.80105504989624D1
C     CO2
      C(14) = Y(14)/4.400995063781738D1
C     C2H2
      C(15) = Y(15)/2.603824067115784D1
C     C2H4
      C(16) = Y(16)/2.805418062210083D1
C     C2H6
      C(17) = Y(17)/3.007012057304382D1
C     CH2CHO
      C(18) = Y(18)/4.304561078548431D1
C     aC3H5
      C(19) = Y(19)/4.107330095767975D1
C     C3H6
      C(20) = Y(20)/4.208127093315125D1
C     C2H3CHO
      C(21) = Y(21)/5.606473112106323D1
C     C4H7
      C(22) = Y(22)/5.510039126873016D1
C     C4H81
      C(23) = Y(23)/5.610836124420166D1
C     C5H9
      C(24) = Y(24)/6.912748157978058D1
C     C5H10
      C(25) = Y(25)/7.013545155525208D1
C     C6H12
      C(26) = Y(26)/8.416254186630249D1
C     C7H14
      C(27) = Y(27)/9.818963217735291D1
C     C8H16
      C(28) = Y(28)/1.122167224884033D2
C     C9H18
      C(29) = Y(29)/1.262438127994537D2
C     PXC9H19
      C(30) = Y(30)/1.272517827749252D2
C     C10H20
      C(31) = Y(31)/1.402709031105042D2
C     C12H24
      C(32) = Y(32)/1.68325083732605D2
C     C12H25O2
      C(33) = Y(33)/2.013318539857864D2
C     OC12H23OOH
      C(34) = Y(34)/2.163232841491699D2
C     N2
      C(35) = Y(35)/2.801339912414551D1
C
      DO K = 1, 35
         C(K) = RHO * C(K)
      ENDDO
C
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE RATT (T, RF, RB, RKLOW)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
      PARAMETER (RU=8.31451D7, SMALL=1.D-200, PATM=1.01325D6)
      DIMENSION RF(*), RB(*), RKLOW(*)
      DIMENSION SMH(52), EG(52)
C
      ALOGT = LOG(T)
      TI = 1D0/T
      TI2 = TI*TI
C
      CALL RDSMH (T, SMH)
      EG(1) = EXP(SMH(1))
      EG(2) = EXP(SMH(2))
      EG(3) = EXP(SMH(3))
      EG(4) = EXP(SMH(4))
      EG(5) = EXP(SMH(5))
      EG(6) = EXP(SMH(6))
      EG(7) = EXP(SMH(7))
      EG(8) = EXP(SMH(8))
      EG(9) = EXP(SMH(9))
      EG(10) = EXP(SMH(10))
      EG(11) = EXP(SMH(11))
      EG(12) = EXP(SMH(12))
      EG(13) = EXP(SMH(13))
      EG(14) = EXP(SMH(14))
      EG(15) = EXP(SMH(15))
      EG(16) = EXP(SMH(16))
      EG(17) = EXP(SMH(17))
      EG(18) = EXP(SMH(18))
      EG(19) = EXP(SMH(19))
      EG(20) = EXP(SMH(20))
      EG(21) = EXP(SMH(21))
      EG(22) = EXP(SMH(22))
      EG(23) = EXP(SMH(23))
      EG(24) = EXP(SMH(24))
      EG(25) = EXP(SMH(25))
      EG(26) = EXP(SMH(26))
      EG(27) = EXP(SMH(27))
      EG(28) = EXP(SMH(28))
      EG(29) = EXP(SMH(29))
      EG(30) = EXP(SMH(30))
      EG(31) = EXP(SMH(31))
      EG(32) = EXP(SMH(32))
      EG(33) = EXP(SMH(33))
      EG(34) = EXP(SMH(34))
      EG(35) = EXP(SMH(35))
      EG(36) = EXP(SMH(36))
      EG(37) = EXP(SMH(37))
      EG(38) = EXP(SMH(38))
      EG(39) = EXP(SMH(39))
      EG(40) = EXP(SMH(40))
      EG(41) = EXP(SMH(41))
      EG(42) = EXP(SMH(42))
      EG(43) = EXP(SMH(43))
      EG(44) = EXP(SMH(44))
      EG(45) = EXP(SMH(45))
      EG(46) = EXP(SMH(46))
      EG(47) = EXP(SMH(47))
      EG(48) = EXP(SMH(48))
      EG(51) = EXP(SMH(51))
      EG(52) = EXP(SMH(52))
      PFAC1 = PATM / (RU*T)
      PFAC2 = PFAC1*PFAC1
      PFAC3 = PFAC2*PFAC1
C
C
C     R1: H + O2 = O + OH
      RF(1) = EXP(3.221148868927627D1 -7.468872617445274D3*TI)
      EQK = EG(3)*EG(4)/EG(2)/EG(9)
      RB(1) = RF(1) / MAX(EQK, SMALL)
C     R2: O + H2 = H + OH
      RF(2) = EXP(1.073400250738888D1 +2.7D0*ALOGT 
     * -3.150136339425897D3*TI)
      EQK = EG(2)*EG(4)/EG(3)/EG(6)
      RB(2) = RF(2) / MAX(EQK, SMALL)
C     R3: OH + H2 = H + H2O
      RF(3) = EXP(1.844439727056968D1 +1.6D0*ALOGT 
     * -1.659749470543394D3*TI)
      EQK = EG(2)*EG(7)/EG(4)/EG(6)
      RB(3) = RF(3) / MAX(EQK, SMALL)
C     R4: 2OH = O + H2O
      RF(4) = EXP(1.058986184880864D1 +2.4D0*ALOGT 
     * +1.061787168720231D3*TI)
      EQK = EG(3)*EG(7)/EG(4)/EG(4)
      RB(4) = RF(4) / MAX(EQK, SMALL)
C     R5: H + O2 = HO2
      RF(5) = EXP(2.926339399964515D1 +4.4D-1*ALOGT)
      EQK = EG(5)/EG(2)/EG(9)/PFAC1
      RB(5) = RF(5) / MAX(EQK, SMALL)
C     R6: H + HO2 = 2OH
      RF(6) = EXP(3.194650722679419D1 -1.484489169537763D2*TI)
      EQK = EG(4)*EG(4)/EG(2)/EG(5)
      RB(6) = RF(6) / MAX(EQK, SMALL)
C     R7: H2 + O2 = H + HO2
      RF(7) = EXP(1.329058600981878D1 +2.433D0*ALOGT 
     * -2.692309815207098D4*TI)
      EQK = EG(2)*EG(5)/EG(6)/EG(9)
      RB(7) = RF(7) / MAX(EQK, SMALL)
C     R8: OH + HO2 = H2O + O2
      RF(8) = EXP(3.09952086719568D1 +2.525694776551521D2*TI)
      EQK = EG(7)*EG(9)/EG(4)/EG(5)
      RB(8) = RF(8) / MAX(EQK, SMALL)
C     R9: H + HO2 = O + H2O
      RF(9) = EXP(2.900978721062765D1 -3.376583839863861D2*TI)
      EQK = EG(3)*EG(7)/EG(2)/EG(5)
      RB(9) = RF(9) / MAX(EQK, SMALL)
C     R10: O + HO2 = OH + O2
      RF(10) = 4.D13
      EQK = EG(4)*EG(9)/EG(3)/EG(5)
      RB(10) = RF(10) / MAX(EQK, SMALL)
C     R11: 2HO2 = H2O2 + O2
      RF(11) = EXP(2.559080028740199D1 +8.20243168253069D2*TI)
      EQK = EG(8)*EG(9)/EG(5)/EG(5)
      RB(11) = RF(11) / MAX(EQK, SMALL)
C     R12: 2HO2 = H2O2 + O2
      RF(12) = EXP(3.353310785188531D1 -6.038600011679037D3*TI)
      EQK = EG(8)*EG(9)/EG(5)/EG(5)
      RB(12) = RF(12) / MAX(EQK, SMALL)
C     R13: H + H2O2 = OH + H2O
      RF(13) = EXP(3.081323295642516D1 -1.997770170530481D3*TI)
      EQK = EG(4)*EG(7)/EG(2)/EG(8)
      RB(13) = RF(13) / MAX(EQK, SMALL)
C     R14: H + H2O2 = HO2 + H2
      RF(14) = EXP(1.561556883000703D1 +2.D0*ALOGT 
     * -2.616726671727583D3*TI)
      EQK = EG(5)*EG(6)/EG(2)/EG(8)
      RB(14) = RF(14) / MAX(EQK, SMALL)
C     R15: O + H2O2 = OH + HO2
      RF(15) = EXP(1.608039378377431D1 +2.D0*ALOGT 
     * -1.997770170530481D3*TI)
      EQK = EG(4)*EG(5)/EG(3)/EG(8)
      RB(15) = RF(15) / MAX(EQK, SMALL)
C     R16: OH + H2O2 = HO2 + H2O
      RF(16) = EXP(2.832416829648849D1 -2.148735170822457D2*TI)
      EQK = EG(5)*EG(7)/EG(4)/EG(8)
      RB(16) = RF(16) / MAX(EQK, SMALL)
C     R17: OH + H2O2 = HO2 + H2O
      RF(17) = EXP(9.538806728516803D1 -7.D0*ALOGT 
     * -1.892094670326098D4*TI)
      EQK = EG(5)*EG(7)/EG(4)/EG(8)
      RB(17) = RF(17) / MAX(EQK, SMALL)
C     R18: 2OH = H2O2
      RF(18) = EXP(3.234055131724088D1 -3.7D-1*ALOGT)
      EQK = EG(8)/EG(4)/EG(4)/PFAC1
      RB(18) = RF(18) / MAX(EQK, SMALL)
C     R19: 2H = H2
      RF(19) = EXP(4.202314503819682D1 -1.D0*ALOGT)
      EQK = EG(6)/EG(2)/EG(2)/PFAC1
      RB(19) = RF(19) / MAX(EQK, SMALL)
C     R20: H + OH = H2O
      RF(20) = EXP(5.213847658679322D1 -2.D0*ALOGT)
      EQK = EG(7)/EG(2)/EG(4)/PFAC1
      RB(20) = RF(20) / MAX(EQK, SMALL)
C     R21: 2O = O2
      RF(21) = EXP(3.932626813769273D1 -1.D0*ALOGT)
      EQK = EG(9)/EG(3)/EG(3)/PFAC1
      RB(21) = RF(21) / MAX(EQK, SMALL)
C     R22: 2H = H2
      RF(22) = EXP(3.903858606524095D1 -6.D-1*ALOGT)
      EQK = EG(6)/EG(2)/EG(2)/PFAC1
      RB(22) = RF(22) / MAX(EQK, SMALL)
C     R23: 2H = H2
      RF(23) = EXP(4.547615992139523D1 -1.25D0*ALOGT)
      EQK = EG(6)/EG(2)/EG(2)/PFAC1
      RB(23) = RF(23) / MAX(EQK, SMALL)
C     R24: 2H = H2
      RF(24) = EXP(4.775644995211934D1 -2.D0*ALOGT)
      EQK = EG(6)/EG(2)/EG(2)/PFAC1
      RB(24) = RF(24) / MAX(EQK, SMALL)
C     R25: H + O = OH
      RF(25) = EXP(4.369021565896671D1 -1.D0*ALOGT)
      EQK = EG(4)/EG(2)/EG(3)/PFAC1
      RB(25) = RF(25) / MAX(EQK, SMALL)
C     R26: OH + CO = H + CO2
      RF(26) = EXP(1.116280045189523D1 +2.053D0*ALOGT 
     * +1.789790721794902D2*TI)
      EQK = EG(2)*EG(18)/EG(4)/EG(17)
      RB(26) = RF(26) / MAX(EQK, SMALL)
C     R27: OH + CO = H + CO2
      RF(27) = EXP(2.938143762162222D1 -6.64D-1*ALOGT 
     * -1.669823868229545D2*TI)
      EQK = EG(2)*EG(18)/EG(4)/EG(17)
      RB(27) = RF(27) / MAX(EQK, SMALL)
C     R28: HO2 + CO = OH + CO2
      RF(28) = EXP(1.196400108433045D1 +2.18D0*ALOGT 
     * -9.0290204129627D3*TI)
      EQK = EG(4)*EG(18)/EG(5)/EG(17)
      RB(28) = RF(28) / MAX(EQK, SMALL)
C     R29: O + CO = CO2
      RF(29) = EXP(2.333480513766778D1 -1.199668535653569D3*TI)
      EQK = EG(18)/EG(3)/EG(17)/PFAC1
      RB(29) = RF(29) / MAX(EQK, SMALL)
C     R30: O2 + CO = O + CO2
      RF(30) = EXP(2.774345654525834D1 -2.400343504642417D4*TI)
      EQK = EG(3)*EG(18)/EG(9)/EG(17)
      RB(30) = RF(30) / MAX(EQK, SMALL)
C     R31: HCO = H + CO
      RF(31) = EXP(3.976988501176528D1 -1.D0*ALOGT 
     * -8.554683349878635D3*TI)
      EQK = EG(2)*EG(17)/EG(14)*PFAC1
      RB(31) = RF(31) / MAX(EQK, SMALL)
C     R32: H + HCO = H2 + CO
      RF(32) = 1.2D14
      EQK = EG(6)*EG(17)/EG(2)/EG(14)
      RB(32) = RF(32) / MAX(EQK, SMALL)
C     R33: O + HCO = OH + CO
      RF(33) = 3.D13
      EQK = EG(4)*EG(17)/EG(3)/EG(14)
      RB(33) = RF(33) / MAX(EQK, SMALL)
C     R34: O + HCO = H + CO2
      RF(34) = 3.D13
      EQK = EG(2)*EG(18)/EG(3)/EG(14)
      RB(34) = RF(34) / MAX(EQK, SMALL)
C     R35: OH + HCO = H2O + CO
      RF(35) = 3.02D13
      EQK = EG(7)*EG(17)/EG(4)/EG(14)
      RB(35) = RF(35) / MAX(EQK, SMALL)
C     R36: O2 + HCO = HO2 + CO
      RF(36) = EXP(2.321150027682709D1 +8.070000000000001D-1*ALOGT 
     * +3.658385173742216D2*TI)
      EQK = EG(5)*EG(17)/EG(9)/EG(14)
      RB(36) = RF(36) / MAX(EQK, SMALL)
C     R37: HCO = H + CO
      RF(37) = EXP(4.225479166155327D1 -1.D0*ALOGT 
     * -8.554683349878635D3*TI)
      EQK = EG(2)*EG(17)/EG(14)*PFAC1
      RB(37) = RF(37) / MAX(EQK, SMALL)
C     R38: H2 + CO = CH2O
      RF(38) = EXP(1.757671067365784D1 +1.5D0*ALOGT 
     * -4.005604674413761D4*TI)
      EQK = EG(15)/EG(6)/EG(17)/PFAC1
      RB(38) = RF(38) / MAX(EQK, SMALL)
C     R39: H + HCO = CH2O
      RF(39) = EXP(2.77171988121696D1 +4.8D-1*ALOGT 
     * +1.308363335863791D2*TI)
      EQK = EG(15)/EG(2)/EG(14)/PFAC1
      RB(39) = RF(39) / MAX(EQK, SMALL)
C     R40: H + CH2 = CH3
      RF(40) = EXP(3.775765221977888D1 -8.D-1*ALOGT)
      EQK = EG(12)/EG(2)/EG(10)/PFAC1
      RB(40) = RF(40) / MAX(EQK, SMALL)
C     R41: O + CH2 = H + HCO
      RF(41) = 8.D13
      EQK = EG(2)*EG(14)/EG(3)/EG(10)
      RB(41) = RF(41) / MAX(EQK, SMALL)
C     R42: OH + CH2 = H + CH2O
      RF(42) = 2.D13
      EQK = EG(2)*EG(15)/EG(4)/EG(10)
      RB(42) = RF(42) / MAX(EQK, SMALL)
C     R43: H2 + CH2 = H + CH3
      RF(43) = EXP(1.312236337740433D1 +2.D0*ALOGT 
     * -3.63825650703662D3*TI)
      EQK = EG(2)*EG(12)/EG(6)/EG(10)
      RB(43) = RF(43) / MAX(EQK, SMALL)
C     R44: O2 + CH2 = OH + HCO
      RF(44) = EXP(2.999187511704657D1 -7.548250014598796D2*TI)
      EQK = EG(4)*EG(14)/EG(9)/EG(10)
      RB(44) = RF(44) / MAX(EQK, SMALL)
C     R45: O2 + CH2 = 2H + CO2
      RF(45) = EXP(2.860180003308677D1 -7.548250014598796D2*TI)
      EQK = EG(2)*EG(2)*EG(18)/EG(9)/EG(10)*PFAC1
      RB(45) = RF(45) / MAX(EQK, SMALL)
C     R46: HO2 + CH2 = OH + CH2O
      RF(46) = 2.D13
      EQK = EG(4)*EG(15)/EG(5)/EG(10)
      RB(46) = RF(46) / MAX(EQK, SMALL)
C     R47: 2CH2 = H2 + C2H2
      RF(47) = 3.2D13
      EQK = EG(6)*EG(19)/EG(10)/EG(10)
      RB(47) = RF(47) / MAX(EQK, SMALL)
C     R48: CH2* = CH2
      RF(48) = EXP(3.033907131703076D1 -3.019300005839518D2*TI)
      EQK = EG(10)/EG(11)
      RB(48) = RF(48) / MAX(EQK, SMALL)
C     R49: O + CH2* = H2 + CO
      RF(49) = 1.5D13
      EQK = EG(6)*EG(17)/EG(3)/EG(11)
      RB(49) = RF(49) / MAX(EQK, SMALL)
C     R50: O + CH2* = H + HCO
      RF(50) = 1.5D13
      EQK = EG(2)*EG(14)/EG(3)/EG(11)
      RB(50) = RF(50) / MAX(EQK, SMALL)
C     R51: OH + CH2* = H + CH2O
      RF(51) = 3.D13
      EQK = EG(2)*EG(15)/EG(4)/EG(11)
      RB(51) = RF(51) / MAX(EQK, SMALL)
C     R52: H2 + CH2* = H + CH3
      RF(52) = 7.D13
      EQK = EG(2)*EG(12)/EG(6)/EG(11)
      RB(52) = RF(52) / MAX(EQK, SMALL)
C     R53: O2 + CH2* = H + OH + CO
      RF(53) = 2.8D13
      EQK = EG(2)*EG(4)*EG(17)/EG(9)/EG(11)*PFAC1
      RB(53) = RF(53) / MAX(EQK, SMALL)
C     R54: O2 + CH2* = H2O + CO
      RF(54) = 1.2D13
      EQK = EG(7)*EG(17)/EG(9)/EG(11)
      RB(54) = RF(54) / MAX(EQK, SMALL)
C     R55: CH2* = CH2
      RF(55) = 3.D13
      EQK = EG(10)/EG(11)
      RB(55) = RF(55) / MAX(EQK, SMALL)
C     R56: CH2* = CH2
      RF(56) = 9.D12
      EQK = EG(10)/EG(11)
      RB(56) = RF(56) / MAX(EQK, SMALL)
C     R57: CH2* = CH2
      RF(57) = 7.D12
      EQK = EG(10)/EG(11)
      RB(57) = RF(57) / MAX(EQK, SMALL)
C     R58: CH2* + CO2 = CH2O + CO
      RF(58) = 1.4D13
      EQK = EG(15)*EG(17)/EG(11)/EG(18)
      RB(58) = RF(58) / MAX(EQK, SMALL)
C     R59: H + CH2O = CH3O
      RF(59) = EXP(2.701483497650473D1 +4.54D-1*ALOGT 
     * -1.308363335863791D3*TI)
      EQK = EG(16)/EG(2)/EG(15)/PFAC1
      RB(59) = RF(59) / MAX(EQK, SMALL)
C     R60: H + CH2O = H2 + HCO
      RF(60) = EXP(2.385876005287556D1 +1.05D0*ALOGT 
     * -1.648034586520737D3*TI)
      EQK = EG(6)*EG(14)/EG(2)/EG(15)
      RB(60) = RF(60) / MAX(EQK, SMALL)
C     R61: O + CH2O = OH + HCO
      RF(61) = EXP(3.129458276205819D1 -1.781387003445316D3*TI)
      EQK = EG(4)*EG(14)/EG(3)/EG(15)
      RB(61) = RF(61) / MAX(EQK, SMALL)
C     R62: OH + CH2O = H2O + HCO
      RF(62) = EXP(2.195582609812426D1 +1.18D0*ALOGT 
     * +2.249378504350441D2*TI)
      EQK = EG(7)*EG(14)/EG(4)/EG(15)
      RB(62) = RF(62) / MAX(EQK, SMALL)
C     R63: O2 + CH2O = HO2 + HCO
      RF(63) = EXP(3.223619130191664D1 -2.012866670559679D4*TI)
      EQK = EG(5)*EG(14)/EG(9)/EG(15)
      RB(63) = RF(63) / MAX(EQK, SMALL)
C     R64: HO2 + CH2O = H2O2 + HCO
      RF(64) = EXP(2.763102111592855D1 -4.025733341119358D3*TI)
      EQK = EG(8)*EG(14)/EG(5)/EG(15)
      RB(64) = RF(64) / MAX(EQK, SMALL)
C     R65: H + CH3 = CH4
      RF(65) = EXP(3.708037838837523D1 -6.3D-1*ALOGT 
     * -1.927319837060892D2*TI)
      EQK = EG(13)/EG(2)/EG(12)/PFAC1
      RB(65) = RF(65) / MAX(EQK, SMALL)
C     R66: O + CH3 = H + CH2O
      RF(66) = 8.43D13
      EQK = EG(2)*EG(15)/EG(3)/EG(12)
      RB(66) = RF(66) / MAX(EQK, SMALL)
C     R67: OH + CH3 = H2O + CH2
      RF(67) = EXP(1.784086224869942D1 +1.6D0*ALOGT 
     * -2.727434338608365D3*TI)
      EQK = EG(7)*EG(10)/EG(4)/EG(12)
      RB(67) = RF(67) / MAX(EQK, SMALL)
C     R68: OH + CH3 = H2O + CH2*
      RF(68) = 2.501D13
      EQK = EG(7)*EG(11)/EG(4)/EG(12)
      RB(68) = RF(68) / MAX(EQK, SMALL)
C     R69: O2 + CH3 = O + CH3O
      RF(69) = EXP(3.105950935782661D1 -1.449264002802969D4*TI)
      EQK = EG(3)*EG(16)/EG(9)/EG(12)
      RB(69) = RF(69) / MAX(EQK, SMALL)
C     R70: O2 + CH3 = OH + CH2O
      RF(70) = EXP(2.430678477540252D1 -4.498757008700882D3*TI)
      EQK = EG(4)*EG(15)/EG(9)/EG(12)
      RB(70) = RF(70) / MAX(EQK, SMALL)
C     R71: HO2 + CH3 = O2 + CH4
      RF(71) = 1.D12
      EQK = EG(9)*EG(13)/EG(5)/EG(12)
      RB(71) = RF(71) / MAX(EQK, SMALL)
C     R72: HO2 + CH3 = OH + CH3O
      RF(72) = 1.34D13
      EQK = EG(4)*EG(16)/EG(5)/EG(12)
      RB(72) = RF(72) / MAX(EQK, SMALL)
C     R73: H2O2 + CH3 = HO2 + CH4
      RF(73) = EXP(1.010642839653282D1 +2.47D0*ALOGT 
     * -2.606662338374784D3*TI)
      EQK = EG(5)*EG(13)/EG(8)/EG(12)
      RB(73) = RF(73) / MAX(EQK, SMALL)
C     R74: CH3 + HCO = CH4 + CO
      RF(74) = 8.48D12
      EQK = EG(13)*EG(17)/EG(12)/EG(14)
      RB(74) = RF(74) / MAX(EQK, SMALL)
C     R75: CH3 + CH2O = CH4 + HCO
      RF(75) = EXP(8.107720061910534D0 +2.81D0*ALOGT 
     * -2.94884967236993D3*TI)
      EQK = EG(13)*EG(14)/EG(12)/EG(15)
      RB(75) = RF(75) / MAX(EQK, SMALL)
C     R76: CH2 + CH3 = H + C2H4
      RF(76) = 4.D13
      EQK = EG(2)*EG(21)/EG(10)/EG(12)
      RB(76) = RF(76) / MAX(EQK, SMALL)
C     R77: CH2* + CH3 = H + C2H4
      RF(77) = EXP(3.011592776571655D1 +2.868335005547542D2*TI)
      EQK = EG(2)*EG(21)/EG(11)/EG(12)
      RB(77) = RF(77) / MAX(EQK, SMALL)
C     R78: 2CH3 = C2H6
      RF(78) = EXP(3.759277757658865D1 -9.7D-1*ALOGT 
     * -3.119943339367502D2*TI)
      EQK = EG(23)/EG(12)/EG(12)/PFAC1
      RB(78) = RF(78) / MAX(EQK, SMALL)
C     R79: 2CH3 = H + C2H5
      RF(79) = EXP(2.923845702569198D1 +1.D-1*ALOGT 
     * -5.334096676983148D3*TI)
      EQK = EG(2)*EG(22)/EG(12)/EG(12)
      RB(79) = RF(79) / MAX(EQK, SMALL)
C     R80: H + CH3O = H2 + CH2O
      RF(80) = 2.D13
      EQK = EG(6)*EG(15)/EG(2)/EG(16)
      RB(80) = RF(80) / MAX(EQK, SMALL)
C     R81: H + CH3O = OH + CH3
      RF(81) = 3.2D13
      EQK = EG(4)*EG(12)/EG(2)/EG(16)
      RB(81) = RF(81) / MAX(EQK, SMALL)
C     R82: H + CH3O = H2O + CH2*
      RF(82) = 1.6D13
      EQK = EG(7)*EG(11)/EG(2)/EG(16)
      RB(82) = RF(82) / MAX(EQK, SMALL)
C     R83: O + CH3O = OH + CH2O
      RF(83) = 1.D13
      EQK = EG(4)*EG(15)/EG(3)/EG(16)
      RB(83) = RF(83) / MAX(EQK, SMALL)
C     R84: OH + CH3O = H2O + CH2O
      RF(84) = 5.D12
      EQK = EG(7)*EG(15)/EG(4)/EG(16)
      RB(84) = RF(84) / MAX(EQK, SMALL)
C     R85: O2 + CH3O = HO2 + CH2O
      RF(85) = EXP(-2.847965319932889D1 +7.6D0*ALOGT 
     * +1.776354836768916D3*TI)
      EQK = EG(5)*EG(15)/EG(9)/EG(16)
      RB(85) = RF(85) / MAX(EQK, SMALL)
C     R86: H + CH4 = H2 + CH3
      RF(86) = EXP(2.030775039298474D1 +1.62D0*ALOGT 
     * -5.454868677216729D3*TI)
      EQK = EG(6)*EG(12)/EG(2)/EG(13)
      RB(86) = RF(86) / MAX(EQK, SMALL)
C     R87: O + CH4 = OH + CH3
      RF(87) = EXP(2.074306846424259D1 +1.5D0*ALOGT 
     * -4.327663341703309D3*TI)
      EQK = EG(4)*EG(12)/EG(3)/EG(13)
      RB(87) = RF(87) / MAX(EQK, SMALL)
C     R88: OH + CH4 = H2O + CH3
      RF(88) = EXP(1.842068074395237D1 +1.6D0*ALOGT 
     * -1.570036003036549D3*TI)
      EQK = EG(7)*EG(12)/EG(4)/EG(13)
      RB(88) = RF(88) / MAX(EQK, SMALL)
C     R89: CH2 + CH4 = 2CH3
      RF(89) = EXP(1.471567190790855D1 +2.D0*ALOGT 
     * -4.161601841382136D3*TI)
      EQK = EG(12)*EG(12)/EG(10)/EG(13)
      RB(89) = RF(89) / MAX(EQK, SMALL)
C     R90: CH2* + CH4 = 2CH3
      RF(90) = EXP(3.040360983816833D1 +2.868335005547542D2*TI)
      EQK = EG(12)*EG(12)/EG(11)/EG(13)
      RB(90) = RF(90) / MAX(EQK, SMALL)
C     R91: C2H3 = H + C2H2
      RF(91) = EXP(1.746876283443506D1 +1.62D0*ALOGT 
     * -1.864327174605727D4*TI)
      EQK = EG(2)*EG(19)/EG(20)*PFAC1
      RB(91) = RF(91) / MAX(EQK, SMALL)
C     R92: O + C2H2 = CH2 + CO
      RF(92) = EXP(1.522160754638034D1 +2.D0*ALOGT 
     * -9.561116685158474D2*TI)
      EQK = EG(10)*EG(17)/EG(3)/EG(19)
      RB(92) = RF(92) / MAX(EQK, SMALL)
C     R93: OH + C2H2 = CH3 + CO
      RF(93) = EXP(-7.635493904311701D0 +4.D0*ALOGT 
     * +1.006433335279839D3*TI)
      EQK = EG(12)*EG(17)/EG(4)/EG(19)
      RB(93) = RF(93) / MAX(EQK, SMALL)
C     R94: HCO + C2H2 = CO + C2H3
      RF(94) = EXP(1.611809565095832D1 +2.D0*ALOGT 
     * -3.019300005839518D3*TI)
      EQK = EG(17)*EG(20)/EG(14)/EG(19)
      RB(94) = RF(94) / MAX(EQK, SMALL)
C     R95: CH3 + C2H2 = aC3H5
      RF(95) = EXP(1.230228267232072D2 -1.282D1*ALOGT 
     * -1.797993153477433D4*TI)
      EQK = EG(25)/EG(12)/EG(19)/PFAC1
      RB(95) = RF(95) / MAX(EQK, SMALL)
C     R96: H + C2H3 = C2H4
      RF(96) = EXP(2.943602581190662D1 +2.7D-1*ALOGT 
     * -1.409006669391775D2*TI)
      EQK = EG(21)/EG(2)/EG(20)/PFAC1
      RB(96) = RF(96) / MAX(EQK, SMALL)
C     R97: H + C2H3 = H2 + C2H2
      RF(97) = 9.D13
      EQK = EG(6)*EG(19)/EG(2)/EG(20)
      RB(97) = RF(97) / MAX(EQK, SMALL)
C     R98: O + C2H3 = CH3 + CO
      RF(98) = 4.8D13
      EQK = EG(12)*EG(17)/EG(3)/EG(20)
      RB(98) = RF(98) / MAX(EQK, SMALL)
C     R99: OH + C2H3 = H2O + C2H2
      RF(99) = 3.011D13
      EQK = EG(7)*EG(19)/EG(4)/EG(20)
      RB(99) = RF(99) / MAX(EQK, SMALL)
C     R100: O2 + C2H3 = HO2 + C2H2
      RF(100) = EXP(1.410818017192709D1 +1.61D0*ALOGT 
     * +1.929332703731452D2*TI)
      EQK = EG(5)*EG(19)/EG(9)/EG(20)
      RB(100) = RF(100) / MAX(EQK, SMALL)
C     R101: O2 + C2H3 = O + CH2CHO
      RF(101) = EXP(2.642704831160261D1 +2.9D-1*ALOGT 
     * -5.535383344039117D0*TI)
      EQK = EG(3)*EG(24)/EG(9)/EG(20)
      RB(101) = RF(101) / MAX(EQK, SMALL)
C     R102: O2 + C2H3 = HCO + CH2O
      RF(102) = EXP(3.836741779139978D1 -1.39D0*ALOGT 
     * -5.082488343163189D2*TI)
      EQK = EG(14)*EG(15)/EG(9)/EG(20)
      RB(102) = RF(102) / MAX(EQK, SMALL)
C     R103: HO2 + C2H3 = OH + CH2CHO
      RF(103) = 1.D13
      EQK = EG(4)*EG(24)/EG(5)/EG(20)
      RB(103) = RF(103) / MAX(EQK, SMALL)
C     R104: H2O2 + C2H3 = HO2 + C2H4
      RF(104) = EXP(2.321647128954911D1 +2.999171339133922D2*TI)
      EQK = EG(5)*EG(21)/EG(8)/EG(20)
      RB(104) = RF(104) / MAX(EQK, SMALL)
C     R105: HCO + C2H3 = CO + C2H4
      RF(105) = 9.033D13
      EQK = EG(17)*EG(21)/EG(14)/EG(20)
      RB(105) = RF(105) / MAX(EQK, SMALL)
C     R106: HCO + C2H3 = C2H3CHO
      RF(106) = 1.8D13
      EQK = EG(28)/EG(14)/EG(20)/PFAC1
      RB(106) = RF(106) / MAX(EQK, SMALL)
C     R107: CH3 + C2H3 = CH4 + C2H2
      RF(107) = 3.92D11
      EQK = EG(13)*EG(19)/EG(12)/EG(20)
      RB(107) = RF(107) / MAX(EQK, SMALL)
C     R108: CH3 + C2H3 = C3H6
      RF(108) = 2.5D13
      EQK = EG(26)/EG(12)/EG(20)/PFAC1
      RB(108) = RF(108) / MAX(EQK, SMALL)
C     R109: CH3 + C2H3 = H + aC3H5
      RF(109) = EXP(5.566750733996526D1 -2.83D0*ALOGT 
     * -9.368887918120025D3*TI)
      EQK = EG(2)*EG(25)/EG(12)/EG(20)
      RB(109) = RF(109) / MAX(EQK, SMALL)
C     R110: 2C2H3 = C2H2 + C2H4
      RF(110) = 9.6D11
      EQK = EG(19)*EG(21)/EG(20)/EG(20)
      RB(110) = RF(110) / MAX(EQK, SMALL)
C     R111: CH2CHO = CH3 + CO
      RF(111) = EXP(9.646011254645141D1 -9.147D0*ALOGT 
     * -2.360086171231223D4*TI)
      EQK = EG(12)*EG(17)/EG(24)*PFAC1
      RB(111) = RF(111) / MAX(EQK, SMALL)
C     R112: H + CH2CHO = CH3 + HCO
      RF(112) = 9.D13
      EQK = EG(12)*EG(14)/EG(2)/EG(24)
      RB(112) = RF(112) / MAX(EQK, SMALL)
C     R113: O2 + CH2CHO = OH + CH2O + CO
      RF(113) = 1.8D10
      EQK = EG(4)*EG(15)*EG(17)/EG(9)/EG(24)*PFAC1
      RB(113) = RF(113) / MAX(EQK, SMALL)
C     R114: H + C2H4 = C2H5
      RF(114) = EXP(2.210329058505271D1 +1.28D0*ALOGT 
     * -6.494665277561096D2*TI)
      EQK = EG(22)/EG(2)/EG(21)/PFAC1
      RB(114) = RF(114) / MAX(EQK, SMALL)
C     R115: H + C2H4 = H2 + C2H3
      RF(115) = EXP(1.774143646856141D1 +1.9D0*ALOGT 
     * -6.51665584593696D3*TI)
      EQK = EG(6)*EG(20)/EG(2)/EG(21)
      RB(115) = RF(115) / MAX(EQK, SMALL)
C     R116: O + C2H4 = OH + C2H3
      RF(116) = EXP(1.653020530178515D1 +1.9D0*ALOGT 
     * -1.8820303369733D3*TI)
      EQK = EG(4)*EG(20)/EG(3)/EG(21)
      RB(116) = RF(116) / MAX(EQK, SMALL)
C     R117: O + C2H4 = CH3 + HCO
      RF(117) = EXP(1.677042083699801D1 +1.83D0*ALOGT 
     * -1.107076668807823D2*TI)
      EQK = EG(12)*EG(14)/EG(3)/EG(21)
      RB(117) = RF(117) / MAX(EQK, SMALL)
C     R118: O + C2H4 = CH2 + CH2O
      RF(118) = EXP(1.285839783156986D1 +1.83D0*ALOGT 
     * -1.107076668807823D2*TI)
      EQK = EG(10)*EG(15)/EG(3)/EG(21)
      RB(118) = RF(118) / MAX(EQK, SMALL)
C     R119: OH + C2H4 = H2O + C2H3
      RF(119) = EXP(1.509644440342634D1 +2.D0*ALOGT 
     * -1.258041669099799D3*TI)
      EQK = EG(7)*EG(20)/EG(4)/EG(21)
      RB(119) = RF(119) / MAX(EQK, SMALL)
C     R120: HCO + C2H4 = CO + C2H5
      RF(120) = EXP(1.611809565095832D1 +2.D0*ALOGT 
     * -4.025733341119358D3*TI)
      EQK = EG(17)*EG(22)/EG(14)/EG(21)
      RB(120) = RF(120) / MAX(EQK, SMALL)
C     R121: CH2 + C2H4 = H + aC3H5
      RF(121) = EXP(3.062675338948254D1 -3.019300005839518D3*TI)
      EQK = EG(2)*EG(25)/EG(10)/EG(21)
      RB(121) = RF(121) / MAX(EQK, SMALL)
C     R122: CH2* + C2H4 = H + aC3H5
      RF(122) = 5.D13
      EQK = EG(2)*EG(25)/EG(11)/EG(21)
      RB(122) = RF(122) / MAX(EQK, SMALL)
C     R123: CH3 + C2H4 = CH4 + C2H3
      RF(123) = EXP(1.233270529646354D1 +2.D0*ALOGT 
     * -4.629593342287261D3*TI)
      EQK = EG(13)*EG(20)/EG(12)/EG(21)
      RB(123) = RF(123) / MAX(EQK, SMALL)
C     R124: nC3H7 = CH3 + C2H4
      RF(124) = EXP(3.219536930739638D1 -1.561126048719317D4*TI)
      EQK = EG(12)*EG(21)/EG(27)*PFAC1
      RB(124) = RF(124) / MAX(EQK, SMALL)
C     R125: O2 + C2H4 = HO2 + C2H3
      RF(125) = EXP(3.137344133697052D1 -3.059557339250712D4*TI)
      EQK = EG(5)*EG(20)/EG(9)/EG(21)
      RB(125) = RF(125) / MAX(EQK, SMALL)
C     R126: C2H3 + C2H4 = C4H7
      RF(126) = EXP(8.95688865694205D1 -8.470000000000001D0*ALOGT 
     * -7.155741013839658D3*TI)
      EQK = EG(29)/EG(20)/EG(21)/PFAC1
      RB(126) = RF(126) / MAX(EQK, SMALL)
C     R127: H + C2H5 = C2H6
      RF(127) = EXP(4.079452643666405D1 -9.9D-1*ALOGT 
     * -7.95082334871073D2*TI)
      EQK = EG(23)/EG(2)/EG(22)/PFAC1
      RB(127) = RF(127) / MAX(EQK, SMALL)
C     R128: H + C2H5 = H2 + C2H4
      RF(128) = 2.D12
      EQK = EG(6)*EG(21)/EG(2)/EG(22)
      RB(128) = RF(128) / MAX(EQK, SMALL)
C     R129: O + C2H5 = CH3 + CH2O
      RF(129) = 1.604D13
      EQK = EG(12)*EG(15)/EG(3)/EG(22)
      RB(129) = RF(129) / MAX(EQK, SMALL)
C     R130: O2 + C2H5 = HO2 + C2H4
      RF(130) = 2.D10
      EQK = EG(5)*EG(21)/EG(9)/EG(22)
      RB(130) = RF(130) / MAX(EQK, SMALL)
C     R131: HO2 + C2H5 = O2 + C2H6
      RF(131) = 3.D11
      EQK = EG(9)*EG(23)/EG(5)/EG(22)
      RB(131) = RF(131) / MAX(EQK, SMALL)
C     R132: HO2 + C2H5 = H2O2 + C2H4
      RF(132) = 3.D11
      EQK = EG(8)*EG(21)/EG(5)/EG(22)
      RB(132) = RF(132) / MAX(EQK, SMALL)
C     R133: HO2 + C2H5 = OH + CH3 + CH2O
      RF(133) = 2.4D13
      EQK = EG(4)*EG(12)*EG(15)/EG(5)/EG(22)*PFAC1
      RB(133) = RF(133) / MAX(EQK, SMALL)
C     R134: H2O2 + C2H5 = HO2 + C2H6
      RF(134) = EXP(2.288658886260695D1 -4.901330342812818D2*TI)
      EQK = EG(5)*EG(23)/EG(8)/EG(22)
      RB(134) = RF(134) / MAX(EQK, SMALL)
C     R135: C2H3 + C2H5 = C4H81
      RF(135) = 1.5D13
      EQK = EG(30)/EG(20)/EG(22)/PFAC1
      RB(135) = RF(135) / MAX(EQK, SMALL)
C     R136: C2H3 + C2H5 = CH3 + aC3H5
      RF(136) = EXP(7.504369952894507D1 -5.22D0*ALOGT 
     * -9.937019535885494D3*TI)
      EQK = EG(12)*EG(25)/EG(20)/EG(22)
      RB(136) = RF(136) / MAX(EQK, SMALL)
C     R137: H + C2H6 = H2 + C2H5
      RF(137) = EXP(1.856044268632752D1 +1.9D0*ALOGT 
     * -3.789221507328595D3*TI)
      EQK = EG(6)*EG(22)/EG(2)/EG(23)
      RB(137) = RF(137) / MAX(EQK, SMALL)
C     R138: O + C2H6 = OH + C2H5
      RF(138) = EXP(1.831309553327243D1 +1.92D0*ALOGT 
     * -2.863302838871143D3*TI)
      EQK = EG(4)*EG(22)/EG(3)/EG(23)
      RB(138) = RF(138) / MAX(EQK, SMALL)
C     R139: OH + C2H6 = H2O + C2H5
      RF(139) = EXP(1.507963728510996D1 +2.12D0*ALOGT 
     * -4.377985008467301D2*TI)
      EQK = EG(7)*EG(22)/EG(4)/EG(23)
      RB(139) = RF(139) / MAX(EQK, SMALL)
C     R140: CH2* + C2H6 = CH3 + C2H5
      RF(140) = EXP(3.131990057004249D1 +2.767691672019558D2*TI)
      EQK = EG(12)*EG(22)/EG(11)/EG(23)
      RB(140) = RF(140) / MAX(EQK, SMALL)
C     R141: CH3 + C2H6 = CH4 + C2H5
      RF(141) = EXP(1.563033530012333D1 +1.74D0*ALOGT 
     * -5.258614176837161D3*TI)
      EQK = EG(13)*EG(22)/EG(12)/EG(23)
      RB(141) = RF(141) / MAX(EQK, SMALL)
C     R142: H + aC3H5 = C3H6
      RF(142) = 2.D14
      EQK = EG(26)/EG(2)/EG(25)/PFAC1
      RB(142) = RF(142) / MAX(EQK, SMALL)
C     R143: O + aC3H5 = H + C2H3CHO
      RF(143) = 6.D13
      EQK = EG(2)*EG(28)/EG(3)/EG(25)
      RB(143) = RF(143) / MAX(EQK, SMALL)
C     R144: OH + aC3H5 = 2H + C2H3CHO
      RF(144) = EXP(7.511780750109878D1 -5.16D0*ALOGT 
     * -1.515990532932022D4*TI)
      EQK = EG(2)*EG(2)*EG(28)/EG(4)/EG(25)*PFAC1
      RB(144) = RF(144) / MAX(EQK, SMALL)
C     R145: O2 + aC3H5 = OH + C2H3CHO
      RF(145) = EXP(3.05324427100113D1 -4.1D-1*ALOGT 
     * -1.150302980558092D4*TI)
      EQK = EG(4)*EG(28)/EG(9)/EG(25)
      RB(145) = RF(145) / MAX(EQK, SMALL)
C     R146: HO2 + aC3H5 = O2 + C3H6
      RF(146) = 2.66D12
      EQK = EG(9)*EG(26)/EG(5)/EG(25)
      RB(146) = RF(146) / MAX(EQK, SMALL)
C     R147: HO2 + aC3H5 = OH + CH2O + C2H3
      RF(147) = 6.6D12
      EQK = EG(4)*EG(15)*EG(20)/EG(5)/EG(25)*PFAC1
      RB(147) = RF(147) / MAX(EQK, SMALL)
C     R148: HCO + aC3H5 = CO + C3H6
      RF(148) = 6.D13
      EQK = EG(17)*EG(26)/EG(14)/EG(25)
      RB(148) = RF(148) / MAX(EQK, SMALL)
C     R149: CH3 + aC3H5 = C4H81
      RF(149) = EXP(3.223619130191664D1 -3.2D-1*ALOGT 
     * +1.319937319219509D2*TI)
      EQK = EG(30)/EG(12)/EG(25)/PFAC1
      RB(149) = RF(149) / MAX(EQK, SMALL)
C     R150: H + C3H6 = nC3H7
      RF(150) = EXP(3.021878515115626D1 -1.640838588173486D3*TI)
      EQK = EG(27)/EG(2)/EG(26)/PFAC1
      RB(150) = RF(150) / MAX(EQK, SMALL)
C     R151: H + C3H6 = CH3 + C2H4
      RF(151) = EXP(5.043372849455479D1 -2.39D0*ALOGT 
     * -5.625962344214302D3*TI)
      EQK = EG(12)*EG(21)/EG(2)/EG(26)
      RB(151) = RF(151) / MAX(EQK, SMALL)
C     R152: H + C3H6 = H2 + aC3H5
      RF(152) = EXP(1.206104687347992D1 +2.5D0*ALOGT 
     * -1.2530095024234D3*TI)
      EQK = EG(6)*EG(25)/EG(2)/EG(26)
      RB(152) = RF(152) / MAX(EQK, SMALL)
C     R153: O + C3H6 = 2H + C2H3CHO
      RF(153) = EXP(1.750439001207821D1 +1.65D0*ALOGT 
     * -1.645518503182537D2*TI)
      EQK = EG(2)*EG(2)*EG(28)/EG(3)/EG(26)*PFAC1
      RB(153) = RF(153) / MAX(EQK, SMALL)
C     R154: O + C3H6 = HCO + C2H5
      RF(154) = EXP(1.737085861945369D1 +1.65D0*ALOGT 
     * +4.891266009460019D2*TI)
      EQK = EG(14)*EG(22)/EG(3)/EG(26)
      RB(154) = RF(154) / MAX(EQK, SMALL)
C     R155: O + C3H6 = OH + aC3H5
      RF(155) = EXP(2.591622268783662D1 +7.D-1*ALOGT 
     * -2.958914005722728D3*TI)
      EQK = EG(4)*EG(25)/EG(3)/EG(26)
      RB(155) = RF(155) / MAX(EQK, SMALL)
C     R156: OH + C3H6 = H2O + aC3H5
      RF(156) = EXP(1.494691266945537D1 +2.D0*ALOGT 
     * +1.499585669566961D2*TI)
      EQK = EG(7)*EG(25)/EG(4)/EG(26)
      RB(156) = RF(156) / MAX(EQK, SMALL)
C     R157: HO2 + C3H6 = H2O2 + aC3H5
      RF(157) = EXP(9.169518377455928D0 +2.6D0*ALOGT 
     * -6.999743846871283D3*TI)
      EQK = EG(8)*EG(25)/EG(5)/EG(26)
      RB(157) = RF(157) / MAX(EQK, SMALL)
C     R158: CH3 + C3H6 = CH4 + aC3H5
      RF(158) = EXP(7.884573603642703D-1 +3.5D0*ALOGT 
     * -2.855754588856544D3*TI)
      EQK = EG(13)*EG(25)/EG(12)/EG(26)
      RB(158) = RF(158) / MAX(EQK, SMALL)
C     R159: H + C2H3CHO = HCO + C2H4
      RF(159) = EXP(2.540539706407063D1 +4.54D-1*ALOGT 
     * -2.928721005664333D3*TI)
      EQK = EG(14)*EG(21)/EG(2)/EG(28)
      RB(159) = RF(159) / MAX(EQK, SMALL)
C     R160: O + C2H3CHO = OH + CO + C2H3
      RF(160) = EXP(3.103221849759071D1 -1.781387003445316D3*TI)
      EQK = EG(4)*EG(17)*EG(20)/EG(3)/EG(28)*PFAC1
      RB(160) = RF(160) / MAX(EQK, SMALL)
C     R161: OH + C2H3CHO = H2O + CO + C2H3
      RF(161) = EXP(2.195582609812426D1 +1.18D0*ALOGT 
     * +2.249378504350441D2*TI)
      EQK = EG(7)*EG(17)*EG(20)/EG(4)/EG(28)*PFAC1
      RB(161) = RF(161) / MAX(EQK, SMALL)
C     R162: H + nC3H7 = CH3 + C2H5
      RF(162) = EXP(5.657037505150728D1 -2.92D0*ALOGT 
     * -6.292724428837196D3*TI)
      EQK = EG(12)*EG(22)/EG(2)/EG(27)
      RB(162) = RF(162) / MAX(EQK, SMALL)
C     R163: H + nC3H7 = H2 + C3H6
      RF(163) = 1.8D12
      EQK = EG(6)*EG(26)/EG(2)/EG(27)
      RB(163) = RF(163) / MAX(EQK, SMALL)
C     R164: O + nC3H7 = CH2O + C2H5
      RF(164) = 9.6D13
      EQK = EG(15)*EG(22)/EG(3)/EG(27)
      RB(164) = RF(164) / MAX(EQK, SMALL)
C     R165: OH + nC3H7 = H2O + C3H6
      RF(165) = 2.4D13
      EQK = EG(7)*EG(26)/EG(4)/EG(27)
      RB(165) = RF(165) / MAX(EQK, SMALL)
C     R166: O2 + nC3H7 = HO2 + C3H6
      RF(166) = 9.D10
      EQK = EG(5)*EG(26)/EG(9)/EG(27)
      RB(166) = RF(166) / MAX(EQK, SMALL)
C     R167: HO2 + nC3H7 = OH + CH2O + C2H5
      RF(167) = 2.4D13
      EQK = EG(4)*EG(15)*EG(22)/EG(5)/EG(27)*PFAC1
      RB(167) = RF(167) / MAX(EQK, SMALL)
C     R168: CH3 + nC3H7 = CH4 + C3H6
      RF(168) = 1.1D13
      EQK = EG(13)*EG(26)/EG(12)/EG(27)
      RB(168) = RF(168) / MAX(EQK, SMALL)
C     R169: H + C4H7 = C4H81
      RF(169) = 3.6D13
      EQK = EG(30)/EG(2)/EG(29)/PFAC1
      RB(169) = RF(169) / MAX(EQK, SMALL)
C     R170: H + C4H7 = CH3 + aC3H5
      RF(170) = EXP(4.904743413343491D1 -2.D0*ALOGT 
     * -5.535383344039116D3*TI)
      EQK = EG(12)*EG(25)/EG(2)/EG(29)
      RB(170) = RF(170) / MAX(EQK, SMALL)
C     R171: HO2 + C4H7 = OH + CH2O + aC3H5
      RF(171) = 2.4D13
      EQK = EG(4)*EG(15)*EG(25)/EG(5)/EG(29)*PFAC1
      RB(171) = RF(171) / MAX(EQK, SMALL)
C     R172: HCO + C4H7 = CO + C4H81
      RF(172) = 6.D13
      EQK = EG(17)*EG(30)/EG(14)/EG(29)
      RB(172) = RF(172) / MAX(EQK, SMALL)
C     R173: H + C4H81 = pC4H9
      RF(173) = EXP(3.021878515115626D1 -1.640838588173486D3*TI)
      EQK = EG(31)/EG(2)/EG(30)/PFAC1
      RB(173) = RF(173) / MAX(EQK, SMALL)
C     R174: H + C4H81 = C2H4 + C2H5
      RF(174) = EXP(5.112687567511474D1 -2.39D0*ALOGT 
     * -5.625962344214302D3*TI)
      EQK = EG(21)*EG(22)/EG(2)/EG(30)
      RB(174) = RF(174) / MAX(EQK, SMALL)
C     R175: H + C4H81 = CH3 + C3H6
      RF(175) = EXP(5.182002285567469D1 -2.39D0*ALOGT 
     * -5.625962344214302D3*TI)
      EQK = EG(12)*EG(26)/EG(2)/EG(30)
      RB(175) = RF(175) / MAX(EQK, SMALL)
C     R176: H + C4H81 = H2 + C4H7
      RF(176) = EXP(1.338472764187182D1 +2.54D0*ALOGT 
     * -3.399731806575298D3*TI)
      EQK = EG(6)*EG(29)/EG(2)/EG(30)
      RB(176) = RF(176) / MAX(EQK, SMALL)
C     R177: O + C4H81 = HCO + nC3H7
      RF(177) = EXP(1.96146032124248D1 +1.45D0*ALOGT 
     * +2.022931003912477D2*TI)
      EQK = EG(14)*EG(27)/EG(3)/EG(30)
      RB(177) = RF(177) / MAX(EQK, SMALL)
C     R178: O + C4H81 = OH + C4H7
      RF(178) = EXP(3.033907131703076D1 -2.898528005605937D3*TI)
      EQK = EG(4)*EG(29)/EG(3)/EG(30)
      RB(178) = RF(178) / MAX(EQK, SMALL)
C     R179: O + C4H81 = OH + C4H7
      RF(179) = EXP(3.088911765395003D1 -2.249378504350441D3*TI)
      EQK = EG(4)*EG(29)/EG(3)/EG(30)
      RB(179) = RF(179) / MAX(EQK, SMALL)
C     R180: OH + C4H81 = H2O + C4H7
      RF(180) = EXP(6.551080335043404D0 +2.66D0*ALOGT 
     * -2.651951838462377D2*TI)
      EQK = EG(7)*EG(29)/EG(4)/EG(30)
      RB(180) = RF(180) / MAX(EQK, SMALL)
C     R181: O2 + C4H81 = HO2 + C4H7
      RF(181) = EXP(3.062675338948254D1 -2.562882488290111D4*TI)
      EQK = EG(5)*EG(29)/EG(9)/EG(30)
      RB(181) = RF(181) / MAX(EQK, SMALL)
C     R182: HO2 + C4H81 = H2O2 + C4H7
      RF(182) = EXP(2.763102111592855D1 -7.216127013956449D3*TI)
      EQK = EG(8)*EG(29)/EG(5)/EG(30)
      RB(182) = RF(182) / MAX(EQK, SMALL)
C     R183: CH3 + C4H81 = CH4 + C4H7
      RF(183) = EXP(-7.985076962177716D-1 +3.65D0*ALOGT 
     * -3.599508823628346D3*TI)
      EQK = EG(13)*EG(29)/EG(12)/EG(30)
      RB(183) = RF(183) / MAX(EQK, SMALL)
C     R184: H + pC4H9 = 2C2H5
      RF(184) = EXP(5.657037505150728D1 -2.92D0*ALOGT 
     * -6.292724428837196D3*TI)
      EQK = EG(22)*EG(22)/EG(2)/EG(31)
      RB(184) = RF(184) / MAX(EQK, SMALL)
C     R185: H + pC4H9 = H2 + C4H81
      RF(185) = 1.8D12
      EQK = EG(6)*EG(30)/EG(2)/EG(31)
      RB(185) = RF(185) / MAX(EQK, SMALL)
C     R186: O + pC4H9 = CH2O + nC3H7
      RF(186) = 9.6D13
      EQK = EG(15)*EG(27)/EG(3)/EG(31)
      RB(186) = RF(186) / MAX(EQK, SMALL)
C     R187: OH + pC4H9 = H2O + C4H81
      RF(187) = 2.4D13
      EQK = EG(7)*EG(30)/EG(4)/EG(31)
      RB(187) = RF(187) / MAX(EQK, SMALL)
C     R188: O2 + pC4H9 = HO2 + C4H81
      RF(188) = 2.7D11
      EQK = EG(5)*EG(30)/EG(9)/EG(31)
      RB(188) = RF(188) / MAX(EQK, SMALL)
C     R189: HO2 + pC4H9 = OH + CH2O + nC3H7
      RF(189) = 2.4D13
      EQK = EG(4)*EG(15)*EG(27)/EG(5)/EG(31)*PFAC1
      RB(189) = RF(189) / MAX(EQK, SMALL)
C     R190: CH3 + pC4H9 = CH4 + C4H81
      RF(190) = 1.1D13
      EQK = EG(13)*EG(30)/EG(12)/EG(31)
      RB(190) = RF(190) / MAX(EQK, SMALL)
C     R191: C5H9 = C2H4 + aC3H5
      RF(191) = EXP(3.084989694079675D1 -1.510612153188287D4*TI)
      RB(191) = 0D0
C     R192: C5H9 = C2H3 + C3H6
      RF(192) = EXP(3.084989694079675D1 -1.510612153188287D4*TI)
      RB(192) = 0D0
C     R193: H + C5H10 = PXC5H11
      RF(193) = EXP(3.021878515115626D1 -1.640838588173486D3*TI)
      EQK = EG(34)/EG(2)/EG(33)/PFAC1
      RB(193) = RF(193) / MAX(EQK, SMALL)
C     R194: H + C5H10 = C2H4 + nC3H7
      RF(194) = EXP(5.043372849455479D1 -2.39D0*ALOGT 
     * -5.625962344214302D3*TI)
      EQK = EG(21)*EG(27)/EG(2)/EG(33)
      RB(194) = RF(194) / MAX(EQK, SMALL)
C     R195: H + C5H10 = C2H5 + C3H6
      RF(195) = EXP(5.112687567511474D1 -2.39D0*ALOGT 
     * -5.625962344214302D3*TI)
      EQK = EG(22)*EG(26)/EG(2)/EG(33)
      RB(195) = RF(195) / MAX(EQK, SMALL)
C     R196: C2H4 + nC3H7 = PXC5H11
      RF(196) = EXP(2.642704831160261D1 -3.673481673771414D3*TI)
      EQK = EG(34)/EG(21)/EG(27)/PFAC1
      RB(196) = RF(196) / MAX(EQK, SMALL)
C     R197: H + C6H12 = PXC6H13
      RF(197) = EXP(3.021878515115626D1 -1.640838588173486D3*TI)
      EQK = EG(36)/EG(2)/EG(35)/PFAC1
      RB(197) = RF(197) / MAX(EQK, SMALL)
C     R198: H + C6H12 = C2H4 + pC4H9
      RF(198) = EXP(5.043372849455479D1 -2.39D0*ALOGT 
     * -5.625962344214302D3*TI)
      EQK = EG(21)*EG(31)/EG(2)/EG(35)
      RB(198) = RF(198) / MAX(EQK, SMALL)
C     R199: H + C6H12 = C3H6 + nC3H7
      RF(199) = EXP(5.112687567511474D1 -2.39D0*ALOGT 
     * -5.625962344214302D3*TI)
      EQK = EG(26)*EG(27)/EG(2)/EG(35)
      RB(199) = RF(199) / MAX(EQK, SMALL)
C     R200: C2H4 + pC4H9 = PXC6H13
      RF(200) = EXP(2.642704831160261D1 -3.673481673771414D3*TI)
      EQK = EG(36)/EG(21)/EG(31)/PFAC1
      RB(200) = RF(200) / MAX(EQK, SMALL)
C     R201: H + C7H14 = PXC7H15
      RF(201) = EXP(3.021878515115626D1 -1.640838588173486D3*TI)
      EQK = EG(38)/EG(2)/EG(37)/PFAC1
      RB(201) = RF(201) / MAX(EQK, SMALL)
C     R202: H + C7H14 = C2H4 + PXC5H11
      RF(202) = EXP(5.043372849455479D1 -2.39D0*ALOGT 
     * -5.625962344214302D3*TI)
      EQK = EG(21)*EG(34)/EG(2)/EG(37)
      RB(202) = RF(202) / MAX(EQK, SMALL)
C     R203: H + C7H14 = C3H6 + pC4H9
      RF(203) = EXP(5.112687567511474D1 -2.39D0*ALOGT 
     * -5.625962344214302D3*TI)
      EQK = EG(26)*EG(31)/EG(2)/EG(37)
      RB(203) = RF(203) / MAX(EQK, SMALL)
C     R204: C2H4 + PXC5H11 = PXC7H15
      RF(204) = EXP(2.642704831160261D1 -3.673481673771414D3*TI)
      EQK = EG(38)/EG(21)/EG(34)/PFAC1
      RB(204) = RF(204) / MAX(EQK, SMALL)
C     R205: H + C8H16 = PXC8H17
      RF(205) = EXP(3.021878515115626D1 -1.640838588173486D3*TI)
      EQK = EG(40)/EG(2)/EG(39)/PFAC1
      RB(205) = RF(205) / MAX(EQK, SMALL)
C     R206: H + C8H16 = C2H4 + PXC6H13
      RF(206) = EXP(5.043372849455479D1 -2.39D0*ALOGT 
     * -5.625962344214302D3*TI)
      EQK = EG(21)*EG(36)/EG(2)/EG(39)
      RB(206) = RF(206) / MAX(EQK, SMALL)
C     R207: H + C8H16 = C3H6 + PXC5H11
      RF(207) = EXP(5.112687567511474D1 -2.39D0*ALOGT 
     * -5.625962344214302D3*TI)
      EQK = EG(26)*EG(34)/EG(2)/EG(39)
      RB(207) = RF(207) / MAX(EQK, SMALL)
C     R208: C2H4 + PXC6H13 = PXC8H17
      RF(208) = EXP(2.642704831160261D1 -3.673481673771414D3*TI)
      EQK = EG(40)/EG(21)/EG(36)/PFAC1
      RB(208) = RF(208) / MAX(EQK, SMALL)
C     R209: H + C9H18 = PXC9H19
      RF(209) = EXP(3.021878515115626D1 -1.640838588173486D3*TI)
      EQK = EG(42)/EG(2)/EG(41)/PFAC1
      RB(209) = RF(209) / MAX(EQK, SMALL)
C     R210: H + C9H18 = C2H4 + PXC7H15
      RF(210) = EXP(5.043372849455479D1 -2.39D0*ALOGT 
     * -5.625962344214302D3*TI)
      EQK = EG(21)*EG(38)/EG(2)/EG(41)
      RB(210) = RF(210) / MAX(EQK, SMALL)
C     R211: H + C9H18 = C3H6 + PXC6H13
      RF(211) = EXP(5.112687567511474D1 -2.39D0*ALOGT 
     * -5.625962344214302D3*TI)
      EQK = EG(26)*EG(36)/EG(2)/EG(41)
      RB(211) = RF(211) / MAX(EQK, SMALL)
C     R212: C2H4 + PXC7H15 = PXC9H19
      RF(212) = EXP(2.642704831160261D1 -3.673481673771414D3*TI)
      EQK = EG(42)/EG(21)/EG(38)/PFAC1
      RB(212) = RF(212) / MAX(EQK, SMALL)
C     R213: H + C10H20 = PXC10H21
      RF(213) = EXP(3.021878515115626D1 -1.640838588173486D3*TI)
      EQK = EG(44)/EG(2)/EG(43)/PFAC1
      RB(213) = RF(213) / MAX(EQK, SMALL)
C     R214: H + C10H20 = C2H4 + PXC8H17
      RF(214) = EXP(5.043372849455479D1 -2.39D0*ALOGT 
     * -5.625962344214302D3*TI)
      EQK = EG(21)*EG(40)/EG(2)/EG(43)
      RB(214) = RF(214) / MAX(EQK, SMALL)
C     R215: H + C10H20 = C3H6 + PXC7H15
      RF(215) = EXP(5.112687567511474D1 -2.39D0*ALOGT 
     * -5.625962344214302D3*TI)
      EQK = EG(26)*EG(38)/EG(2)/EG(43)
      RB(215) = RF(215) / MAX(EQK, SMALL)
C     R216: C2H4 + PXC8H17 = PXC10H21
      RF(216) = EXP(2.642704831160261D1 -3.673481673771414D3*TI)
      EQK = EG(44)/EG(21)/EG(40)/PFAC1
      RB(216) = RF(216) / MAX(EQK, SMALL)
C     R217: C12H24 = C5H9 + PXC7H15
      RF(217) = EXP(3.80941244564001D1 -3.569663043070622D4*TI)
      EQK = EG(32)*EG(38)/EG(48)*PFAC1
      RB(217) = RF(217) / MAX(EQK, SMALL)
C     R218: C2H4 + PXC10H21 = PXC12H25
      RF(218) = EXP(2.642704831160261D1 -3.673481673771414D3*TI)
      EQK = EG(45)/EG(21)/EG(44)/PFAC1
      RB(218) = RF(218) / MAX(EQK, SMALL)
C     R219: PXC12H25 = S3XC12H25
      RF(219) = EXP(2.893121277799503D1 -6.D-1*ALOGT 
     * -7.246320014014844D3*TI)
      EQK = EG(47)/EG(45)
      RB(219) = RF(219) / MAX(EQK, SMALL)
C     R220: C3H6 + PXC9H19 = SXC12H25
      RF(220) = EXP(2.642704831160261D1 -3.673481673771414D3*TI)
      EQK = EG(46)/EG(26)/EG(42)/PFAC1
      RB(220) = RF(220) / MAX(EQK, SMALL)
C     R221: C4H81 + PXC8H17 = SXC12H25
      RF(221) = EXP(2.642704831160261D1 -3.673481673771414D3*TI)
      EQK = EG(46)/EG(30)/EG(40)/PFAC1
      RB(221) = RF(221) / MAX(EQK, SMALL)
C     R222: C5H10 + PXC7H15 = S3XC12H25
      RF(222) = EXP(2.642704831160261D1 -3.673481673771414D3*TI)
      EQK = EG(47)/EG(33)/EG(38)/PFAC1
      RB(222) = RF(222) / MAX(EQK, SMALL)
C     R223: C2H5 + C10H20 = S3XC12H25
      RF(223) = EXP(2.642704831160261D1 -3.673481673771414D3*TI)
      EQK = EG(47)/EG(22)/EG(43)/PFAC1
      RB(223) = RF(223) / MAX(EQK, SMALL)
C     R224: C6H12 + PXC6H13 = S3XC12H25
      RF(224) = EXP(2.642704831160261D1 -3.673481673771414D3*TI)
      EQK = EG(47)/EG(35)/EG(36)/PFAC1
      RB(224) = RF(224) / MAX(EQK, SMALL)
C     R225: nC3H7 + C9H18 = S3XC12H25
      RF(225) = EXP(2.642704831160261D1 -3.673481673771414D3*TI)
      EQK = EG(47)/EG(27)/EG(41)/PFAC1
      RB(225) = RF(225) / MAX(EQK, SMALL)
C     R226: PXC5H11 + C7H14 = S3XC12H25
      RF(226) = EXP(2.642704831160261D1 -3.673481673771414D3*TI)
      EQK = EG(47)/EG(34)/EG(37)/PFAC1
      RB(226) = RF(226) / MAX(EQK, SMALL)
C     R227: pC4H9 + C8H16 = S3XC12H25
      RF(227) = EXP(2.642704831160261D1 -3.673481673771414D3*TI)
      EQK = EG(47)/EG(31)/EG(39)/PFAC1
      RB(227) = RF(227) / MAX(EQK, SMALL)
C     R228: C2H5 + PXC10H21 = NC12H26
      RF(228) = EXP(3.28674630787585D1 -5.D-1*ALOGT)
      EQK = EG(1)/EG(22)/EG(44)/PFAC1
      RB(228) = RF(228) / MAX(EQK, SMALL)
C     R229: nC3H7 + PXC9H19 = NC12H26
      RF(229) = EXP(3.28674630787585D1 -5.D-1*ALOGT)
      EQK = EG(1)/EG(27)/EG(42)/PFAC1
      RB(229) = RF(229) / MAX(EQK, SMALL)
C     R230: pC4H9 + PXC8H17 = NC12H26
      RF(230) = EXP(3.28674630787585D1 -5.D-1*ALOGT)
      EQK = EG(1)/EG(31)/EG(40)/PFAC1
      RB(230) = RF(230) / MAX(EQK, SMALL)
C     R231: PXC5H11 + PXC7H15 = NC12H26
      RF(231) = EXP(3.28674630787585D1 -5.D-1*ALOGT)
      EQK = EG(1)/EG(34)/EG(38)/PFAC1
      RB(231) = RF(231) / MAX(EQK, SMALL)
C     R232: 2PXC6H13 = NC12H26
      RF(232) = EXP(3.28674630787585D1 -5.D-1*ALOGT)
      EQK = EG(1)/EG(36)/EG(36)/PFAC1
      RB(232) = RF(232) / MAX(EQK, SMALL)
C     R233: NC12H26 + H = H2 + PXC12H25
      RF(233) = EXP(1.407787482243177D1 +2.54D0*ALOGT 
     * -3.399731806575298D3*TI)
      EQK = EG(6)*EG(45)/EG(1)/EG(2)
      RB(233) = RF(233) / MAX(EQK, SMALL)
C     R234: NC12H26 + H = H2 + SXC12H25
      RF(234) = EXP(1.477102200299171D1 +2.4D0*ALOGT 
     * -2.249881721018081D3*TI)
      EQK = EG(6)*EG(46)/EG(1)/EG(2)
      RB(234) = RF(234) / MAX(EQK, SMALL)
C     R235: NC12H26 + H = H2 + S3XC12H25
      RF(235) = EXP(1.517648711109987D1 +2.4D0*ALOGT 
     * -2.249881721018081D3*TI)
      EQK = EG(6)*EG(47)/EG(1)/EG(2)
      RB(235) = RF(235) / MAX(EQK, SMALL)
C     R236: NC12H26 + O = OH + PXC12H25
      RF(236) = EXP(1.215477935114262D1 +2.68D0*ALOGT 
     * -1.869953136949942D3*TI)
      EQK = EG(4)*EG(45)/EG(1)/EG(3)
      RB(236) = RF(236) / MAX(EQK, SMALL)
C     R237: NC12H26 + O = OH + SXC12H25
      RF(237) = EXP(1.146373522077946D1 +2.71D0*ALOGT 
     * -1.059774302049671D3*TI)
      EQK = EG(4)*EG(46)/EG(1)/EG(3)
      RB(237) = RF(237) / MAX(EQK, SMALL)
C     R238: NC12H26 + O = OH + S3XC12H25
      RF(238) = EXP(1.186920032888762D1 +2.71D0*ALOGT 
     * -1.059774302049671D3*TI)
      EQK = EG(4)*EG(47)/EG(1)/EG(3)
      RB(238) = RF(238) / MAX(EQK, SMALL)
C     R239: NC12H26 + OH = H2O + PXC12H25
      RF(239) = EXP(8.131530710604253D0 +2.66D0*ALOGT 
     * -2.651951838462377D2*TI)
      EQK = EG(7)*EG(45)/EG(1)/EG(4)
      RB(239) = RF(239) / MAX(EQK, SMALL)
C     R240: NC12H26 + OH = H2O + SXC12H25
      RF(240) = EXP(1.121182037218631D1 +2.39D0*ALOGT 
     * -1.977641503824884D2*TI)
      EQK = EG(7)*EG(46)/EG(1)/EG(4)
      RB(240) = RF(240) / MAX(EQK, SMALL)
C     R241: NC12H26 + OH = H2O + S3XC12H25
      RF(241) = EXP(1.15228757958234D1 +2.39D0*ALOGT 
     * -1.977641503824884D2*TI)
      EQK = EG(7)*EG(47)/EG(1)/EG(4)
      RB(241) = RF(241) / MAX(EQK, SMALL)
C     R242: NC12H26 + O2 = HO2 + PXC12H25
      RF(242) = EXP(3.131990057004249D1 -2.562882488290111D4*TI)
      EQK = EG(5)*EG(45)/EG(1)/EG(9)
      RB(242) = RF(242) / MAX(EQK, SMALL)
C     R243: NC12H26 + O2 = HO2 + SXC12H25
      RF(243) = EXP(3.201304775060243D1 -2.394808121298378D4*TI)
      EQK = EG(5)*EG(46)/EG(1)/EG(9)
      RB(243) = RF(243) / MAX(EQK, SMALL)
C     R244: NC12H26 + O2 = HO2 + S3XC12H25
      RF(244) = EXP(3.241851285871059D1 -2.394808121298378D4*TI)
      EQK = EG(5)*EG(47)/EG(1)/EG(9)
      RB(244) = RF(244) / MAX(EQK, SMALL)
C     R245: NC12H26 + HO2 = H2O2 + PXC12H25
      RF(245) = EXP(1.112136326203106D1 +2.55D0*ALOGT 
     * -8.298042849382275D3*TI)
      EQK = EG(8)*EG(45)/EG(1)/EG(5)
      RB(245) = RF(245) / MAX(EQK, SMALL)
C     R246: NC12H26 + HO2 = H2O2 + SXC12H25
      RF(246) = EXP(1.139639164871428D1 +2.6D0*ALOGT 
     * -6.999743846871283D3*TI)
      EQK = EG(8)*EG(46)/EG(1)/EG(5)
      RB(246) = RF(246) / MAX(EQK, SMALL)
C     R247: NC12H26 + HO2 = H2O2 + S3XC12H25
      RF(247) = EXP(1.139075783099602D1 +2.6D0*ALOGT 
     * -6.999743846871283D3*TI)
      EQK = EG(8)*EG(47)/EG(1)/EG(5)
      RB(247) = RF(247) / MAX(EQK, SMALL)
C     R248: NC12H26 + CH3 = CH4 + PXC12H25
      RF(248) = EXP(5.933268452777344D-1 +3.65D0*ALOGT 
     * -3.599508823628346D3*TI)
      EQK = EG(13)*EG(45)/EG(1)/EG(12)
      RB(248) = RF(248) / MAX(EQK, SMALL)
C     R249: NC12H26 + CH3 = CH4 + SXC12H25
      RF(249) = EXP(1.791759469228055D0 +3.46D0*ALOGT 
     * -2.75762733866676D3*TI)
      EQK = EG(13)*EG(46)/EG(1)/EG(12)
      RB(249) = RF(249) / MAX(EQK, SMALL)
C     R250: NC12H26 + CH3 = CH4 + S3XC12H25
      RF(250) = EXP(2.19722457733622D0 +3.46D0*ALOGT 
     * -2.75762733866676D3*TI)
      EQK = EG(13)*EG(47)/EG(1)/EG(12)
      RB(250) = RF(250) / MAX(EQK, SMALL)
C     R251: O2 + PXC12H25 = C12H25O2
      RF(251) = 5.D13
      RB(251) = 0D0
C     R252: C12H25O2 = O2 + PXC12H25
      RF(252) = EXP(3.094520712060107D1 -1.37881366933338D4*TI)
      RB(252) = 0D0
C     R253: O2 + SXC12H25 = C12H25O2
      RF(253) = 5.D13
      RB(253) = 0D0
C     R254: C12H25O2 = O2 + SXC12H25
      RF(254) = EXP(3.094520712060107D1 -1.37881366933338D4*TI)
      RB(254) = 0D0
C     R255: O2 + S3XC12H25 = C12H25O2
      RF(255) = 5.D13
      RB(255) = 0D0
C     R256: C12H25O2 = O2 + S3XC12H25
      RF(256) = EXP(3.094520712060107D1 -1.37881366933338D4*TI)
      RB(256) = 0D0
C     R257: C12H25O2 = C12OOH
      RF(257) = EXP(2.804313076675538D1 -9.561116685158473D3*TI)
      RB(257) = 0D0
C     R258: C12OOH = C12H25O2
      RF(258) = EXP(2.53284360229345D1 -5.786991677859076D3*TI)
      RB(258) = 0D0
C     R259: O2 + PXC12H25 = HO2 + C12H24
      RF(259) = EXP(2.658119899142987D1 -3.019300005839518D3*TI)
      RB(259) = 0D0
C     R260: HO2 + C12H24 = O2 + PXC12H25
      RF(260) = EXP(2.647900805053332D1 -9.812725018978434D3*TI)
      RB(260) = 0D0
C     R261: O2 + SXC12H25 = HO2 + C12H24
      RF(261) = EXP(2.658119899142987D1 -3.019300005839518D3*TI)
      RB(261) = 0D0
C     R262: HO2 + C12H24 = O2 + SXC12H25
      RF(262) = EXP(2.647900805053332D1 -9.812725018978434D3*TI)
      RB(262) = 0D0
C     R263: O2 + S3XC12H25 = HO2 + C12H24
      RF(263) = EXP(2.658119899142987D1 -3.019300005839518D3*TI)
      RB(263) = 0D0
C     R264: HO2 + C12H24 = O2 + S3XC12H25
      RF(264) = EXP(2.647900805053332D1 -9.812725018978434D3*TI)
      RB(264) = 0D0
C     R265: O2 + C12OOH = O2C12H24OOH
      RF(265) = 4.6D10
      RB(265) = 0D0
C     R266: O2C12H24OOH = O2 + C12OOH
      RF(266) = EXP(3.085388896206629D1 -1.37881366933338D4*TI)
      RB(266) = 0D0
C     R267: O2C12H24OOH = OH + OC12H23OOH
      RF(267) = EXP(2.521190220667855D1 -8.554683349878635D3*TI)
      EQK = EG(4)*EG(52)/EG(51)*PFAC1
      RB(267) = RF(267) / MAX(EQK, SMALL)
C     R268: OC12H23OOH = OH + 3C2H4 + C2H5 + 2CH2CHO
      RF(268) = EXP(3.51265630598128D1 -2.116780912427322D4*TI)
      RB(268) = 0D0
C
      RKLOW(1) = EXP(4.559410099735222D1 -1.4D0*ALOGT)
      RKLOW(2) = EXP(3.984208130296976D1 -5.84D-1*ALOGT 
     * +1.153875818898336D3*TI)
      RKLOW(3) = EXP(5.542160680152843D1 -2.79D0*ALOGT 
     * -2.108981054078904D3*TI)
      RKLOW(4) = EXP(6.379313832844233D1 -3.42D0*ALOGT 
     * -4.244632591542722D4*TI)
      RKLOW(5) = EXP(5.556214682430743D1 -2.57D0*ALOGT 
     * -7.170837513868855D2*TI)
      RKLOW(6) = EXP(6.333294832064492D1 -3.14D0*ALOGT 
     * -6.189565011971013D2*TI)
      RKLOW(7) = EXP(6.986601015018565D1 -4.8D0*ALOGT 
     * -2.797884672077953D3*TI)
      RKLOW(8) = EXP(7.689235621931073D1 -4.76D0*ALOGT 
     * -1.227848669041404D3*TI)
      RKLOW(9) = EXP(1.15700234196288D2 -9.67D0*ALOGT 
     * -3.1300076727203D3*TI)
      RKLOW(10) = EXP(6.311175598946197D1 -3.4D0*ALOGT 
     * -1.801451258417455D4*TI)
      RKLOW(11) = EXP(6.941402502644259D1 -3.86D0*ALOGT 
     * -1.670679336564533D3*TI)
      RKLOW(12) = EXP(1.350015492208952D2 -1.194D1*ALOGT 
     * -4.916326199508487D3*TI)
      RKLOW(13) = EXP(4.299728058997825D1 -3.800594204017257D2*TI)
      RKLOW(14) = EXP(9.509412345149228D1 -7.08D0*ALOGT 
     * -3.364003423172863D3*TI)
      RKLOW(15) = EXP(1.293830201385977D2 -1.179D1*ALOGT 
     * -4.521150150410858D3*TI)
      RKLOW(16) = EXP(1.384402845218764D2 -1.2D1*ALOGT 
     * -3.003096429141513D3*TI)
      RKLOW(17) = EXP(1.3951864295364D2 -1.281D1*ALOGT 
     * -3.145104172749498D3*TI)
      RKLOW(18) = EXP(8.933241371888575D1 -6.66D0*ALOGT 
     * -3.522516673479438D3*TI)
      RKLOW(19) = EXP(1.11626024542475D2 -9.32D0*ALOGT 
     * -2.935564752344236D3*TI)
      RKLOW(20) = EXP(8.933241371888575D1 -6.66D0*ALOGT 
     * -3.522516673479438D3*TI)
      RKLOW(21) = EXP(8.933241371888575D1 -6.66D0*ALOGT 
     * -3.522516673479438D3*TI)
      RKLOW(22) = EXP(8.933241371888575D1 -6.66D0*ALOGT 
     * -3.522516673479438D3*TI)
      RKLOW(23) = EXP(8.933241371888575D1 -6.66D0*ALOGT 
     * -3.522516673479438D3*TI)
      RKLOW(24) = EXP(8.933241371888575D1 -6.66D0*ALOGT 
     * -3.522516673479438D3*TI)
      RKLOW(25) = EXP(8.933241371888575D1 -6.66D0*ALOGT 
     * -3.522516673479438D3*TI)
      RKLOW(26) = EXP(8.933241371888575D1 -6.66D0*ALOGT 
     * -3.522516673479438D3*TI)
C
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE RDSMH  (T, SMH)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
      DIMENSION SMH(*), TN(5)
C
      TLOG = LOG(T)
      TI = 1D0/T
C
      TN(1) = TLOG - 1D0
      TN(2) = T
      TN(3) = TN(2)*T
      TN(4) = TN(3)*T
      TN(5) = TN(4)*T
C NC12H26
      IF (T .GT. 1391) THEN
      SMH(1) = -1.72670922D2 +5.48843465D4*TI 
     *         +3.85095037D1*TN(1) +2.81775024D-2*TN(2) 
     *         -3.191553333333333D-6*TN(3) +2.46687385D-10*TN(4) 
     *         -8.5622075D-15*TN(5) 
      ELSE
      SMH(1) = 5.00994626D1 +4.00654253D4*TI 
     *         -2.62181594D0*TN(1) +7.36188555D-2*TN(2) 
     *         -1.573283785D-5*TN(3) +2.562010566666667D-9*TN(4) 
     *         -2.01801115D-13*TN(5) 
      ENDIF
C H
      IF (T .GT. 1000) THEN
      SMH(2) = -4.46682914D-1 -2.54736599D4*TI 
     *         +2.50000001D0*TN(1) -1.154214865D-11*TN(2) 
     *         +2.692699133333334D-15*TN(3) -3.945960291666667D-19*TN(4) 
     *         +2.490986785D-23*TN(5) 
      ELSE
      SMH(2) = -4.46682853D-1 -2.54736599D4*TI 
     *         +2.5D0*TN(1) +3.526664095D-13*TN(2) 
     *         -3.326532733333333D-16*TN(3) +1.917346933333333D-19*TN(4) 
     *         -4.63866166D-23*TN(5) 
      ENDIF
C O
      IF (T .GT. 1000) THEN
      SMH(3) = 4.78433864D0 -2.92175791D4*TI 
     *         +2.56942078D0*TN(1) -4.298705685D-5*TN(2) 
     *         +6.991409816666667D-9*TN(3) -8.348149916666666D-13*TN(4) 
     *         +6.141684549999999D-17*TN(5) 
      ELSE
      SMH(3) = 2.05193346D0 -2.91222592D4*TI 
     *         +3.1682671D0*TN(1) -1.63965942D-3*TN(2) 
     *         +1.107177326666667D-6*TN(3) -5.106721866666666D-10*TN(4) 
     *         +1.056329855D-13*TN(5) 
      ENDIF
C OH
      IF (T .GT. 1000) THEN
      SMH(4) = 5.70164073D0 -3.71885774D3*TI 
     *         +2.86472886D0*TN(1) +5.2825224D-4*TN(2) 
     *         -4.318045966666667D-8*TN(3) +2.54348895D-12*TN(4) 
     *         -6.6597938D-17*TN(5) 
      ELSE
      SMH(4) = -6.9043296D-1 -3.38153812D3*TI 
     *         +4.12530561D0*TN(1) -1.612724695D-3*TN(2) 
     *         +1.087941151666667D-6*TN(3) -4.832113691666666D-10*TN(4) 
     *         +1.031186895D-13*TN(5) 
      ENDIF
C HO2
      IF (T .GT. 1000) THEN
      SMH(5) = 3.78510215D0 -1.11856713D2*TI 
     *         +4.0172109D0*TN(1) +1.119910065D-3*TN(2) 
     *         -1.056096916666667D-7*TN(3) +9.520530833333334D-12*TN(4) 
     *         -5.39542675D-16*TN(5) 
      ELSE
      SMH(5) = 3.71666245D0 -2.9480804D2*TI 
     *         +4.30179801D0*TN(1) -2.374560255D-3*TN(2) 
     *         +3.526381516666666D-6*TN(3) -2.02303245D-9*TN(4) 
     *         +4.646125620000001D-13*TN(5) 
      ENDIF
C H2
      IF (T .GT. 1000) THEN
      SMH(6) = -3.20502331D0 +9.50158922D2*TI 
     *         +3.3372792D0*TN(1) -2.470123655D-5*TN(2) 
     *         +8.324279633333333D-8*TN(3) -1.496386616666667D-11*TN(4) 
     *         +1.00127688D-15*TN(5) 
      ELSE
      SMH(6) = 6.830102380000001D-1 +9.17935173D2*TI 
     *         +2.34433112D0*TN(1) +3.990260375D-3*TN(2) 
     *         -3.2463585D-6*TN(3) +1.67976745D-9*TN(4) 
     *         -3.688058805D-13*TN(5) 
      ENDIF
C H2O
      IF (T .GT. 1000) THEN
      SMH(7) = 4.9667701D0 +3.00042971D4*TI 
     *         +3.03399249D0*TN(1) +1.08845902D-3*TN(2) 
     *         -2.734541966666666D-8*TN(3) -8.08683225D-12*TN(4) 
     *         +8.4100496D-16*TN(5) 
      ELSE
      SMH(7) = -8.49032208D-1 +3.02937267D4*TI 
     *         +4.19864056D0*TN(1) -1.01821705D-3*TN(2) 
     *         +1.086733685D-6*TN(3) -4.57330885D-10*TN(4) 
     *         +8.85989085D-14*TN(5) 
      ENDIF
C H2O2
      IF (T .GT. 1000) THEN
      SMH(8) = 2.91615662D0 +1.78617877D4*TI 
     *         +4.16500285D0*TN(1) +2.45415847D-3*TN(2) 
     *         -3.168987083333333D-7*TN(3) +3.09321655D-11*TN(4) 
     *         -1.439541525D-15*TN(5) 
      ELSE
      SMH(8) = 3.43505074D0 +1.77025821D4*TI 
     *         +4.27611269D0*TN(1) -2.714112085D-4*TN(2) 
     *         +2.78892835D-6*TN(3) -1.798090108333333D-9*TN(4) 
     *         +4.312271815D-13*TN(5) 
      ENDIF
C O2
      IF (T .GT. 1000) THEN
      SMH(9) = 5.45323129D0 +1.08845772D3*TI 
     *         +3.28253784D0*TN(1) +7.4154377D-4*TN(2) 
     *         -1.263277781666667D-7*TN(3) +1.745587958333333D-11*TN(4) 
     *         -1.08358897D-15*TN(5) 
      ELSE
      SMH(9) = 3.65767573D0 +1.06394356D3*TI 
     *         +3.78245636D0*TN(1) -1.49836708D-3*TN(2) 
     *         +1.641217001666667D-6*TN(3) -8.067745908333334D-10*TN(4) 
     *         +1.621864185D-13*TN(5) 
      ENDIF
C CH2
      IF (T .GT. 1000) THEN
      SMH(10) = 6.17119324D0 -4.6263604D4*TI 
     *         +2.87410113D0*TN(1) +1.82819646D-3*TN(2) 
     *         -2.348243283333333D-7*TN(3) +2.168162908333333D-11*TN(4) 
     *         -9.38637835D-16*TN(5) 
      ELSE
      SMH(10) = 1.56253185D0 -4.60040401D4*TI 
     *         +3.76267867D0*TN(1) +4.844360715D-4*TN(2) 
     *         +4.658164016666667D-7*TN(3) -3.209092941666667D-10*TN(4) 
     *         +8.43708595D-14*TN(5) 
      ENDIF
C CH2*
      IF (T .GT. 1000) THEN
      SMH(11) = 8.62650169D0 -5.09259997D4*TI 
     *         +2.29203842D0*TN(1) +2.327943185D-3*TN(2) 
     *         -3.353199116666667D-7*TN(3) +3.48255D-11*TN(4) 
     *         -1.698581825D-15*TN(5) 
      ELSE
      SMH(11) = -7.69118967D-1 -5.04968163D4*TI 
     *         +4.19860411D0*TN(1) -1.183307095D-3*TN(2) 
     *         +1.372160366666667D-6*TN(3) -5.573466508333334D-10*TN(4) 
     *         +9.71573685D-14*TN(5) 
      ENDIF
C CH3
      IF (T .GT. 1000) THEN
      SMH(12) = 8.48007179D0 -1.67755843D4*TI 
     *         +2.28571772D0*TN(1) +3.619950185D-3*TN(2) 
     *         -4.978572466666667D-7*TN(3) +4.9640387D-11*TN(4) 
     *         -2.33577197D-15*TN(5) 
      ELSE
      SMH(12) = 1.60456433D0 -1.64449988D4*TI 
     *         +3.6735904D0*TN(1) +1.005475875D-3*TN(2) 
     *         +9.550364266666668D-7*TN(3) -5.725978541666666D-10*TN(4) 
     *         +1.27192867D-13*TN(5) 
      ENDIF
C CH4
      IF (T .GT. 1000) THEN
      SMH(13) = 1.8437318D1 +9.468344590000001D3*TI 
     *         +7.4851495D-2*TN(1) +6.69547335D-3*TN(2) 
     *         -9.554763483333333D-7*TN(3) +1.019104458333333D-10*TN(4) 
     *         -5.0907615D-15*TN(5) 
      ELSE
      SMH(13) = -4.64130376D0 +1.02466476D4*TI 
     *         +5.14987613D0*TN(1) -6.8354894D-3*TN(2) 
     *         +8.19667665D-6*TN(3) -4.039525216666667D-9*TN(4) 
     *         +8.3346978D-13*TN(5) 
      ENDIF
C HCO
      IF (T .GT. 1000) THEN
      SMH(14) = 9.79834492D0 -4.01191815D3*TI 
     *         +2.77217438D0*TN(1) +2.47847763D-3*TN(2) 
     *         -4.140760216666667D-7*TN(3) +4.909681483333334D-11*TN(4) 
     *         -2.667543555D-15*TN(5) 
      ELSE
      SMH(14) = 3.39437243D0 -3.83956496D3*TI 
     *         +4.22118584D0*TN(1) -1.62196266D-3*TN(2) 
     *         +2.296657433333333D-6*TN(3) -1.109534108333333D-9*TN(4) 
     *         +2.168844325D-13*TN(5) 
      ENDIF
C CH2O
      IF (T .GT. 1000) THEN
      SMH(15) = 1.3656323D1 +1.39958323D4*TI 
     *         +1.76069008D0*TN(1) +4.60000041D-3*TN(2) 
     *         -7.370980216666666D-7*TN(3) +8.386767666666666D-11*TN(4) 
     *         -4.419278200000001D-15*TN(5) 
      ELSE
      SMH(15) = 6.028129D-1 +1.43089567D4*TI 
     *         +4.79372315D0*TN(1) -4.954166845D-3*TN(2) 
     *         +6.220333466666666D-6*TN(3) -3.160710508333333D-9*TN(4) 
     *         +6.5886326D-13*TN(5) 
      ENDIF
C CH3O
      IF (T .GT. 1000) THEN
      SMH(16) = -1.96680028D0 -3.7811194D2*TI 
     *         +4.75779238D0*TN(1) +3.72071237D-3*TN(2) 
     *         -4.495086266666667D-7*TN(3) +3.6507542D-11*TN(4) 
     *         -1.31768549D-15*TN(5) 
      ELSE
      SMH(16) = 6.57240864D0 -1.2956976D3*TI 
     *         +3.71180502D0*TN(1) -1.40231653D-3*TN(2) 
     *         +6.275849516666667D-6*TN(3) -3.942267408333333D-9*TN(4) 
     *         +9.329421000000001D-13*TN(5) 
      ENDIF
C CO
      IF (T .GT. 1000) THEN
      SMH(17) = 7.81868772D0 +1.41518724D4*TI 
     *         +2.71518561D0*TN(1) +1.031263715D-3*TN(2) 
     *         -1.664709618333334D-7*TN(3) +1.9171084D-11*TN(4) 
     *         -1.01823858D-15*TN(5) 
      ELSE
      SMH(17) = 3.50840928D0 +1.4344086D4*TI 
     *         +3.57953347D0*TN(1) -3.0517684D-4*TN(2) 
     *         +1.69469055D-7*TN(3) +7.558382366666667D-11*TN(4) 
     *         -4.522122495D-14*TN(5) 
      ENDIF
C CO2
      IF (T .GT. 1000) THEN
      SMH(18) = 2.27163806D0 +4.8759166D4*TI 
     *         +3.85746029D0*TN(1) +2.20718513D-3*TN(2) 
     *         -3.691356733333334D-7*TN(3) +4.362418233333334D-11*TN(4) 
     *         -2.36042082D-15*TN(5) 
      ELSE
      SMH(18) = 9.90105222D0 +4.83719697D4*TI 
     *         +2.35677352D0*TN(1) +4.492298385D-3*TN(2) 
     *         -1.187260448333333D-6*TN(3) +2.049325183333333D-10*TN(4) 
     *         -7.184977399999999D-15*TN(5) 
      ENDIF
C C2H2
      IF (T .GT. 1000) THEN
      SMH(19) = -1.23028121D0 -2.59359992D4*TI 
     *         +4.14756964D0*TN(1) +2.98083332D-3*TN(2) 
     *         -3.9549142D-7*TN(3) +3.895101425D-11*TN(4) 
     *         -1.806176065D-15*TN(5) 
      ELSE
      SMH(19) = 1.39397051D1 -2.64289807D4*TI 
     *         +8.08681094D-1*TN(1) +1.168078145D-2*TN(2) 
     *         -5.91953025D-6*TN(3) +2.334603641666667D-9*TN(4) 
     *         -4.25036487D-13*TN(5) 
      ENDIF
C C2H3
      IF (T .GT. 1000) THEN
      SMH(20) = 7.78732378D0 -3.46128739D4*TI 
     *         +3.016724D0*TN(1) +5.1651146D-3*TN(2) 
     *         -7.801372483333333D-7*TN(3) +8.480274D-11*TN(4) 
     *         -4.313035205D-15*TN(5) 
      ELSE
      SMH(20) = 8.51054025D0 -3.48598468D4*TI 
     *         +3.21246645D0*TN(1) +7.5739581D-4*TN(2) 
     *         +4.320156866666666D-6*TN(3) -2.980482058333333D-9*TN(4) 
     *         +7.35754365D-13*TN(5) 
      ENDIF
C C2H4
      IF (T .GT. 1000) THEN
      SMH(21) = 1.03053693D1 -4.93988614D3*TI 
     *         +2.03611116D0*TN(1) +7.32270755D-3*TN(2) 
     *         -1.118463191666667D-6*TN(3) +1.226857691666667D-10*TN(4) 
     *         -6.28530305D-15*TN(5) 
      ELSE
      SMH(21) = 4.09733096D0 -5.08977593D3*TI 
     *         +3.95920148D0*TN(1) -3.785261235D-3*TN(2) 
     *         +9.516504866666667D-6*TN(3) -5.763239608333333D-9*TN(4) 
     *         +1.349421865D-12*TN(5) 
      ENDIF
C C2H5
      IF (T .GT. 1000) THEN
      SMH(22) = 1.34624343D1 -1.285752D4*TI 
     *         +1.95465642D0*TN(1) +8.698636100000001D-3*TN(2) 
     *         -1.330344446666667D-6*TN(3) +1.460147408333333D-10*TN(4) 
     *         -7.4820788D-15*TN(5) 
      ELSE
      SMH(22) = 4.70720924D0 -1.28416265D4*TI 
     *         +4.30646568D0*TN(1) -2.09329446D-3*TN(2) 
     *         +8.28571345D-6*TN(3) -4.992721716666666D-9*TN(4) 
     *         +1.15254502D-12*TN(5) 
      ENDIF
C C2H6
      IF (T .GT. 1000) THEN
      SMH(23) = 1.51156107D1 +1.14263932D4*TI 
     *         +1.0718815D0*TN(1) +1.084263385D-2*TN(2) 
     *         -1.67093445D-6*TN(3) +1.845100008333333D-10*TN(4) 
     *         -9.5001445D-15*TN(5) 
      ELSE
      SMH(23) = 2.66682316D0 +1.15222055D4*TI 
     *         +4.29142492D0*TN(1) -2.75077135D-3*TN(2) 
     *         +9.990638133333334D-6*TN(3) -5.903885708333334D-9*TN(4) 
     *         +1.343428855D-12*TN(5) 
      ENDIF
C CH2CHO
      IF (T .GT. 1000) THEN
      SMH(24) = -5.0320879D0 +9.695D2*TI 
     *         +5.9756699D0*TN(1) +4.0652957D-3*TN(2) 
     *         -4.5727075D-7*TN(3) +3.391920083333333D-11*TN(4) 
     *         -1.08800855D-15*TN(5) 
      ELSE
      SMH(24) = 9.571453500000001D0 -6.2D1*TI 
     *         +3.4090624D0*TN(1) +5.369287D-3*TN(2) 
     *         +3.1524875D-7*TN(3) -5.965485916666667D-10*TN(4) 
     *         +1.43369255D-13*TN(5) 
      ENDIF
C aC3H5
      IF (T .GT. 1000) THEN
      SMH(25) = -1.124305D1 -1.7482449D4*TI 
     *         +6.5007877D0*TN(1) +7.1623655D-3*TN(2) 
     *         -9.463605333333332D-7*TN(3) +9.234000833333333D-11*TN(4) 
     *         -4.518194349999999D-15*TN(5) 
      ELSE
      SMH(25) = 1.7173214D1 -1.9245629D4*TI 
     *         +1.3631835D0*TN(1) +9.906910499999999D-3*TN(2) 
     *         +2.082843333333334D-6*TN(3) -2.779629583333333D-9*TN(4) 
     *         +7.9232855D-13*TN(5) 
      ENDIF
C C3H6
      IF (T .GT. 1000) THEN
      SMH(26) = -1.331335D1 +9.235703D2*TI 
     *         +6.732257D0*TN(1) +7.45417D-3*TN(2) 
     *         -8.249831666666666D-7*TN(3) +6.010018333333334D-11*TN(4) 
     *         -1.883102D-15*TN(5) 
      ELSE
      SMH(26) = 1.614534D1 -1.074826D3*TI 
     *         +1.493307D0*TN(1) +1.046259D-2*TN(2) 
     *         +7.47799D-7*TN(3) -1.39076D-9*TN(4) 
     *         +3.579073D-13*TN(5) 
      ENDIF
C nC3H7
      IF (T .GT. 1000) THEN
      SMH(27) = -1.5515297D1 -7.9762236D3*TI 
     *         +7.7097479D0*TN(1) +8.015742500000001D-3*TN(2) 
     *         -8.786706333333332D-7*TN(3) +6.324029333333334D-11*TN(4) 
     *         -1.94313595D-15*TN(5) 
      ELSE
      SMH(27) = 2.1136034D1 -1.0312346D4*TI 
     *         +1.0491173D0*TN(1) +1.30044865D-2*TN(2) 
     *         +3.923752666666667D-7*TN(3) -1.632927666666667D-9*TN(4) 
     *         +4.68601035D-13*TN(5) 
      ENDIF
C C2H3CHO
      IF (T .GT. 1000) THEN
      SMH(28) = -4.8588004D0 +1.0784054D4*TI 
     *         +5.8111868D0*TN(1) +8.557128000000001D-3*TN(2) 
     *         -1.247236016666667D-6*TN(3) +1.187687416666667D-10*TN(4) 
     *         -4.58734205D-15*TN(5) 
      ELSE
      SMH(28) = 1.9498077D1 +9.335734399999999D3*TI 
     *         +1.2713498D0*TN(1) +1.3115527D-2*TN(2) 
     *         -1.548538416666667D-6*TN(3) -3.986439333333333D-10*TN(4) 
     *         +1.67402715D-13*TN(5) 
      ENDIF
C C4H7
      IF (T .GT. 1000) THEN
      SMH(29) = -8.889308D0 -2.0955008D4*TI 
     *         +7.0134835D0*TN(1) +1.1317279D-2*TN(2) 
     *         -1.5424245D-6*TN(3) +1.400660583333333D-10*TN(4) 
     *         -5.2043085D-15*TN(5) 
      ELSE
      SMH(29) = 2.3437878D1 -2.2653328D4*TI 
     *         +7.4449432D-1*TN(1) +1.98394285D-2*TN(2) 
     *         -3.816347666666667D-6*TN(3) +1.779414416666667D-10*TN(4) 
     *         +1.15481875D-13*TN(5) 
      ENDIF
C C4H81
      IF (T .GT. 1000) THEN
      SMH(30) = 1.5543201D1 +2.1397231D3*TI 
     *         +2.0535841D0*TN(1) +1.71752535D-2*TN(2) 
     *         -2.6471995D-6*TN(3) +2.757471833333334D-10*TN(4) 
     *         -1.26805225D-14*TN(5) 
      ELSE
      SMH(30) = 2.1062469D1 +1.7904004D3*TI 
     *         +1.181138D0*TN(1) +1.542669D-2*TN(2) 
     *         +8.477541166666667D-7*TN(3) -2.054574D-9*TN(4) 
     *         +5.555096499999999D-13*TN(5) 
      ENDIF
C pC4H9
      IF (T .GT. 1000) THEN
      SMH(31) = -1.7891747D1 -4.964405800000001D3*TI 
     *         +8.6822395D0*TN(1) +1.18455355D-2*TN(2) 
     *         -1.265814416666667D-6*TN(3) +5.535594666666666D-11*TN(4) 
     *         +2.7422568D-15*TN(5) 
      ELSE
      SMH(31) = 2.2169268D1 -7.322104D3*TI 
     *         +1.2087042D0*TN(1) +1.91487485D-2*TN(2) 
     *         -1.211008483333333D-6*TN(3) -1.28571225D-9*TN(4) 
     *         +4.34297175D-13*TN(5) 
      ENDIF
C C5H9
      IF (T .GT. 1000) THEN
      SMH(32) = -3.3125885D1 +1.7218359D3*TI 
     *         +1.013864D1*TN(1) +1.1357069D-2*TN(2) 
     *         -1.298507716666667D-6*TN(3) +9.897101666666667D-11*TN(4) 
     *         -3.2966224D-15*TN(5) 
      ELSE
      SMH(32) = 3.6459244D1 -2.8121887D3*TI 
     *         -2.4190111D0*TN(1) +2.02151945D-2*TN(2) 
     *         +1.130038983333333D-6*TN(3) -2.810395166666667D-9*TN(4) 
     *         +7.558356500000001D-13*TN(5) 
      ENDIF
C C5H10
      IF (T .GT. 1392) THEN
      SMH(33) = -5.23683936D1 +1.00898205D4*TI 
     *         +1.45851539D1*TN(1) +1.120362355D-2*TN(2) 
     *         -1.272246708333333D-6*TN(3) +9.849080500000001D-11*TN(4) 
     *         -3.421925695D-15*TN(5) 
      ELSE
      SMH(33) = 3.2273979D1 +4.46546666D3*TI 
     *         -1.06223481D0*TN(1) +2.87109147D-2*TN(2) 
     *         -6.241448166666667D-6*TN(3) +1.061374908333333D-9*TN(4) 
     *         -8.980489449999999D-14*TN(5) 
      ENDIF
C PXC5H11
      IF (T .GT. 1390) THEN
      SMH(34) = -5.44829293D1 +9.80712307D2*TI 
     *         +1.52977446D1*TN(1) +1.19867655D-2*TN(2) 
     *         -1.363988246666667D-6*TN(3) +1.057358966666667D-10*TN(4) 
     *         -3.677045275D-15*TN(5) 
      ELSE
      SMH(34) = 2.87238666D1 -4.7161146D3*TI 
     *         +5.24384081D-2*TN(1) +2.80398479D-2*TN(2) 
     *         -5.525763383333334D-6*TN(3) +8.146114841666667D-10*TN(4) 
     *         -5.700483D-14*TN(5) 
      ENDIF
C C6H12
      IF (T .GT. 1392) THEN
      SMH(35) = -6.838188510000001D1 +1.4206286D4*TI 
     *         +1.78337529D1*TN(1) +1.33688829D-2*TN(2) 
     *         -1.516727955D-6*TN(3) +1.173498066666667D-10*TN(4) 
     *         -4.07562122D-15*TN(5) 
      ELSE
      SMH(35) = 3.53120691D1 +7.34368617D3*TI 
     *         -1.35275205D0*TN(1) +3.49327713D-2*TN(2) 
     *         -7.656800366666667D-6*TN(3) +1.308061191666667D-9*TN(4) 
     *         -1.106480875D-13*TN(5) 
      ENDIF
C PXC6H13
      IF (T .GT. 1390) THEN
      SMH(36) = -7.04490943D1 +5.09299041D3*TI 
     *         +1.8538547D1*TN(1) +1.41553981D-2*TN(2) 
     *         -1.60884541D-6*TN(3) +1.246229875D-10*TN(4) 
     *         -4.33168032D-15*TN(5) 
      ELSE
      SMH(36) = 3.16075093D1 -1.83280393D3*TI 
     *         -2.04871465D-1*TN(1) +3.41900636D-2*TN(2) 
     *         -6.9074652D-6*TN(3) +1.05129835D-9*TN(4) 
     *         -7.656002899999999D-14*TN(5) 
      ENDIF
C C7H14
      IF (T .GT. 1392) THEN
      SMH(37) = -8.44391108D1 +1.83260065D4*TI 
     *         +2.10898039D1*TN(1) +1.55303939D-2*TN(2) 
     *         -1.76074655D-6*TN(3) +1.361714833333333D-10*TN(4) 
     *         -4.727991095D-15*TN(5) 
      ELSE
      SMH(37) = 3.85068032D1 +1.02168601D4*TI 
     *         -1.67720549D0*TN(1) +4.123058005D-2*TN(2) 
     *         -9.1084018D-6*TN(3) +1.565519191666667D-9*TN(4) 
     *         -1.328689915D-13*TN(5) 
      ENDIF
C PXC7H15
      IF (T .GT. 1390) THEN
      SMH(38) = -8.649543110000001D1 +9.20938221D3*TI 
     *         +2.17940709D1*TN(1) +1.631401215D-2*TN(2) 
     *         -1.852304066666667D-6*TN(3) +1.4338929D-10*TN(4) 
     *         -4.981834995D-15*TN(5) 
      ELSE
      SMH(38) = 3.46564011D1 +1.04590223D3*TI 
     *         -4.99570406D-1*TN(1) +4.044132335D-2*TN(2) 
     *         -8.342212566666667D-6*TN(3) +1.304577566666667D-9*TN(4) 
     *         -9.83081135D-14*TN(5) 
      ENDIF
C C8H16
      IF (T .GT. 1392) THEN
      SMH(39) = -1.00537716D2 +2.24485674D4*TI 
     *         +2.43540125D1*TN(1) +1.76833231D-2*TN(2) 
     *         -2.003473133333334D-6*TN(3) +1.548792108333334D-10*TN(4) 
     *         -5.3761131D-15*TN(5) 
      ELSE
      SMH(39) = 4.11878981D1 +1.31074559D4*TI 
     *         -1.89226915D0*TN(1) +4.730331785D-2*TN(2) 
     *         -1.045642535D-5*TN(3) +1.792985908333333D-9*TN(4) 
     *         -1.513593415D-13*TN(5) 
      ENDIF
C PXC8H17
      IF (T .GT. 1390) THEN
      SMH(40) = -1.02557384D2 +1.33300535D4*TI 
     *         +2.50510356D1*TN(1) +1.84740081D-2*TN(2) 
     *         -2.096087733333333D-6*TN(3) +1.621903408333333D-10*TN(4) 
     *         -5.6334449D-15*TN(5) 
      ELSE
      SMH(40) = 3.76130631D1 +3.92689511D3*TI 
     *         -7.72759438D-1*TN(1) +4.662748525D-2*TN(2) 
     *         -9.740787416666667D-6*TN(3) +1.54641845D-9*TN(4) 
     *         -1.185637415D-13*TN(5) 
      ENDIF
C C9H18
      IF (T .GT. 1392) THEN
      SMH(41) = -1.16618623D2 +2.65709061D4*TI 
     *         +2.76142176D1*TN(1) +1.984126435D-2*TN(2) 
     *         -2.246990766666667D-6*TN(3) +1.7365871D-10*TN(4) 
     *         -6.0269647D-15*TN(5) 
      ELSE
      SMH(41) = 4.41245128D1 +1.59890847D4*TI 
     *         -2.16108263D0*TN(1) +5.34791485D-2*TN(2) 
     *         -1.184955406666667D-5*TN(3) +2.033092308333333D-9*TN(4) 
     *         -1.713857735D-13*TN(5) 
      ENDIF
C PXC9H19
      IF (T .GT. 1390) THEN
      SMH(42) = -1.16837897D2 +1.7451603D4*TI 
     *         +2.83097514D1*TN(1) +2.06328672D-2*TN(2) 
     *         -2.339721483333333D-6*TN(3) +1.809790591666667D-10*TN(4) 
     *         -6.284615349999999D-15*TN(5) 
      ELSE
      SMH(42) = 4.23518992D1 +6.80818512D3*TI 
     *         -1.04387292D0*TN(1) +5.28086415D-2*TN(2) 
     *         -1.113666618333333D-5*TN(3) +1.787384716666667D-9*TN(4) 
     *         -1.387021375D-13*TN(5) 
      ENDIF
C C10H20
      IF (T .GT. 1392) THEN
      SMH(43) = -1.32705172D2 +3.06937307D4*TI 
     *         +3.08753903D1*TN(1) +2.19985763D-2*TN(2) 
     *         -2.4904255D-6*TN(3) +1.924313983333333D-10*TN(4) 
     *         -6.67757385D-15*TN(5) 
      ELSE
      SMH(43) = 4.70571383D1 +1.88708365D4*TI 
     *         -2.42901688D0*TN(1) +5.9652799D-2*TN(2) 
     *         -1.324148375D-5*TN(3) +2.272804966666667D-9*TN(4) 
     *         -1.913591865D-13*TN(5) 
      ENDIF
C PXC10H21
      IF (T .GT. 1390) THEN
      SMH(44) = -1.34708986D2 +2.15737832D4*TI 
     *         +3.1569716D1*TN(1) +2.279092015D-2*TN(2) 
     *         -2.583249416666667D-6*TN(3) +1.997591108333333D-10*TN(4) 
     *         -6.93547795D-15*TN(5) 
      ELSE
      SMH(44) = 4.35010452D1 +9.689675499999999D3*TI 
     *         -1.31358348D0*TN(1) +5.89864065D-2*TN(2) 
     *         -1.253071798333333D-5*TN(3) +2.027759216666667D-9*TN(4) 
     *         -1.58761426D-13*TN(5) 
      ENDIF
C PXC12H25
      IF (T .GT. 1390) THEN
      SMH(45) = -1.66882734D2 +2.98194375D4*TI 
     *         +3.80921885D1*TN(1) +2.71053924D-2*TN(2) 
     *         -3.07009195D-6*TN(3) +2.373018108333333D-10*TN(4) 
     *         -8.236587400000001D-15*TN(5) 
      ELSE
      SMH(45) = 4.93702421D1 +1.54530435D4*TI 
     *         -1.85028741D0*TN(1) +7.1335354D-2*TN(2) 
     *         -1.531527591666667D-5*TN(3) +2.5073616D-9*TN(4) 
     *         -1.9872715D-13*TN(5) 
      ENDIF
C SXC12H25
      IF (T .GT. 1385) THEN
      SMH(46) = -1.65805933D2 +3.12144988D4*TI 
     *         +3.79688268D1*TN(1) +2.69359732D-2*TN(2) 
     *         -3.036187716666667D-6*TN(3) +2.339787525D-10*TN(4) 
     *         -8.105421000000001D-15*TN(5) 
      ELSE
      SMH(46) = 4.83521895D1 +1.67660539D4*TI 
     *         -1.36787089D0*TN(1) +6.867767399999999D-2*TN(2) 
     *         -1.373460263333333D-5*TN(3) +1.970179683333333D-9*TN(4) 
     *         -1.23717966D-13*TN(5) 
      ENDIF
C S3XC12H25
      IF (T .GT. 1385) THEN
      SMH(47) = -1.65805933D2 +3.12144988D4*TI 
     *         +3.79688268D1*TN(1) +2.69359732D-2*TN(2) 
     *         -3.036187716666667D-6*TN(3) +2.339787525D-10*TN(4) 
     *         -8.105421000000001D-15*TN(5) 
      ELSE
      SMH(47) = 4.83521895D1 +1.67660539D4*TI 
     *         -1.36787089D0*TN(1) +6.867767399999999D-2*TN(2) 
     *         -1.373460263333333D-5*TN(3) +1.970179683333333D-9*TN(4) 
     *         -1.23717966D-13*TN(5) 
      ENDIF
C C12H24
      IF (T .GT. 1391) THEN
      SMH(48) = -1.64892663D2 +3.89405962D4*TI 
     *         +3.74002111D1*TN(1) +2.631153765D-2*TN(2) 
     *         -2.977071983333333D-6*TN(3) +2.299582191666667D-10*TN(4) 
     *         -7.97812495D-15*TN(5) 
      ELSE
      SMH(48) = 5.2915887D1 +2.46345299D4*TI 
     *         -2.96342681D0*TN(1) +7.199618000000001D-2*TN(2) 
     *         -1.602306691666667D-5*TN(3) +2.751453941666667D-9*TN(4) 
     *         -2.31199095D-13*TN(5) 
      ENDIF
C O2C12H24OOH
      IF (T .GT. 1000) THEN
      SMH(51) = -1.3775D2 +5.12675D4*TI 
     *         +3.50907D1*TN(1) +2.55295D-2*TN(2) 
     *         -2.572416666666666D-6*TN(3) +1.871891666666667D-10*TN(4) 
     *         -6.44505D-15*TN(5) 
      ELSE
      SMH(51) = 4.13429D1 +4.16875D4*TI 
     *         +4.81972D-1*TN(1) +7.251000000000001D-2*TN(2) 
     *         -1.665513333333333D-5*TN(3) +2.170183333333333D-9*TN(4) 
     *         +5.967900000000001D-14*TN(5) 
      ENDIF
C OC12H23OOH
      IF (T .GT. 1000) THEN
      SMH(52) = -7.77662D1 +7.18258D4*TI 
     *         +2.36731D1*TN(1) +3.08196D-2*TN(2) 
     *         -3.497266666666667D-6*TN(3) +2.776383333333333D-10*TN(4) 
     *         -1.01795D-14*TN(5) 
      ELSE
      SMH(52) = 6.84155D0 +6.653610000000001D4*TI 
     *         +8.80733D0*TN(1) +3.253115D-2*TN(2) 
     *         +1.15843D-5*TN(3) -1.057541666666667D-8*TN(4) 
     *         +2.554955D-12*TN(5) 
      ENDIF
C
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE RATX (T, C, RF, RB, RKLOW)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
      PARAMETER (SMALL = 1D-200)
      DIMENSION C(*), RF(*), RB(*), RKLOW(*)
C
      ALOGT = LOG(T)
      CTOT = 0.0
      DO K = 1, 35
         CTOT = CTOT + C(K)
      ENDDO
C
C     R1: H + O2 = O + OH
      RF(1) = RF(1)*C(2)*C(9)
      RB(1) = RB(1)*C(3)*C(4)
C     R2: O + H2 = H + OH
      RF(2) = RF(2)*C(3)*C(6)
      RB(2) = RB(2)*C(2)*C(4)
C     R3: OH + H2 = H + H2O
      RF(3) = RF(3)*C(4)*C(6)
      RB(3) = RB(3)*C(2)*C(7)
C     R4: 2OH = O + H2O
      RF(4) = RF(4)*C(4)*C(4)
      RB(4) = RB(4)*C(3)*C(7)
C     R5: H + O2 = HO2
      CTB = CTOT-1.5D-1*C(9)+1.089D1*C(7)
     * +9.000000000000008D-2*C(13)+1.18D0*C(14)
      PR = RKLOW(1) * CTB / RF(5)
      PCOR = PR / (1.0 + PR)
      PRLOG = LOG10(MAX(PR,SMALL))
      FC = (PRLOG -1.983099029051326D-1)
     *     /(1.132308094493256D0 -0.14D0*(PRLOG -1.983099029051326D-1))
      FC = EXP(-6.931471805599453D-1 /(1.0D0 + FC*FC))
      PCOR = FC * PCOR
      RF(5) = RF(5) * PCOR
      RB(5) = RB(5) * PCOR
      RF(5) = RF(5)*C(2)*C(9)
      RB(5) = RB(5)*C(5)
C     R6: H + HO2 = 2OH
      RF(6) = RF(6)*C(2)*C(5)
      RB(6) = RB(6)*C(4)*C(4)
C     R7: H2 + O2 = H + HO2
      RF(7) = RF(7)*C(6)*C(9)
      RB(7) = RB(7)*C(2)*C(5)
C     R8: OH + HO2 = H2O + O2
      RF(8) = RF(8)*C(4)*C(5)
      RB(8) = RB(8)*C(7)*C(9)
C     R9: H + HO2 = O + H2O
      RF(9) = RF(9)*C(2)*C(5)
      RB(9) = RB(9)*C(3)*C(7)
C     R10: O + HO2 = OH + O2
      RF(10) = RF(10)*C(3)*C(5)
      RB(10) = RB(10)*C(4)*C(9)
C     R11: 2HO2 = H2O2 + O2
      RF(11) = RF(11)*C(5)*C(5)
      RB(11) = RB(11)*C(8)*C(9)
C     R12: 2HO2 = H2O2 + O2
      RF(12) = RF(12)*C(5)*C(5)
      RB(12) = RB(12)*C(8)*C(9)
C     R13: H + H2O2 = OH + H2O
      RF(13) = RF(13)*C(2)*C(8)
      RB(13) = RB(13)*C(4)*C(7)
C     R14: H + H2O2 = HO2 + H2
      RF(14) = RF(14)*C(2)*C(8)
      RB(14) = RB(14)*C(5)*C(6)
C     R15: O + H2O2 = OH + HO2
      RF(15) = RF(15)*C(3)*C(8)
      RB(15) = RB(15)*C(4)*C(5)
C     R16: OH + H2O2 = HO2 + H2O
      RF(16) = RF(16)*C(4)*C(8)
      RB(16) = RB(16)*C(5)*C(7)
C     R17: OH + H2O2 = HO2 + H2O
      RF(17) = RF(17)*C(4)*C(8)
      RB(17) = RB(17)*C(5)*C(7)
C     R18: 2OH = H2O2
      CTB = CTOT+C(6)+5.D0*C(7)
     * +7.5D-1*C(13)+2.6D0*C(14)
      PR = RKLOW(2) * CTB / RF(18)
      PCOR = PR / (1.0 + PR)
      PRLOG = LOG10(MAX(PR,SMALL))
      FCENT = 2.654D-1*EXP(-T/9.4D1)
     *      + 7.346D-1*EXP(-T/1.756D3)
     *     + EXP(-5.182D3/T)
      FCLOG = LOG10(MAX(FCENT,SMALL))
      XN    = 0.75 - 1.27*FCLOG
      CPRLOG= PRLOG - (0.4 + 0.67*FCLOG)
      FLOG = FCLOG/(1.0 + (CPRLOG/(XN-0.14*CPRLOG))**2)
      FC = 10.0**FLOG
      PCOR = FC * PCOR
      RF(18) = RF(18) * PCOR
      RB(18) = RB(18) * PCOR
      RF(18) = RF(18)*C(4)*C(4)
      RB(18) = RB(18)*C(8)
C     R19: 2H = H2
      CTB = CTOT-C(6)-C(7)
     * -C(14)
      RF(19) = RF(19)*CTB*C(2)*C(2)
      RB(19) = RB(19)*CTB*C(6)
C     R20: H + OH = H2O
      CTB = CTOT+C(6)+5.3D0*C(7)
     * +7.5D-1*C(13)+2.6D0*C(14)
      RF(20) = RF(20)*CTB*C(2)*C(4)
      RB(20) = RB(20)*CTB*C(7)
C     R21: 2O = O2
      CTB = CTOT+1.4D0*C(6)+1.44D1*C(7)
     * +7.5D-1*C(13)+2.6D0*C(14)
      RF(21) = RF(21)*CTB*C(3)*C(3)
      RB(21) = RB(21)*CTB*C(9)
C     R22: 2H = H2
      RF(22) = RF(22)*C(2)*C(2)*C(6)
      RB(22) = RB(22)*C(6)*C(6)
C     R23: 2H = H2
      RF(23) = RF(23)*C(2)*C(2)*C(7)
      RB(23) = RB(23)*C(6)*C(7)
C     R24: 2H = H2
      RF(24) = RF(24)*C(2)*C(2)*C(14)
      RB(24) = RB(24)*C(6)*C(14)
C     R25: H + O = OH
      CTB = CTOT+C(6)+1.1D1*C(7)
     * +7.5D-1*C(13)+2.6D0*C(14)
      RF(25) = RF(25)*CTB*C(2)*C(3)
      RB(25) = RB(25)*CTB*C(4)
C     R26: OH + CO = H + CO2
      RF(26) = RF(26)*C(4)*C(13)
      RB(26) = RB(26)*C(2)*C(14)
C     R27: OH + CO = H + CO2
      RF(27) = RF(27)*C(4)*C(13)
      RB(27) = RB(27)*C(2)*C(14)
C     R28: HO2 + CO = OH + CO2
      RF(28) = RF(28)*C(5)*C(13)
      RB(28) = RB(28)*C(4)*C(14)
C     R29: O + CO = CO2
      CTB = CTOT+C(6)+1.1D1*C(7)
     * +7.5D-1*C(13)+2.6D0*C(14)
      PR = RKLOW(3) * CTB / RF(29)
      PCOR = PR / (1.0 + PR)
      RF(29) = RF(29) * PCOR
      RB(29) = RB(29) * PCOR
      RF(29) = RF(29)*C(3)*C(13)
      RB(29) = RB(29)*C(14)
C     R30: O2 + CO = O + CO2
      RF(30) = RF(30)*C(9)*C(13)
      RB(30) = RB(30)*C(3)*C(14)
C     R31: HCO = H + CO
      CTB = CTOT+C(6)-C(7)
     * +7.5D-1*C(13)+2.6D0*C(14)
      RF(31) = RF(31)*CTB
      RB(31) = RB(31)*CTB*C(2)*C(13)
C     R32: H + HCO = H2 + CO
      RF(32) = RF(32)*C(2)
      RB(32) = RB(32)*C(6)*C(13)
C     R33: O + HCO = OH + CO
      RF(33) = RF(33)*C(3)
      RB(33) = RB(33)*C(4)*C(13)
C     R34: O + HCO = H + CO2
      RF(34) = RF(34)*C(3)
      RB(34) = RB(34)*C(2)*C(14)
C     R35: OH + HCO = H2O + CO
      RF(35) = RF(35)*C(4)
      RB(35) = RB(35)*C(7)*C(13)
C     R36: O2 + HCO = HO2 + CO
      RF(36) = RF(36)*C(9)
      RB(36) = RB(36)*C(5)*C(13)
C     R37: HCO = H + CO
      RF(37) = RF(37)*C(7)
      RB(37) = RB(37)*C(2)*C(7)*C(13)
C     R38: H2 + CO = CH2O
      CTB = CTOT+C(6)+5.D0*C(7)
     * +C(11)+5.D-1*C(13)
     * +C(14)+2.D0*C(17)
      PR = RKLOW(4) * CTB / RF(38)
      PCOR = PR / (1.0 + PR)
      PRLOG = LOG10(MAX(PR,SMALL))
      FCENT = 6.799999999999995D-2*EXP(-T/1.97D2)
     *      + 9.320000000000001D-1*EXP(-T/1.54D3)
     *     + EXP(-1.03D4/T)
      FCLOG = LOG10(MAX(FCENT,SMALL))
      XN    = 0.75 - 1.27*FCLOG
      CPRLOG= PRLOG - (0.4 + 0.67*FCLOG)
      FLOG = FCLOG/(1.0 + (CPRLOG/(XN-0.14*CPRLOG))**2)
      FC = 10.0**FLOG
      PCOR = FC * PCOR
      RF(38) = RF(38) * PCOR
      RB(38) = RB(38) * PCOR
      RF(38) = RF(38)*C(6)*C(13)
      RB(38) = RB(38)*C(12)
C     R39: H + HCO = CH2O
      CTB = CTOT+C(6)+5.D0*C(7)
     * +C(11)+5.D-1*C(13)
     * +C(14)+2.D0*C(17)
      PR = RKLOW(5) * CTB / RF(39)
      PCOR = PR / (1.0 + PR)
      PRLOG = LOG10(MAX(PR,SMALL))
      FCENT = 2.176D-1*EXP(-T/2.71D2)
     *      + 7.824D-1*EXP(-T/2.755D3)
     *     + EXP(-6.57D3/T)
      FCLOG = LOG10(MAX(FCENT,SMALL))
      XN    = 0.75 - 1.27*FCLOG
      CPRLOG= PRLOG - (0.4 + 0.67*FCLOG)
      FLOG = FCLOG/(1.0 + (CPRLOG/(XN-0.14*CPRLOG))**2)
      FC = 10.0**FLOG
      PCOR = FC * PCOR
      RF(39) = RF(39) * PCOR
      RB(39) = RB(39) * PCOR
      RF(39) = RF(39)*C(2)
      RB(39) = RB(39)*C(12)
C     R40: H + CH2 = CH3
      CTB = CTOT+C(6)+5.D0*C(7)
     * +C(11)+5.D-1*C(13)
     * +C(14)+2.D0*C(17)
      PR = RKLOW(6) * CTB / RF(40)
      PCOR = PR / (1.0 + PR)
      PRLOG = LOG10(MAX(PR,SMALL))
      FCENT = 3.2D-1*EXP(-T/7.8D1)
     *      + 6.800000000000001D-1*EXP(-T/1.995D3)
     *     + EXP(-5.59D3/T)
      FCLOG = LOG10(MAX(FCENT,SMALL))
      XN    = 0.75 - 1.27*FCLOG
      CPRLOG= PRLOG - (0.4 + 0.67*FCLOG)
      FLOG = FCLOG/(1.0 + (CPRLOG/(XN-0.14*CPRLOG))**2)
      FC = 10.0**FLOG
      PCOR = FC * PCOR
      RF(40) = RF(40) * PCOR
      RB(40) = RB(40) * PCOR
      RF(40) = RF(40)*C(2)
      RB(40) = RB(40)*C(10)
C     R41: O + CH2 = H + HCO
      RF(41) = RF(41)*C(3)
      RB(41) = RB(41)*C(2)
C     R42: OH + CH2 = H + CH2O
      RF(42) = RF(42)*C(4)
      RB(42) = RB(42)*C(2)*C(12)
C     R43: H2 + CH2 = H + CH3
      RF(43) = RF(43)*C(6)
      RB(43) = RB(43)*C(2)*C(10)
C     R44: O2 + CH2 = OH + HCO
      RF(44) = RF(44)*C(9)
      RB(44) = RB(44)*C(4)
C     R45: O2 + CH2 = 2H + CO2
      RF(45) = RF(45)*C(9)
      RB(45) = RB(45)*C(2)*C(2)*C(14)
C     R46: HO2 + CH2 = OH + CH2O
      RF(46) = RF(46)*C(5)
      RB(46) = RB(46)*C(4)*C(12)
C     R47: 2CH2 = H2 + C2H2
      RB(47) = RB(47)*C(6)*C(15)
C     R48: CH2* = CH2
      RF(48) = RF(48)*C(35)
      RB(48) = RB(48)*C(35)
C     R49: O + CH2* = H2 + CO
      RF(49) = RF(49)*C(3)
      RB(49) = RB(49)*C(6)*C(13)
C     R50: O + CH2* = H + HCO
      RF(50) = RF(50)*C(3)
      RB(50) = RB(50)*C(2)
C     R51: OH + CH2* = H + CH2O
      RF(51) = RF(51)*C(4)
      RB(51) = RB(51)*C(2)*C(12)
C     R52: H2 + CH2* = H + CH3
      RF(52) = RF(52)*C(6)
      RB(52) = RB(52)*C(2)*C(10)
C     R53: O2 + CH2* = H + OH + CO
      RF(53) = RF(53)*C(9)
      RB(53) = RB(53)*C(2)*C(4)*C(13)
C     R54: O2 + CH2* = H2O + CO
      RF(54) = RF(54)*C(9)
      RB(54) = RB(54)*C(7)*C(13)
C     R55: CH2* = CH2
      RF(55) = RF(55)*C(7)
      RB(55) = RB(55)*C(7)
C     R56: CH2* = CH2
      RF(56) = RF(56)*C(13)
      RB(56) = RB(56)*C(13)
C     R57: CH2* = CH2
      RF(57) = RF(57)*C(14)
      RB(57) = RB(57)*C(14)
C     R58: CH2* + CO2 = CH2O + CO
      RF(58) = RF(58)*C(14)
      RB(58) = RB(58)*C(12)*C(13)
C     R59: H + CH2O = CH3O
      CTB = CTOT+C(6)+5.D0*C(7)
     * +C(11)+5.D-1*C(13)
     * +C(14)+2.D0*C(17)
      PR = RKLOW(7) * CTB / RF(59)
      PCOR = PR / (1.0 + PR)
      PRLOG = LOG10(MAX(PR,SMALL))
      FCENT = 2.42D-1*EXP(-T/9.4D1)
     *      + 7.58D-1*EXP(-T/1.555D3)
     *     + EXP(-4.2D3/T)
      FCLOG = LOG10(MAX(FCENT,SMALL))
      XN    = 0.75 - 1.27*FCLOG
      CPRLOG= PRLOG - (0.4 + 0.67*FCLOG)
      FLOG = FCLOG/(1.0 + (CPRLOG/(XN-0.14*CPRLOG))**2)
      FC = 10.0**FLOG
      PCOR = FC * PCOR
      RF(59) = RF(59) * PCOR
      RB(59) = RB(59) * PCOR
      RF(59) = RF(59)*C(2)*C(12)
C     R60: H + CH2O = H2 + HCO
      RF(60) = RF(60)*C(2)*C(12)
      RB(60) = RB(60)*C(6)
C     R61: O + CH2O = OH + HCO
      RF(61) = RF(61)*C(3)*C(12)
      RB(61) = RB(61)*C(4)
C     R62: OH + CH2O = H2O + HCO
      RF(62) = RF(62)*C(4)*C(12)
      RB(62) = RB(62)*C(7)
C     R63: O2 + CH2O = HO2 + HCO
      RF(63) = RF(63)*C(9)*C(12)
      RB(63) = RB(63)*C(5)
C     R64: HO2 + CH2O = H2O2 + HCO
      RF(64) = RF(64)*C(5)*C(12)
      RB(64) = RB(64)*C(8)
C     R65: H + CH3 = CH4
      CTB = CTOT+C(6)+5.D0*C(7)
     * +C(11)+5.D-1*C(13)
     * +C(14)+2.D0*C(17)
      PR = RKLOW(8) * CTB / RF(65)
      PCOR = PR / (1.0 + PR)
      PRLOG = LOG10(MAX(PR,SMALL))
      FCENT = 2.17D-1*EXP(-T/7.4D1)
     *      + 7.83D-1*EXP(-T/2.941D3)
     *     + EXP(-6.964D3/T)
      FCLOG = LOG10(MAX(FCENT,SMALL))
      XN    = 0.75 - 1.27*FCLOG
      CPRLOG= PRLOG - (0.4 + 0.67*FCLOG)
      FLOG = FCLOG/(1.0 + (CPRLOG/(XN-0.14*CPRLOG))**2)
      FC = 10.0**FLOG
      PCOR = FC * PCOR
      RF(65) = RF(65) * PCOR
      RB(65) = RB(65) * PCOR
      RF(65) = RF(65)*C(2)*C(10)
      RB(65) = RB(65)*C(11)
C     R66: O + CH3 = H + CH2O
      RF(66) = RF(66)*C(3)*C(10)
      RB(66) = RB(66)*C(2)*C(12)
C     R67: OH + CH3 = H2O + CH2
      RF(67) = RF(67)*C(4)*C(10)
      RB(67) = RB(67)*C(7)
C     R68: OH + CH3 = H2O + CH2*
      RF(68) = RF(68)*C(4)*C(10)
      RB(68) = RB(68)*C(7)
C     R69: O2 + CH3 = O + CH3O
      RF(69) = RF(69)*C(9)*C(10)
      RB(69) = RB(69)*C(3)
C     R70: O2 + CH3 = OH + CH2O
      RF(70) = RF(70)*C(9)*C(10)
      RB(70) = RB(70)*C(4)*C(12)
C     R71: HO2 + CH3 = O2 + CH4
      RF(71) = RF(71)*C(5)*C(10)
      RB(71) = RB(71)*C(9)*C(11)
C     R72: HO2 + CH3 = OH + CH3O
      RF(72) = RF(72)*C(5)*C(10)
      RB(72) = RB(72)*C(4)
C     R73: H2O2 + CH3 = HO2 + CH4
      RF(73) = RF(73)*C(8)*C(10)
      RB(73) = RB(73)*C(5)*C(11)
C     R74: CH3 + HCO = CH4 + CO
      RF(74) = RF(74)*C(10)
      RB(74) = RB(74)*C(11)*C(13)
C     R75: CH3 + CH2O = CH4 + HCO
      RF(75) = RF(75)*C(10)*C(12)
      RB(75) = RB(75)*C(11)
C     R76: CH2 + CH3 = H + C2H4
      RF(76) = RF(76)*C(10)
      RB(76) = RB(76)*C(2)*C(16)
C     R77: CH2* + CH3 = H + C2H4
      RF(77) = RF(77)*C(10)
      RB(77) = RB(77)*C(2)*C(16)
C     R78: 2CH3 = C2H6
      CTB = CTOT+C(6)+5.D0*C(7)
     * +C(11)+5.D-1*C(13)
     * +C(14)+2.D0*C(17)
      PR = RKLOW(9) * CTB / RF(78)
      PCOR = PR / (1.0 + PR)
      PRLOG = LOG10(MAX(PR,SMALL))
      FCENT = 4.675D-1*EXP(-T/1.51D2)
     *      + 5.325D-1*EXP(-T/1.038D3)
     *     + EXP(-4.97D3/T)
      FCLOG = LOG10(MAX(FCENT,SMALL))
      XN    = 0.75 - 1.27*FCLOG
      CPRLOG= PRLOG - (0.4 + 0.67*FCLOG)
      FLOG = FCLOG/(1.0 + (CPRLOG/(XN-0.14*CPRLOG))**2)
      FC = 10.0**FLOG
      PCOR = FC * PCOR
      RF(78) = RF(78) * PCOR
      RB(78) = RB(78) * PCOR
      RF(78) = RF(78)*C(10)*C(10)
      RB(78) = RB(78)*C(17)
C     R79: 2CH3 = H + C2H5
      RF(79) = RF(79)*C(10)*C(10)
      RB(79) = RB(79)*C(2)
C     R80: H + CH3O = H2 + CH2O
      RF(80) = RF(80)*C(2)
      RB(80) = RB(80)*C(6)*C(12)
C     R81: H + CH3O = OH + CH3
      RF(81) = RF(81)*C(2)
      RB(81) = RB(81)*C(4)*C(10)
C     R82: H + CH3O = H2O + CH2*
      RF(82) = RF(82)*C(2)
      RB(82) = RB(82)*C(7)
C     R83: O + CH3O = OH + CH2O
      RF(83) = RF(83)*C(3)
      RB(83) = RB(83)*C(4)*C(12)
C     R84: OH + CH3O = H2O + CH2O
      RF(84) = RF(84)*C(4)
      RB(84) = RB(84)*C(7)*C(12)
C     R85: O2 + CH3O = HO2 + CH2O
      RF(85) = RF(85)*C(9)
      RB(85) = RB(85)*C(5)*C(12)
C     R86: H + CH4 = H2 + CH3
      RF(86) = RF(86)*C(2)*C(11)
      RB(86) = RB(86)*C(6)*C(10)
C     R87: O + CH4 = OH + CH3
      RF(87) = RF(87)*C(3)*C(11)
      RB(87) = RB(87)*C(4)*C(10)
C     R88: OH + CH4 = H2O + CH3
      RF(88) = RF(88)*C(4)*C(11)
      RB(88) = RB(88)*C(7)*C(10)
C     R89: CH2 + CH4 = 2CH3
      RF(89) = RF(89)*C(11)
      RB(89) = RB(89)*C(10)*C(10)
C     R90: CH2* + CH4 = 2CH3
      RF(90) = RF(90)*C(11)
      RB(90) = RB(90)*C(10)*C(10)
C     R91: C2H3 = H + C2H2
      CTB = CTOT+C(6)+5.D0*C(7)
     * +C(11)+5.D-1*C(13)
     * +C(14)+2.D0*C(17)
     * +2.D0*C(15)+2.D0*C(16)
      PR = RKLOW(10) * CTB / RF(91)
      PCOR = PR / (1.0 + PR)
      PRLOG = LOG10(MAX(PR,SMALL))
      FCENT = -9.816D-1*EXP(-T/5.3837D3)
     *      + 1.9816D0*EXP(-T/4.2932D0)
     *     + EXP(7.95D-2/T)
      FCLOG = LOG10(MAX(FCENT,SMALL))
      XN    = 0.75 - 1.27*FCLOG
      CPRLOG= PRLOG - (0.4 + 0.67*FCLOG)
      FLOG = FCLOG/(1.0 + (CPRLOG/(XN-0.14*CPRLOG))**2)
      FC = 10.0**FLOG
      PCOR = FC * PCOR
      RF(91) = RF(91) * PCOR
      RB(91) = RB(91) * PCOR
      RB(91) = RB(91)*C(2)*C(15)
C     R92: O + C2H2 = CH2 + CO
      RF(92) = RF(92)*C(3)*C(15)
      RB(92) = RB(92)*C(13)
C     R93: OH + C2H2 = CH3 + CO
      RF(93) = RF(93)*C(4)*C(15)
      RB(93) = RB(93)*C(10)*C(13)
C     R94: HCO + C2H2 = CO + C2H3
      RF(94) = RF(94)*C(15)
      RB(94) = RB(94)*C(13)
C     R95: CH3 + C2H2 = aC3H5
      RF(95) = RF(95)*C(10)*C(15)
      RB(95) = RB(95)*C(19)
C     R96: H + C2H3 = C2H4
      CTB = CTOT+C(6)+5.D0*C(7)
     * +C(11)+5.D-1*C(13)
     * +C(14)+2.D0*C(17)
     * +2.D0*C(15)+2.D0*C(16)
      PR = RKLOW(11) * CTB / RF(96)
      PCOR = PR / (1.0 + PR)
      PRLOG = LOG10(MAX(PR,SMALL))
      FCENT = 2.18D-1*EXP(-T/2.075D2)
     *      + 7.82D-1*EXP(-T/2.663D3)
     *     + EXP(-6.095D3/T)
      FCLOG = LOG10(MAX(FCENT,SMALL))
      XN    = 0.75 - 1.27*FCLOG
      CPRLOG= PRLOG - (0.4 + 0.67*FCLOG)
      FLOG = FCLOG/(1.0 + (CPRLOG/(XN-0.14*CPRLOG))**2)
      FC = 10.0**FLOG
      PCOR = FC * PCOR
      RF(96) = RF(96) * PCOR
      RB(96) = RB(96) * PCOR
      RF(96) = RF(96)*C(2)
      RB(96) = RB(96)*C(16)
C     R97: H + C2H3 = H2 + C2H2
      RF(97) = RF(97)*C(2)
      RB(97) = RB(97)*C(6)*C(15)
C     R98: O + C2H3 = CH3 + CO
      RF(98) = RF(98)*C(3)
      RB(98) = RB(98)*C(10)*C(13)
C     R99: OH + C2H3 = H2O + C2H2
      RF(99) = RF(99)*C(4)
      RB(99) = RB(99)*C(7)*C(15)
C     R100: O2 + C2H3 = HO2 + C2H2
      RF(100) = RF(100)*C(9)
      RB(100) = RB(100)*C(5)*C(15)
C     R101: O2 + C2H3 = O + CH2CHO
      RF(101) = RF(101)*C(9)
      RB(101) = RB(101)*C(3)*C(18)
C     R102: O2 + C2H3 = HCO + CH2O
      RF(102) = RF(102)*C(9)
      RB(102) = RB(102)*C(12)
C     R103: HO2 + C2H3 = OH + CH2CHO
      RF(103) = RF(103)*C(5)
      RB(103) = RB(103)*C(4)*C(18)
C     R104: H2O2 + C2H3 = HO2 + C2H4
      RF(104) = RF(104)*C(8)
      RB(104) = RB(104)*C(5)*C(16)
C     R105: HCO + C2H3 = CO + C2H4
      RB(105) = RB(105)*C(13)*C(16)
C     R106: HCO + C2H3 = C2H3CHO
      RB(106) = RB(106)*C(21)
C     R107: CH3 + C2H3 = CH4 + C2H2
      RF(107) = RF(107)*C(10)
      RB(107) = RB(107)*C(11)*C(15)
C     R108: CH3 + C2H3 = C3H6
      CTB = CTOT+C(6)+5.D0*C(7)
     * +C(11)+5.D-1*C(13)
     * +C(14)+2.D0*C(17)
     * +2.D0*C(16)
      PR = RKLOW(12) * CTB / RF(108)
      PCOR = PR / (1.0 + PR)
      PRLOG = LOG10(MAX(PR,SMALL))
      FCENT = 8.25D-1*EXP(-T/1.3406D3)
     *      + 1.75D-1*EXP(-T/6.D4)
     *     + EXP(-1.01398D4/T)
      FCLOG = LOG10(MAX(FCENT,SMALL))
      XN    = 0.75 - 1.27*FCLOG
      CPRLOG= PRLOG - (0.4 + 0.67*FCLOG)
      FLOG = FCLOG/(1.0 + (CPRLOG/(XN-0.14*CPRLOG))**2)
      FC = 10.0**FLOG
      PCOR = FC * PCOR
      RF(108) = RF(108) * PCOR
      RB(108) = RB(108) * PCOR
      RF(108) = RF(108)*C(10)
      RB(108) = RB(108)*C(20)
C     R109: CH3 + C2H3 = H + aC3H5
      RF(109) = RF(109)*C(10)
      RB(109) = RB(109)*C(2)*C(19)
C     R110: 2C2H3 = C2H2 + C2H4
      RB(110) = RB(110)*C(15)*C(16)
C     R111: CH2CHO = CH3 + CO
      RF(111) = RF(111)*C(18)
      RB(111) = RB(111)*C(10)*C(13)
C     R112: H + CH2CHO = CH3 + HCO
      RF(112) = RF(112)*C(2)*C(18)
      RB(112) = RB(112)*C(10)
C     R113: O2 + CH2CHO = OH + CH2O + CO
      RF(113) = RF(113)*C(9)*C(18)
      RB(113) = RB(113)*C(4)*C(12)*C(13)
C     R114: H + C2H4 = C2H5
      CTB = CTOT+C(6)+5.D0*C(7)
     * +C(11)+5.D-1*C(13)
     * +C(14)+2.D0*C(17)
      PR = RKLOW(13) * CTB / RF(114)
      PCOR = PR / (1.0 + PR)
      PRLOG = LOG10(MAX(PR,SMALL))
      FCENT = 2.4D-1*EXP(-T/4.D1)
     *      + 7.6D-1*EXP(-T/1.025D3)
      FCLOG = LOG10(MAX(FCENT,SMALL))
      XN    = 0.75 - 1.27*FCLOG
      CPRLOG= PRLOG - (0.4 + 0.67*FCLOG)
      FLOG = FCLOG/(1.0 + (CPRLOG/(XN-0.14*CPRLOG))**2)
      FC = 10.0**FLOG
      PCOR = FC * PCOR
      RF(114) = RF(114) * PCOR
      RB(114) = RB(114) * PCOR
      RF(114) = RF(114)*C(2)*C(16)
C     R115: H + C2H4 = H2 + C2H3
      RF(115) = RF(115)*C(2)*C(16)
      RB(115) = RB(115)*C(6)
C     R116: O + C2H4 = OH + C2H3
      RF(116) = RF(116)*C(3)*C(16)
      RB(116) = RB(116)*C(4)
C     R117: O + C2H4 = CH3 + HCO
      RF(117) = RF(117)*C(3)*C(16)
      RB(117) = RB(117)*C(10)
C     R118: O + C2H4 = CH2 + CH2O
      RF(118) = RF(118)*C(3)*C(16)
      RB(118) = RB(118)*C(12)
C     R119: OH + C2H4 = H2O + C2H3
      RF(119) = RF(119)*C(4)*C(16)
      RB(119) = RB(119)*C(7)
C     R120: HCO + C2H4 = CO + C2H5
      RF(120) = RF(120)*C(16)
      RB(120) = RB(120)*C(13)
C     R121: CH2 + C2H4 = H + aC3H5
      RF(121) = RF(121)*C(16)
      RB(121) = RB(121)*C(2)*C(19)
C     R122: CH2* + C2H4 = H + aC3H5
      RF(122) = RF(122)*C(16)
      RB(122) = RB(122)*C(2)*C(19)
C     R123: CH3 + C2H4 = CH4 + C2H3
      RF(123) = RF(123)*C(10)*C(16)
      RB(123) = RB(123)*C(11)
C     R124: nC3H7 = CH3 + C2H4
      RB(124) = RB(124)*C(10)*C(16)
C     R125: O2 + C2H4 = HO2 + C2H3
      RF(125) = RF(125)*C(9)*C(16)
      RB(125) = RB(125)*C(5)
C     R126: C2H3 + C2H4 = C4H7
      RF(126) = RF(126)*C(16)
      RB(126) = RB(126)*C(22)
C     R127: H + C2H5 = C2H6
      CTB = CTOT+C(6)+5.D0*C(7)
     * +C(11)+5.D-1*C(13)
     * +C(14)+2.D0*C(17)
      PR = RKLOW(14) * CTB / RF(127)
      PCOR = PR / (1.0 + PR)
      PRLOG = LOG10(MAX(PR,SMALL))
      FCENT = 1.578000000000001D-1*EXP(-T/1.25D2)
     *      + 8.422D-1*EXP(-T/2.219D3)
     *     + EXP(-6.882D3/T)
      FCLOG = LOG10(MAX(FCENT,SMALL))
      XN    = 0.75 - 1.27*FCLOG
      CPRLOG= PRLOG - (0.4 + 0.67*FCLOG)
      FLOG = FCLOG/(1.0 + (CPRLOG/(XN-0.14*CPRLOG))**2)
      FC = 10.0**FLOG
      PCOR = FC * PCOR
      RF(127) = RF(127) * PCOR
      RB(127) = RB(127) * PCOR
      RF(127) = RF(127)*C(2)
      RB(127) = RB(127)*C(17)
C     R128: H + C2H5 = H2 + C2H4
      RF(128) = RF(128)*C(2)
      RB(128) = RB(128)*C(6)*C(16)
C     R129: O + C2H5 = CH3 + CH2O
      RF(129) = RF(129)*C(3)
      RB(129) = RB(129)*C(10)*C(12)
C     R130: O2 + C2H5 = HO2 + C2H4
      RF(130) = RF(130)*C(9)
      RB(130) = RB(130)*C(5)*C(16)
C     R131: HO2 + C2H5 = O2 + C2H6
      RF(131) = RF(131)*C(5)
      RB(131) = RB(131)*C(9)*C(17)
C     R132: HO2 + C2H5 = H2O2 + C2H4
      RF(132) = RF(132)*C(5)
      RB(132) = RB(132)*C(8)*C(16)
C     R133: HO2 + C2H5 = OH + CH3 + CH2O
      RF(133) = RF(133)*C(5)
      RB(133) = RB(133)*C(4)*C(10)*C(12)
C     R134: H2O2 + C2H5 = HO2 + C2H6
      RF(134) = RF(134)*C(8)
      RB(134) = RB(134)*C(5)*C(17)
C     R135: C2H3 + C2H5 = C4H81
      CTB = CTOT+C(6)+5.D0*C(7)
     * +C(11)+5.D-1*C(13)
     * +C(14)+2.D0*C(17)
      PR = RKLOW(15) * CTB / RF(135)
      PCOR = PR / (1.0 + PR)
      PRLOG = LOG10(MAX(PR,SMALL))
      FCENT = 8.020000000000001D-1*EXP(-T/2.2779D3)
     *      + 1.98D-1*EXP(-T/6.D4)
     *     + EXP(-5.7232D3/T)
      FCLOG = LOG10(MAX(FCENT,SMALL))
      XN    = 0.75 - 1.27*FCLOG
      CPRLOG= PRLOG - (0.4 + 0.67*FCLOG)
      FLOG = FCLOG/(1.0 + (CPRLOG/(XN-0.14*CPRLOG))**2)
      FC = 10.0**FLOG
      PCOR = FC * PCOR
      RF(135) = RF(135) * PCOR
      RB(135) = RB(135) * PCOR
      RB(135) = RB(135)*C(23)
C     R136: C2H3 + C2H5 = CH3 + aC3H5
      RB(136) = RB(136)*C(10)*C(19)
C     R137: H + C2H6 = H2 + C2H5
      RF(137) = RF(137)*C(2)*C(17)
      RB(137) = RB(137)*C(6)
C     R138: O + C2H6 = OH + C2H5
      RF(138) = RF(138)*C(3)*C(17)
      RB(138) = RB(138)*C(4)
C     R139: OH + C2H6 = H2O + C2H5
      RF(139) = RF(139)*C(4)*C(17)
      RB(139) = RB(139)*C(7)
C     R140: CH2* + C2H6 = CH3 + C2H5
      RF(140) = RF(140)*C(17)
      RB(140) = RB(140)*C(10)
C     R141: CH3 + C2H6 = CH4 + C2H5
      RF(141) = RF(141)*C(10)*C(17)
      RB(141) = RB(141)*C(11)
C     R142: H + aC3H5 = C3H6
      CTB = CTOT+C(6)+5.D0*C(7)
     * +C(11)+5.D-1*C(13)
     * +C(14)+2.D0*C(17)
      PR = RKLOW(16) * CTB / RF(142)
      PCOR = PR / (1.0 + PR)
      PRLOG = LOG10(MAX(PR,SMALL))
      FCENT = 9.8D-1*EXP(-T/1.0966D3)
     *      + 2.D-2*EXP(-T/1.0966D3)
     *     + EXP(-6.8595D3/T)
      FCLOG = LOG10(MAX(FCENT,SMALL))
      XN    = 0.75 - 1.27*FCLOG
      CPRLOG= PRLOG - (0.4 + 0.67*FCLOG)
      FLOG = FCLOG/(1.0 + (CPRLOG/(XN-0.14*CPRLOG))**2)
      FC = 10.0**FLOG
      PCOR = FC * PCOR
      RF(142) = RF(142) * PCOR
      RB(142) = RB(142) * PCOR
      RF(142) = RF(142)*C(2)*C(19)
      RB(142) = RB(142)*C(20)
C     R143: O + aC3H5 = H + C2H3CHO
      RF(143) = RF(143)*C(3)*C(19)
      RB(143) = RB(143)*C(2)*C(21)
C     R144: OH + aC3H5 = 2H + C2H3CHO
      RF(144) = RF(144)*C(4)*C(19)
      RB(144) = RB(144)*C(2)*C(2)*C(21)
C     R145: O2 + aC3H5 = OH + C2H3CHO
      RF(145) = RF(145)*C(9)*C(19)
      RB(145) = RB(145)*C(4)*C(21)
C     R146: HO2 + aC3H5 = O2 + C3H6
      RF(146) = RF(146)*C(5)*C(19)
      RB(146) = RB(146)*C(9)*C(20)
C     R147: HO2 + aC3H5 = OH + CH2O + C2H3
      RF(147) = RF(147)*C(5)*C(19)
      RB(147) = RB(147)*C(4)*C(12)
C     R148: HCO + aC3H5 = CO + C3H6
      RF(148) = RF(148)*C(19)
      RB(148) = RB(148)*C(13)*C(20)
C     R149: CH3 + aC3H5 = C4H81
      CTB = CTOT+C(6)+5.D0*C(7)
     * +C(11)+5.D-1*C(13)
     * +C(14)+2.D0*C(17)
      PR = RKLOW(17) * CTB / RF(149)
      PCOR = PR / (1.0 + PR)
      PRLOG = LOG10(MAX(PR,SMALL))
      FCENT = 8.96D-1*EXP(-T/1.606D3)
     *      + 1.04D-1*EXP(-T/6.D4)
     *     + EXP(-6.1184D3/T)
      FCLOG = LOG10(MAX(FCENT,SMALL))
      XN    = 0.75 - 1.27*FCLOG
      CPRLOG= PRLOG - (0.4 + 0.67*FCLOG)
      FLOG = FCLOG/(1.0 + (CPRLOG/(XN-0.14*CPRLOG))**2)
      FC = 10.0**FLOG
      PCOR = FC * PCOR
      RF(149) = RF(149) * PCOR
      RB(149) = RB(149) * PCOR
      RF(149) = RF(149)*C(10)*C(19)
      RB(149) = RB(149)*C(23)
C     R150: H + C3H6 = nC3H7
      CTB = CTOT+C(6)+5.D0*C(7)
     * +C(11)+5.D-1*C(13)
     * +C(14)+2.D0*C(17)
      PR = RKLOW(18) * CTB / RF(150)
      PCOR = PR / (1.0 + PR)
      PRLOG = LOG10(MAX(PR,SMALL))
      FCENT = 0.D0*EXP(-T/1.D3)
     *      + 1.D0*EXP(-T/1.31D3)
     *     + EXP(-4.8097D4/T)
      FCLOG = LOG10(MAX(FCENT,SMALL))
      XN    = 0.75 - 1.27*FCLOG
      CPRLOG= PRLOG - (0.4 + 0.67*FCLOG)
      FLOG = FCLOG/(1.0 + (CPRLOG/(XN-0.14*CPRLOG))**2)
      FC = 10.0**FLOG
      PCOR = FC * PCOR
      RF(150) = RF(150) * PCOR
      RB(150) = RB(150) * PCOR
      RF(150) = RF(150)*C(2)*C(20)
C     R151: H + C3H6 = CH3 + C2H4
      RF(151) = RF(151)*C(2)*C(20)
      RB(151) = RB(151)*C(10)*C(16)
C     R152: H + C3H6 = H2 + aC3H5
      RF(152) = RF(152)*C(2)*C(20)
      RB(152) = RB(152)*C(6)*C(19)
C     R153: O + C3H6 = 2H + C2H3CHO
      RF(153) = RF(153)*C(3)*C(20)
      RB(153) = RB(153)*C(2)*C(2)*C(21)
C     R154: O + C3H6 = HCO + C2H5
      RF(154) = RF(154)*C(3)*C(20)
C     R155: O + C3H6 = OH + aC3H5
      RF(155) = RF(155)*C(3)*C(20)
      RB(155) = RB(155)*C(4)*C(19)
C     R156: OH + C3H6 = H2O + aC3H5
      RF(156) = RF(156)*C(4)*C(20)
      RB(156) = RB(156)*C(7)*C(19)
C     R157: HO2 + C3H6 = H2O2 + aC3H5
      RF(157) = RF(157)*C(5)*C(20)
      RB(157) = RB(157)*C(8)*C(19)
C     R158: CH3 + C3H6 = CH4 + aC3H5
      RF(158) = RF(158)*C(10)*C(20)
      RB(158) = RB(158)*C(11)*C(19)
C     R159: H + C2H3CHO = HCO + C2H4
      RF(159) = RF(159)*C(2)*C(21)
      RB(159) = RB(159)*C(16)
C     R160: O + C2H3CHO = OH + CO + C2H3
      RF(160) = RF(160)*C(3)*C(21)
      RB(160) = RB(160)*C(4)*C(13)
C     R161: OH + C2H3CHO = H2O + CO + C2H3
      RF(161) = RF(161)*C(4)*C(21)
      RB(161) = RB(161)*C(7)*C(13)
C     R162: H + nC3H7 = CH3 + C2H5
      RF(162) = RF(162)*C(2)
      RB(162) = RB(162)*C(10)
C     R163: H + nC3H7 = H2 + C3H6
      RF(163) = RF(163)*C(2)
      RB(163) = RB(163)*C(6)*C(20)
C     R164: O + nC3H7 = CH2O + C2H5
      RF(164) = RF(164)*C(3)
      RB(164) = RB(164)*C(12)
C     R165: OH + nC3H7 = H2O + C3H6
      RF(165) = RF(165)*C(4)
      RB(165) = RB(165)*C(7)*C(20)
C     R166: O2 + nC3H7 = HO2 + C3H6
      RF(166) = RF(166)*C(9)
      RB(166) = RB(166)*C(5)*C(20)
C     R167: HO2 + nC3H7 = OH + CH2O + C2H5
      RF(167) = RF(167)*C(5)
      RB(167) = RB(167)*C(4)*C(12)
C     R168: CH3 + nC3H7 = CH4 + C3H6
      RF(168) = RF(168)*C(10)
      RB(168) = RB(168)*C(11)*C(20)
C     R169: H + C4H7 = C4H81
      CTB = CTOT+C(6)+5.D0*C(7)
     * +C(11)+5.D-1*C(13)
     * +C(14)+2.D0*C(17)
      PR = RKLOW(19) * CTB / RF(169)
      PCOR = PR / (1.0 + PR)
      PRLOG = LOG10(MAX(PR,SMALL))
      FCENT = 5.02D-1*EXP(-T/1.314D3)
     *      + 4.98D-1*EXP(-T/1.314D3)
     *     + EXP(-5.D4/T)
      FCLOG = LOG10(MAX(FCENT,SMALL))
      XN    = 0.75 - 1.27*FCLOG
      CPRLOG= PRLOG - (0.4 + 0.67*FCLOG)
      FLOG = FCLOG/(1.0 + (CPRLOG/(XN-0.14*CPRLOG))**2)
      FC = 10.0**FLOG
      PCOR = FC * PCOR
      RF(169) = RF(169) * PCOR
      RB(169) = RB(169) * PCOR
      RF(169) = RF(169)*C(2)*C(22)
      RB(169) = RB(169)*C(23)
C     R170: H + C4H7 = CH3 + aC3H5
      RF(170) = RF(170)*C(2)*C(22)
      RB(170) = RB(170)*C(10)*C(19)
C     R171: HO2 + C4H7 = OH + CH2O + aC3H5
      RF(171) = RF(171)*C(5)*C(22)
      RB(171) = RB(171)*C(4)*C(12)*C(19)
C     R172: HCO + C4H7 = CO + C4H81
      RF(172) = RF(172)*C(22)
      RB(172) = RB(172)*C(13)*C(23)
C     R173: H + C4H81 = pC4H9
      CTB = CTOT+C(6)+5.D0*C(7)
     * +C(11)+5.D-1*C(13)
     * +C(14)+2.D0*C(17)
      PR = RKLOW(20) * CTB / RF(173)
      PCOR = PR / (1.0 + PR)
      PRLOG = LOG10(MAX(PR,SMALL))
      FCENT = 0.D0*EXP(-T/1.D3)
     *      + 1.D0*EXP(-T/1.31D3)
     *     + EXP(-4.8097D4/T)
      FCLOG = LOG10(MAX(FCENT,SMALL))
      XN    = 0.75 - 1.27*FCLOG
      CPRLOG= PRLOG - (0.4 + 0.67*FCLOG)
      FLOG = FCLOG/(1.0 + (CPRLOG/(XN-0.14*CPRLOG))**2)
      FC = 10.0**FLOG
      PCOR = FC * PCOR
      RF(173) = RF(173) * PCOR
      RB(173) = RB(173) * PCOR
      RF(173) = RF(173)*C(2)*C(23)
C     R174: H + C4H81 = C2H4 + C2H5
      RF(174) = RF(174)*C(2)*C(23)
      RB(174) = RB(174)*C(16)
C     R175: H + C4H81 = CH3 + C3H6
      RF(175) = RF(175)*C(2)*C(23)
      RB(175) = RB(175)*C(10)*C(20)
C     R176: H + C4H81 = H2 + C4H7
      RF(176) = RF(176)*C(2)*C(23)
      RB(176) = RB(176)*C(6)*C(22)
C     R177: O + C4H81 = HCO + nC3H7
      RF(177) = RF(177)*C(3)*C(23)
C     R178: O + C4H81 = OH + C4H7
      RF(178) = RF(178)*C(3)*C(23)
      RB(178) = RB(178)*C(4)*C(22)
C     R179: O + C4H81 = OH + C4H7
      RF(179) = RF(179)*C(3)*C(23)
      RB(179) = RB(179)*C(4)*C(22)
C     R180: OH + C4H81 = H2O + C4H7
      RF(180) = RF(180)*C(4)*C(23)
      RB(180) = RB(180)*C(7)*C(22)
C     R181: O2 + C4H81 = HO2 + C4H7
      RF(181) = RF(181)*C(9)*C(23)
      RB(181) = RB(181)*C(5)*C(22)
C     R182: HO2 + C4H81 = H2O2 + C4H7
      RF(182) = RF(182)*C(5)*C(23)
      RB(182) = RB(182)*C(8)*C(22)
C     R183: CH3 + C4H81 = CH4 + C4H7
      RF(183) = RF(183)*C(10)*C(23)
      RB(183) = RB(183)*C(11)*C(22)
C     R184: H + pC4H9 = 2C2H5
      RF(184) = RF(184)*C(2)
C     R185: H + pC4H9 = H2 + C4H81
      RF(185) = RF(185)*C(2)
      RB(185) = RB(185)*C(6)*C(23)
C     R186: O + pC4H9 = CH2O + nC3H7
      RF(186) = RF(186)*C(3)
      RB(186) = RB(186)*C(12)
C     R187: OH + pC4H9 = H2O + C4H81
      RF(187) = RF(187)*C(4)
      RB(187) = RB(187)*C(7)*C(23)
C     R188: O2 + pC4H9 = HO2 + C4H81
      RF(188) = RF(188)*C(9)
      RB(188) = RB(188)*C(5)*C(23)
C     R189: HO2 + pC4H9 = OH + CH2O + nC3H7
      RF(189) = RF(189)*C(5)
      RB(189) = RB(189)*C(4)*C(12)
C     R190: CH3 + pC4H9 = CH4 + C4H81
      RF(190) = RF(190)*C(10)
      RB(190) = RB(190)*C(11)*C(23)
C     R191: C5H9 = C2H4 + aC3H5
      RF(191) = RF(191)*C(24)
C     R192: C5H9 = C2H3 + C3H6
      RF(192) = RF(192)*C(24)
C     R193: H + C5H10 = PXC5H11
      CTB = CTOT+C(6)+5.D0*C(7)
     * +C(11)+5.D-1*C(13)
     * +C(14)+2.D0*C(17)
      PR = RKLOW(21) * CTB / RF(193)
      PCOR = PR / (1.0 + PR)
      PRLOG = LOG10(MAX(PR,SMALL))
      FCENT = 0.D0*EXP(-T/1.D3)
     *      + 1.D0*EXP(-T/1.31D3)
     *     + EXP(-4.8097D4/T)
      FCLOG = LOG10(MAX(FCENT,SMALL))
      XN    = 0.75 - 1.27*FCLOG
      CPRLOG= PRLOG - (0.4 + 0.67*FCLOG)
      FLOG = FCLOG/(1.0 + (CPRLOG/(XN-0.14*CPRLOG))**2)
      FC = 10.0**FLOG
      PCOR = FC * PCOR
      RF(193) = RF(193) * PCOR
      RB(193) = RB(193) * PCOR
      RF(193) = RF(193)*C(2)*C(25)
C     R194: H + C5H10 = C2H4 + nC3H7
      RF(194) = RF(194)*C(2)*C(25)
      RB(194) = RB(194)*C(16)
C     R195: H + C5H10 = C2H5 + C3H6
      RF(195) = RF(195)*C(2)*C(25)
      RB(195) = RB(195)*C(20)
C     R196: C2H4 + nC3H7 = PXC5H11
      RF(196) = RF(196)*C(16)
C     R197: H + C6H12 = PXC6H13
      CTB = CTOT+C(6)+5.D0*C(7)
     * +C(11)+5.D-1*C(13)
     * +C(14)+2.D0*C(17)
      PR = RKLOW(22) * CTB / RF(197)
      PCOR = PR / (1.0 + PR)
      PRLOG = LOG10(MAX(PR,SMALL))
      FCENT = 0.D0*EXP(-T/1.D3)
     *      + 1.D0*EXP(-T/1.31D3)
     *     + EXP(-4.8097D4/T)
      FCLOG = LOG10(MAX(FCENT,SMALL))
      XN    = 0.75 - 1.27*FCLOG
      CPRLOG= PRLOG - (0.4 + 0.67*FCLOG)
      FLOG = FCLOG/(1.0 + (CPRLOG/(XN-0.14*CPRLOG))**2)
      FC = 10.0**FLOG
      PCOR = FC * PCOR
      RF(197) = RF(197) * PCOR
      RB(197) = RB(197) * PCOR
      RF(197) = RF(197)*C(2)*C(26)
C     R198: H + C6H12 = C2H4 + pC4H9
      RF(198) = RF(198)*C(2)*C(26)
      RB(198) = RB(198)*C(16)
C     R199: H + C6H12 = C3H6 + nC3H7
      RF(199) = RF(199)*C(2)*C(26)
      RB(199) = RB(199)*C(20)
C     R200: C2H4 + pC4H9 = PXC6H13
      RF(200) = RF(200)*C(16)
C     R201: H + C7H14 = PXC7H15
      CTB = CTOT+C(6)+5.D0*C(7)
     * +C(11)+5.D-1*C(13)
     * +C(14)+2.D0*C(17)
      PR = RKLOW(23) * CTB / RF(201)
      PCOR = PR / (1.0 + PR)
      PRLOG = LOG10(MAX(PR,SMALL))
      FCENT = 0.D0*EXP(-T/1.D3)
     *      + 1.D0*EXP(-T/1.31D3)
     *     + EXP(-4.8097D4/T)
      FCLOG = LOG10(MAX(FCENT,SMALL))
      XN    = 0.75 - 1.27*FCLOG
      CPRLOG= PRLOG - (0.4 + 0.67*FCLOG)
      FLOG = FCLOG/(1.0 + (CPRLOG/(XN-0.14*CPRLOG))**2)
      FC = 10.0**FLOG
      PCOR = FC * PCOR
      RF(201) = RF(201) * PCOR
      RB(201) = RB(201) * PCOR
      RF(201) = RF(201)*C(2)*C(27)
C     R202: H + C7H14 = C2H4 + PXC5H11
      RF(202) = RF(202)*C(2)*C(27)
      RB(202) = RB(202)*C(16)
C     R203: H + C7H14 = C3H6 + pC4H9
      RF(203) = RF(203)*C(2)*C(27)
      RB(203) = RB(203)*C(20)
C     R204: C2H4 + PXC5H11 = PXC7H15
      RF(204) = RF(204)*C(16)
C     R205: H + C8H16 = PXC8H17
      CTB = CTOT+C(6)+5.D0*C(7)
     * +C(11)+5.D-1*C(13)
     * +C(14)+2.D0*C(17)
      PR = RKLOW(24) * CTB / RF(205)
      PCOR = PR / (1.0 + PR)
      PRLOG = LOG10(MAX(PR,SMALL))
      FCENT = 0.D0*EXP(-T/1.D3)
     *      + 1.D0*EXP(-T/1.31D3)
     *     + EXP(-4.8097D4/T)
      FCLOG = LOG10(MAX(FCENT,SMALL))
      XN    = 0.75 - 1.27*FCLOG
      CPRLOG= PRLOG - (0.4 + 0.67*FCLOG)
      FLOG = FCLOG/(1.0 + (CPRLOG/(XN-0.14*CPRLOG))**2)
      FC = 10.0**FLOG
      PCOR = FC * PCOR
      RF(205) = RF(205) * PCOR
      RB(205) = RB(205) * PCOR
      RF(205) = RF(205)*C(2)*C(28)
C     R206: H + C8H16 = C2H4 + PXC6H13
      RF(206) = RF(206)*C(2)*C(28)
      RB(206) = RB(206)*C(16)
C     R207: H + C8H16 = C3H6 + PXC5H11
      RF(207) = RF(207)*C(2)*C(28)
      RB(207) = RB(207)*C(20)
C     R208: C2H4 + PXC6H13 = PXC8H17
      RF(208) = RF(208)*C(16)
C     R209: H + C9H18 = PXC9H19
      CTB = CTOT+C(6)+5.D0*C(7)
     * +C(11)+5.D-1*C(13)
     * +C(14)+2.D0*C(17)
      PR = RKLOW(25) * CTB / RF(209)
      PCOR = PR / (1.0 + PR)
      PRLOG = LOG10(MAX(PR,SMALL))
      FCENT = 0.D0*EXP(-T/1.D3)
     *      + 1.D0*EXP(-T/1.31D3)
     *     + EXP(-4.8097D4/T)
      FCLOG = LOG10(MAX(FCENT,SMALL))
      XN    = 0.75 - 1.27*FCLOG
      CPRLOG= PRLOG - (0.4 + 0.67*FCLOG)
      FLOG = FCLOG/(1.0 + (CPRLOG/(XN-0.14*CPRLOG))**2)
      FC = 10.0**FLOG
      PCOR = FC * PCOR
      RF(209) = RF(209) * PCOR
      RB(209) = RB(209) * PCOR
      RF(209) = RF(209)*C(2)*C(29)
      RB(209) = RB(209)*C(30)
C     R210: H + C9H18 = C2H4 + PXC7H15
      RF(210) = RF(210)*C(2)*C(29)
      RB(210) = RB(210)*C(16)
C     R211: H + C9H18 = C3H6 + PXC6H13
      RF(211) = RF(211)*C(2)*C(29)
      RB(211) = RB(211)*C(20)
C     R212: C2H4 + PXC7H15 = PXC9H19
      RF(212) = RF(212)*C(16)
      RB(212) = RB(212)*C(30)
C     R213: H + C10H20 = PXC10H21
      CTB = CTOT+C(6)+5.D0*C(7)
     * +C(11)+5.D-1*C(13)
     * +C(14)+2.D0*C(17)
      PR = RKLOW(26) * CTB / RF(213)
      PCOR = PR / (1.0 + PR)
      PRLOG = LOG10(MAX(PR,SMALL))
      FCENT = 0.D0*EXP(-T/1.D3)
     *      + 1.D0*EXP(-T/1.31D3)
     *     + EXP(-4.8097D4/T)
      FCLOG = LOG10(MAX(FCENT,SMALL))
      XN    = 0.75 - 1.27*FCLOG
      CPRLOG= PRLOG - (0.4 + 0.67*FCLOG)
      FLOG = FCLOG/(1.0 + (CPRLOG/(XN-0.14*CPRLOG))**2)
      FC = 10.0**FLOG
      PCOR = FC * PCOR
      RF(213) = RF(213) * PCOR
      RB(213) = RB(213) * PCOR
      RF(213) = RF(213)*C(2)*C(31)
C     R214: H + C10H20 = C2H4 + PXC8H17
      RF(214) = RF(214)*C(2)*C(31)
      RB(214) = RB(214)*C(16)
C     R215: H + C10H20 = C3H6 + PXC7H15
      RF(215) = RF(215)*C(2)*C(31)
      RB(215) = RB(215)*C(20)
C     R216: C2H4 + PXC8H17 = PXC10H21
      RF(216) = RF(216)*C(16)
C     R217: C12H24 = C5H9 + PXC7H15
      RF(217) = RF(217)*C(32)
      RB(217) = RB(217)*C(24)
C     R218: C2H4 + PXC10H21 = PXC12H25
      RF(218) = RF(218)*C(16)
C     R219: PXC12H25 = S3XC12H25
C     R220: C3H6 + PXC9H19 = SXC12H25
      RF(220) = RF(220)*C(20)*C(30)
C     R221: C4H81 + PXC8H17 = SXC12H25
      RF(221) = RF(221)*C(23)
C     R222: C5H10 + PXC7H15 = S3XC12H25
      RF(222) = RF(222)*C(25)
C     R223: C2H5 + C10H20 = S3XC12H25
      RF(223) = RF(223)*C(31)
C     R224: C6H12 + PXC6H13 = S3XC12H25
      RF(224) = RF(224)*C(26)
C     R225: nC3H7 + C9H18 = S3XC12H25
      RF(225) = RF(225)*C(29)
C     R226: PXC5H11 + C7H14 = S3XC12H25
      RF(226) = RF(226)*C(27)
C     R227: pC4H9 + C8H16 = S3XC12H25
      RF(227) = RF(227)*C(28)
C     R228: C2H5 + PXC10H21 = NC12H26
      RB(228) = RB(228)*C(1)
C     R229: nC3H7 + PXC9H19 = NC12H26
      RF(229) = RF(229)*C(30)
      RB(229) = RB(229)*C(1)
C     R230: pC4H9 + PXC8H17 = NC12H26
      RB(230) = RB(230)*C(1)
C     R231: PXC5H11 + PXC7H15 = NC12H26
      RB(231) = RB(231)*C(1)
C     R232: 2PXC6H13 = NC12H26
      RB(232) = RB(232)*C(1)
C     R233: NC12H26 + H = H2 + PXC12H25
      RF(233) = RF(233)*C(1)*C(2)
      RB(233) = RB(233)*C(6)
C     R234: NC12H26 + H = H2 + SXC12H25
      RF(234) = RF(234)*C(1)*C(2)
      RB(234) = RB(234)*C(6)
C     R235: NC12H26 + H = H2 + S3XC12H25
      RF(235) = RF(235)*C(1)*C(2)
      RB(235) = RB(235)*C(6)
C     R236: NC12H26 + O = OH + PXC12H25
      RF(236) = RF(236)*C(1)*C(3)
      RB(236) = RB(236)*C(4)
C     R237: NC12H26 + O = OH + SXC12H25
      RF(237) = RF(237)*C(1)*C(3)
      RB(237) = RB(237)*C(4)
C     R238: NC12H26 + O = OH + S3XC12H25
      RF(238) = RF(238)*C(1)*C(3)
      RB(238) = RB(238)*C(4)
C     R239: NC12H26 + OH = H2O + PXC12H25
      RF(239) = RF(239)*C(1)*C(4)
      RB(239) = RB(239)*C(7)
C     R240: NC12H26 + OH = H2O + SXC12H25
      RF(240) = RF(240)*C(1)*C(4)
      RB(240) = RB(240)*C(7)
C     R241: NC12H26 + OH = H2O + S3XC12H25
      RF(241) = RF(241)*C(1)*C(4)
      RB(241) = RB(241)*C(7)
C     R242: NC12H26 + O2 = HO2 + PXC12H25
      RF(242) = RF(242)*C(1)*C(9)
      RB(242) = RB(242)*C(5)
C     R243: NC12H26 + O2 = HO2 + SXC12H25
      RF(243) = RF(243)*C(1)*C(9)
      RB(243) = RB(243)*C(5)
C     R244: NC12H26 + O2 = HO2 + S3XC12H25
      RF(244) = RF(244)*C(1)*C(9)
      RB(244) = RB(244)*C(5)
C     R245: NC12H26 + HO2 = H2O2 + PXC12H25
      RF(245) = RF(245)*C(1)*C(5)
      RB(245) = RB(245)*C(8)
C     R246: NC12H26 + HO2 = H2O2 + SXC12H25
      RF(246) = RF(246)*C(1)*C(5)
      RB(246) = RB(246)*C(8)
C     R247: NC12H26 + HO2 = H2O2 + S3XC12H25
      RF(247) = RF(247)*C(1)*C(5)
      RB(247) = RB(247)*C(8)
C     R248: NC12H26 + CH3 = CH4 + PXC12H25
      RF(248) = RF(248)*C(1)*C(10)
      RB(248) = RB(248)*C(11)
C     R249: NC12H26 + CH3 = CH4 + SXC12H25
      RF(249) = RF(249)*C(1)*C(10)
      RB(249) = RB(249)*C(11)
C     R250: NC12H26 + CH3 = CH4 + S3XC12H25
      RF(250) = RF(250)*C(1)*C(10)
      RB(250) = RB(250)*C(11)
C     R251: O2 + PXC12H25 = C12H25O2
      RF(251) = RF(251)*C(9)
C     R252: C12H25O2 = O2 + PXC12H25
      RF(252) = RF(252)*C(33)
C     R253: O2 + SXC12H25 = C12H25O2
      RF(253) = RF(253)*C(9)
C     R254: C12H25O2 = O2 + SXC12H25
      RF(254) = RF(254)*C(33)
C     R255: O2 + S3XC12H25 = C12H25O2
      RF(255) = RF(255)*C(9)
C     R256: C12H25O2 = O2 + S3XC12H25
      RF(256) = RF(256)*C(33)
C     R257: C12H25O2 = C12OOH
      RF(257) = RF(257)*C(33)
C     R258: C12OOH = C12H25O2
C     R259: O2 + PXC12H25 = HO2 + C12H24
      RF(259) = RF(259)*C(9)
C     R260: HO2 + C12H24 = O2 + PXC12H25
      RF(260) = RF(260)*C(5)*C(32)
C     R261: O2 + SXC12H25 = HO2 + C12H24
      RF(261) = RF(261)*C(9)
C     R262: HO2 + C12H24 = O2 + SXC12H25
      RF(262) = RF(262)*C(5)*C(32)
C     R263: O2 + S3XC12H25 = HO2 + C12H24
      RF(263) = RF(263)*C(9)
C     R264: HO2 + C12H24 = O2 + S3XC12H25
      RF(264) = RF(264)*C(5)*C(32)
C     R265: O2 + C12OOH = O2C12H24OOH
      RF(265) = RF(265)*C(9)
C     R266: O2C12H24OOH = O2 + C12OOH
C     R267: O2C12H24OOH = OH + OC12H23OOH
      RB(267) = RB(267)*C(4)*C(34)
C     R268: OC12H23OOH = OH + 3C2H4 + C2H5 + 2CH2CHO
      RF(268) = RF(268)*C(34)
C
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE RDOT(RF, RB, WDOT)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
      DIMENSION RF(*), RB(*), WDOT(*)
C
      DO K = 1, 35
         WDOT(K) = 0D0
      ENDDO
C
      ROP = RF(1)-RB(1)
      WDOT(2) = WDOT(2) -ROP
      WDOT(3) = WDOT(3) +ROP
      WDOT(4) = WDOT(4) +ROP
      WDOT(9) = WDOT(9) -ROP
      ROP = RF(2)-RB(2)
      WDOT(2) = WDOT(2) +ROP
      WDOT(3) = WDOT(3) -ROP
      WDOT(4) = WDOT(4) +ROP
      WDOT(6) = WDOT(6) -ROP
      ROP = RF(3)-RB(3)
      WDOT(2) = WDOT(2) +ROP
      WDOT(4) = WDOT(4) -ROP
      WDOT(6) = WDOT(6) -ROP
      WDOT(7) = WDOT(7) +ROP
      ROP = RF(4)-RB(4)
      WDOT(3) = WDOT(3) +ROP
      WDOT(4) = WDOT(4) -2*ROP
      WDOT(7) = WDOT(7) +ROP
      ROP = RF(5)-RB(5)
      WDOT(2) = WDOT(2) -ROP
      WDOT(5) = WDOT(5) +ROP
      WDOT(9) = WDOT(9) -ROP
      ROP = RF(6)-RB(6)
      WDOT(2) = WDOT(2) -ROP
      WDOT(4) = WDOT(4) +2*ROP
      WDOT(5) = WDOT(5) -ROP
      ROP = RF(7)-RB(7)
      WDOT(2) = WDOT(2) +ROP
      WDOT(5) = WDOT(5) +ROP
      WDOT(6) = WDOT(6) -ROP
      WDOT(9) = WDOT(9) -ROP
      ROP = RF(8)-RB(8)
      WDOT(4) = WDOT(4) -ROP
      WDOT(5) = WDOT(5) -ROP
      WDOT(7) = WDOT(7) +ROP
      WDOT(9) = WDOT(9) +ROP
      ROP = RF(9)-RB(9)
      WDOT(2) = WDOT(2) -ROP
      WDOT(3) = WDOT(3) +ROP
      WDOT(5) = WDOT(5) -ROP
      WDOT(7) = WDOT(7) +ROP
      ROP = RF(10)-RB(10)
      WDOT(3) = WDOT(3) -ROP
      WDOT(4) = WDOT(4) +ROP
      WDOT(5) = WDOT(5) -ROP
      WDOT(9) = WDOT(9) +ROP
      ROP = RF(11)-RB(11)
      WDOT(5) = WDOT(5) -2*ROP
      WDOT(8) = WDOT(8) +ROP
      WDOT(9) = WDOT(9) +ROP
      ROP = RF(12)-RB(12)
      WDOT(5) = WDOT(5) -2*ROP
      WDOT(8) = WDOT(8) +ROP
      WDOT(9) = WDOT(9) +ROP
      ROP = RF(13)-RB(13)
      WDOT(2) = WDOT(2) -ROP
      WDOT(4) = WDOT(4) +ROP
      WDOT(7) = WDOT(7) +ROP
      WDOT(8) = WDOT(8) -ROP
      ROP = RF(14)-RB(14)
      WDOT(2) = WDOT(2) -ROP
      WDOT(5) = WDOT(5) +ROP
      WDOT(6) = WDOT(6) +ROP
      WDOT(8) = WDOT(8) -ROP
      ROP = RF(15)-RB(15)
      WDOT(3) = WDOT(3) -ROP
      WDOT(4) = WDOT(4) +ROP
      WDOT(5) = WDOT(5) +ROP
      WDOT(8) = WDOT(8) -ROP
      ROP = RF(16)-RB(16)
      WDOT(4) = WDOT(4) -ROP
      WDOT(5) = WDOT(5) +ROP
      WDOT(7) = WDOT(7) +ROP
      WDOT(8) = WDOT(8) -ROP
      ROP = RF(17)-RB(17)
      WDOT(4) = WDOT(4) -ROP
      WDOT(5) = WDOT(5) +ROP
      WDOT(7) = WDOT(7) +ROP
      WDOT(8) = WDOT(8) -ROP
      ROP = RF(18)-RB(18)
      WDOT(4) = WDOT(4) -2*ROP
      WDOT(8) = WDOT(8) +ROP
      ROP = RF(19)-RB(19)
      WDOT(2) = WDOT(2) -2*ROP
      WDOT(6) = WDOT(6) +ROP
      ROP = RF(20)-RB(20)
      WDOT(2) = WDOT(2) -ROP
      WDOT(4) = WDOT(4) -ROP
      WDOT(7) = WDOT(7) +ROP
      ROP = RF(21)-RB(21)
      WDOT(3) = WDOT(3) -2*ROP
      WDOT(9) = WDOT(9) +ROP
      ROP = RF(22)-RB(22)
      WDOT(2) = WDOT(2) -2*ROP
      WDOT(6) = WDOT(6) +ROP
      ROP = RF(23)-RB(23)
      WDOT(2) = WDOT(2) -2*ROP
      WDOT(6) = WDOT(6) +ROP
      ROP = RF(24)-RB(24)
      WDOT(2) = WDOT(2) -2*ROP
      WDOT(6) = WDOT(6) +ROP
      ROP = RF(25)-RB(25)
      WDOT(2) = WDOT(2) -ROP
      WDOT(3) = WDOT(3) -ROP
      WDOT(4) = WDOT(4) +ROP
      ROP = RF(26)-RB(26)
      WDOT(2) = WDOT(2) +ROP
      WDOT(4) = WDOT(4) -ROP
      WDOT(13) = WDOT(13) -ROP
      WDOT(14) = WDOT(14) +ROP
      ROP = RF(27)-RB(27)
      WDOT(2) = WDOT(2) +ROP
      WDOT(4) = WDOT(4) -ROP
      WDOT(13) = WDOT(13) -ROP
      WDOT(14) = WDOT(14) +ROP
      ROP = RF(28)-RB(28)
      WDOT(4) = WDOT(4) +ROP
      WDOT(5) = WDOT(5) -ROP
      WDOT(13) = WDOT(13) -ROP
      WDOT(14) = WDOT(14) +ROP
      ROP = RF(29)-RB(29)
      WDOT(3) = WDOT(3) -ROP
      WDOT(13) = WDOT(13) -ROP
      WDOT(14) = WDOT(14) +ROP
      ROP = RF(30)-RB(30)
      WDOT(3) = WDOT(3) +ROP
      WDOT(9) = WDOT(9) -ROP
      WDOT(13) = WDOT(13) -ROP
      WDOT(14) = WDOT(14) +ROP
      ROP = RF(31)-RB(31)
      WDOT(2) = WDOT(2) +ROP
      WDOT(13) = WDOT(13) +ROP
      ROP = RF(32)-RB(32)
      WDOT(2) = WDOT(2) -ROP
      WDOT(6) = WDOT(6) +ROP
      WDOT(13) = WDOT(13) +ROP
      ROP = RF(33)-RB(33)
      WDOT(3) = WDOT(3) -ROP
      WDOT(4) = WDOT(4) +ROP
      WDOT(13) = WDOT(13) +ROP
      ROP = RF(34)-RB(34)
      WDOT(2) = WDOT(2) +ROP
      WDOT(3) = WDOT(3) -ROP
      WDOT(14) = WDOT(14) +ROP
      ROP = RF(35)-RB(35)
      WDOT(4) = WDOT(4) -ROP
      WDOT(7) = WDOT(7) +ROP
      WDOT(13) = WDOT(13) +ROP
      ROP = RF(36)-RB(36)
      WDOT(5) = WDOT(5) +ROP
      WDOT(9) = WDOT(9) -ROP
      WDOT(13) = WDOT(13) +ROP
      ROP = RF(37)-RB(37)
      WDOT(2) = WDOT(2) +ROP
      WDOT(13) = WDOT(13) +ROP
      ROP = RF(38)-RB(38)
      WDOT(6) = WDOT(6) -ROP
      WDOT(12) = WDOT(12) +ROP
      WDOT(13) = WDOT(13) -ROP
      ROP = RF(39)-RB(39)
      WDOT(2) = WDOT(2) -ROP
      WDOT(12) = WDOT(12) +ROP
      ROP = RF(40)-RB(40)
      WDOT(2) = WDOT(2) -ROP
      WDOT(10) = WDOT(10) +ROP
      ROP = RF(41)-RB(41)
      WDOT(2) = WDOT(2) +ROP
      WDOT(3) = WDOT(3) -ROP
      ROP = RF(42)-RB(42)
      WDOT(2) = WDOT(2) +ROP
      WDOT(4) = WDOT(4) -ROP
      WDOT(12) = WDOT(12) +ROP
      ROP = RF(43)-RB(43)
      WDOT(2) = WDOT(2) +ROP
      WDOT(6) = WDOT(6) -ROP
      WDOT(10) = WDOT(10) +ROP
      ROP = RF(44)-RB(44)
      WDOT(4) = WDOT(4) +ROP
      WDOT(9) = WDOT(9) -ROP
      ROP = RF(45)-RB(45)
      WDOT(2) = WDOT(2) +2*ROP
      WDOT(9) = WDOT(9) -ROP
      WDOT(14) = WDOT(14) +ROP
      ROP = RF(46)-RB(46)
      WDOT(4) = WDOT(4) +ROP
      WDOT(5) = WDOT(5) -ROP
      WDOT(12) = WDOT(12) +ROP
      ROP = RF(47)-RB(47)
      WDOT(6) = WDOT(6) +ROP
      WDOT(15) = WDOT(15) +ROP
      ROP = RF(48)-RB(48)
      ROP = RF(49)-RB(49)
      WDOT(3) = WDOT(3) -ROP
      WDOT(6) = WDOT(6) +ROP
      WDOT(13) = WDOT(13) +ROP
      ROP = RF(50)-RB(50)
      WDOT(2) = WDOT(2) +ROP
      WDOT(3) = WDOT(3) -ROP
      ROP = RF(51)-RB(51)
      WDOT(2) = WDOT(2) +ROP
      WDOT(4) = WDOT(4) -ROP
      WDOT(12) = WDOT(12) +ROP
      ROP = RF(52)-RB(52)
      WDOT(2) = WDOT(2) +ROP
      WDOT(6) = WDOT(6) -ROP
      WDOT(10) = WDOT(10) +ROP
      ROP = RF(53)-RB(53)
      WDOT(2) = WDOT(2) +ROP
      WDOT(4) = WDOT(4) +ROP
      WDOT(9) = WDOT(9) -ROP
      WDOT(13) = WDOT(13) +ROP
      ROP = RF(54)-RB(54)
      WDOT(7) = WDOT(7) +ROP
      WDOT(9) = WDOT(9) -ROP
      WDOT(13) = WDOT(13) +ROP
      ROP = RF(55)-RB(55)
      ROP = RF(56)-RB(56)
      ROP = RF(57)-RB(57)
      ROP = RF(58)-RB(58)
      WDOT(12) = WDOT(12) +ROP
      WDOT(13) = WDOT(13) +ROP
      WDOT(14) = WDOT(14) -ROP
      ROP = RF(59)-RB(59)
      WDOT(2) = WDOT(2) -ROP
      WDOT(12) = WDOT(12) -ROP
      ROP = RF(60)-RB(60)
      WDOT(2) = WDOT(2) -ROP
      WDOT(6) = WDOT(6) +ROP
      WDOT(12) = WDOT(12) -ROP
      ROP = RF(61)-RB(61)
      WDOT(3) = WDOT(3) -ROP
      WDOT(4) = WDOT(4) +ROP
      WDOT(12) = WDOT(12) -ROP
      ROP = RF(62)-RB(62)
      WDOT(4) = WDOT(4) -ROP
      WDOT(7) = WDOT(7) +ROP
      WDOT(12) = WDOT(12) -ROP
      ROP = RF(63)-RB(63)
      WDOT(5) = WDOT(5) +ROP
      WDOT(9) = WDOT(9) -ROP
      WDOT(12) = WDOT(12) -ROP
      ROP = RF(64)-RB(64)
      WDOT(5) = WDOT(5) -ROP
      WDOT(8) = WDOT(8) +ROP
      WDOT(12) = WDOT(12) -ROP
      ROP = RF(65)-RB(65)
      WDOT(2) = WDOT(2) -ROP
      WDOT(10) = WDOT(10) -ROP
      WDOT(11) = WDOT(11) +ROP
      ROP = RF(66)-RB(66)
      WDOT(2) = WDOT(2) +ROP
      WDOT(3) = WDOT(3) -ROP
      WDOT(10) = WDOT(10) -ROP
      WDOT(12) = WDOT(12) +ROP
      ROP = RF(67)-RB(67)
      WDOT(4) = WDOT(4) -ROP
      WDOT(7) = WDOT(7) +ROP
      WDOT(10) = WDOT(10) -ROP
      ROP = RF(68)-RB(68)
      WDOT(4) = WDOT(4) -ROP
      WDOT(7) = WDOT(7) +ROP
      WDOT(10) = WDOT(10) -ROP
      ROP = RF(69)-RB(69)
      WDOT(3) = WDOT(3) +ROP
      WDOT(9) = WDOT(9) -ROP
      WDOT(10) = WDOT(10) -ROP
      ROP = RF(70)-RB(70)
      WDOT(4) = WDOT(4) +ROP
      WDOT(9) = WDOT(9) -ROP
      WDOT(10) = WDOT(10) -ROP
      WDOT(12) = WDOT(12) +ROP
      ROP = RF(71)-RB(71)
      WDOT(5) = WDOT(5) -ROP
      WDOT(9) = WDOT(9) +ROP
      WDOT(10) = WDOT(10) -ROP
      WDOT(11) = WDOT(11) +ROP
      ROP = RF(72)-RB(72)
      WDOT(4) = WDOT(4) +ROP
      WDOT(5) = WDOT(5) -ROP
      WDOT(10) = WDOT(10) -ROP
      ROP = RF(73)-RB(73)
      WDOT(5) = WDOT(5) +ROP
      WDOT(8) = WDOT(8) -ROP
      WDOT(10) = WDOT(10) -ROP
      WDOT(11) = WDOT(11) +ROP
      ROP = RF(74)-RB(74)
      WDOT(10) = WDOT(10) -ROP
      WDOT(11) = WDOT(11) +ROP
      WDOT(13) = WDOT(13) +ROP
      ROP = RF(75)-RB(75)
      WDOT(10) = WDOT(10) -ROP
      WDOT(11) = WDOT(11) +ROP
      WDOT(12) = WDOT(12) -ROP
      ROP = RF(76)-RB(76)
      WDOT(2) = WDOT(2) +ROP
      WDOT(10) = WDOT(10) -ROP
      WDOT(16) = WDOT(16) +ROP
      ROP = RF(77)-RB(77)
      WDOT(2) = WDOT(2) +ROP
      WDOT(10) = WDOT(10) -ROP
      WDOT(16) = WDOT(16) +ROP
      ROP = RF(78)-RB(78)
      WDOT(10) = WDOT(10) -2*ROP
      WDOT(17) = WDOT(17) +ROP
      ROP = RF(79)-RB(79)
      WDOT(2) = WDOT(2) +ROP
      WDOT(10) = WDOT(10) -2*ROP
      ROP = RF(80)-RB(80)
      WDOT(2) = WDOT(2) -ROP
      WDOT(6) = WDOT(6) +ROP
      WDOT(12) = WDOT(12) +ROP
      ROP = RF(81)-RB(81)
      WDOT(2) = WDOT(2) -ROP
      WDOT(4) = WDOT(4) +ROP
      WDOT(10) = WDOT(10) +ROP
      ROP = RF(82)-RB(82)
      WDOT(2) = WDOT(2) -ROP
      WDOT(7) = WDOT(7) +ROP
      ROP = RF(83)-RB(83)
      WDOT(3) = WDOT(3) -ROP
      WDOT(4) = WDOT(4) +ROP
      WDOT(12) = WDOT(12) +ROP
      ROP = RF(84)-RB(84)
      WDOT(4) = WDOT(4) -ROP
      WDOT(7) = WDOT(7) +ROP
      WDOT(12) = WDOT(12) +ROP
      ROP = RF(85)-RB(85)
      WDOT(5) = WDOT(5) +ROP
      WDOT(9) = WDOT(9) -ROP
      WDOT(12) = WDOT(12) +ROP
      ROP = RF(86)-RB(86)
      WDOT(2) = WDOT(2) -ROP
      WDOT(6) = WDOT(6) +ROP
      WDOT(10) = WDOT(10) +ROP
      WDOT(11) = WDOT(11) -ROP
      ROP = RF(87)-RB(87)
      WDOT(3) = WDOT(3) -ROP
      WDOT(4) = WDOT(4) +ROP
      WDOT(10) = WDOT(10) +ROP
      WDOT(11) = WDOT(11) -ROP
      ROP = RF(88)-RB(88)
      WDOT(4) = WDOT(4) -ROP
      WDOT(7) = WDOT(7) +ROP
      WDOT(10) = WDOT(10) +ROP
      WDOT(11) = WDOT(11) -ROP
      ROP = RF(89)-RB(89)
      WDOT(10) = WDOT(10) +2*ROP
      WDOT(11) = WDOT(11) -ROP
      ROP = RF(90)-RB(90)
      WDOT(10) = WDOT(10) +2*ROP
      WDOT(11) = WDOT(11) -ROP
      ROP = RF(91)-RB(91)
      WDOT(2) = WDOT(2) +ROP
      WDOT(15) = WDOT(15) +ROP
      ROP = RF(92)-RB(92)
      WDOT(3) = WDOT(3) -ROP
      WDOT(13) = WDOT(13) +ROP
      WDOT(15) = WDOT(15) -ROP
      ROP = RF(93)-RB(93)
      WDOT(4) = WDOT(4) -ROP
      WDOT(10) = WDOT(10) +ROP
      WDOT(13) = WDOT(13) +ROP
      WDOT(15) = WDOT(15) -ROP
      ROP = RF(94)-RB(94)
      WDOT(13) = WDOT(13) +ROP
      WDOT(15) = WDOT(15) -ROP
      ROP = RF(95)-RB(95)
      WDOT(10) = WDOT(10) -ROP
      WDOT(15) = WDOT(15) -ROP
      WDOT(19) = WDOT(19) +ROP
      ROP = RF(96)-RB(96)
      WDOT(2) = WDOT(2) -ROP
      WDOT(16) = WDOT(16) +ROP
      ROP = RF(97)-RB(97)
      WDOT(2) = WDOT(2) -ROP
      WDOT(6) = WDOT(6) +ROP
      WDOT(15) = WDOT(15) +ROP
      ROP = RF(98)-RB(98)
      WDOT(3) = WDOT(3) -ROP
      WDOT(10) = WDOT(10) +ROP
      WDOT(13) = WDOT(13) +ROP
      ROP = RF(99)-RB(99)
      WDOT(4) = WDOT(4) -ROP
      WDOT(7) = WDOT(7) +ROP
      WDOT(15) = WDOT(15) +ROP
      ROP = RF(100)-RB(100)
      WDOT(5) = WDOT(5) +ROP
      WDOT(9) = WDOT(9) -ROP
      WDOT(15) = WDOT(15) +ROP
      ROP = RF(101)-RB(101)
      WDOT(3) = WDOT(3) +ROP
      WDOT(9) = WDOT(9) -ROP
      WDOT(18) = WDOT(18) +ROP
      ROP = RF(102)-RB(102)
      WDOT(9) = WDOT(9) -ROP
      WDOT(12) = WDOT(12) +ROP
      ROP = RF(103)-RB(103)
      WDOT(4) = WDOT(4) +ROP
      WDOT(5) = WDOT(5) -ROP
      WDOT(18) = WDOT(18) +ROP
      ROP = RF(104)-RB(104)
      WDOT(5) = WDOT(5) +ROP
      WDOT(8) = WDOT(8) -ROP
      WDOT(16) = WDOT(16) +ROP
      ROP = RF(105)-RB(105)
      WDOT(13) = WDOT(13) +ROP
      WDOT(16) = WDOT(16) +ROP
      ROP = RF(106)-RB(106)
      WDOT(21) = WDOT(21) +ROP
      ROP = RF(107)-RB(107)
      WDOT(10) = WDOT(10) -ROP
      WDOT(11) = WDOT(11) +ROP
      WDOT(15) = WDOT(15) +ROP
      ROP = RF(108)-RB(108)
      WDOT(10) = WDOT(10) -ROP
      WDOT(20) = WDOT(20) +ROP
      ROP = RF(109)-RB(109)
      WDOT(2) = WDOT(2) +ROP
      WDOT(10) = WDOT(10) -ROP
      WDOT(19) = WDOT(19) +ROP
      ROP = RF(110)-RB(110)
      WDOT(15) = WDOT(15) +ROP
      WDOT(16) = WDOT(16) +ROP
      ROP = RF(111)-RB(111)
      WDOT(10) = WDOT(10) +ROP
      WDOT(13) = WDOT(13) +ROP
      WDOT(18) = WDOT(18) -ROP
      ROP = RF(112)-RB(112)
      WDOT(2) = WDOT(2) -ROP
      WDOT(10) = WDOT(10) +ROP
      WDOT(18) = WDOT(18) -ROP
      ROP = RF(113)-RB(113)
      WDOT(4) = WDOT(4) +ROP
      WDOT(9) = WDOT(9) -ROP
      WDOT(12) = WDOT(12) +ROP
      WDOT(13) = WDOT(13) +ROP
      WDOT(18) = WDOT(18) -ROP
      ROP = RF(114)-RB(114)
      WDOT(2) = WDOT(2) -ROP
      WDOT(16) = WDOT(16) -ROP
      ROP = RF(115)-RB(115)
      WDOT(2) = WDOT(2) -ROP
      WDOT(6) = WDOT(6) +ROP
      WDOT(16) = WDOT(16) -ROP
      ROP = RF(116)-RB(116)
      WDOT(3) = WDOT(3) -ROP
      WDOT(4) = WDOT(4) +ROP
      WDOT(16) = WDOT(16) -ROP
      ROP = RF(117)-RB(117)
      WDOT(3) = WDOT(3) -ROP
      WDOT(10) = WDOT(10) +ROP
      WDOT(16) = WDOT(16) -ROP
      ROP = RF(118)-RB(118)
      WDOT(3) = WDOT(3) -ROP
      WDOT(12) = WDOT(12) +ROP
      WDOT(16) = WDOT(16) -ROP
      ROP = RF(119)-RB(119)
      WDOT(4) = WDOT(4) -ROP
      WDOT(7) = WDOT(7) +ROP
      WDOT(16) = WDOT(16) -ROP
      ROP = RF(120)-RB(120)
      WDOT(13) = WDOT(13) +ROP
      WDOT(16) = WDOT(16) -ROP
      ROP = RF(121)-RB(121)
      WDOT(2) = WDOT(2) +ROP
      WDOT(16) = WDOT(16) -ROP
      WDOT(19) = WDOT(19) +ROP
      ROP = RF(122)-RB(122)
      WDOT(2) = WDOT(2) +ROP
      WDOT(16) = WDOT(16) -ROP
      WDOT(19) = WDOT(19) +ROP
      ROP = RF(123)-RB(123)
      WDOT(10) = WDOT(10) -ROP
      WDOT(11) = WDOT(11) +ROP
      WDOT(16) = WDOT(16) -ROP
      ROP = RF(124)-RB(124)
      WDOT(10) = WDOT(10) +ROP
      WDOT(16) = WDOT(16) +ROP
      ROP = RF(125)-RB(125)
      WDOT(5) = WDOT(5) +ROP
      WDOT(9) = WDOT(9) -ROP
      WDOT(16) = WDOT(16) -ROP
      ROP = RF(126)-RB(126)
      WDOT(16) = WDOT(16) -ROP
      WDOT(22) = WDOT(22) +ROP
      ROP = RF(127)-RB(127)
      WDOT(2) = WDOT(2) -ROP
      WDOT(17) = WDOT(17) +ROP
      ROP = RF(128)-RB(128)
      WDOT(2) = WDOT(2) -ROP
      WDOT(6) = WDOT(6) +ROP
      WDOT(16) = WDOT(16) +ROP
      ROP = RF(129)-RB(129)
      WDOT(3) = WDOT(3) -ROP
      WDOT(10) = WDOT(10) +ROP
      WDOT(12) = WDOT(12) +ROP
      ROP = RF(130)-RB(130)
      WDOT(5) = WDOT(5) +ROP
      WDOT(9) = WDOT(9) -ROP
      WDOT(16) = WDOT(16) +ROP
      ROP = RF(131)-RB(131)
      WDOT(5) = WDOT(5) -ROP
      WDOT(9) = WDOT(9) +ROP
      WDOT(17) = WDOT(17) +ROP
      ROP = RF(132)-RB(132)
      WDOT(5) = WDOT(5) -ROP
      WDOT(8) = WDOT(8) +ROP
      WDOT(16) = WDOT(16) +ROP
      ROP = RF(133)-RB(133)
      WDOT(4) = WDOT(4) +ROP
      WDOT(5) = WDOT(5) -ROP
      WDOT(10) = WDOT(10) +ROP
      WDOT(12) = WDOT(12) +ROP
      ROP = RF(134)-RB(134)
      WDOT(5) = WDOT(5) +ROP
      WDOT(8) = WDOT(8) -ROP
      WDOT(17) = WDOT(17) +ROP
      ROP = RF(135)-RB(135)
      WDOT(23) = WDOT(23) +ROP
      ROP = RF(136)-RB(136)
      WDOT(10) = WDOT(10) +ROP
      WDOT(19) = WDOT(19) +ROP
      ROP = RF(137)-RB(137)
      WDOT(2) = WDOT(2) -ROP
      WDOT(6) = WDOT(6) +ROP
      WDOT(17) = WDOT(17) -ROP
      ROP = RF(138)-RB(138)
      WDOT(3) = WDOT(3) -ROP
      WDOT(4) = WDOT(4) +ROP
      WDOT(17) = WDOT(17) -ROP
      ROP = RF(139)-RB(139)
      WDOT(4) = WDOT(4) -ROP
      WDOT(7) = WDOT(7) +ROP
      WDOT(17) = WDOT(17) -ROP
      ROP = RF(140)-RB(140)
      WDOT(10) = WDOT(10) +ROP
      WDOT(17) = WDOT(17) -ROP
      ROP = RF(141)-RB(141)
      WDOT(10) = WDOT(10) -ROP
      WDOT(11) = WDOT(11) +ROP
      WDOT(17) = WDOT(17) -ROP
      ROP = RF(142)-RB(142)
      WDOT(2) = WDOT(2) -ROP
      WDOT(19) = WDOT(19) -ROP
      WDOT(20) = WDOT(20) +ROP
      ROP = RF(143)-RB(143)
      WDOT(2) = WDOT(2) +ROP
      WDOT(3) = WDOT(3) -ROP
      WDOT(19) = WDOT(19) -ROP
      WDOT(21) = WDOT(21) +ROP
      ROP = RF(144)-RB(144)
      WDOT(2) = WDOT(2) +2*ROP
      WDOT(4) = WDOT(4) -ROP
      WDOT(19) = WDOT(19) -ROP
      WDOT(21) = WDOT(21) +ROP
      ROP = RF(145)-RB(145)
      WDOT(4) = WDOT(4) +ROP
      WDOT(9) = WDOT(9) -ROP
      WDOT(19) = WDOT(19) -ROP
      WDOT(21) = WDOT(21) +ROP
      ROP = RF(146)-RB(146)
      WDOT(5) = WDOT(5) -ROP
      WDOT(9) = WDOT(9) +ROP
      WDOT(19) = WDOT(19) -ROP
      WDOT(20) = WDOT(20) +ROP
      ROP = RF(147)-RB(147)
      WDOT(4) = WDOT(4) +ROP
      WDOT(5) = WDOT(5) -ROP
      WDOT(12) = WDOT(12) +ROP
      WDOT(19) = WDOT(19) -ROP
      ROP = RF(148)-RB(148)
      WDOT(13) = WDOT(13) +ROP
      WDOT(19) = WDOT(19) -ROP
      WDOT(20) = WDOT(20) +ROP
      ROP = RF(149)-RB(149)
      WDOT(10) = WDOT(10) -ROP
      WDOT(19) = WDOT(19) -ROP
      WDOT(23) = WDOT(23) +ROP
      ROP = RF(150)-RB(150)
      WDOT(2) = WDOT(2) -ROP
      WDOT(20) = WDOT(20) -ROP
      ROP = RF(151)-RB(151)
      WDOT(2) = WDOT(2) -ROP
      WDOT(10) = WDOT(10) +ROP
      WDOT(16) = WDOT(16) +ROP
      WDOT(20) = WDOT(20) -ROP
      ROP = RF(152)-RB(152)
      WDOT(2) = WDOT(2) -ROP
      WDOT(6) = WDOT(6) +ROP
      WDOT(19) = WDOT(19) +ROP
      WDOT(20) = WDOT(20) -ROP
      ROP = RF(153)-RB(153)
      WDOT(2) = WDOT(2) +2*ROP
      WDOT(3) = WDOT(3) -ROP
      WDOT(20) = WDOT(20) -ROP
      WDOT(21) = WDOT(21) +ROP
      ROP = RF(154)-RB(154)
      WDOT(3) = WDOT(3) -ROP
      WDOT(20) = WDOT(20) -ROP
      ROP = RF(155)-RB(155)
      WDOT(3) = WDOT(3) -ROP
      WDOT(4) = WDOT(4) +ROP
      WDOT(19) = WDOT(19) +ROP
      WDOT(20) = WDOT(20) -ROP
      ROP = RF(156)-RB(156)
      WDOT(4) = WDOT(4) -ROP
      WDOT(7) = WDOT(7) +ROP
      WDOT(19) = WDOT(19) +ROP
      WDOT(20) = WDOT(20) -ROP
      ROP = RF(157)-RB(157)
      WDOT(5) = WDOT(5) -ROP
      WDOT(8) = WDOT(8) +ROP
      WDOT(19) = WDOT(19) +ROP
      WDOT(20) = WDOT(20) -ROP
      ROP = RF(158)-RB(158)
      WDOT(10) = WDOT(10) -ROP
      WDOT(11) = WDOT(11) +ROP
      WDOT(19) = WDOT(19) +ROP
      WDOT(20) = WDOT(20) -ROP
      ROP = RF(159)-RB(159)
      WDOT(2) = WDOT(2) -ROP
      WDOT(16) = WDOT(16) +ROP
      WDOT(21) = WDOT(21) -ROP
      ROP = RF(160)-RB(160)
      WDOT(3) = WDOT(3) -ROP
      WDOT(4) = WDOT(4) +ROP
      WDOT(13) = WDOT(13) +ROP
      WDOT(21) = WDOT(21) -ROP
      ROP = RF(161)-RB(161)
      WDOT(4) = WDOT(4) -ROP
      WDOT(7) = WDOT(7) +ROP
      WDOT(13) = WDOT(13) +ROP
      WDOT(21) = WDOT(21) -ROP
      ROP = RF(162)-RB(162)
      WDOT(2) = WDOT(2) -ROP
      WDOT(10) = WDOT(10) +ROP
      ROP = RF(163)-RB(163)
      WDOT(2) = WDOT(2) -ROP
      WDOT(6) = WDOT(6) +ROP
      WDOT(20) = WDOT(20) +ROP
      ROP = RF(164)-RB(164)
      WDOT(3) = WDOT(3) -ROP
      WDOT(12) = WDOT(12) +ROP
      ROP = RF(165)-RB(165)
      WDOT(4) = WDOT(4) -ROP
      WDOT(7) = WDOT(7) +ROP
      WDOT(20) = WDOT(20) +ROP
      ROP = RF(166)-RB(166)
      WDOT(5) = WDOT(5) +ROP
      WDOT(9) = WDOT(9) -ROP
      WDOT(20) = WDOT(20) +ROP
      ROP = RF(167)-RB(167)
      WDOT(4) = WDOT(4) +ROP
      WDOT(5) = WDOT(5) -ROP
      WDOT(12) = WDOT(12) +ROP
      ROP = RF(168)-RB(168)
      WDOT(10) = WDOT(10) -ROP
      WDOT(11) = WDOT(11) +ROP
      WDOT(20) = WDOT(20) +ROP
      ROP = RF(169)-RB(169)
      WDOT(2) = WDOT(2) -ROP
      WDOT(22) = WDOT(22) -ROP
      WDOT(23) = WDOT(23) +ROP
      ROP = RF(170)-RB(170)
      WDOT(2) = WDOT(2) -ROP
      WDOT(10) = WDOT(10) +ROP
      WDOT(19) = WDOT(19) +ROP
      WDOT(22) = WDOT(22) -ROP
      ROP = RF(171)-RB(171)
      WDOT(4) = WDOT(4) +ROP
      WDOT(5) = WDOT(5) -ROP
      WDOT(12) = WDOT(12) +ROP
      WDOT(19) = WDOT(19) +ROP
      WDOT(22) = WDOT(22) -ROP
      ROP = RF(172)-RB(172)
      WDOT(13) = WDOT(13) +ROP
      WDOT(22) = WDOT(22) -ROP
      WDOT(23) = WDOT(23) +ROP
      ROP = RF(173)-RB(173)
      WDOT(2) = WDOT(2) -ROP
      WDOT(23) = WDOT(23) -ROP
      ROP = RF(174)-RB(174)
      WDOT(2) = WDOT(2) -ROP
      WDOT(16) = WDOT(16) +ROP
      WDOT(23) = WDOT(23) -ROP
      ROP = RF(175)-RB(175)
      WDOT(2) = WDOT(2) -ROP
      WDOT(10) = WDOT(10) +ROP
      WDOT(20) = WDOT(20) +ROP
      WDOT(23) = WDOT(23) -ROP
      ROP = RF(176)-RB(176)
      WDOT(2) = WDOT(2) -ROP
      WDOT(6) = WDOT(6) +ROP
      WDOT(22) = WDOT(22) +ROP
      WDOT(23) = WDOT(23) -ROP
      ROP = RF(177)-RB(177)
      WDOT(3) = WDOT(3) -ROP
      WDOT(23) = WDOT(23) -ROP
      ROP = RF(178)-RB(178)
      WDOT(3) = WDOT(3) -ROP
      WDOT(4) = WDOT(4) +ROP
      WDOT(22) = WDOT(22) +ROP
      WDOT(23) = WDOT(23) -ROP
      ROP = RF(179)-RB(179)
      WDOT(3) = WDOT(3) -ROP
      WDOT(4) = WDOT(4) +ROP
      WDOT(22) = WDOT(22) +ROP
      WDOT(23) = WDOT(23) -ROP
      ROP = RF(180)-RB(180)
      WDOT(4) = WDOT(4) -ROP
      WDOT(7) = WDOT(7) +ROP
      WDOT(22) = WDOT(22) +ROP
      WDOT(23) = WDOT(23) -ROP
      ROP = RF(181)-RB(181)
      WDOT(5) = WDOT(5) +ROP
      WDOT(9) = WDOT(9) -ROP
      WDOT(22) = WDOT(22) +ROP
      WDOT(23) = WDOT(23) -ROP
      ROP = RF(182)-RB(182)
      WDOT(5) = WDOT(5) -ROP
      WDOT(8) = WDOT(8) +ROP
      WDOT(22) = WDOT(22) +ROP
      WDOT(23) = WDOT(23) -ROP
      ROP = RF(183)-RB(183)
      WDOT(10) = WDOT(10) -ROP
      WDOT(11) = WDOT(11) +ROP
      WDOT(22) = WDOT(22) +ROP
      WDOT(23) = WDOT(23) -ROP
      ROP = RF(184)-RB(184)
      WDOT(2) = WDOT(2) -ROP
      ROP = RF(185)-RB(185)
      WDOT(2) = WDOT(2) -ROP
      WDOT(6) = WDOT(6) +ROP
      WDOT(23) = WDOT(23) +ROP
      ROP = RF(186)-RB(186)
      WDOT(3) = WDOT(3) -ROP
      WDOT(12) = WDOT(12) +ROP
      ROP = RF(187)-RB(187)
      WDOT(4) = WDOT(4) -ROP
      WDOT(7) = WDOT(7) +ROP
      WDOT(23) = WDOT(23) +ROP
      ROP = RF(188)-RB(188)
      WDOT(5) = WDOT(5) +ROP
      WDOT(9) = WDOT(9) -ROP
      WDOT(23) = WDOT(23) +ROP
      ROP = RF(189)-RB(189)
      WDOT(4) = WDOT(4) +ROP
      WDOT(5) = WDOT(5) -ROP
      WDOT(12) = WDOT(12) +ROP
      ROP = RF(190)-RB(190)
      WDOT(10) = WDOT(10) -ROP
      WDOT(11) = WDOT(11) +ROP
      WDOT(23) = WDOT(23) +ROP
      ROP = RF(191)
      WDOT(16) = WDOT(16) +ROP
      WDOT(19) = WDOT(19) +ROP
      WDOT(24) = WDOT(24) -ROP
      ROP = RF(192)
      WDOT(20) = WDOT(20) +ROP
      WDOT(24) = WDOT(24) -ROP
      ROP = RF(193)-RB(193)
      WDOT(2) = WDOT(2) -ROP
      WDOT(25) = WDOT(25) -ROP
      ROP = RF(194)-RB(194)
      WDOT(2) = WDOT(2) -ROP
      WDOT(16) = WDOT(16) +ROP
      WDOT(25) = WDOT(25) -ROP
      ROP = RF(195)-RB(195)
      WDOT(2) = WDOT(2) -ROP
      WDOT(20) = WDOT(20) +ROP
      WDOT(25) = WDOT(25) -ROP
      ROP = RF(196)-RB(196)
      WDOT(16) = WDOT(16) -ROP
      ROP = RF(197)-RB(197)
      WDOT(2) = WDOT(2) -ROP
      WDOT(26) = WDOT(26) -ROP
      ROP = RF(198)-RB(198)
      WDOT(2) = WDOT(2) -ROP
      WDOT(16) = WDOT(16) +ROP
      WDOT(26) = WDOT(26) -ROP
      ROP = RF(199)-RB(199)
      WDOT(2) = WDOT(2) -ROP
      WDOT(20) = WDOT(20) +ROP
      WDOT(26) = WDOT(26) -ROP
      ROP = RF(200)-RB(200)
      WDOT(16) = WDOT(16) -ROP
      ROP = RF(201)-RB(201)
      WDOT(2) = WDOT(2) -ROP
      WDOT(27) = WDOT(27) -ROP
      ROP = RF(202)-RB(202)
      WDOT(2) = WDOT(2) -ROP
      WDOT(16) = WDOT(16) +ROP
      WDOT(27) = WDOT(27) -ROP
      ROP = RF(203)-RB(203)
      WDOT(2) = WDOT(2) -ROP
      WDOT(20) = WDOT(20) +ROP
      WDOT(27) = WDOT(27) -ROP
      ROP = RF(204)-RB(204)
      WDOT(16) = WDOT(16) -ROP
      ROP = RF(205)-RB(205)
      WDOT(2) = WDOT(2) -ROP
      WDOT(28) = WDOT(28) -ROP
      ROP = RF(206)-RB(206)
      WDOT(2) = WDOT(2) -ROP
      WDOT(16) = WDOT(16) +ROP
      WDOT(28) = WDOT(28) -ROP
      ROP = RF(207)-RB(207)
      WDOT(2) = WDOT(2) -ROP
      WDOT(20) = WDOT(20) +ROP
      WDOT(28) = WDOT(28) -ROP
      ROP = RF(208)-RB(208)
      WDOT(16) = WDOT(16) -ROP
      ROP = RF(209)-RB(209)
      WDOT(2) = WDOT(2) -ROP
      WDOT(29) = WDOT(29) -ROP
      WDOT(30) = WDOT(30) +ROP
      ROP = RF(210)-RB(210)
      WDOT(2) = WDOT(2) -ROP
      WDOT(16) = WDOT(16) +ROP
      WDOT(29) = WDOT(29) -ROP
      ROP = RF(211)-RB(211)
      WDOT(2) = WDOT(2) -ROP
      WDOT(20) = WDOT(20) +ROP
      WDOT(29) = WDOT(29) -ROP
      ROP = RF(212)-RB(212)
      WDOT(16) = WDOT(16) -ROP
      WDOT(30) = WDOT(30) +ROP
      ROP = RF(213)-RB(213)
      WDOT(2) = WDOT(2) -ROP
      WDOT(31) = WDOT(31) -ROP
      ROP = RF(214)-RB(214)
      WDOT(2) = WDOT(2) -ROP
      WDOT(16) = WDOT(16) +ROP
      WDOT(31) = WDOT(31) -ROP
      ROP = RF(215)-RB(215)
      WDOT(2) = WDOT(2) -ROP
      WDOT(20) = WDOT(20) +ROP
      WDOT(31) = WDOT(31) -ROP
      ROP = RF(216)-RB(216)
      WDOT(16) = WDOT(16) -ROP
      ROP = RF(217)-RB(217)
      WDOT(24) = WDOT(24) +ROP
      WDOT(32) = WDOT(32) -ROP
      ROP = RF(218)-RB(218)
      WDOT(16) = WDOT(16) -ROP
      ROP = RF(219)-RB(219)
      ROP = RF(220)-RB(220)
      WDOT(20) = WDOT(20) -ROP
      WDOT(30) = WDOT(30) -ROP
      ROP = RF(221)-RB(221)
      WDOT(23) = WDOT(23) -ROP
      ROP = RF(222)-RB(222)
      WDOT(25) = WDOT(25) -ROP
      ROP = RF(223)-RB(223)
      WDOT(31) = WDOT(31) -ROP
      ROP = RF(224)-RB(224)
      WDOT(26) = WDOT(26) -ROP
      ROP = RF(225)-RB(225)
      WDOT(29) = WDOT(29) -ROP
      ROP = RF(226)-RB(226)
      WDOT(27) = WDOT(27) -ROP
      ROP = RF(227)-RB(227)
      WDOT(28) = WDOT(28) -ROP
      ROP = RF(228)-RB(228)
      WDOT(1) = WDOT(1) +ROP
      ROP = RF(229)-RB(229)
      WDOT(1) = WDOT(1) +ROP
      WDOT(30) = WDOT(30) -ROP
      ROP = RF(230)-RB(230)
      WDOT(1) = WDOT(1) +ROP
      ROP = RF(231)-RB(231)
      WDOT(1) = WDOT(1) +ROP
      ROP = RF(232)-RB(232)
      WDOT(1) = WDOT(1) +ROP
      ROP = RF(233)-RB(233)
      WDOT(1) = WDOT(1) -ROP
      WDOT(2) = WDOT(2) -ROP
      WDOT(6) = WDOT(6) +ROP
      ROP = RF(234)-RB(234)
      WDOT(1) = WDOT(1) -ROP
      WDOT(2) = WDOT(2) -ROP
      WDOT(6) = WDOT(6) +ROP
      ROP = RF(235)-RB(235)
      WDOT(1) = WDOT(1) -ROP
      WDOT(2) = WDOT(2) -ROP
      WDOT(6) = WDOT(6) +ROP
      ROP = RF(236)-RB(236)
      WDOT(1) = WDOT(1) -ROP
      WDOT(3) = WDOT(3) -ROP
      WDOT(4) = WDOT(4) +ROP
      ROP = RF(237)-RB(237)
      WDOT(1) = WDOT(1) -ROP
      WDOT(3) = WDOT(3) -ROP
      WDOT(4) = WDOT(4) +ROP
      ROP = RF(238)-RB(238)
      WDOT(1) = WDOT(1) -ROP
      WDOT(3) = WDOT(3) -ROP
      WDOT(4) = WDOT(4) +ROP
      ROP = RF(239)-RB(239)
      WDOT(1) = WDOT(1) -ROP
      WDOT(4) = WDOT(4) -ROP
      WDOT(7) = WDOT(7) +ROP
      ROP = RF(240)-RB(240)
      WDOT(1) = WDOT(1) -ROP
      WDOT(4) = WDOT(4) -ROP
      WDOT(7) = WDOT(7) +ROP
      ROP = RF(241)-RB(241)
      WDOT(1) = WDOT(1) -ROP
      WDOT(4) = WDOT(4) -ROP
      WDOT(7) = WDOT(7) +ROP
      ROP = RF(242)-RB(242)
      WDOT(1) = WDOT(1) -ROP
      WDOT(5) = WDOT(5) +ROP
      WDOT(9) = WDOT(9) -ROP
      ROP = RF(243)-RB(243)
      WDOT(1) = WDOT(1) -ROP
      WDOT(5) = WDOT(5) +ROP
      WDOT(9) = WDOT(9) -ROP
      ROP = RF(244)-RB(244)
      WDOT(1) = WDOT(1) -ROP
      WDOT(5) = WDOT(5) +ROP
      WDOT(9) = WDOT(9) -ROP
      ROP = RF(245)-RB(245)
      WDOT(1) = WDOT(1) -ROP
      WDOT(5) = WDOT(5) -ROP
      WDOT(8) = WDOT(8) +ROP
      ROP = RF(246)-RB(246)
      WDOT(1) = WDOT(1) -ROP
      WDOT(5) = WDOT(5) -ROP
      WDOT(8) = WDOT(8) +ROP
      ROP = RF(247)-RB(247)
      WDOT(1) = WDOT(1) -ROP
      WDOT(5) = WDOT(5) -ROP
      WDOT(8) = WDOT(8) +ROP
      ROP = RF(248)-RB(248)
      WDOT(1) = WDOT(1) -ROP
      WDOT(10) = WDOT(10) -ROP
      WDOT(11) = WDOT(11) +ROP
      ROP = RF(249)-RB(249)
      WDOT(1) = WDOT(1) -ROP
      WDOT(10) = WDOT(10) -ROP
      WDOT(11) = WDOT(11) +ROP
      ROP = RF(250)-RB(250)
      WDOT(1) = WDOT(1) -ROP
      WDOT(10) = WDOT(10) -ROP
      WDOT(11) = WDOT(11) +ROP
      ROP = RF(251)
      WDOT(9) = WDOT(9) -ROP
      WDOT(33) = WDOT(33) +ROP
      ROP = RF(252)
      WDOT(9) = WDOT(9) +ROP
      WDOT(33) = WDOT(33) -ROP
      ROP = RF(253)
      WDOT(9) = WDOT(9) -ROP
      WDOT(33) = WDOT(33) +ROP
      ROP = RF(254)
      WDOT(9) = WDOT(9) +ROP
      WDOT(33) = WDOT(33) -ROP
      ROP = RF(255)
      WDOT(9) = WDOT(9) -ROP
      WDOT(33) = WDOT(33) +ROP
      ROP = RF(256)
      WDOT(9) = WDOT(9) +ROP
      WDOT(33) = WDOT(33) -ROP
      ROP = RF(257)
      WDOT(33) = WDOT(33) -ROP
      ROP = RF(258)
      WDOT(33) = WDOT(33) +ROP
      ROP = RF(259)
      WDOT(5) = WDOT(5) +ROP
      WDOT(9) = WDOT(9) -ROP
      WDOT(32) = WDOT(32) +ROP
      ROP = RF(260)
      WDOT(5) = WDOT(5) -ROP
      WDOT(9) = WDOT(9) +ROP
      WDOT(32) = WDOT(32) -ROP
      ROP = RF(261)
      WDOT(5) = WDOT(5) +ROP
      WDOT(9) = WDOT(9) -ROP
      WDOT(32) = WDOT(32) +ROP
      ROP = RF(262)
      WDOT(5) = WDOT(5) -ROP
      WDOT(9) = WDOT(9) +ROP
      WDOT(32) = WDOT(32) -ROP
      ROP = RF(263)
      WDOT(5) = WDOT(5) +ROP
      WDOT(9) = WDOT(9) -ROP
      WDOT(32) = WDOT(32) +ROP
      ROP = RF(264)
      WDOT(5) = WDOT(5) -ROP
      WDOT(9) = WDOT(9) +ROP
      WDOT(32) = WDOT(32) -ROP
      ROP = RF(265)
      WDOT(9) = WDOT(9) -ROP
      ROP = RF(266)
      WDOT(9) = WDOT(9) +ROP
      ROP = RF(267)-RB(267)
      WDOT(4) = WDOT(4) +ROP
      WDOT(34) = WDOT(34) +ROP
      ROP = RF(268)
      WDOT(4) = WDOT(4) +ROP
      WDOT(16) = WDOT(16) +3*ROP
      WDOT(18) = WDOT(18) +2*ROP
      WDOT(34) = WDOT(34) -ROP
C
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE QSSA(RF, RB, XQ)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
      DIMENSION RF(*), RB(*), XQ(*)
      DOUBLE PRECISION SMALL
      PARAMETER (SMALL = 1.0D-50)
C
      RF(47) = 0.D0
      RF(105) = 0.D0
      RF(106) = 0.D0
      RF(110) = 0.D0
      RF(135) = 0.D0
      RF(136) = 0.D0
      RB(154) = 0.D0
      RB(177) = 0.D0
      RB(184) = 0.D0
      RF(228) = 0.D0
      RF(230) = 0.D0
      RF(231) = 0.D0
      RF(232) = 0.D0
C
C     CH2
      DEN = +RF( 40) +RF( 41) +RF( 42) +RF( 43) +RF( 44) 
     *  +RF( 45) +RF( 46) +RF( 76) +RF( 89) +RF(121) +RB( 48) 
     *  +RB( 55) +RB( 56) +RB( 57) +RB( 67) +RB( 92) +RB(118) 
      A1_0 = ( +RB( 40) +RB( 42) +RB( 43) +RB( 45) +RB( 46) 
     *  +2*RB( 47) +RF( 67) +RB( 76) +RB( 89) +RF( 92) +RF(118) 
     *  +RB(121) )/DMAX1(DEN,SMALL)
      A1_2 = ( +RF( 48) +RF( 55) +RF( 56) +RF( 57) )/DMAX1(DEN,SMALL)
      A1_3 = ( +RB( 41) +RB( 44) )/DMAX1(DEN,SMALL)
C     CH2*
      DEN = +RF( 48) +RF( 49) +RF( 50) +RF( 51) +RF( 52) 
     *  +RF( 53) +RF( 54) +RF( 55) +RF( 56) +RF( 57) +RF( 58) 
     *  +RF( 77) +RF( 90) +RF(122) +RF(140) +RB( 68) +RB( 82) 
      A2_0 = ( +RB( 49) +RB( 51) +RB( 52) +RB( 53) +RB( 54) 
     *  +RB( 58) +RF( 68) +RB( 77) +RB( 90) +RB(122) )/DMAX1(DEN,SMALL)
      A2_1 = ( +RB( 48) +RB( 55) +RB( 56) +RB( 57) )/DMAX1(DEN,SMALL)
      A2_3 = ( +RB( 50) )/DMAX1(DEN,SMALL)
      A2_4 = ( +RF( 82) )/DMAX1(DEN,SMALL)
      A2_6 = ( +RB(140) )/DMAX1(DEN,SMALL)
C     HCO
      DEN = +RF( 31) +RF( 32) +RF( 33) +RF( 34) +RF( 35) 
     *  +RF( 36) +RF( 37) +RF( 39) +RF( 74) +RF( 94) +RF(120) 
     *  +RF(148) +RF(172) +RB( 41) +RB( 44) +RB( 50) +RB( 60) 
     *  +RB( 61) +RB( 62) +RB( 63) +RB( 64) +RB( 75) +RB(102) 
     *  +RB(112) +RB(117) +RB(159) 
      A3_0 = ( +RB( 31) +RB( 32) +RB( 33) +RB( 34) +RB( 35) 
     *  +RB( 36) +RB( 37) +RB( 39) +RF( 60) +RF( 61) +RF( 62) 
     *  +RF( 63) +RF( 64) +RB( 74) +RF( 75) +RB(105) +RB(106) 
     *  +RF(112) +RF(117) +RB(148) +RF(154) +RF(159) +RB(172) 
     *  +RF(177) )/DMAX1(DEN,SMALL)
      A3_1 = ( +RF( 41) +RF( 44) )/DMAX1(DEN,SMALL)
      A3_2 = ( +RF( 50) )/DMAX1(DEN,SMALL)
      A3_5 = ( +RB( 94) +RF(102) )/DMAX1(DEN,SMALL)
      A3_6 = ( +RB(120) )/DMAX1(DEN,SMALL)
C     CH3O
      DEN = +RF( 80) +RF( 81) +RF( 82) +RF( 83) +RF( 84) 
     *  +RF( 85) +RB( 59) +RB( 69) +RB( 72) 
      A4_0 = ( +RF( 59) +RF( 69) +RF( 72) +RB( 80) +RB( 81) 
     *  +RB( 83) +RB( 84) +RB( 85) )/DMAX1(DEN,SMALL)
      A4_2 = ( +RB( 82) )/DMAX1(DEN,SMALL)
C     C2H3
      DEN = +RF( 91) +RF( 96) +RF( 97) +RF( 98) +RF( 99) 
     *  +RF(100) +RF(101) +RF(102) +RF(103) +RF(104) +RF(107) 
     *  +RF(108) +RF(109) +RF(126) +RB( 94) +RB(115) +RB(116) 
     *  +RB(119) +RB(123) +RB(125) +RB(147) +RB(160) +RB(161) 
     *  +RB(192) 
      A5_0 = ( +RB( 91) +RB( 96) +RB( 97) +RB( 98) +RB( 99) 
     *  +RB(100) +RB(101) +RB(103) +RB(104) +RB(105) +RB(106) 
     *  +RB(107) +RB(108) +RB(109) +2*RB(110) +RF(115) +RF(116) 
     *  +RF(119) +RF(123) +RF(125) +RB(126) +RB(135) +RB(136) 
     *  +RF(147) +RF(160) +RF(161) +RF(192) )/DMAX1(DEN,SMALL)
      A5_3 = ( +RF( 94) +RB(102) )/DMAX1(DEN,SMALL)
C     C2H5
      DEN = +RF(127) +RF(128) +RF(129) +RF(130) +RF(131) 
     *  +RF(132) +RF(133) +RF(134) +RF(223) +RB( 79) +RB(114) 
     *  +RB(120) +RB(137) +RB(138) +RB(139) +RB(140) +RB(141) 
     *  +RB(162) +RB(164) +RB(167) +RB(174) +RB(195) +RB(268) 
      A6_0 = ( +RF( 79) +RF(114) +RB(127) +RB(128) +RB(129) 
     *  +RB(130) +RB(131) +RB(132) +RB(133) +RB(134) +RB(135) 
     *  +RB(136) +RF(137) +RF(138) +RF(139) +RF(141) +RF(154) 
     *  +RF(174) +RF(195) +RB(228) +RF(268) )/DMAX1(DEN,SMALL)
      A6_2 = ( +RF(140) )/DMAX1(DEN,SMALL)
      A6_3 = ( +RF(120) )/DMAX1(DEN,SMALL)
      A6_7 = ( +RF(162) +RF(164) +RF(167) )/DMAX1(DEN,SMALL)
      A6_8 = ( +2*RF(184) )/DMAX1(DEN,SMALL)
      A6_16 = ( +RB(223) )/DMAX1(DEN,SMALL)
C     nC3H7
      DEN = +RF(124) +RF(162) +RF(163) +RF(164) +RF(165) 
     *  +RF(166) +RF(167) +RF(168) +RF(196) +RF(225) +RF(229) 
     *  +RB(150) +RB(186) +RB(189) +RB(194) +RB(199) 
      A7_0 = ( +RB(124) +RF(150) +RB(163) +RB(165) +RB(166) 
     *  +RB(168) +RF(177) +RF(194) +RF(199) +RB(229) )/DMAX1(DEN,SMALL)
      A7_6 = ( +RB(162) +RB(164) +RB(167) )/DMAX1(DEN,SMALL)
      A7_8 = ( +RF(186) +RF(189) )/DMAX1(DEN,SMALL)
      A7_9 = ( +RB(196) )/DMAX1(DEN,SMALL)
      A7_16 = ( +RB(225) )/DMAX1(DEN,SMALL)
C     pC4H9
      DEN = +RF(184) +RF(185) +RF(186) +RF(187) +RF(188) 
     *  +RF(189) +RF(190) +RF(200) +RF(227) +RB(173) +RB(198) 
     *  +RB(203) 
      A8_0 = ( +RF(173) +RB(185) +RB(187) +RB(188) +RB(190) 
     *  +RF(198) +RF(203) +RB(230) )/DMAX1(DEN,SMALL)
      A8_7 = ( +RB(186) +RB(189) )/DMAX1(DEN,SMALL)
      A8_10 = ( +RB(200) )/DMAX1(DEN,SMALL)
      A8_16 = ( +RB(227) )/DMAX1(DEN,SMALL)
C     PXC5H11
      DEN = +RF(204) +RF(226) +RB(193) +RB(196) +RB(202) 
     *  +RB(207) 
      A9_0 = ( +RF(193) +RF(202) +RF(207) +RB(231) )/DMAX1(DEN,SMALL)
      A9_7 = ( +RF(196) )/DMAX1(DEN,SMALL)
      A9_11 = ( +RB(204) )/DMAX1(DEN,SMALL)
      A9_16 = ( +RB(226) )/DMAX1(DEN,SMALL)
C     PXC6H13
      DEN = +RF(208) +RF(224) +RB(197) +RB(200) +RB(206) 
     *  +RB(211) 
      A10_0 = ( +RF(197) +RF(206) +RF(211) +2*RB(232) )/DMAX1(DEN,SMALL)
      A10_8 = ( +RF(200) )/DMAX1(DEN,SMALL)
      A10_12 = ( +RB(208) )/DMAX1(DEN,SMALL)
      A10_16 = ( +RB(224) )/DMAX1(DEN,SMALL)
C     PXC7H15
      DEN = +RF(212) +RF(222) +RB(201) +RB(204) +RB(210) 
     *  +RB(215) +RB(217) 
      A11_0 = ( +RF(201) +RF(210) +RB(212) +RF(215) +RF(217) 
     *  +RB(231) )/DMAX1(DEN,SMALL)
      A11_9 = ( +RF(204) )/DMAX1(DEN,SMALL)
      A11_16 = ( +RB(222) )/DMAX1(DEN,SMALL)
C     PXC8H17
      DEN = +RF(216) +RF(221) +RB(205) +RB(208) +RB(214) 
      A12_0 = ( +RF(205) +RF(214) +RB(230) )/DMAX1(DEN,SMALL)
      A12_10 = ( +RF(208) )/DMAX1(DEN,SMALL)
      A12_13 = ( +RB(216) )/DMAX1(DEN,SMALL)
      A12_15 = ( +RB(221) )/DMAX1(DEN,SMALL)
C     PXC10H21
      DEN = +RF(218) +RB(213) +RB(216) 
      A13_0 = ( +RF(213) +RB(228) )/DMAX1(DEN,SMALL)
      A13_12 = ( +RF(216) )/DMAX1(DEN,SMALL)
      A13_14 = ( +RB(218) )/DMAX1(DEN,SMALL)
C     PXC12H25
      DEN = +RF(219) +RF(251) +RF(259) +RB(218) +RB(233) 
     *  +RB(236) +RB(239) +RB(242) +RB(245) +RB(248) +RB(252) 
     *  +RB(260) 
      A14_0 = ( +RF(233) +RF(236) +RF(239) +RF(242) +RF(245) 
     *  +RF(248) +RB(251) +RF(252) +RB(259) +RF(260) )/DMAX1(DEN,SMALL)
      A14_13 = ( +RF(218) )/DMAX1(DEN,SMALL)
      A14_16 = ( +RB(219) )/DMAX1(DEN,SMALL)
C     SXC12H25
      DEN = +RF(253) +RF(261) +RB(220) +RB(221) +RB(234) 
     *  +RB(237) +RB(240) +RB(243) +RB(246) +RB(249) +RB(254) 
     *  +RB(262) 
      A15_0 = ( +RF(220) +RF(234) +RF(237) +RF(240) +RF(243) 
     *  +RF(246) +RF(249) +RB(253) +RF(254) +RB(261) +RF(262) )
     *  /DMAX1(DEN,SMALL)
      A15_12 = ( +RF(221) )/DMAX1(DEN,SMALL)
C     S3XC12H25
      DEN = +RF(255) +RF(263) +RB(219) +RB(222) +RB(223) 
     *  +RB(224) +RB(225) +RB(226) +RB(227) +RB(235) +RB(238) 
     *  +RB(241) +RB(244) +RB(247) +RB(250) +RB(256) +RB(264) 
      A16_0 = ( +RF(235) +RF(238) +RF(241) +RF(244) +RF(247) 
     *  +RF(250) +RB(255) +RF(256) +RB(263) +RF(264) )/DMAX1(DEN,SMALL)
      A16_6 = ( +RF(223) )/DMAX1(DEN,SMALL)
      A16_7 = ( +RF(225) )/DMAX1(DEN,SMALL)
      A16_8 = ( +RF(227) )/DMAX1(DEN,SMALL)
      A16_9 = ( +RF(226) )/DMAX1(DEN,SMALL)
      A16_10 = ( +RF(224) )/DMAX1(DEN,SMALL)
      A16_11 = ( +RF(222) )/DMAX1(DEN,SMALL)
      A16_14 = ( +RF(219) )/DMAX1(DEN,SMALL)
C     C12OOH
      DEN = +RF(258) +RF(265) +RB(257) +RB(266) 
      A17_0 = ( +RF(257) +RB(258) )/DMAX1(DEN,SMALL)
      A17_18 = ( +RB(265) +RF(266) )/DMAX1(DEN,SMALL)
C     O2C12H24OOH
      DEN = +RF(266) +RF(267) +RB(265) 
      A18_0 = ( +RB(267) )/DMAX1(DEN,SMALL)
      A18_17 = ( +RF(265) +RB(266) )/DMAX1(DEN,SMALL)
C
      A12_0 = A12_0 + A12_15*A15_0
      DEN = 1 -A12_15*A15_12
      A12_0 = A12_0/DMAX1(DEN,SMALL)
      A12_10 = A12_10/DMAX1(DEN,SMALL)
      A12_13 = A12_13/DMAX1(DEN,SMALL)
      A2_0 = A2_0 + A2_4*A4_0
      DEN = 1 -A2_4*A4_2
      A2_0 = A2_0/DMAX1(DEN,SMALL)
      A2_6 = A2_6/DMAX1(DEN,SMALL)
      A2_3 = A2_3/DMAX1(DEN,SMALL)
      A2_1 = A2_1/DMAX1(DEN,SMALL)
      A3_0 = A3_0 + A3_5*A5_0
      DEN = 1 -A3_5*A5_3
      A3_0 = A3_0/DMAX1(DEN,SMALL)
      A3_6 = A3_6/DMAX1(DEN,SMALL)
      A3_2 = A3_2/DMAX1(DEN,SMALL)
      A3_1 = A3_1/DMAX1(DEN,SMALL)
      A12_0 = A12_0 + A12_13*A13_0
      A12_14 = A12_13*A13_14
      DEN = 1 -A12_13*A13_12
      A12_0 = A12_0/DMAX1(DEN,SMALL)
      A12_10 = A12_10/DMAX1(DEN,SMALL)
      A12_14 = A12_14/DMAX1(DEN,SMALL)
      A14_0 = A14_0 + A14_13*A13_0
      A14_12 = A14_13*A13_12
      DEN = 1 -A14_13*A13_14
      A14_0 = A14_0/DMAX1(DEN,SMALL)
      A14_16 = A14_16/DMAX1(DEN,SMALL)
      A14_12 = A14_12/DMAX1(DEN,SMALL)
      A16_0 = A16_0 + A16_14*A14_0
      A16_12 = A16_14*A14_12
      DEN = 1 -A16_14*A14_16
      A16_0 = A16_0/DMAX1(DEN,SMALL)
      A16_6 = A16_6/DMAX1(DEN,SMALL)
      A16_7 = A16_7/DMAX1(DEN,SMALL)
      A16_9 = A16_9/DMAX1(DEN,SMALL)
      A16_8 = A16_8/DMAX1(DEN,SMALL)
      A16_12 = A16_12/DMAX1(DEN,SMALL)
      A16_10 = A16_10/DMAX1(DEN,SMALL)
      A16_11 = A16_11/DMAX1(DEN,SMALL)
      A12_0 = A12_0 + A12_14*A14_0
      A12_16 = A12_14*A14_16
      DEN = 1 -A12_14*A14_12
      A12_0 = A12_0/DMAX1(DEN,SMALL)
      A12_16 = A12_16/DMAX1(DEN,SMALL)
      A12_10 = A12_10/DMAX1(DEN,SMALL)
      A16_0 = A16_0 + A16_11*A11_0
      A16_9 = A16_9 + A16_11*A11_9
      DEN = 1 -A16_11*A11_16
      A16_0 = A16_0/DMAX1(DEN,SMALL)
      A16_6 = A16_6/DMAX1(DEN,SMALL)
      A16_7 = A16_7/DMAX1(DEN,SMALL)
      A16_9 = A16_9/DMAX1(DEN,SMALL)
      A16_8 = A16_8/DMAX1(DEN,SMALL)
      A16_12 = A16_12/DMAX1(DEN,SMALL)
      A16_10 = A16_10/DMAX1(DEN,SMALL)
      A9_0 = A9_0 + A9_11*A11_0
      A9_16 = A9_16 + A9_11*A11_16
      DEN = 1 -A9_11*A11_9
      A9_0 = A9_0/DMAX1(DEN,SMALL)
      A9_16 = A9_16/DMAX1(DEN,SMALL)
      A9_7 = A9_7/DMAX1(DEN,SMALL)
      A16_0 = A16_0 + A16_10*A10_0
      A16_8 = A16_8 + A16_10*A10_8
      A16_12 = A16_12 + A16_10*A10_12
      DEN = 1 -A16_10*A10_16
      A16_0 = A16_0/DMAX1(DEN,SMALL)
      A16_6 = A16_6/DMAX1(DEN,SMALL)
      A16_7 = A16_7/DMAX1(DEN,SMALL)
      A16_9 = A16_9/DMAX1(DEN,SMALL)
      A16_8 = A16_8/DMAX1(DEN,SMALL)
      A16_12 = A16_12/DMAX1(DEN,SMALL)
      A8_0 = A8_0 + A8_10*A10_0
      A8_16 = A8_16 + A8_10*A10_16
      A8_12 = A8_10*A10_12
      DEN = 1 -A8_10*A10_8
      A8_0 = A8_0/DMAX1(DEN,SMALL)
      A8_16 = A8_16/DMAX1(DEN,SMALL)
      A8_7 = A8_7/DMAX1(DEN,SMALL)
      A8_12 = A8_12/DMAX1(DEN,SMALL)
      A12_0 = A12_0 + A12_10*A10_0
      A12_16 = A12_16 + A12_10*A10_16
      A12_8 = A12_10*A10_8
      DEN = 1 -A12_10*A10_12
      A12_0 = A12_0/DMAX1(DEN,SMALL)
      A12_16 = A12_16/DMAX1(DEN,SMALL)
      A12_8 = A12_8/DMAX1(DEN,SMALL)
      A16_0 = A16_0 + A16_12*A12_0
      A16_8 = A16_8 + A16_12*A12_8
      DEN = 1 -A16_12*A12_16
      A16_0 = A16_0/DMAX1(DEN,SMALL)
      A16_6 = A16_6/DMAX1(DEN,SMALL)
      A16_7 = A16_7/DMAX1(DEN,SMALL)
      A16_9 = A16_9/DMAX1(DEN,SMALL)
      A16_8 = A16_8/DMAX1(DEN,SMALL)
      A8_0 = A8_0 + A8_12*A12_0
      A8_16 = A8_16 + A8_12*A12_16
      DEN = 1 -A8_12*A12_8
      A8_0 = A8_0/DMAX1(DEN,SMALL)
      A8_16 = A8_16/DMAX1(DEN,SMALL)
      A8_7 = A8_7/DMAX1(DEN,SMALL)
      A2_0 = A2_0 + A2_1*A1_0
      A2_3 = A2_3 + A2_1*A1_3
      DEN = 1 -A2_1*A1_2
      A2_0 = A2_0/DMAX1(DEN,SMALL)
      A2_6 = A2_6/DMAX1(DEN,SMALL)
      A2_3 = A2_3/DMAX1(DEN,SMALL)
      A3_0 = A3_0 + A3_1*A1_0
      A3_2 = A3_2 + A3_1*A1_2
      DEN = 1 -A3_1*A1_3
      A3_0 = A3_0/DMAX1(DEN,SMALL)
      A3_6 = A3_6/DMAX1(DEN,SMALL)
      A3_2 = A3_2/DMAX1(DEN,SMALL)
      A16_0 = A16_0 + A16_8*A8_0
      A16_7 = A16_7 + A16_8*A8_7
      DEN = 1 -A16_8*A8_16
      A16_0 = A16_0/DMAX1(DEN,SMALL)
      A16_6 = A16_6/DMAX1(DEN,SMALL)
      A16_7 = A16_7/DMAX1(DEN,SMALL)
      A16_9 = A16_9/DMAX1(DEN,SMALL)
      A6_0 = A6_0 + A6_8*A8_0
      A6_16 = A6_16 + A6_8*A8_16
      A6_7 = A6_7 + A6_8*A8_7
      A7_0 = A7_0 + A7_8*A8_0
      A7_16 = A7_16 + A7_8*A8_16
      DEN = 1 -A7_8*A8_7
      A7_0 = A7_0/DMAX1(DEN,SMALL)
      A7_16 = A7_16/DMAX1(DEN,SMALL)
      A7_6 = A7_6/DMAX1(DEN,SMALL)
      A7_9 = A7_9/DMAX1(DEN,SMALL)
      A16_0 = A16_0 + A16_9*A9_0
      A16_7 = A16_7 + A16_9*A9_7
      DEN = 1 -A16_9*A9_16
      A16_0 = A16_0/DMAX1(DEN,SMALL)
      A16_6 = A16_6/DMAX1(DEN,SMALL)
      A16_7 = A16_7/DMAX1(DEN,SMALL)
      A7_0 = A7_0 + A7_9*A9_0
      A7_16 = A7_16 + A7_9*A9_16
      DEN = 1 -A7_9*A9_7
      A7_0 = A7_0/DMAX1(DEN,SMALL)
      A7_16 = A7_16/DMAX1(DEN,SMALL)
      A7_6 = A7_6/DMAX1(DEN,SMALL)
      A16_0 = A16_0 + A16_7*A7_0
      A16_6 = A16_6 + A16_7*A7_6
      DEN = 1 -A16_7*A7_16
      A16_0 = A16_0/DMAX1(DEN,SMALL)
      A16_6 = A16_6/DMAX1(DEN,SMALL)
      A6_0 = A6_0 + A6_7*A7_0
      A6_16 = A6_16 + A6_7*A7_16
      DEN = 1 -A6_7*A7_6
      A6_0 = A6_0/DMAX1(DEN,SMALL)
      A6_16 = A6_16/DMAX1(DEN,SMALL)
      A6_2 = A6_2/DMAX1(DEN,SMALL)
      A6_3 = A6_3/DMAX1(DEN,SMALL)
      A6_0 = A6_0 + A6_3*A3_0
      A6_2 = A6_2 + A6_3*A3_2
      DEN = 1 -A6_3*A3_6
      A6_0 = A6_0/DMAX1(DEN,SMALL)
      A6_16 = A6_16/DMAX1(DEN,SMALL)
      A6_2 = A6_2/DMAX1(DEN,SMALL)
      A2_0 = A2_0 + A2_3*A3_0
      A2_6 = A2_6 + A2_3*A3_6
      DEN = 1 -A2_3*A3_2
      A2_0 = A2_0/DMAX1(DEN,SMALL)
      A2_6 = A2_6/DMAX1(DEN,SMALL)
      A6_0 = A6_0 + A6_2*A2_0
      DEN = 1 -A6_2*A2_6
      A6_0 = A6_0/DMAX1(DEN,SMALL)
      A6_16 = A6_16/DMAX1(DEN,SMALL)
      A16_0 = A16_0 + A16_6*A6_0
      DEN = 1 -A16_6*A6_16
      A16_0 = A16_0/DMAX1(DEN,SMALL)
      XQ(16) = A16_0
      XQ(6) = A6_0 +A6_16*XQ(16)
      XQ(2) = A2_0 +A2_6*XQ(6)
      XQ(3) = A3_0 +A3_6*XQ(6) +A3_2*XQ(2)
      XQ(7) = A7_0 +A7_16*XQ(16) +A7_6*XQ(6)
      XQ(9) = A9_0 +A9_16*XQ(16) +A9_7*XQ(7)
      XQ(8) = A8_0 +A8_16*XQ(16) +A8_7*XQ(7)
      XQ(1) = A1_0 +A1_2*XQ(2) +A1_3*XQ(3)
      XQ(12) = A12_0 +A12_16*XQ(16) +A12_8*XQ(8)
      XQ(10) = A10_0 +A10_16*XQ(16) +A10_8*XQ(8) +A10_12*XQ(12)
      XQ(11) = A11_0 +A11_16*XQ(16) +A11_9*XQ(9)
      XQ(14) = A14_0 +A14_16*XQ(16) +A14_12*XQ(12)
      XQ(13) = A13_0 +A13_12*XQ(12) +A13_14*XQ(14)
      XQ(5) = A5_0 +A5_3*XQ(3)
      XQ(4) = A4_0 +A4_2*XQ(2)
      XQ(15) = A15_0 +A15_12*XQ(12)
      A17_0 = A17_0 + A17_18*A18_0
      DEN = 1 -A17_18*A18_17
      A17_0 = A17_0/DMAX1(DEN,SMALL)
      XQ(17) = A17_0
      XQ(18) = A18_0 +A18_17*XQ(17)
C
      RF( 31) = RF( 31)*XQ( 3)
      RF( 32) = RF( 32)*XQ( 3)
      RF( 33) = RF( 33)*XQ( 3)
      RF( 34) = RF( 34)*XQ( 3)
      RF( 35) = RF( 35)*XQ( 3)
      RF( 36) = RF( 36)*XQ( 3)
      RF( 37) = RF( 37)*XQ( 3)
      RF( 39) = RF( 39)*XQ( 3)
      RF( 40) = RF( 40)*XQ( 1)
      RF( 41) = RF( 41)*XQ( 1)
      RB( 41) = RB( 41)*XQ( 3)
      RF( 42) = RF( 42)*XQ( 1)
      RF( 43) = RF( 43)*XQ( 1)
      RF( 44) = RF( 44)*XQ( 1)
      RB( 44) = RB( 44)*XQ( 3)
      RF( 45) = RF( 45)*XQ( 1)
      RF( 46) = RF( 46)*XQ( 1)
      RF( 48) = RF( 48)*XQ( 2)
      RB( 48) = RB( 48)*XQ( 1)
      RF( 49) = RF( 49)*XQ( 2)
      RF( 50) = RF( 50)*XQ( 2)
      RB( 50) = RB( 50)*XQ( 3)
      RF( 51) = RF( 51)*XQ( 2)
      RF( 52) = RF( 52)*XQ( 2)
      RF( 53) = RF( 53)*XQ( 2)
      RF( 54) = RF( 54)*XQ( 2)
      RF( 55) = RF( 55)*XQ( 2)
      RB( 55) = RB( 55)*XQ( 1)
      RF( 56) = RF( 56)*XQ( 2)
      RB( 56) = RB( 56)*XQ( 1)
      RF( 57) = RF( 57)*XQ( 2)
      RB( 57) = RB( 57)*XQ( 1)
      RF( 58) = RF( 58)*XQ( 2)
      RB( 59) = RB( 59)*XQ( 4)
      RB( 60) = RB( 60)*XQ( 3)
      RB( 61) = RB( 61)*XQ( 3)
      RB( 62) = RB( 62)*XQ( 3)
      RB( 63) = RB( 63)*XQ( 3)
      RB( 64) = RB( 64)*XQ( 3)
      RB( 67) = RB( 67)*XQ( 1)
      RB( 68) = RB( 68)*XQ( 2)
      RB( 69) = RB( 69)*XQ( 4)
      RB( 72) = RB( 72)*XQ( 4)
      RF( 74) = RF( 74)*XQ( 3)
      RB( 75) = RB( 75)*XQ( 3)
      RF( 76) = RF( 76)*XQ( 1)
      RF( 77) = RF( 77)*XQ( 2)
      RB( 79) = RB( 79)*XQ( 6)
      RF( 80) = RF( 80)*XQ( 4)
      RF( 81) = RF( 81)*XQ( 4)
      RF( 82) = RF( 82)*XQ( 4)
      RB( 82) = RB( 82)*XQ( 2)
      RF( 83) = RF( 83)*XQ( 4)
      RF( 84) = RF( 84)*XQ( 4)
      RF( 85) = RF( 85)*XQ( 4)
      RF( 89) = RF( 89)*XQ( 1)
      RF( 90) = RF( 90)*XQ( 2)
      RF( 91) = RF( 91)*XQ( 5)
      RB( 92) = RB( 92)*XQ( 1)
      RF( 94) = RF( 94)*XQ( 3)
      RB( 94) = RB( 94)*XQ( 5)
      RF( 96) = RF( 96)*XQ( 5)
      RF( 97) = RF( 97)*XQ( 5)
      RF( 98) = RF( 98)*XQ( 5)
      RF( 99) = RF( 99)*XQ( 5)
      RF(100) = RF(100)*XQ( 5)
      RF(101) = RF(101)*XQ( 5)
      RF(102) = RF(102)*XQ( 5)
      RB(102) = RB(102)*XQ( 3)
      RF(103) = RF(103)*XQ( 5)
      RF(104) = RF(104)*XQ( 5)
      RF(107) = RF(107)*XQ( 5)
      RF(108) = RF(108)*XQ( 5)
      RF(109) = RF(109)*XQ( 5)
      RB(112) = RB(112)*XQ( 3)
      RB(114) = RB(114)*XQ( 6)
      RB(115) = RB(115)*XQ( 5)
      RB(116) = RB(116)*XQ( 5)
      RB(117) = RB(117)*XQ( 3)
      RB(118) = RB(118)*XQ( 1)
      RB(119) = RB(119)*XQ( 5)
      RF(120) = RF(120)*XQ( 3)
      RB(120) = RB(120)*XQ( 6)
      RF(121) = RF(121)*XQ( 1)
      RF(122) = RF(122)*XQ( 2)
      RB(123) = RB(123)*XQ( 5)
      RF(124) = RF(124)*XQ( 7)
      RB(125) = RB(125)*XQ( 5)
      RF(126) = RF(126)*XQ( 5)
      RF(127) = RF(127)*XQ( 6)
      RF(128) = RF(128)*XQ( 6)
      RF(129) = RF(129)*XQ( 6)
      RF(130) = RF(130)*XQ( 6)
      RF(131) = RF(131)*XQ( 6)
      RF(132) = RF(132)*XQ( 6)
      RF(133) = RF(133)*XQ( 6)
      RF(134) = RF(134)*XQ( 6)
      RB(137) = RB(137)*XQ( 6)
      RB(138) = RB(138)*XQ( 6)
      RB(139) = RB(139)*XQ( 6)
      RF(140) = RF(140)*XQ( 2)
      RB(140) = RB(140)*XQ( 6)
      RB(141) = RB(141)*XQ( 6)
      RB(147) = RB(147)*XQ( 5)
      RF(148) = RF(148)*XQ( 3)
      RB(150) = RB(150)*XQ( 7)
      RB(159) = RB(159)*XQ( 3)
      RB(160) = RB(160)*XQ( 5)
      RB(161) = RB(161)*XQ( 5)
      RF(162) = RF(162)*XQ( 7)
      RB(162) = RB(162)*XQ( 6)
      RF(163) = RF(163)*XQ( 7)
      RF(164) = RF(164)*XQ( 7)
      RB(164) = RB(164)*XQ( 6)
      RF(165) = RF(165)*XQ( 7)
      RF(166) = RF(166)*XQ( 7)
      RF(167) = RF(167)*XQ( 7)
      RB(167) = RB(167)*XQ( 6)
      RF(168) = RF(168)*XQ( 7)
      RF(172) = RF(172)*XQ( 3)
      RB(173) = RB(173)*XQ( 8)
      RB(174) = RB(174)*XQ( 6)
      RF(184) = RF(184)*XQ( 8)
      RF(185) = RF(185)*XQ( 8)
      RF(186) = RF(186)*XQ( 8)
      RB(186) = RB(186)*XQ( 7)
      RF(187) = RF(187)*XQ( 8)
      RF(188) = RF(188)*XQ( 8)
      RF(189) = RF(189)*XQ( 8)
      RB(189) = RB(189)*XQ( 7)
      RF(190) = RF(190)*XQ( 8)
      RB(192) = RB(192)*XQ( 5)
      RB(193) = RB(193)*XQ( 9)
      RB(194) = RB(194)*XQ( 7)
      RB(195) = RB(195)*XQ( 6)
      RF(196) = RF(196)*XQ( 7)
      RB(196) = RB(196)*XQ( 9)
      RB(197) = RB(197)*XQ(10)
      RB(198) = RB(198)*XQ( 8)
      RB(199) = RB(199)*XQ( 7)
      RF(200) = RF(200)*XQ( 8)
      RB(200) = RB(200)*XQ(10)
      RB(201) = RB(201)*XQ(11)
      RB(202) = RB(202)*XQ( 9)
      RB(203) = RB(203)*XQ( 8)
      RF(204) = RF(204)*XQ( 9)
      RB(204) = RB(204)*XQ(11)
      RB(205) = RB(205)*XQ(12)
      RB(206) = RB(206)*XQ(10)
      RB(207) = RB(207)*XQ( 9)
      RF(208) = RF(208)*XQ(10)
      RB(208) = RB(208)*XQ(12)
      RB(210) = RB(210)*XQ(11)
      RB(211) = RB(211)*XQ(10)
      RF(212) = RF(212)*XQ(11)
      RB(213) = RB(213)*XQ(13)
      RB(214) = RB(214)*XQ(12)
      RB(215) = RB(215)*XQ(11)
      RF(216) = RF(216)*XQ(12)
      RB(216) = RB(216)*XQ(13)
      RB(217) = RB(217)*XQ(11)
      RF(218) = RF(218)*XQ(13)
      RB(218) = RB(218)*XQ(14)
      RF(219) = RF(219)*XQ(14)
      RB(219) = RB(219)*XQ(16)
      RB(220) = RB(220)*XQ(15)
      RF(221) = RF(221)*XQ(12)
      RB(221) = RB(221)*XQ(15)
      RF(222) = RF(222)*XQ(11)
      RB(222) = RB(222)*XQ(16)
      RF(223) = RF(223)*XQ( 6)
      RB(223) = RB(223)*XQ(16)
      RF(224) = RF(224)*XQ(10)
      RB(224) = RB(224)*XQ(16)
      RF(225) = RF(225)*XQ( 7)
      RB(225) = RB(225)*XQ(16)
      RF(226) = RF(226)*XQ( 9)
      RB(226) = RB(226)*XQ(16)
      RF(227) = RF(227)*XQ( 8)
      RB(227) = RB(227)*XQ(16)
      RF(229) = RF(229)*XQ( 7)
      RB(233) = RB(233)*XQ(14)
      RB(234) = RB(234)*XQ(15)
      RB(235) = RB(235)*XQ(16)
      RB(236) = RB(236)*XQ(14)
      RB(237) = RB(237)*XQ(15)
      RB(238) = RB(238)*XQ(16)
      RB(239) = RB(239)*XQ(14)
      RB(240) = RB(240)*XQ(15)
      RB(241) = RB(241)*XQ(16)
      RB(242) = RB(242)*XQ(14)
      RB(243) = RB(243)*XQ(15)
      RB(244) = RB(244)*XQ(16)
      RB(245) = RB(245)*XQ(14)
      RB(246) = RB(246)*XQ(15)
      RB(247) = RB(247)*XQ(16)
      RB(248) = RB(248)*XQ(14)
      RB(249) = RB(249)*XQ(15)
      RB(250) = RB(250)*XQ(16)
      RF(251) = RF(251)*XQ(14)
      RB(252) = RB(252)*XQ(14)
      RF(253) = RF(253)*XQ(15)
      RB(254) = RB(254)*XQ(15)
      RF(255) = RF(255)*XQ(16)
      RB(256) = RB(256)*XQ(16)
      RB(257) = RB(257)*XQ(17)
      RF(258) = RF(258)*XQ(17)
      RF(259) = RF(259)*XQ(14)
      RB(260) = RB(260)*XQ(14)
      RF(261) = RF(261)*XQ(15)
      RB(262) = RB(262)*XQ(15)
      RF(263) = RF(263)*XQ(16)
      RB(264) = RB(264)*XQ(16)
      RF(265) = RF(265)*XQ(17)
      RB(265) = RB(265)*XQ(18)
      RF(266) = RF(266)*XQ(18)
      RB(266) = RB(266)*XQ(17)
      RF(267) = RF(267)*XQ(18)
      RB(268) = RB(268)*XQ( 6)
C
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE STIF(RF, RB, DIFF, DT, C)
        IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
      DIMENSION RF(*), RB(*), DIFF(*), C(*)
      TC = 1./DT
C
C     NC12H26
      DDOT=+RB(228)+RB(229)+RB(230)+RB(231)+RB(232)+RF(233)+RF(234)
     *+RF(235)+RF(236)+RF(237)+RF(238)+RF(239)+RF(240)+RF(241)+RF(242)
     *+RF(243)+RF(244)+RF(245)+RF(246)+RF(247)+RF(248)+RF(249)+RF(250)
      TINV=DDOT/C(1)
      IF (TINV .GT. TC) THEN
      CDOT=+RF(228)+RF(229)+RF(230)+RF(231)+RF(232)+RB(233)+RB(234)
     *+RB(235)+RB(236)+RB(237)+RB(238)+RB(239)+RB(240)+RB(241)+RB(242)
     *+RB(243)+RB(244)+RB(245)+RB(246)+RB(247)+RB(248)+RB(249)+RB(250)
      C0=C(1)*(CDOT+DIFF(1)/170.34102368)/DDOT
      C0=C(1)*(CDOT+DIFF(1)/170.34102368+(C(1)-C0)/DT)/DDOT
      R=C0/C(1)
      RB(228)=RB(228)*R
      RB(229)=RB(229)*R
      RB(230)=RB(230)*R
      RB(231)=RB(231)*R
      RB(232)=RB(232)*R
      RF(233)=RF(233)*R
      RF(234)=RF(234)*R
      RF(235)=RF(235)*R
      RF(236)=RF(236)*R
      RF(237)=RF(237)*R
      RF(238)=RF(238)*R
      RF(239)=RF(239)*R
      RF(240)=RF(240)*R
      RF(241)=RF(241)*R
      RF(242)=RF(242)*R
      RF(243)=RF(243)*R
      RF(244)=RF(244)*R
      RF(245)=RF(245)*R
      RF(246)=RF(246)*R
      RF(247)=RF(247)*R
      RF(248)=RF(248)*R
      RF(249)=RF(249)*R
      RF(250)=RF(250)*R
      ENDIF
C
C     H
      DDOT=+RF(1)+RB(2)+RB(3)+RF(5)+RF(6)+RB(7)+RF(9)+RF(13)+RF(14)
     *+2*RF(19)+RF(20)+2*RF(22)+2*RF(23)+2*RF(24)+RF(25)+RB(26)+RB(27)
     *+RB(31)+RF(32)+RB(34)+RB(37)+RF(39)+RF(40)+RB(41)+RB(42)+RB(43)
     *+2*RB(45)+RB(50)+RB(51)+RB(52)+RB(53)+RF(59)+RF(60)+RF(65)+RB(66)
     *+RB(76)+RB(77)+RB(79)+RF(80)+RF(81)+RF(82)+RF(86)+RB(91)+RF(96)
     *+RF(97)+RB(109)+RF(112)+RF(114)+RF(115)+RB(121)+RB(122)+RF(127)
     *+RF(128)+RF(137)+RF(142)+RB(143)+2*RB(144)+RF(150)+RF(151)+RF(152)
     *+2*RB(153)+RF(159)+RF(162)+RF(163)+RF(169)+RF(170)+RF(173)+RF(174)
     *+RF(175)+RF(176)+RF(184)+RF(185)+RF(193)+RF(194)+RF(195)+RF(197)
     *+RF(198)+RF(199)+RF(201)+RF(202)+RF(203)+RF(205)+RF(206)+RF(207)
     *+RF(209)+RF(210)+RF(211)+RF(213)+RF(214)+RF(215)+RF(233)+RF(234)
     *+RF(235)
      TINV=DDOT/C(2)
      IF (TINV .GT. TC) THEN
      CDOT=+RB(1)+RF(2)+RF(3)+RB(5)+RB(6)+RF(7)+RB(9)+RB(13)+RB(14)
     *+2*RB(19)+RB(20)+2*RB(22)+2*RB(23)+2*RB(24)+RB(25)+RF(26)+RF(27)
     *+RF(31)+RB(32)+RF(34)+RF(37)+RB(39)+RB(40)+RF(41)+RF(42)+RF(43)
     *+2*RF(45)+RF(50)+RF(51)+RF(52)+RF(53)+RB(59)+RB(60)+RB(65)+RF(66)
     *+RF(76)+RF(77)+RF(79)+RB(80)+RB(81)+RB(82)+RB(86)+RF(91)+RB(96)
     *+RB(97)+RF(109)+RB(112)+RB(114)+RB(115)+RF(121)+RF(122)+RB(127)
     *+RB(128)+RB(137)+RB(142)+RF(143)+2*RF(144)+RB(150)+RB(151)+RB(152)
     *+2*RF(153)+RB(159)+RB(162)+RB(163)+RB(169)+RB(170)+RB(173)+RB(174)
     *+RB(175)+RB(176)+RB(184)+RB(185)+RB(193)+RB(194)+RB(195)+RB(197)
     *+RB(198)+RB(199)+RB(201)+RB(202)+RB(203)+RB(205)+RB(206)+RB(207)
     *+RB(209)+RB(210)+RB(211)+RB(213)+RB(214)+RB(215)+RB(233)+RB(234)
     *+RB(235)
      C0=C(2)*(CDOT+DIFF(2)/1.00796998)/DDOT
      C0=C(2)*(CDOT+DIFF(2)/1.00796998+(C(2)-C0)/DT)/DDOT
      R=C0/C(2)
      RF(1)=RF(1)*R
      RB(2)=RB(2)*R
      RB(3)=RB(3)*R
      RF(5)=RF(5)*R
      RF(6)=RF(6)*R
      RB(7)=RB(7)*R
      RF(9)=RF(9)*R
      RF(13)=RF(13)*R
      RF(14)=RF(14)*R
      RF(19)=RF(19)*R
      RF(19)=RF(19)*R
      RF(20)=RF(20)*R
      RF(22)=RF(22)*R
      RF(22)=RF(22)*R
      RF(23)=RF(23)*R
      RF(23)=RF(23)*R
      RF(24)=RF(24)*R
      RF(24)=RF(24)*R
      RF(25)=RF(25)*R
      RB(26)=RB(26)*R
      RB(27)=RB(27)*R
      RB(31)=RB(31)*R
      RF(32)=RF(32)*R
      RB(34)=RB(34)*R
      RB(37)=RB(37)*R
      RF(39)=RF(39)*R
      RF(40)=RF(40)*R
      RB(41)=RB(41)*R
      RB(42)=RB(42)*R
      RB(43)=RB(43)*R
      RB(45)=RB(45)*R
      RB(45)=RB(45)*R
      RB(50)=RB(50)*R
      RB(51)=RB(51)*R
      RB(52)=RB(52)*R
      RB(53)=RB(53)*R
      RF(59)=RF(59)*R
      RF(60)=RF(60)*R
      RF(65)=RF(65)*R
      RB(66)=RB(66)*R
      RB(76)=RB(76)*R
      RB(77)=RB(77)*R
      RB(79)=RB(79)*R
      RF(80)=RF(80)*R
      RF(81)=RF(81)*R
      RF(82)=RF(82)*R
      RF(86)=RF(86)*R
      RB(91)=RB(91)*R
      RF(96)=RF(96)*R
      RF(97)=RF(97)*R
      RB(109)=RB(109)*R
      RF(112)=RF(112)*R
      RF(114)=RF(114)*R
      RF(115)=RF(115)*R
      RB(121)=RB(121)*R
      RB(122)=RB(122)*R
      RF(127)=RF(127)*R
      RF(128)=RF(128)*R
      RF(137)=RF(137)*R
      RF(142)=RF(142)*R
      RB(143)=RB(143)*R
      RB(144)=RB(144)*R
      RB(144)=RB(144)*R
      RF(150)=RF(150)*R
      RF(151)=RF(151)*R
      RF(152)=RF(152)*R
      RB(153)=RB(153)*R
      RB(153)=RB(153)*R
      RF(159)=RF(159)*R
      RF(162)=RF(162)*R
      RF(163)=RF(163)*R
      RF(169)=RF(169)*R
      RF(170)=RF(170)*R
      RF(173)=RF(173)*R
      RF(174)=RF(174)*R
      RF(175)=RF(175)*R
      RF(176)=RF(176)*R
      RF(184)=RF(184)*R
      RF(185)=RF(185)*R
      RF(193)=RF(193)*R
      RF(194)=RF(194)*R
      RF(195)=RF(195)*R
      RF(197)=RF(197)*R
      RF(198)=RF(198)*R
      RF(199)=RF(199)*R
      RF(201)=RF(201)*R
      RF(202)=RF(202)*R
      RF(203)=RF(203)*R
      RF(205)=RF(205)*R
      RF(206)=RF(206)*R
      RF(207)=RF(207)*R
      RF(209)=RF(209)*R
      RF(210)=RF(210)*R
      RF(211)=RF(211)*R
      RF(213)=RF(213)*R
      RF(214)=RF(214)*R
      RF(215)=RF(215)*R
      RF(233)=RF(233)*R
      RF(234)=RF(234)*R
      RF(235)=RF(235)*R
      ENDIF
C
C     O
      DDOT=+RB(1)+RF(2)+RB(4)+RB(9)+RF(10)+RF(15)+2*RF(21)+RF(25)+RF(29)
     *+RB(30)+RF(33)+RF(34)+RF(41)+RF(49)+RF(50)+RF(61)+RF(66)+RB(69)
     *+RF(83)+RF(87)+RF(92)+RF(98)+RB(101)+RF(116)+RF(117)+RF(118)
     *+RF(129)+RF(138)+RF(143)+RF(153)+RF(154)+RF(155)+RF(160)+RF(164)
     *+RF(177)+RF(178)+RF(179)+RF(186)+RF(236)+RF(237)+RF(238)
      TINV=DDOT/C(3)
      IF (TINV .GT. TC) THEN
      CDOT=+RF(1)+RB(2)+RF(4)+RF(9)+RB(10)+RB(15)+2*RB(21)+RB(25)+RB(29)
     *+RF(30)+RB(33)+RB(34)+RB(41)+RB(49)+RB(50)+RB(61)+RB(66)+RF(69)
     *+RB(83)+RB(87)+RB(92)+RB(98)+RF(101)+RB(116)+RB(117)+RB(118)
     *+RB(129)+RB(138)+RB(143)+RB(153)+RB(154)+RB(155)+RB(160)+RB(164)
     *+RB(177)+RB(178)+RB(179)+RB(186)+RB(236)+RB(237)+RB(238)
      C0=C(3)*(CDOT+DIFF(3)/15.99940014)/DDOT
      C0=C(3)*(CDOT+DIFF(3)/15.99940014+(C(3)-C0)/DT)/DDOT
      R=C0/C(3)
      RB(1)=RB(1)*R
      RF(2)=RF(2)*R
      RB(4)=RB(4)*R
      RB(9)=RB(9)*R
      RF(10)=RF(10)*R
      RF(15)=RF(15)*R
      RF(21)=RF(21)*R
      RF(21)=RF(21)*R
      RF(25)=RF(25)*R
      RF(29)=RF(29)*R
      RB(30)=RB(30)*R
      RF(33)=RF(33)*R
      RF(34)=RF(34)*R
      RF(41)=RF(41)*R
      RF(49)=RF(49)*R
      RF(50)=RF(50)*R
      RF(61)=RF(61)*R
      RF(66)=RF(66)*R
      RB(69)=RB(69)*R
      RF(83)=RF(83)*R
      RF(87)=RF(87)*R
      RF(92)=RF(92)*R
      RF(98)=RF(98)*R
      RB(101)=RB(101)*R
      RF(116)=RF(116)*R
      RF(117)=RF(117)*R
      RF(118)=RF(118)*R
      RF(129)=RF(129)*R
      RF(138)=RF(138)*R
      RF(143)=RF(143)*R
      RF(153)=RF(153)*R
      RF(154)=RF(154)*R
      RF(155)=RF(155)*R
      RF(160)=RF(160)*R
      RF(164)=RF(164)*R
      RF(177)=RF(177)*R
      RF(178)=RF(178)*R
      RF(179)=RF(179)*R
      RF(186)=RF(186)*R
      RF(236)=RF(236)*R
      RF(237)=RF(237)*R
      RF(238)=RF(238)*R
      ENDIF
C
C     OH
      DDOT=+RB(1)+RB(2)+RF(3)+2*RF(4)+2*RB(6)+RF(8)+RB(10)+RB(13)+RB(15)
     *+RF(16)+RF(17)+2*RF(18)+RF(20)+RB(25)+RF(26)+RF(27)+RB(28)+RB(33)
     *+RF(35)+RF(42)+RB(44)+RB(46)+RF(51)+RB(53)+RB(61)+RF(62)+RF(67)
     *+RF(68)+RB(70)+RB(72)+RB(81)+RB(83)+RF(84)+RB(87)+RF(88)+RF(93)
     *+RF(99)+RB(103)+RB(113)+RB(116)+RF(119)+RB(133)+RB(138)+RF(139)
     *+RF(144)+RB(145)+RB(147)+RB(155)+RF(156)+RB(160)+RF(161)+RF(165)
     *+RB(167)+RB(171)+RB(178)+RB(179)+RF(180)+RF(187)+RB(189)+RB(236)
     *+RB(237)+RB(238)+RF(239)+RF(240)+RF(241)+RB(267)
      TINV=DDOT/C(4)
      IF (TINV .GT. TC) THEN
      CDOT=+RF(1)+RF(2)+RB(3)+2*RB(4)+2*RF(6)+RB(8)+RF(10)+RF(13)+RF(15)
     *+RB(16)+RB(17)+2*RB(18)+RB(20)+RF(25)+RB(26)+RB(27)+RF(28)+RF(33)
     *+RB(35)+RB(42)+RF(44)+RF(46)+RB(51)+RF(53)+RF(61)+RB(62)+RB(67)
     *+RB(68)+RF(70)+RF(72)+RF(81)+RF(83)+RB(84)+RF(87)+RB(88)+RB(93)
     *+RB(99)+RF(103)+RF(113)+RF(116)+RB(119)+RF(133)+RF(138)+RB(139)
     *+RB(144)+RF(145)+RF(147)+RF(155)+RB(156)+RF(160)+RB(161)+RB(165)
     *+RF(167)+RF(171)+RF(178)+RF(179)+RB(180)+RB(187)+RF(189)+RF(236)
     *+RF(237)+RF(238)+RB(239)+RB(240)+RB(241)+RF(267)+RF(268)
      C0=C(4)*(CDOT+DIFF(4)/17.00737011)/DDOT
      C0=C(4)*(CDOT+DIFF(4)/17.00737011+(C(4)-C0)/DT)/DDOT
      R=C0/C(4)
      RB(1)=RB(1)*R
      RB(2)=RB(2)*R
      RF(3)=RF(3)*R
      RF(4)=RF(4)*R
      RF(4)=RF(4)*R
      RB(6)=RB(6)*R
      RB(6)=RB(6)*R
      RF(8)=RF(8)*R
      RB(10)=RB(10)*R
      RB(13)=RB(13)*R
      RB(15)=RB(15)*R
      RF(16)=RF(16)*R
      RF(17)=RF(17)*R
      RF(18)=RF(18)*R
      RF(18)=RF(18)*R
      RF(20)=RF(20)*R
      RB(25)=RB(25)*R
      RF(26)=RF(26)*R
      RF(27)=RF(27)*R
      RB(28)=RB(28)*R
      RB(33)=RB(33)*R
      RF(35)=RF(35)*R
      RF(42)=RF(42)*R
      RB(44)=RB(44)*R
      RB(46)=RB(46)*R
      RF(51)=RF(51)*R
      RB(53)=RB(53)*R
      RB(61)=RB(61)*R
      RF(62)=RF(62)*R
      RF(67)=RF(67)*R
      RF(68)=RF(68)*R
      RB(70)=RB(70)*R
      RB(72)=RB(72)*R
      RB(81)=RB(81)*R
      RB(83)=RB(83)*R
      RF(84)=RF(84)*R
      RB(87)=RB(87)*R
      RF(88)=RF(88)*R
      RF(93)=RF(93)*R
      RF(99)=RF(99)*R
      RB(103)=RB(103)*R
      RB(113)=RB(113)*R
      RB(116)=RB(116)*R
      RF(119)=RF(119)*R
      RB(133)=RB(133)*R
      RB(138)=RB(138)*R
      RF(139)=RF(139)*R
      RF(144)=RF(144)*R
      RB(145)=RB(145)*R
      RB(147)=RB(147)*R
      RB(155)=RB(155)*R
      RF(156)=RF(156)*R
      RB(160)=RB(160)*R
      RF(161)=RF(161)*R
      RF(165)=RF(165)*R
      RB(167)=RB(167)*R
      RB(171)=RB(171)*R
      RB(178)=RB(178)*R
      RB(179)=RB(179)*R
      RF(180)=RF(180)*R
      RF(187)=RF(187)*R
      RB(189)=RB(189)*R
      RB(236)=RB(236)*R
      RB(237)=RB(237)*R
      RB(238)=RB(238)*R
      RF(239)=RF(239)*R
      RF(240)=RF(240)*R
      RF(241)=RF(241)*R
      RB(267)=RB(267)*R
      ENDIF
C
C     HO2
      DDOT=+RB(5)+RF(6)+RB(7)+RF(8)+RF(9)+RF(10)+2*RF(11)+2*RF(12)
     *+RB(14)+RB(15)+RB(16)+RB(17)+RF(28)+RB(36)+RF(46)+RB(63)+RF(64)
     *+RF(71)+RF(72)+RB(73)+RB(85)+RB(100)+RF(103)+RB(104)+RB(125)
     *+RB(130)+RF(131)+RF(132)+RF(133)+RB(134)+RF(146)+RF(147)+RF(157)
     *+RB(166)+RF(167)+RF(171)+RB(181)+RF(182)+RB(188)+RF(189)+RB(242)
     *+RB(243)+RB(244)+RF(245)+RF(246)+RF(247)+RF(260)+RF(262)+RF(264)
      TINV=DDOT/C(5)
      IF (TINV .GT. TC) THEN
      CDOT=+RF(5)+RB(6)+RF(7)+RB(8)+RB(9)+RB(10)+2*RB(11)+2*RB(12)
     *+RF(14)+RF(15)+RF(16)+RF(17)+RB(28)+RF(36)+RB(46)+RF(63)+RB(64)
     *+RB(71)+RB(72)+RF(73)+RF(85)+RF(100)+RB(103)+RF(104)+RF(125)
     *+RF(130)+RB(131)+RB(132)+RB(133)+RF(134)+RB(146)+RB(147)+RB(157)
     *+RF(166)+RB(167)+RB(171)+RF(181)+RB(182)+RF(188)+RB(189)+RF(242)
     *+RF(243)+RF(244)+RB(245)+RB(246)+RB(247)+RF(259)+RF(261)+RF(263)
      C0=C(5)*(CDOT+DIFF(5)/33.00677025)/DDOT
      C0=C(5)*(CDOT+DIFF(5)/33.00677025+(C(5)-C0)/DT)/DDOT
      R=C0/C(5)
      RB(5)=RB(5)*R
      RF(6)=RF(6)*R
      RB(7)=RB(7)*R
      RF(8)=RF(8)*R
      RF(9)=RF(9)*R
      RF(10)=RF(10)*R
      RF(11)=RF(11)*R
      RF(11)=RF(11)*R
      RF(12)=RF(12)*R
      RF(12)=RF(12)*R
      RB(14)=RB(14)*R
      RB(15)=RB(15)*R
      RB(16)=RB(16)*R
      RB(17)=RB(17)*R
      RF(28)=RF(28)*R
      RB(36)=RB(36)*R
      RF(46)=RF(46)*R
      RB(63)=RB(63)*R
      RF(64)=RF(64)*R
      RF(71)=RF(71)*R
      RF(72)=RF(72)*R
      RB(73)=RB(73)*R
      RB(85)=RB(85)*R
      RB(100)=RB(100)*R
      RF(103)=RF(103)*R
      RB(104)=RB(104)*R
      RB(125)=RB(125)*R
      RB(130)=RB(130)*R
      RF(131)=RF(131)*R
      RF(132)=RF(132)*R
      RF(133)=RF(133)*R
      RB(134)=RB(134)*R
      RF(146)=RF(146)*R
      RF(147)=RF(147)*R
      RF(157)=RF(157)*R
      RB(166)=RB(166)*R
      RF(167)=RF(167)*R
      RF(171)=RF(171)*R
      RB(181)=RB(181)*R
      RF(182)=RF(182)*R
      RB(188)=RB(188)*R
      RF(189)=RF(189)*R
      RB(242)=RB(242)*R
      RB(243)=RB(243)*R
      RB(244)=RB(244)*R
      RF(245)=RF(245)*R
      RF(246)=RF(246)*R
      RF(247)=RF(247)*R
      RF(260)=RF(260)*R
      RF(262)=RF(262)*R
      RF(264)=RF(264)*R
      ENDIF
C
C     H2
      DDOT=+RF(2)+RF(3)+RF(7)+RB(14)+RB(19)+RF(22)+2*RB(22)+RB(23)
     *+RB(24)+RB(32)+RF(38)+RF(43)+RB(47)+RB(49)+RF(52)+RB(60)+RB(80)
     *+RB(86)+RB(97)+RB(115)+RB(128)+RB(137)+RB(152)+RB(163)+RB(176)
     *+RB(185)+RB(233)+RB(234)+RB(235)
      TINV=DDOT/C(6)
      IF (TINV .GT. TC) THEN
      CDOT=+RB(2)+RB(3)+RB(7)+RF(14)+RF(19)+RB(22)+2*RF(22)+RF(23)
     *+RF(24)+RF(32)+RB(38)+RB(43)+RF(47)+RF(49)+RB(52)+RF(60)+RF(80)
     *+RF(86)+RF(97)+RF(115)+RF(128)+RF(137)+RF(152)+RF(163)+RF(176)
     *+RF(185)+RF(233)+RF(234)+RF(235)
      C0=C(6)*(CDOT+DIFF(6)/2.01593995)/DDOT
      C0=C(6)*(CDOT+DIFF(6)/2.01593995+(C(6)-C0)/DT)/DDOT
      R=C0/C(6)
      RF(2)=RF(2)*R
      RF(3)=RF(3)*R
      RF(7)=RF(7)*R
      RB(14)=RB(14)*R
      RB(19)=RB(19)*R
      RF(22)=RF(22)*R
      RB(22)=RB(22)*R
      RB(22)=RB(22)*R
      RB(23)=RB(23)*R
      RB(24)=RB(24)*R
      RB(32)=RB(32)*R
      RF(38)=RF(38)*R
      RF(43)=RF(43)*R
      RB(47)=RB(47)*R
      RB(49)=RB(49)*R
      RF(52)=RF(52)*R
      RB(60)=RB(60)*R
      RB(80)=RB(80)*R
      RB(86)=RB(86)*R
      RB(97)=RB(97)*R
      RB(115)=RB(115)*R
      RB(128)=RB(128)*R
      RB(137)=RB(137)*R
      RB(152)=RB(152)*R
      RB(163)=RB(163)*R
      RB(176)=RB(176)*R
      RB(185)=RB(185)*R
      RB(233)=RB(233)*R
      RB(234)=RB(234)*R
      RB(235)=RB(235)*R
      ENDIF
C
C     H2O2
      DDOT=+RB(11)+RB(12)+RF(13)+RF(14)+RF(15)+RF(16)+RF(17)+RB(18)
     *+RB(64)+RF(73)+RF(104)+RB(132)+RF(134)+RB(157)+RB(182)+RB(245)
     *+RB(246)+RB(247)
      TINV=DDOT/C(8)
      IF (TINV .GT. TC) THEN
      CDOT=+RF(11)+RF(12)+RB(13)+RB(14)+RB(15)+RB(16)+RB(17)+RF(18)
     *+RF(64)+RB(73)+RB(104)+RF(132)+RB(134)+RF(157)+RF(182)+RF(245)
     *+RF(246)+RF(247)
      C0=C(8)*(CDOT+DIFF(8)/34.01474023)/DDOT
      C0=C(8)*(CDOT+DIFF(8)/34.01474023+(C(8)-C0)/DT)/DDOT
      R=C0/C(8)
      RB(11)=RB(11)*R
      RB(12)=RB(12)*R
      RF(13)=RF(13)*R
      RF(14)=RF(14)*R
      RF(15)=RF(15)*R
      RF(16)=RF(16)*R
      RF(17)=RF(17)*R
      RB(18)=RB(18)*R
      RB(64)=RB(64)*R
      RF(73)=RF(73)*R
      RF(104)=RF(104)*R
      RB(132)=RB(132)*R
      RF(134)=RF(134)*R
      RB(157)=RB(157)*R
      RB(182)=RB(182)*R
      RB(245)=RB(245)*R
      RB(246)=RB(246)*R
      RB(247)=RB(247)*R
      ENDIF
C
C     CH3
      DDOT=+RB(40)+RB(43)+RB(52)+RF(65)+RF(66)+RF(67)+RF(68)+RF(69)
     *+RF(70)+RF(71)+RF(72)+RF(73)+RF(74)+RF(75)+RF(76)+RF(77)+2*RF(78)
     *+2*RF(79)+RB(81)+RB(86)+RB(87)+RB(88)+2*RB(89)+2*RB(90)+RB(93)
     *+RF(95)+RB(98)+RF(107)+RF(108)+RF(109)+RB(111)+RB(112)+RB(117)
     *+RF(123)+RB(124)+RB(129)+RB(133)+RB(136)+RB(140)+RF(141)+RF(149)
     *+RB(151)+RF(158)+RB(162)+RF(168)+RB(170)+RB(175)+RF(183)+RF(190)
     *+RF(248)+RF(249)+RF(250)
      TINV=DDOT/C(10)
      IF (TINV .GT. TC) THEN
      CDOT=+RF(40)+RF(43)+RF(52)+RB(65)+RB(66)+RB(67)+RB(68)+RB(69)
     *+RB(70)+RB(71)+RB(72)+RB(73)+RB(74)+RB(75)+RB(76)+RB(77)+2*RB(78)
     *+2*RB(79)+RF(81)+RF(86)+RF(87)+RF(88)+2*RF(89)+2*RF(90)+RF(93)
     *+RB(95)+RF(98)+RB(107)+RB(108)+RB(109)+RF(111)+RF(112)+RF(117)
     *+RB(123)+RF(124)+RF(129)+RF(133)+RF(136)+RF(140)+RB(141)+RB(149)
     *+RF(151)+RB(158)+RF(162)+RB(168)+RF(170)+RF(175)+RB(183)+RB(190)
     *+RB(248)+RB(249)+RB(250)
      C0=C(10)*(CDOT+DIFF(10)/15.03506029)/DDOT
      C0=C(10)*(CDOT+DIFF(10)/15.03506029+(C(10)-C0)/DT)/DDOT
      R=C0/C(10)
      RB(40)=RB(40)*R
      RB(43)=RB(43)*R
      RB(52)=RB(52)*R
      RF(65)=RF(65)*R
      RF(66)=RF(66)*R
      RF(67)=RF(67)*R
      RF(68)=RF(68)*R
      RF(69)=RF(69)*R
      RF(70)=RF(70)*R
      RF(71)=RF(71)*R
      RF(72)=RF(72)*R
      RF(73)=RF(73)*R
      RF(74)=RF(74)*R
      RF(75)=RF(75)*R
      RF(76)=RF(76)*R
      RF(77)=RF(77)*R
      RF(78)=RF(78)*R
      RF(78)=RF(78)*R
      RF(79)=RF(79)*R
      RF(79)=RF(79)*R
      RB(81)=RB(81)*R
      RB(86)=RB(86)*R
      RB(87)=RB(87)*R
      RB(88)=RB(88)*R
      RB(89)=RB(89)*R
      RB(89)=RB(89)*R
      RB(90)=RB(90)*R
      RB(90)=RB(90)*R
      RB(93)=RB(93)*R
      RF(95)=RF(95)*R
      RB(98)=RB(98)*R
      RF(107)=RF(107)*R
      RF(108)=RF(108)*R
      RF(109)=RF(109)*R
      RB(111)=RB(111)*R
      RB(112)=RB(112)*R
      RB(117)=RB(117)*R
      RF(123)=RF(123)*R
      RB(124)=RB(124)*R
      RB(129)=RB(129)*R
      RB(133)=RB(133)*R
      RB(136)=RB(136)*R
      RB(140)=RB(140)*R
      RF(141)=RF(141)*R
      RF(149)=RF(149)*R
      RB(151)=RB(151)*R
      RF(158)=RF(158)*R
      RB(162)=RB(162)*R
      RF(168)=RF(168)*R
      RB(170)=RB(170)*R
      RB(175)=RB(175)*R
      RF(183)=RF(183)*R
      RF(190)=RF(190)*R
      RF(248)=RF(248)*R
      RF(249)=RF(249)*R
      RF(250)=RF(250)*R
      ENDIF
C
C     CH4
      DDOT=+RB(65)+RB(71)+RB(73)+RB(74)+RB(75)+RF(86)+RF(87)+RF(88)
     *+RF(89)+RF(90)+RB(107)+RB(123)+RB(141)+RB(158)+RB(168)+RB(183)
     *+RB(190)+RB(248)+RB(249)+RB(250)
      TINV=DDOT/C(11)
      IF (TINV .GT. TC) THEN
      CDOT=+RF(65)+RF(71)+RF(73)+RF(74)+RF(75)+RB(86)+RB(87)+RB(88)
     *+RB(89)+RB(90)+RF(107)+RF(123)+RF(141)+RF(158)+RF(168)+RF(183)
     *+RF(190)+RF(248)+RF(249)+RF(250)
      C0=C(11)*(CDOT+DIFF(11)/16.04303026)/DDOT
      C0=C(11)*(CDOT+DIFF(11)/16.04303026+(C(11)-C0)/DT)/DDOT
      R=C0/C(11)
      RB(65)=RB(65)*R
      RB(71)=RB(71)*R
      RB(73)=RB(73)*R
      RB(74)=RB(74)*R
      RB(75)=RB(75)*R
      RF(86)=RF(86)*R
      RF(87)=RF(87)*R
      RF(88)=RF(88)*R
      RF(89)=RF(89)*R
      RF(90)=RF(90)*R
      RB(107)=RB(107)*R
      RB(123)=RB(123)*R
      RB(141)=RB(141)*R
      RB(158)=RB(158)*R
      RB(168)=RB(168)*R
      RB(183)=RB(183)*R
      RB(190)=RB(190)*R
      RB(248)=RB(248)*R
      RB(249)=RB(249)*R
      RB(250)=RB(250)*R
      ENDIF
C
C     CH2O
      DDOT=+RB(38)+RB(39)+RB(42)+RB(46)+RB(51)+RB(58)+RF(59)+RF(60)
     *+RF(61)+RF(62)+RF(63)+RF(64)+RB(66)+RB(70)+RF(75)+RB(80)+RB(83)
     *+RB(84)+RB(85)+RB(102)+RB(113)+RB(118)+RB(129)+RB(133)+RB(147)
     *+RB(164)+RB(167)+RB(171)+RB(186)+RB(189)
      TINV=DDOT/C(12)
      IF (TINV .GT. TC) THEN
      CDOT=+RF(38)+RF(39)+RF(42)+RF(46)+RF(51)+RF(58)+RB(59)+RB(60)
     *+RB(61)+RB(62)+RB(63)+RB(64)+RF(66)+RF(70)+RB(75)+RF(80)+RF(83)
     *+RF(84)+RF(85)+RF(102)+RF(113)+RF(118)+RF(129)+RF(133)+RF(147)
     *+RF(164)+RF(167)+RF(171)+RF(186)+RF(189)
      C0=C(12)*(CDOT+DIFF(12)/30.02649045)/DDOT
      C0=C(12)*(CDOT+DIFF(12)/30.02649045+(C(12)-C0)/DT)/DDOT
      R=C0/C(12)
      RB(38)=RB(38)*R
      RB(39)=RB(39)*R
      RB(42)=RB(42)*R
      RB(46)=RB(46)*R
      RB(51)=RB(51)*R
      RB(58)=RB(58)*R
      RF(59)=RF(59)*R
      RF(60)=RF(60)*R
      RF(61)=RF(61)*R
      RF(62)=RF(62)*R
      RF(63)=RF(63)*R
      RF(64)=RF(64)*R
      RB(66)=RB(66)*R
      RB(70)=RB(70)*R
      RF(75)=RF(75)*R
      RB(80)=RB(80)*R
      RB(83)=RB(83)*R
      RB(84)=RB(84)*R
      RB(85)=RB(85)*R
      RB(102)=RB(102)*R
      RB(113)=RB(113)*R
      RB(118)=RB(118)*R
      RB(129)=RB(129)*R
      RB(133)=RB(133)*R
      RB(147)=RB(147)*R
      RB(164)=RB(164)*R
      RB(167)=RB(167)*R
      RB(171)=RB(171)*R
      RB(186)=RB(186)*R
      RB(189)=RB(189)*R
      ENDIF
C
C     C2H4
      DDOT=+RB(76)+RB(77)+RB(96)+RB(104)+RB(105)+RB(110)+RF(114)+RF(115)
     *+RF(116)+RF(117)+RF(118)+RF(119)+RF(120)+RF(121)+RF(122)+RF(123)
     *+RB(124)+RF(125)+RF(126)+RB(128)+RB(130)+RB(132)+RB(151)+RB(159)
     *+RB(174)+RB(194)+RF(196)+RB(198)+RF(200)+RB(202)+RF(204)+RB(206)
     *+RF(208)+RB(210)+RF(212)+RB(214)+RF(216)+RF(218)
      TINV=DDOT/C(16)
      IF (TINV .GT. TC) THEN
      CDOT=+RF(76)+RF(77)+RF(96)+RF(104)+RF(105)+RF(110)+RB(114)+RB(115)
     *+RB(116)+RB(117)+RB(118)+RB(119)+RB(120)+RB(121)+RB(122)+RB(123)
     *+RF(124)+RB(125)+RB(126)+RF(128)+RF(130)+RF(132)+RF(151)+RF(159)
     *+RF(174)+RF(191)+RF(194)+RB(196)+RF(198)+RB(200)+RF(202)+RB(204)
     *+RF(206)+RB(208)+RF(210)+RB(212)+RF(214)+RB(216)+RB(218)+3*RF(268)
      C0=C(16)*(CDOT+DIFF(16)/28.05418062)/DDOT
      C0=C(16)*(CDOT+DIFF(16)/28.05418062+(C(16)-C0)/DT)/DDOT
      R=C0/C(16)
      RB(76)=RB(76)*R
      RB(77)=RB(77)*R
      RB(96)=RB(96)*R
      RB(104)=RB(104)*R
      RB(105)=RB(105)*R
      RB(110)=RB(110)*R
      RF(114)=RF(114)*R
      RF(115)=RF(115)*R
      RF(116)=RF(116)*R
      RF(117)=RF(117)*R
      RF(118)=RF(118)*R
      RF(119)=RF(119)*R
      RF(120)=RF(120)*R
      RF(121)=RF(121)*R
      RF(122)=RF(122)*R
      RF(123)=RF(123)*R
      RB(124)=RB(124)*R
      RF(125)=RF(125)*R
      RF(126)=RF(126)*R
      RB(128)=RB(128)*R
      RB(130)=RB(130)*R
      RB(132)=RB(132)*R
      RB(151)=RB(151)*R
      RB(159)=RB(159)*R
      RB(174)=RB(174)*R
      RB(194)=RB(194)*R
      RF(196)=RF(196)*R
      RB(198)=RB(198)*R
      RF(200)=RF(200)*R
      RB(202)=RB(202)*R
      RF(204)=RF(204)*R
      RB(206)=RB(206)*R
      RF(208)=RF(208)*R
      RB(210)=RB(210)*R
      RF(212)=RF(212)*R
      RB(214)=RB(214)*R
      RF(216)=RF(216)*R
      RF(218)=RF(218)*R
      ENDIF
C
C     C2H6
      DDOT=+RB(78)+RB(127)+RB(131)+RB(134)+RF(137)+RF(138)+RF(139)
     *+RF(140)+RF(141)
      TINV=DDOT/C(17)
      IF (TINV .GT. TC) THEN
      CDOT=+RF(78)+RF(127)+RF(131)+RF(134)+RB(137)+RB(138)+RB(139)
     *+RB(140)+RB(141)
      C0=C(17)*(CDOT+DIFF(17)/30.07012057)/DDOT
      C0=C(17)*(CDOT+DIFF(17)/30.07012057+(C(17)-C0)/DT)/DDOT
      R=C0/C(17)
      RB(78)=RB(78)*R
      RB(127)=RB(127)*R
      RB(131)=RB(131)*R
      RB(134)=RB(134)*R
      RF(137)=RF(137)*R
      RF(138)=RF(138)*R
      RF(139)=RF(139)*R
      RF(140)=RF(140)*R
      RF(141)=RF(141)*R
      ENDIF
C
C     CH2CHO
      DDOT=+RB(101)+RB(103)+RF(111)+RF(112)+RF(113)
      TINV=DDOT/C(18)
      IF (TINV .GT. TC) THEN
      CDOT=+RF(101)+RF(103)+RB(111)+RB(112)+RB(113)+2*RF(268)
      C0=C(18)*(CDOT+DIFF(18)/43.04561079)/DDOT
      C0=C(18)*(CDOT+DIFF(18)/43.04561079+(C(18)-C0)/DT)/DDOT
      R=C0/C(18)
      RB(101)=RB(101)*R
      RB(103)=RB(103)*R
      RF(111)=RF(111)*R
      RF(112)=RF(112)*R
      RF(113)=RF(113)*R
      ENDIF
C
C     aC3H5
      DDOT=+RB(95)+RB(109)+RB(121)+RB(122)+RB(136)+RF(142)+RF(143)
     *+RF(144)+RF(145)+RF(146)+RF(147)+RF(148)+RF(149)+RB(152)+RB(155)
     *+RB(156)+RB(157)+RB(158)+RB(170)+RB(171)
      TINV=DDOT/C(19)
      IF (TINV .GT. TC) THEN
      CDOT=+RF(95)+RF(109)+RF(121)+RF(122)+RF(136)+RB(142)+RB(143)
     *+RB(144)+RB(145)+RB(146)+RB(147)+RB(148)+RB(149)+RF(152)+RF(155)
     *+RF(156)+RF(157)+RF(158)+RF(170)+RF(171)+RF(191)
      C0=C(19)*(CDOT+DIFF(19)/41.07330096)/DDOT
      C0=C(19)*(CDOT+DIFF(19)/41.07330096+(C(19)-C0)/DT)/DDOT
      R=C0/C(19)
      RB(95)=RB(95)*R
      RB(109)=RB(109)*R
      RB(121)=RB(121)*R
      RB(122)=RB(122)*R
      RB(136)=RB(136)*R
      RF(142)=RF(142)*R
      RF(143)=RF(143)*R
      RF(144)=RF(144)*R
      RF(145)=RF(145)*R
      RF(146)=RF(146)*R
      RF(147)=RF(147)*R
      RF(148)=RF(148)*R
      RF(149)=RF(149)*R
      RB(152)=RB(152)*R
      RB(155)=RB(155)*R
      RB(156)=RB(156)*R
      RB(157)=RB(157)*R
      RB(158)=RB(158)*R
      RB(170)=RB(170)*R
      RB(171)=RB(171)*R
      ENDIF
C
C     C3H6
      DDOT=+RB(108)+RB(142)+RB(146)+RB(148)+RF(150)+RF(151)+RF(152)
     *+RF(153)+RF(154)+RF(155)+RF(156)+RF(157)+RF(158)+RB(163)+RB(165)
     *+RB(166)+RB(168)+RB(175)+RB(195)+RB(199)+RB(203)+RB(207)+RB(211)
     *+RB(215)+RF(220)
      TINV=DDOT/C(20)
      IF (TINV .GT. TC) THEN
      CDOT=+RF(108)+RF(142)+RF(146)+RF(148)+RB(150)+RB(151)+RB(152)
     *+RB(153)+RB(154)+RB(155)+RB(156)+RB(157)+RB(158)+RF(163)+RF(165)
     *+RF(166)+RF(168)+RF(175)+RF(192)+RF(195)+RF(199)+RF(203)+RF(207)
     *+RF(211)+RF(215)+RB(220)
      C0=C(20)*(CDOT+DIFF(20)/42.08127093)/DDOT
      C0=C(20)*(CDOT+DIFF(20)/42.08127093+(C(20)-C0)/DT)/DDOT
      R=C0/C(20)
      RB(108)=RB(108)*R
      RB(142)=RB(142)*R
      RB(146)=RB(146)*R
      RB(148)=RB(148)*R
      RF(150)=RF(150)*R
      RF(151)=RF(151)*R
      RF(152)=RF(152)*R
      RF(153)=RF(153)*R
      RF(154)=RF(154)*R
      RF(155)=RF(155)*R
      RF(156)=RF(156)*R
      RF(157)=RF(157)*R
      RF(158)=RF(158)*R
      RB(163)=RB(163)*R
      RB(165)=RB(165)*R
      RB(166)=RB(166)*R
      RB(168)=RB(168)*R
      RB(175)=RB(175)*R
      RB(195)=RB(195)*R
      RB(199)=RB(199)*R
      RB(203)=RB(203)*R
      RB(207)=RB(207)*R
      RB(211)=RB(211)*R
      RB(215)=RB(215)*R
      RF(220)=RF(220)*R
      ENDIF
C
C     C2H3CHO
      DDOT=+RB(106)+RB(143)+RB(144)+RB(145)+RB(153)+RF(159)+RF(160)
     *+RF(161)
      TINV=DDOT/C(21)
      IF (TINV .GT. TC) THEN
      CDOT=+RF(106)+RF(143)+RF(144)+RF(145)+RF(153)+RB(159)+RB(160)
     *+RB(161)
      C0=C(21)*(CDOT+DIFF(21)/56.06473112)/DDOT
      C0=C(21)*(CDOT+DIFF(21)/56.06473112+(C(21)-C0)/DT)/DDOT
      R=C0/C(21)
      RB(106)=RB(106)*R
      RB(143)=RB(143)*R
      RB(144)=RB(144)*R
      RB(145)=RB(145)*R
      RB(153)=RB(153)*R
      RF(159)=RF(159)*R
      RF(160)=RF(160)*R
      RF(161)=RF(161)*R
      ENDIF
C
C     C4H81
      DDOT=+RB(135)+RB(149)+RB(169)+RB(172)+RF(173)+RF(174)+RF(175)
     *+RF(176)+RF(177)+RF(178)+RF(179)+RF(180)+RF(181)+RF(182)+RF(183)
     *+RB(185)+RB(187)+RB(188)+RB(190)+RF(221)
      TINV=DDOT/C(23)
      IF (TINV .GT. TC) THEN
      CDOT=+RF(135)+RF(149)+RF(169)+RF(172)+RB(173)+RB(174)+RB(175)
     *+RB(176)+RB(177)+RB(178)+RB(179)+RB(180)+RB(181)+RB(182)+RB(183)
     *+RF(185)+RF(187)+RF(188)+RF(190)+RB(221)
      C0=C(23)*(CDOT+DIFF(23)/56.10836124)/DDOT
      C0=C(23)*(CDOT+DIFF(23)/56.10836124+(C(23)-C0)/DT)/DDOT
      R=C0/C(23)
      RB(135)=RB(135)*R
      RB(149)=RB(149)*R
      RB(169)=RB(169)*R
      RB(172)=RB(172)*R
      RF(173)=RF(173)*R
      RF(174)=RF(174)*R
      RF(175)=RF(175)*R
      RF(176)=RF(176)*R
      RF(177)=RF(177)*R
      RF(178)=RF(178)*R
      RF(179)=RF(179)*R
      RF(180)=RF(180)*R
      RF(181)=RF(181)*R
      RF(182)=RF(182)*R
      RF(183)=RF(183)*R
      RB(185)=RB(185)*R
      RB(187)=RB(187)*R
      RB(188)=RB(188)*R
      RB(190)=RB(190)*R
      RF(221)=RF(221)*R
      ENDIF
C
C     C5H9
      DDOT=+RF(191)+RF(192)+RB(217)
      TINV=DDOT/C(24)
      IF (TINV .GT. TC) THEN
      CDOT=+RF(217)
      C0=C(24)*(CDOT+DIFF(24)/69.12748158)/DDOT
      C0=C(24)*(CDOT+DIFF(24)/69.12748158+(C(24)-C0)/DT)/DDOT
      R=C0/C(24)
      RF(191)=RF(191)*R
      RF(192)=RF(192)*R
      RB(217)=RB(217)*R
      ENDIF
C
C     PXC9H19
      DDOT=+RB(209)+RB(212)+RF(220)+RF(229)
      TINV=DDOT/C(30)
      IF (TINV .GT. TC) THEN
      CDOT=+RF(209)+RF(212)+RB(220)+RB(229)
      C0=C(30)*(CDOT+DIFF(30)/127.25178277)/DDOT
      C0=C(30)*(CDOT+DIFF(30)/127.25178277+(C(30)-C0)/DT)/DDOT
      R=C0/C(30)
      RB(209)=RB(209)*R
      RB(212)=RB(212)*R
      RF(220)=RF(220)*R
      RF(229)=RF(229)*R
      ENDIF
C
C     C12H24
      DDOT=+RF(217)+RF(260)+RF(262)+RF(264)
      TINV=DDOT/C(32)
      IF (TINV .GT. TC) THEN
      CDOT=+RB(217)+RF(259)+RF(261)+RF(263)
      C0=C(32)*(CDOT+DIFF(32)/168.32508373)/DDOT
      C0=C(32)*(CDOT+DIFF(32)/168.32508373+(C(32)-C0)/DT)/DDOT
      R=C0/C(32)
      RF(217)=RF(217)*R
      RF(260)=RF(260)*R
      RF(262)=RF(262)*R
      RF(264)=RF(264)*R
      ENDIF
C
C     C12H25O2
      DDOT=+RF(252)+RF(254)+RF(256)+RF(257)
      TINV=DDOT/C(33)
      IF (TINV .GT. TC) THEN
      CDOT=+RF(251)+RF(253)+RF(255)+RF(258)
      C0=C(33)*(CDOT+DIFF(33)/201.33185399)/DDOT
      C0=C(33)*(CDOT+DIFF(33)/201.33185399+(C(33)-C0)/DT)/DDOT
      R=C0/C(33)
      RF(252)=RF(252)*R
      RF(254)=RF(254)*R
      RF(256)=RF(256)*R
      RF(257)=RF(257)*R
      ENDIF
C
C     OC12H23OOH
      DDOT=+RB(267)+RF(268)
      TINV=DDOT/C(34)
      IF (TINV .GT. TC) THEN
      CDOT=+RF(267)
      C0=C(34)*(CDOT+DIFF(34)/216.32328415)/DDOT
      C0=C(34)*(CDOT+DIFF(34)/216.32328415+(C(34)-C0)/DT)/DDOT
      R=C0/C(34)
      RB(267)=RB(267)*R
      RF(268)=RF(268)*R
      ENDIF
      END
