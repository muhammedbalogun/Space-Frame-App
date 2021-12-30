Imports System

Module Module1

    Sub Main()
        Dim M, NJ, NR, NRJ, J, NBI, K, IA, I As Integer
        Dim E, G, CX, CY, CZ, CXZ, XCL, YCL, ZCL, Xp, Yp, Zp, XPS, YPS, ZPS, YPG, ZPG, SQ, COSA, SINA As Double
        Dim NDJ As Integer = 6
        Dim MD As Integer = 2 * NDJ
        Dim NB As Integer = 0
        Console.WriteLine("KING FAHD UNIVERSITY OF PETROLEUM & MINERALS")
        Console.WriteLine("Civil & Environmental Engineering Department")
        Console.WriteLine("CE 511 - Advance Structural Analysis and Vibrations")
        Console.WriteLine("")
        Console.WriteLine("Course Instructor:")
        Console.WriteLine("Dr Saheed Adekunle")
        Console.WriteLine("")
        Console.WriteLine("Done By:")
        Console.WriteLine("Saeed Mohammed Al-Houri  ID 201431280")
        Console.WriteLine("Muhammed Olawale Balogun  ID 202110150")
        Console.WriteLine("")
        Console.WriteLine("")
        Console.WriteLine("###---- SPACE FRAME PROGRAM ----###")
        Console.WriteLine("")
        Console.WriteLine("")
        Console.WriteLine("#-- STRUCTURAL PARAMETERS --#")
        Console.WriteLine("")
        Console.WriteLine("Enter the number of members")
        M = Console.ReadLine()
        Dim JJ(M - 1), JK(M - 1) As Integer
        Dim AX(M - 1), XI(M - 1), YI(M - 1), ZI(M - 1), EL(M - 1), R11(M - 1), R12(M - 1), R13(M - 1), R21(M - 1), R22(M - 1), R23(M - 1), R31(M - 1), R32(M - 1), R33(M - 1) As Double
        Console.WriteLine("Enter the number of joints")
        NJ = Console.ReadLine()
        Dim ND As Integer
        ND = NDJ * NJ
        Dim JRL(ND - 1), ID(ND - 1) As Integer
        Dim x(NJ - 1), y(NJ - 1), z(NJ - 1) As Double
        Console.WriteLine("Enter the number of support restraints")
        NR = Console.ReadLine()
        Dim N As Integer = ND - NR
        Console.WriteLine("Enter the number of restrained joints")
        NRJ = Console.ReadLine()
        Console.WriteLine("Enter the modulus of elasticity E")
        E = Console.ReadLine()
        Console.WriteLine("Enter the shear modulus of elasticity G")
        G = Console.ReadLine()
        Console.WriteLine("")

        Console.WriteLine("")
        Console.WriteLine("#-- JOINT COORDINATES --#")
        Console.WriteLine("")
        K = 0
        Do While K <= NJ - 1
            J = K + 1
            Console.WriteLine("JOINT " & J)
            Console.WriteLine("")
            Console.WriteLine(" X Coordinates for joint " & J)
            x(K) = Console.ReadLine()
            Console.WriteLine(" Y Coordinates for joint " & J)
            y(K) = Console.ReadLine()
            Console.WriteLine(" Z Coordinates for joint " & J)
            z(K) = Console.ReadLine()
            Console.WriteLine("")
            K = K + 1
        Loop
        Console.WriteLine("")
        Console.WriteLine("#-- MEMBER INFORMATION --#")
        Console.WriteLine("")
        J = 0
        Do While J <= M - 1
            K = J + 1
            Console.WriteLine("MEMBER " & K)
            Console.WriteLine("")
            Console.WriteLine(" J joint for member " & K)
            JJ(J) = Console.ReadLine()
            Console.WriteLine(" K joint for member " & K)
            JK(J) = Console.ReadLine()
            Console.WriteLine(" Cross-sectional area of member " & K)
            AX(J) = Console.ReadLine()
            Console.WriteLine(" Torsion constant Ix of member " & K)
            XI(J) = Console.ReadLine()
            Console.WriteLine(" Moment of inertia about the ym axis Iy of member " & K)
            YI(J) = Console.ReadLine()
            Console.WriteLine(" Moment of inertia about the zm axis Iz of member " & K)
            ZI(J) = Console.ReadLine()
            Console.WriteLine(" The angle identifier (1 or 0) of member " & K)
            IA = Console.ReadLine()
            NBI = NDJ * (Math.Abs(JK(J) - JJ(J)) + 1)
            If NBI > NB Then NB = NBI
            XCL = x(JK(J) - 1) - x(JJ(J) - 1)
            YCL = y(JK(J) - 1) - y(JJ(J) - 1)
            ZCL = z(JK(J) - 1) - z(JJ(J) - 1)
            EL(J) = Math.Sqrt(XCL * XCL + YCL * YCL + ZCL * ZCL)
            CX = XCL / EL(J)
            CY = YCL / EL(J)
            CZ = ZCL / EL(J)
            CXZ = Math.Sqrt(CX * CX + CZ * CZ)
            If IA = 0 Then
                GoTo LANE6
            Else
                Console.WriteLine(" The P point Coordinates of member " & K)
                Console.WriteLine("")
                Console.WriteLine(" X Coordinates of point P for member " & K)
                Xp = Console.ReadLine()
                Console.WriteLine(" Y Coordinates of point P for member " & K)
                Yp = Console.ReadLine()
                Console.WriteLine(" Z Coordinates of point P for member " & K)
                Zp = Console.ReadLine()
                XPS = Xp - x(JJ(J) - 1)
                YPS = Yp - y(JJ(J) - 1)
                ZPS = Zp - z(JJ(J) - 1)
            End If
LANE6:
            If CXZ > 0.001 Then
                GoTo LANE7
            Else
                R11(J) = 0.0
                R12(J) = CY
                R13(J) = 0.0
                R21(J) = -CY
                R22(J) = 0.0
                R23(J) = 0.0
                R31(J) = 0.0
                R32(J) = 0.0
                R33(J) = 1.0
            End If
            If IA = 0 Then
                GoTo LANE1
            Else
                SQ = Math.Sqrt(XPS * XPS + ZPS * ZPS)
                COSA = -XPS * CY / SQ
                SINA = ZPS / SQ
                R21(J) = -CY * COSA
                R23(J) = SINA
                R31(J) = CY * SINA
                R33(J) = COSA
                GoTo LANE1
            End If
LANE7:
            R11(J) = CX
            R12(J) = CY
            R13(J) = CZ
            R21(J) = -CX * CY / CXZ
            R22(J) = CXZ
            R23(J) = -CY * CZ / CXZ
            R31(J) = -CZ / CXZ
            R32(J) = 0.0
            R33(J) = CX / CXZ
            If IA = 0 Then
                GoTo LANE1
            Else
                YPG = R21(J) * XPS + R22(J) * YPS + R23(J) * ZPS
                ZPG = R31(J) * XPS + R32(J) * YPS + R33(J) * ZPS
                SQ = Math.Sqrt(YPG * YPG + ZPG * ZPG)
                COSA = YPG / SQ
                SINA = ZPG / SQ
                R21(J) = (-CX * CY * COSA - CZ * SINA) / CXZ
                R22(J) = CXZ * COSA
                R23(J) = (-CY * CZ * COSA + CX * SINA) / CXZ
                R31(J) = (CX * CY * SINA - CZ * COSA) / CXZ
                R32(J) = -CXZ * SINA
                R33(J) = (CY * CZ * SINA + CX * COSA) / CXZ
            End If
LANE1:
            Console.WriteLine("")
            J = J + 1
        Loop
        'IF YOU WANT TO OUTPUT THE LENGTH AND ROTATION MATRICES
        Console.WriteLine("")
        Console.WriteLine("#-- THE JOINT RESTRAINTS --#")
        Console.WriteLine("")
        Console.WriteLine("If there is a support, put 1, otherwise 0")
        Console.WriteLine("")

        J = 0
        Do While J <= NRJ - 1
            Console.WriteLine(" THE RESTRAINT JOINT # ")
            K = Console.ReadLine()
            Console.WriteLine("")
            Console.WriteLine(" The translation in the x direction for joint " & K)
            JRL(6 * (K - 1)) = Console.ReadLine()
            Console.WriteLine(" The translation in the y direction for joint " & K)
            JRL(6 * (K - 1) + 1) = Console.ReadLine()
            Console.WriteLine(" The translation in the z direction for joint " & K)
            JRL(6 * (K - 1) + 2) = Console.ReadLine()
            Console.WriteLine(" The rotation in the x sense for joint " & K)
            JRL(6 * (K - 1) + 3) = Console.ReadLine()
            Console.WriteLine(" The rotation in the y sense for joint " & K)
            JRL(6 * (K - 1) + 4) = Console.ReadLine()
            Console.WriteLine(" The rotation in the z sense for joint " & K)
            JRL(6 * (K - 1) + 5) = Console.ReadLine()
            Console.WriteLine("")
            J = J + 1
        Loop
        Dim N1 As Integer = 0
        J = 0
        Do While J <= ND - 1
            N1 = N1 + JRL(J)
            If JRL(J) > 0 Then
                GoTo LANE4
            Else
                ID(J) = J - N1
                GoTo LANE5
            End If
LANE4:
            ID(J) = N - 1 + N1
LANE5:
            J = J + 1
        Loop

        'START OF STIFF6
        Dim K1, K2, K3, J1, J2, J3, J4, J5, J6, I1, I2, I3, IR, IC, ITEM, IM(11) As Integer
        Dim SCM1A, SCM1B, SCM2Y, SCM3Y, SCM4Y, SCM2Z, SCM3Z, SCM4Z, SM(11, 11), SMRT(11, 11), SMS(11, 11), SFF(N - 1, NB - 1) As Double
        I = 0
        Do While I <= M - 1
            SCM1A = E * AX(I) / EL(I)
            SCM1B = G * XI(I) / EL(I)
            SCM2Y = 4.0 * E * YI(I) / EL(I)
            SCM3Y = 1.5 * SCM2Y / EL(I)
            SCM4Y = 2.0 * SCM3Y / EL(I)
            SCM2Z = 4.0 * E * ZI(I) / EL(I)
            SCM3Z = 1.5 * SCM2Z / EL(I)
            SCM4Z = 2.0 * SCM3Z / EL(I)
            SM(0, 0) = SCM1A
            SM(1, 1) = SCM4Z
            SM(1, 7) = -SCM4Z
            SM(2, 2) = SCM4Y
            SM(2, 8) = -SCM4Y
            SM(3, 3) = SCM1B
            SM(4, 4) = SCM2Y
            SM(4, 10) = SCM2Y / 2.0
            SM(5, 7) = -SCM3Z
            SM(6, 6) = SCM1A
            SM(7, 11) = -SCM3Z
            SM(8, 10) = SCM3Y
            SM(10, 10) = SCM2Y

            SM(0, 6) = -SCM1A
            SM(1, 5) = SCM3Z
            SM(1, 11) = SCM3Z
            SM(2, 4) = -SCM3Y
            SM(2, 10) = -SCM3Y
            SM(3, 9) = -SCM1B
            SM(4, 8) = SCM3Y
            SM(5, 5) = SCM2Z
            SM(5, 11) = SCM2Z / 2.0
            SM(7, 7) = SCM4Z
            SM(8, 8) = SCM4Y
            SM(9, 9) = SCM1B
            SM(11, 11) = SCM2Z
            J = 0
            Do While J <= 10
                K = J + 1
                Do While K <= 11
                    SM(K, J) = SM(J, K)
                    K = K + 1
                Loop
                J = J + 1
            Loop
            K = 0
            Do While K <= 3
                K1 = 3 * K
                K2 = 3 * K + 1
                K3 = 3 * K + 2
                J = 0
                Do While J <= 11
                    SMRT(J, K1) = SM(J, K1) * R11(I) + SM(J, K2) * R21(I) + SM(J, K3) * R31(I)
                    SMRT(J, K2) = SM(J, K1) * R12(I) + SM(J, K2) * R22(I) + SM(J, K3) * R32(I)
                    SMRT(J, K3) = SM(J, K1) * R13(I) + SM(J, K2) * R23(I) + SM(J, K3) * R33(I)
                    J = J + 1
                Loop
                K = K + 1
            Loop
            J = 0
            Do While J <= 3
                J1 = 3 * J
                J2 = 3 * J + 1
                J3 = 3 * J + 2
                K = J1
                Do While K <= 11
                    SMS(J1, K) = R11(I) * SMRT(J1, K) + R21(I) * SMRT(J2, K) + R31(I) * SMRT(J3, K)
                    SMS(J2, K) = R12(I) * SMRT(J1, K) + R22(I) * SMRT(J2, K) + R32(I) * SMRT(J3, K)
                    SMS(J3, K) = R13(I) * SMRT(J1, K) + R23(I) * SMRT(J2, K) + R33(I) * SMRT(J3, K)
                    K = K + 1
                Loop
                J = J + 1
            Loop
            IM(0) = 6 * (JJ(I) - 1)
            IM(1) = 6 * (JJ(I) - 1) + 1
            IM(2) = 6 * (JJ(I) - 1) + 2
            IM(3) = 6 * (JJ(I) - 1) + 3
            IM(4) = 6 * (JJ(I) - 1) + 4
            IM(5) = 6 * (JJ(I) - 1) + 5
            IM(6) = 6 * (JK(I) - 1)
            IM(7) = 6 * (JK(I) - 1) + 1
            IM(8) = 6 * (JK(I) - 1) + 2
            IM(9) = 6 * (JK(I) - 1) + 3
            IM(10) = 6 * (JK(I) - 1) + 4
            IM(11) = 6 * (JK(I) - 1) + 5
            J = 0
            Do While J <= MD - 1
                I1 = IM(J)
                If JRL(I1) > 0 Then GoTo LANE8
                K = J
                Do While K <= MD - 1
                    I2 = IM(K)
                    If JRL(I2) > 0 Then GoTo LANE3
                    IR = ID(I1)
                    IC = ID(I2)
                    If IR < IC Then GoTo LANE2
                    ITEM = IR
                    IR = IC
                    IC = ITEM
LANE2:
                    IC = IC - IR
                    SFF(IR, IC) = SFF(IR, IC) + SMS(J, K)
LANE3:
                    K = K + 1
                Loop
LANE8:
                J = J + 1
            Loop
            I = I + 1
        Loop

        'START OF BANFAC

        Dim SUM, TEMP As Double
        If SFF(0, 0) <= 0.0 Then
            Console.WriteLine("SFF IS NOT POSITIVE DEFINITE ")
            GoTo LANELAST
        End If
        J = 1
        Do While J <= N - 1
            J1 = J - 1
            J2 = J - NB + 1
            If J2 < 0 Then J2 = 0   'CHECK
            If J1 = 0 Then GoTo LANE11
            I = 1
            Do While I <= J1
                I1 = I - 1
                If I1 < J2 Then GoTo LANE10
                SUM = SFF(I, J - I)
                K = J2
                Do While K <= I1
                    SUM = SUM - SFF(K, I - K) * SFF(K, J - K)
                    K = K + 1
                Loop
                SFF(I, J - I) = SUM
LANE10:
                I = I + 1
            Loop
LANE11:
            SUM = SFF(J, 0)
            K = J2
            Do While K <= J1
                TEMP = SFF(K, J - K) / SFF(K, 0)
                SUM = SUM - TEMP * SFF(K, J - K)
                SFF(K, J - K) = TEMP
                K = K + 1
            Loop
            If SUM <= 0.0 Then
                Console.WriteLine("SFF IS NOT POSITIVE DEFINITE ")
                GoTo LANELAST
            End If
            SFF(J, 0) = SUM
            J = J + 1
        Loop

        ' LOAD COMBINATIONS

        Dim LN As Integer = 0
        LN = LN + 1

        'START OF LDATA6

        Dim NLJ, NLM, LML(M - 1), JR As Integer
        Dim AJ(ND - 1), AML(11, M - 1), AC(ND - 1), AE(ND - 1), DF(N - 1) As Double
LANE14:
        Console.WriteLine("")
        Console.WriteLine("#-- LOADING INFORMATION --#")
        Console.WriteLine("")
        Console.WriteLine("Enter the number of joint loads")
        NLJ = Console.ReadLine()
        Console.WriteLine("Enter the number of member loads")
        NLM = Console.ReadLine()
        Console.WriteLine("")
        Console.WriteLine("ACTIONS AT THE JOINTS")
        Console.WriteLine("")
        If NLJ = 0 Then
            Console.WriteLine("NO ACTIONS AT THE JOINTS")
            Console.WriteLine("")
            GoTo LANE12
        End If
        J = 0
        Do While J <= NLJ - 1
            Console.WriteLine(" THE LOADED JOINT ")
            K = Console.ReadLine()
            Console.WriteLine("")
            Console.WriteLine(" The force in the x direction at joint " & K)
            AJ(6 * (K - 1)) = Console.ReadLine()
            Console.WriteLine(" The force in the y direction at joint " & K)
            AJ(6 * (K - 1) + 1) = Console.ReadLine()
            Console.WriteLine(" The force in the z direction at joint " & K)
            AJ(6 * (K - 1) + 2) = Console.ReadLine()
            Console.WriteLine(" The moment in the x sense at joint " & K)
            AJ(6 * (K - 1) + 3) = Console.ReadLine()
            Console.WriteLine(" The moment in the y sense at joint " & K)
            AJ(6 * (K - 1) + 4) = Console.ReadLine()
            Console.WriteLine(" The moment in the z sense at joint " & K)
            AJ(6 * (K - 1) + 5) = Console.ReadLine()
            Console.WriteLine("")
            J = J + 1
        Loop
LANE12:
        Console.WriteLine("ACTIONS AT THE ENDS OF RESTRAINED MEMBERS DUE TO THE LOADS")
        Console.WriteLine("")
        If NLM = 0 Then
            Console.WriteLine("NO ACTIONS AT MEMBERS")
            Console.WriteLine("")
            Console.WriteLine("")
            GoTo LANE13
        End If
        J = 0
        Do While J <= NLM - 1
            Console.WriteLine(" THE LOADED MEMBER #")
            I = Console.ReadLine()
            Console.WriteLine("")
            Console.WriteLine("AT J JOINT")
            Console.WriteLine("")
            Console.WriteLine(" The force in the x directio of member " & I)
            AML(0, I - 1) = Console.ReadLine()
            Console.WriteLine(" The force in the y direction of member " & I)
            AML(1, I - 1) = Console.ReadLine()
            Console.WriteLine(" The force in the z direction of member " & I)
            AML(2, I - 1) = Console.ReadLine()
            Console.WriteLine(" The moment in the x sense of member " & I)
            AML(3, I - 1) = Console.ReadLine()
            Console.WriteLine(" The moment in the y sense of member " & I)
            AML(4, I - 1) = Console.ReadLine()
            Console.WriteLine(" The moment in the z sense of member " & I)
            AML(5, I - 1) = Console.ReadLine()
            Console.WriteLine("")
            Console.WriteLine("AT K JOINT")
            Console.WriteLine("")
            Console.WriteLine(" The force in the x direction at K joint of member " & I)
            AML(6, I - 1) = Console.ReadLine()
            Console.WriteLine(" The force in the y direction at K joint of member " & I)
            AML(7, I - 1) = Console.ReadLine()
            Console.WriteLine(" The force in the z direction at K joint of member " & I)
            AML(8, I - 1) = Console.ReadLine()
            Console.WriteLine(" The moment in the x sense at K joint of member " & I)
            AML(9, I - 1) = Console.ReadLine()
            Console.WriteLine(" The moment in the y sense at K joint of member " & I)
            AML(10, I - 1) = Console.ReadLine()
            Console.WriteLine(" The moment in the z sense at K joint of member " & I)
            AML(11, I - 1) = Console.ReadLine()
            Console.WriteLine("")
            LML(I - 1) = 1
            J = J + 1
        Loop
LANE13:
        'START OF LOADS6

        If NLM = 0 Then
            GoTo LANE16
        End If
        I = 0
        Do While I <= M - 1
            If LML(I) = 0 Then
                GoTo LANE15
            End If
            IM(0) = 6 * (JJ(I) - 1)
            IM(1) = 6 * (JJ(I) - 1) + 1
            IM(2) = 6 * (JJ(I) - 1) + 2
            IM(3) = 6 * (JJ(I) - 1) + 3
            IM(4) = 6 * (JJ(I) - 1) + 4
            IM(5) = 6 * (JJ(I) - 1) + 5
            IM(6) = 6 * (JK(I) - 1)
            IM(7) = 6 * (JK(I) - 1) + 1
            IM(8) = 6 * (JK(I) - 1) + 2
            IM(9) = 6 * (JK(I) - 1) + 3
            IM(10) = 6 * (JK(I) - 1) + 4
            IM(11) = 6 * (JK(I) - 1) + 5
            J = 0
            Do While J <= 3
                J1 = 3 * J
                J2 = 3 * J + 1
                J3 = 3 * J + 2
                I1 = IM(J1)
                I2 = IM(J2)
                I3 = IM(J3)

                AE(I1) = AE(I1) - R11(I) * AML(J1, I) - R21(I) * AML(J2, I) - R31(I) * AML(J3, I)
                AE(I2) = AE(I2) - R12(I) * AML(J1, I) - R22(I) * AML(J2, I) - R32(I) * AML(J3, I)
                AE(I3) = AE(I3) - R13(I) * AML(J1, I) - R23(I) * AML(J2, I) - R33(I) * AML(J3, I)

                J = J + 1
            Loop
LANE15:
            I = I + 1
        Loop
LANE16:
        J = 0
        Do While J <= ND - 1
            JR = ID(J)
            AC(JR) = AJ(J) + AE(J)
            J = J + 1
        Loop

        ' BANSOL

        I = 0
        Do While I <= N - 1
            J = I - NB + 1
            If I <= NB - 1 Then J = 0 'CHECK
            SUM = AC(I)
            K1 = I - 1
            If J > K1 Then GoTo LANE17
            K = J
            Do While K <= K1
                SUM = SUM - SFF(K, I - K) * DF(K)
                K = K + 1
            Loop
LANE17:
            DF(I) = SUM
            I = I + 1
        Loop
        I = 0
        Do While I <= N - 1
            DF(I) = DF(I) / SFF(I, 0)
            I = I + 1
        Loop
        I1 = 0
        Do While I1 <= N - 1
            I = N - I1 - 1
            J = I + NB - 1
            If J > N - 1 Then J = N - 1 'CHECK
            SUM = DF(I)
            K2 = I + 1
            If K2 > J Then GoTo LANE18
            K = K2
            Do While K <= J
                SUM = SUM - SFF(I, K - I) * DF(K)
                K = K + 1
            Loop
LANE18:
            DF(I) = SUM
            I1 = I1 + 1
        Loop

        ' RESUL6

        Dim JE As Integer
        Dim DJ(ND - 1), AMD(MD - 1), AM(11), AR(ND - 1) As Double
        J = N 'CHECK
        K = 0
        Do While K <= ND - 1
            JE = ND - K - 1
            If JRL(JE) = 0 Then GoTo LANE19
            DJ(JE) = 0.0
            GoTo LANE20
LANE19:
            J = J - 1
            DJ(JE) = DF(J)
LANE20:
            K = K + 1
        Loop
        Console.WriteLine("#-- THE RESULTS --#")
        Console.WriteLine("")
        Console.WriteLine("")
        Console.WriteLine("THE JOINT DISPLACEMENTS")
        Console.WriteLine("")
        J = 0
        Do While J <= NJ - 1
            K = J + 1
            Console.WriteLine("THE JOINT DISPLACEMENTS AT JOINT " & (J + 1))
            Console.WriteLine("")
            Console.WriteLine(" THE TRANSLATION DISPLACEMENT IN THE X DIRECTION   " & DJ(6 * J))
            Console.WriteLine(" THE TRANSLATION DISPLACEMENT IN THE Y DIRECTION   " & DJ(6 * J + 1))
            Console.WriteLine(" THE TRANSLATION DISPLACEMENT IN THE Z DIRECTION   " & DJ(6 * J + 2))
            Console.WriteLine(" THE ROTATION DISPLACEMENT IN THE X SENSE   " & DJ(6 * J + 3))
            Console.WriteLine(" THE ROTATION DISPLACEMENT IN THE Y SENSE   " & DJ(6 * J + 4))
            Console.WriteLine(" THE ROTATION DISPLACEMENT IN THE Z SENSE   " & DJ(6 * J + 5))
            Console.WriteLine("")
            J = J + 1
        Loop
        Console.WriteLine("")
        Console.WriteLine("THE MEMBER END-ACTIONS")
        Console.WriteLine("")
        I = 0
        Do While I <= M - 1
            SCM1A = E * AX(I) / EL(I)
            SCM1B = G * XI(I) / EL(I)
            SCM2Y = 4.0 * E * YI(I) / EL(I)
            SCM3Y = 1.5 * SCM2Y / EL(I)
            SCM4Y = 2.0 * SCM3Y / EL(I)
            SCM2Z = 4.0 * E * ZI(I) / EL(I)
            SCM3Z = 1.5 * SCM2Z / EL(I)
            SCM4Z = 2.0 * SCM3Z / EL(I)
            SM(0, 0) = SCM1A
            SM(1, 1) = SCM4Z
            SM(1, 7) = -SCM4Z
            SM(2, 2) = SCM4Y
            SM(2, 8) = -SCM4Y
            SM(3, 3) = SCM1B
            SM(4, 4) = SCM2Y
            SM(4, 10) = SCM2Y / 2.0
            SM(5, 7) = -SCM3Z
            SM(6, 6) = SCM1A
            SM(7, 11) = -SCM3Z
            SM(8, 10) = SCM3Y
            SM(10, 10) = SCM2Y

            SM(0, 6) = -SCM1A
            SM(1, 5) = SCM3Z
            SM(1, 11) = SCM3Z
            SM(2, 4) = -SCM3Y
            SM(2, 10) = -SCM3Y
            SM(3, 9) = -SCM1B
            SM(4, 8) = SCM3Y
            SM(5, 5) = SCM2Z
            SM(5, 11) = SCM2Z / 2.0
            SM(7, 7) = SCM4Z
            SM(8, 8) = SCM4Y
            SM(9, 9) = SCM1B
            SM(11, 11) = SCM2Z
            J = 0
            Do While J <= 10
                K = J + 1
                Do While K <= 11
                    SM(K, J) = SM(J, K)
                    K = K + 1
                Loop
                J = J + 1
            Loop
            K = 0
            Do While K <= 3
                K1 = 3 * K
                K2 = 3 * K + 1
                K3 = 3 * K + 2
                J = 0
                Do While J <= 11
                    SMRT(J, K1) = SM(J, K1) * R11(I) + SM(J, K2) * R21(I) + SM(J, K3) * R31(I)
                    SMRT(J, K2) = SM(J, K1) * R12(I) + SM(J, K2) * R22(I) + SM(J, K3) * R32(I)
                    SMRT(J, K3) = SM(J, K1) * R13(I) + SM(J, K2) * R23(I) + SM(J, K3) * R33(I)
                    J = J + 1
                Loop
                K = K + 1
            Loop
            IM(0) = 6 * (JJ(I) - 1)
            IM(1) = 6 * (JJ(I) - 1) + 1
            IM(2) = 6 * (JJ(I) - 1) + 2
            IM(3) = 6 * (JJ(I) - 1) + 3
            IM(4) = 6 * (JJ(I) - 1) + 4
            IM(5) = 6 * (JJ(I) - 1) + 5
            IM(6) = 6 * (JK(I) - 1)
            IM(7) = 6 * (JK(I) - 1) + 1
            IM(8) = 6 * (JK(I) - 1) + 2
            IM(9) = 6 * (JK(I) - 1) + 3
            IM(10) = 6 * (JK(I) - 1) + 4
            IM(11) = 6 * (JK(I) - 1) + 5
            J = 0
            Do While J <= MD - 1
                AMD(J) = 0.0
                K = 0
                Do While K <= MD - 1
                    I1 = IM(K)
                    AMD(J) = AMD(J) + SMRT(J, K) * DJ(I1)
                    K = K + 1
                Loop
                AM(J) = AML(J, I) + AMD(J)
                J = J + 1
            Loop
            Console.WriteLine("THE MEMBER END-ACTIONS FOR MEMBER " & (I + 1))
            Console.WriteLine("")
            Console.WriteLine("AT J JOINT")
            Console.WriteLine(" THE FORCE IN THE X DIRECTION   " & AM(0))
            Console.WriteLine(" THE FORCE IN THE Y DIRECTION   " & AM(1))
            Console.WriteLine(" THE FORCE IN THE Z DIRECTION   " & AM(2))
            Console.WriteLine(" THE MOMENT IN THE X SENSE   " & AM(3))
            Console.WriteLine(" THE MOMENT IN THE Y SENSE   " & AM(4))
            Console.WriteLine(" THE MOMENT IN THE Z SENSE   " & AM(5))
            Console.WriteLine("")
            Console.WriteLine("AT K JOINT")
            Console.WriteLine(" THE FORCE IN THE X DIRECTION   " & AM(6))
            Console.WriteLine(" THE FORCE IN THE Y DIRECTION   " & AM(7))
            Console.WriteLine(" THE FORCE IN THE Z DIRECTION   " & AM(8))
            Console.WriteLine(" THE MOMENT IN THE X SENSE   " & AM(9))
            Console.WriteLine(" THE MOMENT IN THE Y SENSE   " & AM(10))
            Console.WriteLine(" THE MOMENT IN THE Z SENSE   " & AM(11))
            Console.WriteLine("")
            J = 0
            Do While J <= 3
                J1 = 3 * J
                J2 = 3 * J + 1
                J3 = 3 * J + 2
                I1 = IM(J1)
                I2 = IM(J2)
                I3 = IM(J3)
                If JRL(I1) = 1 Then
                    AR(I1) = AR(I1) + R11(I) * AMD(J1) + R21(I) * AMD(J2) + R31(I) * AMD(J3)
                End If
                If JRL(I2) = 1 Then
                    AR(I2) = AR(I2) + R12(I) * AMD(J1) + R22(I) * AMD(J2) + R32(I) * AMD(J3)
                End If
                If JRL(I3) = 1 Then
                    AR(I3) = AR(I3) + R13(I) * AMD(J1) + R23(I) * AMD(J2) + R33(I) * AMD(J3)
                End If
                J = J + 1
            Loop
            I = I + 1
        Loop
        J = 0
        Do While J <= ND - 1
            If JRL(J) = 0 Then GoTo LANE22
            AR(J) = AR(J) - AJ(J) - AE(J)
LANE22:
            J = J + 1
        Loop
        Console.WriteLine("")
        Console.WriteLine("THE SUPPORT REACTIONS")
        Console.WriteLine("")
        J = 0
        Do While J <= NJ - 1
            J1 = 6 * J
            J2 = 6 * J + 1
            J3 = 6 * J + 2
            J4 = 6 * J + 3
            J5 = 6 * J + 4
            J6 = 6 * J + 5
            N1 = JRL(J1) + JRL(J2) + JRL(J3) + JRL(J4) + JRL(J5) + JRL(J6)
            If N1 = 0 Then GoTo LANE23
            Console.WriteLine("THE SUPPORT REACTIONS FOR JOINT " & (J + 1))
            Console.WriteLine("")
            Console.WriteLine(" THE FORCE IN THE X DIRECTION   " & AR(J1))
            Console.WriteLine(" THE FORCE IN THE Y DIRECTION   " & AR(J2))
            Console.WriteLine(" THE FORCE IN THE Z DIRECTION   " & AR(J3))
            Console.WriteLine(" THE MOMENT IN THE X SENSE   " & AR(J4))
            Console.WriteLine(" THE MOMENT IN THE Y SENSE   " & AR(J5))
            Console.WriteLine(" THE MOMENT IN THE Z SENSE   " & AR(J6))
            Console.WriteLine("")
LANE23:
            J = J + 1
        Loop
LANELAST:
        Console.ReadLine()
        Console.ReadLine()
    End Sub

End Module
