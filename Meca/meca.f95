MODULE typeMeca

IMPLICIT NONE

	INTEGER, PARAMETER :: RKIND = SELECTED_REAL_KIND(10, 100)

    TYPE LONGUEUR
		REAL(KIND = RKIND) :: Ap, Dp, h, l, Ep, d
	END TYPE LONGUEUR
	
	TYPE ANGLE
		REAL(KIND = RKIND) :: theta, phi, alpha
	END TYPE ANGLE
	
	TYPE VARIABLE
		REAL(KIND = RKIND) :: min, max, d
	END TYPE VARIABLE

	TYPE OPTIMAL
		REAL(KIND = RKIND) :: k, lA, lD, mom
	END TYPE OPTIMAL

	REAL(KIND = RKIND), PARAMETER :: PI = 3.14159265359_RKIND
	
END MODULE typeMeca





SUBROUTINE PHI (ang, l)

USE typeMeca

IMPLICIT NONE
	
	TYPE(ANGLE) :: ang
	TYPE(LONGUEUR) :: l
	
	
	ang%phi = 2.0*ASIN(SQRT(((l%h - l%Dp*COS(ang%alpha))**2 + (l%Dp*SIN(ang%alpha) - l%Ap)**2)/((l%Ap-l%Dp)**2 + l%h**2)))
	RETURN
	
END SUBROUTINE PHI



SUBROUTINE THETA (ang, l)

USE typeMeca

IMPLICIT NONE

	
	
	TYPE(ANGLE) :: ang
	TYPE(LONGUEUR) :: l
	
	REAL(KIND = RKIND) :: part
	
	part = ACOS((l%Ap - SIN(ang%alpha)*l%Dp)/SQRT((SIN(ang%alpha)*l%Dp - l%Ap)**2 + (l%h - COS(ang%alpha)*l%Dp)**2))
	ang%theta = 3.0*PI/2.0 - ang%phi/2 - part
	
	
	RETURN
	
END SUBROUTINE THETA






PROGRAM Meca

USE typeMeca

IMPLICIT NONE
	
	
	TYPE(ANGLE) :: ang
	TYPE(LONGUEUR) :: l
	TYPE(VARIABLE) :: alpha, Ap, Dp, kv
	TYPE(OPTIMAL) :: opti

	REAL(KIND = RKIND) :: fd, moment, m , g, momMax, k
	INTEGER :: i,ios
	
	OPEN(10, FILE = 'meca.input')
	
	READ(10,*) l%h
	READ(10,*) l%Ap
	READ(10,*) l%Dp
	READ(10,*) l%l
	
	CLOSE(10)
	
	m = 4
	g = 9.81
	l%Ep = SQRT(l%h**2 + l%l**2)
	l%d = 0.5D00*SQRT((l%Dp - l%Ap)**2 + l%h**2)
	
	alpha%min = ATAN(l%l/l%h)
	alpha%max = PI/2.0 + 10*alpha%d
	alpha%d = 0.001

	Ap%min = 0.13
	Ap%max = 0.3
	Ap%d = 0.01
	l%Ap = Ap%min

	Dp%min = 0.13
	Dp%max = 0.35
	Dp%d = 0.01
	l%Dp = 0.15
	
	kv%min = 0.05
	kv%max = 1.1
	kv%d = 0.01
	k = kv%min

	opti%k = -1
	opti%lD = -1
	opti%lA = -1
	opti%mom = 50
	
	!OPEN(11, FILE = 'ang_calc.dat')
	
	
	
	
	!OPEN(12, FILE = 'moment.dat')!+'/'+CHAR(l%Dp,10)+'/'+CHAR(k,10) +
	DO WHILE (k < kv%max)
		l%Dp = Dp%min
		DO WHILE (l%Dp <= Dp%max)
			l%Ap = Ap%min
			DO WHILE(l%Ap <= Ap%max)
				ang%alpha = alpha%min
				!CALL PHI(ang, l)
				!CALL THETA(ang, l)
				!IF (l%Dp > l%Ap) THEN
				momMax = 0
				DO WHILE(ang%alpha < alpha%max)
					CALL PHI(ang, l)
					CALL THETA(ang, l)
					
					!WRITE(11, '(3(F16.10, 1X))') ang%alpha, ang%theta, ang%phi
					
					moment = l%d*m*g*0.5_RKIND*(l%Ep/l%Dp)*SIN(ang%alpha)*SIN(ang%phi)/COS(-ang%alpha+ang%theta+ang%phi) - k*(ang%phi-PI)
					!WRITE (12,  '(2(F16.10, 1X))') ang%alpha, moment
					!PRINT*, ABS(moment)
					IF (ABS(moment) > momMax) THEN
						momMax = ABS(moment)
					END IF
					ang%alpha = ang%alpha + alpha%d
				END DO
				!PRINT*, momMax
				IF (momMax < opti%mom) THEN
					opti%mom = momMax
					opti%lD = l%Dp
					opti%lA = l%Ap
					opti%k = k
				END IF
				WRITE(12, *)
				WRITE(12, *)
				!END IF
				l%Ap = l%Ap + Ap%d
			END DO
			l%Dp = l%Dp + Dp%d
		END DO
		k = k + kv%d
	END DO
	!CLOSE (11)
	!CLOSE (12)

	PRINT*, 'Moment le plus faible trouve : ', opti%mom
	PRINT*, 'k (entre ', kv%min, ' et ', kv%max, ') : ', opti%k
	PRINT*, 'lD (entre ', Dp%min, ' et ', Dp%max, ') : ', opti%lD
	PRINT*, 'lA (entre ', Ap%min, ' et ', Ap%max, ') : ', opti%lA
	
	OPEN(13, FILE = 'BestMoment.dat')

	l%Dp = opti%lD
	l%Ap = opti%lA
	k = opti%k
	ang%alpha = alpha%min
	DO WHILE(ang%alpha < alpha%max)
		CALL PHI(ang, l)
		CALL THETA(ang, l)
		moment = l%d*m*g*0.5_RKIND*(l%Ep/l%Dp)*SIN(ang%alpha)*SIN(ang%phi)/COS(-ang%alpha+ang%theta+ang%phi) - k*(ang%phi-PI)
		ang%alpha = ang%alpha + alpha%d
		WRITE (13,  '(2(F16.10, 1X))') ang%alpha, moment
	END DO
	CLOSE(13)


END PROGRAM Meca
