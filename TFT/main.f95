PROGRAM calculSerre

IMPLICIT NONE
	
	INTEGER, PARAMETER :: IKIND = SELECTED_INT_KIND(22)
	INTEGER, PARAMETER :: RKIND = SELECTED_REAL_KIND(12,200)
	
	REAL(KIND = RKIND), PARAMETER :: rho = 1.2			!masse volumique de l'air (SI)
	REAL(KIND = RKIND), PARAMETER :: c = 1004				!Capacite calorifique massique de l'air (SI)
	REAL(KIND = RKIND), PARAMETER :: h = 15					!Coef de convection de l'air (SI)
	REAL(KIND = RKIND), PARAMETER :: lambda_a = 0.025	!conductivité thermique de l'air (SI)
	REAL(KIND = RKIND), PARAMETER :: lambda_p = 0.20	!conductivité thermique du plexiglas (SI)
	REAL(KIND = RKIND), PARAMETER :: lambda_b = 0.24	!conductivité thermique du bois (SI)
	INTEGER(KIND = IKIND), PARAMETER :: temps_max = 173200	!temps maximum de la simulation (s)
	INTEGER(KIND = IKIND), PARAMETER :: journee = 86600		!temps d'une journée (s)
	INTEGER(KIND = IKIND), PARAMETER :: dt = 1		!pas de temps, le changer est très risqué si plusieurs journées

	REAL(KIND = RKIND) :: l_mince, l_epais, l_p, l_a, v, p
	
	REAL(KIND = RKIND) :: alpha
	
	INTEGER(KIND = IKIND) :: i

	LOGICAL :: rad_on		!indique si le radiateur est allumé
	
	TYPE :: MATERIAUX
		REAL(KIND = RKIND) :: p
		REAL(KIND = RKIND) :: pb
		REAL(KIND = RKIND) :: b_mince
		REAL(KIND = RKIND) :: b_epais
	END TYPE MATERIAUX

	TYPE :: TEMPERATURE
		REAL(KIND = RKIND) :: ext	!temperature de la piece (K)
		REAL(KIND = RKIND) :: int	!temperature de la serre (K)
		REAL(KIND = RKIND) :: max	!temperature jusqu'à laquelle chauffer (K)
		REAL(KIND = RKIND) :: min	!temperature à partir de laquelle il faut chauffer (K)
		REAL(KIND = RKIND) :: ini	!temperature initiale pour 1 equa diff, variable de calcul
	END TYPE TEMPERATURE
	
	INTEGER(KIND = IKIND) :: temps, temps_on, temps_off, temps_ini, temps_abs, temps_rad_on

	TYPE(TEMPERATURE) :: t
	TYPE(MATERIAUX) :: s, r
	
	OPEN(15, FILE = 'calculSerre.input')
	
	READ(15,*) s%p			!surface avec plexiglas/air/plexiglas
	READ(15,*) s%pb			!surface avec plexiglas/bois/plexiglas
	READ(15,*) s%b_mince	!surface bois mince
	READ(15,*) l_mince		!epaisseur bois mince
	READ(15,*) s%b_epais	!surface de bois epais
	READ(15,*) l_epais		!epaisseur bois epais
	READ(15,*) v			!volume interieur serre
	READ(15,*) l_p			!epaisseur plexiglas
	READ(15,*) l_a			!epaisseur couche d air
	READ(15,*)
	READ(15,*) temps_on		!heures de debut du chauffage (s)
	READ(15,*) temps_off	!heures de fin du chauffage (s)
	READ(15,*) t%ext		!temperature exterieure (K)
	READ(15,*) t%min		!limite fixée de temperature minimum (K)
	READ(15,*) t%max		!limite fixée de temperature maximum (K)
	READ(15,*) p			!puissance du radiateur (W)
	
	CLOSE(15)
	
	r%p = (2/h + 2*l_p/lambda_p + l_a/lambda_a) / (s%p)			!resistances thermiques des différentes couches
	r%pb = (2/h + 2*l_p/lambda_p + l_a/lambda_b) / (s%pb)
	r%b_mince = (2/h + l_mince/lambda_b) / s%b_mince
	r%b_epais = (2/h + l_epais/lambda_b) / s%b_epais
	alpha = 1.0/(rho*c*v*(r%p + r%pb + r%b_mince + r%b_epais))	!variable de calcul
	
	OPEN(10, FILE='temperature.dat')		!fichier de résultats
	
	temps = 0
	temps_ini = 0		!sert à réajuster le referentiel de temps pour chaque equa diff
	t%ini = t%ext		!idem temperature
	temps_abs = 0		!permet de garder le compte de la vraie durée
	temps_rad_on = 0	!durée d'allumage du radiateur pour chaque jour

	DO WHILE (temps_abs < temps_max)		!limite la durée de calcul
		IF ((temps >= temps_on) .AND. (temps < temps_off)) THEN				!si en periode d'allumage du radiateur

			t%int = (t%ini - t%ext)*EXP(-alpha*(temps-temps_ini)) + t%ext		!premier calcul de t
			IF (rad_on .EQV. .TRUE.) THEN
				t%int = t%int + (t%ext + p/(alpha*rho*c*v))*(1.0 - EXP(-alpha*(temps-temps_ini)))		!ajout de la partie liée à la production
				temps_rad_on = temps_rad_on + dt		!radiateur allumé
			END IF

			IF ((rad_on .EQV. .TRUE.) .AND. (t%int > t%max)) THEN	!verifie si on doit eteindre le radiateur
				rad_on = .FALSE.
				t%ini = t%int
				temps_ini = temps
			ELSE IF ((rad_on .EQV. .FALSE.) .AND. (t%int < t%min)) THEN		!verifie si on doit allumer le radiateur
				rad_on = .TRUE.
				t%ini = t%int
				temps_ini = temps
			END IF

			WRITE(10, '(F10.6, 1X, F10.2)') REAL(temps)/3600.0, t%int-272.15		!écris dans le fichier temps (heures) et temperature de la serre (°C)
			temps = temps + dt			!incrementation du temps
			temps_abs = temps_abs + dt

		ELSE
			t%int = (t%ini - t%ext)*EXP(-alpha*(temps-temps_ini)) + t%ext			!periode ou le chauffage est forcement eteint
			WRITE(10, '(F10.6, 1X, F10.2)') REAL(temps)/3600.0, t%int-272.15
			temps = temps + dt
			temps_abs = temps_abs + dt
		END IF

		IF (temps == temps_on) THEN		!gere la transition vers la phase allumée/eteinte du radiateur
			rad_on = .TRUE.
			temps_ini = temps
			t%ini = t%int
		ELSE IF (temps == temps_off) THEN
			rad_on = .FALSE.
			t%ini = t%int
			temps_ini = temps
		END IF

		IF (temps >= journee) THEN		!gere le changement de journée
			temps_ini = temps_ini - temps
			temps = 0
			PRINT*, temps_rad_on
			temps_rad_on = 0
			WRITE(10, *)
			WRITE(10, *)
		END IF

	END DO
	
	CLOSE(10)
	
END PROGRAM calculSerre
