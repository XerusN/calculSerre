MODULE BASE
IMPLICIT NONE
    
    INTEGER, PARAMETER :: RKIND = SELECTED_REAL_KIND(10, 100)
    INTEGER, PARAMETER :: IKIND = SELECTED_INT_KIND(2)
    INTEGER, PARAMETER :: NMAX_X = 178, NMAX_Y = 128, NMAX_Z = 227
    INTEGER(KIND = IKIND), PARAMETER :: BOIS = 1, PLEXI = 3, AIR = 2, AIREXT = 0, HEATER = 4
    REAL(KIND = RKIND), PARAMETER :: H = 25, LBOIS = 0.24, LPLEXI = 0.20, LAIR = 0.025, KELVIN = 273.15, PROD = 50
    REAL(KIND = RKIND), PARAMETER :: TF = 50, DT = 0.1, TINIT = 20 + KELVIN, NAFF = 500
    REAL(KIND = RKIND), PARAMETER :: DV = 10D-9, DA = 10D-6, DL = 10D-3
    
    TYPE :: COORD
        INTEGER :: x, y, z
    END TYPE COORD
    
    
END MODULE BASE





SUBROUTINE READFILE(serreMat, nx, ny, nz)
USE BASE
IMPLICIT NONE
    
    INTEGER(KIND = IKIND), DIMENSION(NMAX_X, NMAX_Y, NMAX_Z) :: serreMat
    INTEGER :: nx, ny, nz
    
    INTEGER :: j, k
    
    OPEN(10, FILE = '~/Documents/VSCode/hugeFile/serre/serrePythonOutput.txt')
    
    READ(10, *) nx, ny, nz
    
    READ(10,*)
    
    READ(10,*)      !We could read the int which represents materials here
    
    READ(10,*)
    
    serreMat(:, :, :) = 0
    
    DO k = 1, nz
        DO j = 1, ny
            READ(10, *) serreMat(1:nx, j, k)
        END DO
        READ(10, *)
    END DO
    
    CLOSE(10)
    
END SUBROUTINE READFILE





SUBROUTINE INITCOEFF(serreCoeff, serreMat, NX, NY, NZ)
USE BASE
IMPLICIT NONE
    
    REAL(KIND = RKIND), DIMENSION(NMAX_X, NMAX_Y, NMAX_Z) :: serreCoeff
    INTEGER(KIND = IKIND), DIMENSION(NMAX_X, NMAX_Y, NMAX_Z) :: serreMat
    INTEGER, INTENT(IN) :: NX, NY, NZ
    
    INTEGER :: i, j, k
    
    DO i = 1, NX
        DO j = 1, NY
            DO k = 1, NZ
                SELECT CASE (serreMat(i, j, k))
                CASE (BOIS)
                    serreCoeff(i, j, k) = - LBOIS * DV / DL
                CASE (PLEXI)
                    serreCoeff(i, j, k) = - LPLEXI * DV / DL
                CASE (AIR)
                    serreCoeff(i, j, k) = - H*DA
                CASE (AIREXT)
                    serreCoeff(i, j, k) = - H*DA
                CASE (HEATER)
                    serreCoeff(i, j, k) = - H*DA
                END SELECT
            END DO
        END DO
    END DO
    
    
    
    
    
END SUBROUTINE INITCOEFF





SUBROUTINE CALCULIT(serreCoeff, serreT, serreMat, TEXT, NX, NY, NZ)
USE BASE
IMPLICIT NONE
    
    REAL(KIND = RKIND), DIMENSION(NMAX_X, NMAX_Y, NMAX_Z) :: serreCoeff, serreT
    INTEGER(KIND = IKIND), DIMENSION(NMAX_X, NMAX_Y, NMAX_Z) :: serreMat
    REAL(KIND = RKIND), INTENT(IN) :: TEXT
    INTEGER, INTENT(IN) :: NX, NY, NZ
    
    REAL(KIND = RKIND), DIMENSION(NMAX_X, NMAX_Y, NMAX_Z) :: tempT
    TYPE(COORD), DIMENSION(6) :: lookAround
    INTEGER :: i, j, k, l
    REAL(KIND = RKIND) :: t, coeff
    
    lookAround(1)%x = -1
    lookAround(1)%y = 0
    lookAround(1)%z = 0
    
    lookAround(2)%x = 1
    lookAround(2)%y = 0
    lookAround(2)%z = 0
    
    lookAround(3)%x = 0
    lookAround(3)%y = -1
    lookAround(3)%z = 0
    
    lookAround(4)%x = 0
    lookAround(4)%y = 1
    lookAround(4)%z = 0
    
    lookAround(5)%x = 0
    lookAround(5)%y = 0
    lookAround(5)%z = -1
    
    lookAround(6)%x = 0
    lookAround(6)%y = 0
    lookAround(6)%z = 1
    
    
    
    
    
    tempT(:, :, :) = serreT(:, :, :)
    
    DO i = 1, NX
        
        
        DO j = 1, NY
            
            
            DO k = 1, NZ
                IF (serreMat(i, j, k) == AIREXT) THEN
                    serreT(i, j, k) = TEXT
                ELSE
                    t = 0
                    DO l = 1, SIZE(lookAround)
                        coeff = serreCoeff(i + lookAround(l)%x, j + lookAround(l)%y, k + lookAround(l)%z)
                        t = t + coeff * (t - tempT(i + lookAround(l)%x, j + lookAround(l)%y, k + lookAround(l)%z))
                    END DO
                    
                    IF (serreMat(i, j, k) == HEATER) THEN
                        t = t + PROD
                    END IF
                    t = t*DT + tempT(i, j, k)
                    serreT(i, j, k) = t
                END IF
                
            END DO
            
            
        END DO
        
        
    END DO
    
    
    
END SUBROUTINE CALCULIT





PROGRAM MAIN
USE BASE
IMPLICIT NONE

    INTEGER(KIND = IKIND), DIMENSION(NMAX_X, NMAX_Y, NMAX_Z) :: serreMat
    REAL(KIND = RKIND), DIMENSION(NMAX_X, NMAX_Y, NMAX_Z) :: serreCoeff, serreT
    INTEGER :: nx, ny, nz, nIt, i, dNa, k, j
    REAL(KIND = RKIND) :: tExt = 20 + KELVIN
    
    
    
    
    CALL READFILE(serreMat, nx, ny, nz)
    
    CALL INITCOEFF(serreCoeff, serreMat, nx, ny, nz)
    
    serreT(:, :, :) = TINIT
    
    nIt = TF/DT + 1 
    dNa = nIt / NAFF + 1
    
    OPEN(10, FILE = '~/Documents/VSCode/hugeFile/serre/serreFortranOutput.txt')
    
    k = 8
    
    DO i = 1, nIt
        PRINT*, i
        CALL CALCULIT(serreCoeff, serreT, serreMat, tExt, nx, ny, nz)
        
        IF (MODULO(i, dNa) == 0) THEN
            DO j = 1, nz
                WRITE(10, *) serreT(:, k, j)
            END DO
            WRITE(10, *)
            WRITE(10, *)
        END IF
    END DO
    
    DO j = 1, ny
        WRITE(10, *) serreT(:, k, j)
    END DO
    WRITE(10, *)
    WRITE(10, *)
    
    CLOSE(10)
    
END PROGRAM

