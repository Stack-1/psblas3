C
C             Parallel Sparse BLAS  version 2.2
C   (C) Copyright 2006/2007/2008
C                      Salvatore Filippone    University of Rome Tor Vergata
C                      Alfredo Buttari        University of Rome Tor Vergata
C
C Redistribution and use in source and binary forms, with or without
C modification, are permitted provided that the following conditions
C are met:
C   1. Redistributions of source code must retain the above copyright
C      notice, this list of conditions and the following disclaimer.
C   2. Redistributions in binary form must reproduce the above copyright
C      notice, this list of conditions, and the following disclaimer in the
C      documentation and/or other materials provided with the distribution.
C   3. The name of the PSBLAS group or the names of its contributors may
C      not be used to endorse or promote products derived from this
C      software without specific written permission.
C
C THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
C ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
C TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
C PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE PSBLAS GROUP OR ITS CONTRIBUTORS
C BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
C CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
C SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
C INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
C CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
C ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
C POSSIBILITY OF SUCH DAMAGE.
C
C 
C     SUBROUTINE DGELP(TRANS,M,N,P,B,LDB,WORK,LWORK,IERROR)
C
C     Purpose
C     =======
C     Computing
C                 B <-- op(P) B
C     where op(P) is permutation matrix P or its transpose.
C     Notice: when B is a diagonal matrix and it is stored as
C     a vector whose entry i is B(i,i), this routine computes
C                 B <-- P B P(-1)
C     where P(-1) is the inverse of P.
C
C     Parameters
C     ==========
C
C     TRANS    - CHARACTER*1
C             On entry TRANS specifies whetherif the routine operates with
C             matrix P or with the transpose of P as follows:
C                TRANS = 'N'         ->  use matrix P
C                TRANS = 'T' OR 'C'  ->  use the transpose of matrix P
C             Unchanged on exit.
C
C     M        - INTEGER
C             On entry: number of rows of matrix B.
C             Unchanged on exit.
C
C     N        - INTEGER
C             On entry: number of columns of matrix B.
C             Unchanged on exit.
C
C     P        - INTEGER array of dimension (M)
C             On entry P specifies the row permutation of matrix B
C             (P(1) == 0 if no permutation).
C             Unchanged on exit.
C
C     B        - DOUBLE PRECISION array of dimension (LDB,*)
C             On entry: dense matrix.
C             On exit: permuted matrix.
C
C     LDB      - INTEGER
C             On entry: leading dimension of B.
C             Unchanged on exit.
C
C     WORK     - DOUBLE PRECISION array of dimension (LWORK)
C             On entry: work area.
C             On exit INT(WORK(1)) contains the minimum value
C             for LWORK satisfying DGELP memory requirements.
C
C     LWORK    - INTEGER
C             On entry LWORK specifies the dimension of WORK
C             LWORK should be set as follows:
C                LWORK = M
C             Unchanged on exit.
C
C     IERROR   - INTEGER
C             On exit IERROR contains the value of error flag as follows:
C             IERROR = 0   no error
C             IERROR > 0   warning
C             IERROR < 0   fatal error
C
C
      SUBROUTINE DGELP(TRANS,M,N,P,B,LDB,WORK,LWORK,IERROR)
      IMPLICIT NONE                                                    
C     .. Scalar Arguments ..
      INTEGER           LDB, M, N, LWORK, IERROR
      CHARACTER         TRANS
C     .. Array Arguments ..
      DOUBLE PRECISION  B(LDB,*), WORK(*)
      INTEGER           P(*)
C     .. Local Scalars ..
      INTEGER           I, J, ERR_ACT
      logical           istran,isnotran
C     .. Local Arrays..
      INTEGER           INT_VAL
C     .. Intrinsic Functions ..
      INTRINSIC         DBLE
      logical     psb_lsame
      external    psb_lsame

      character*20      name
c
C     .. Executable Statements ..
C
C
C     Check on M, N, LDB, LWORK, TRANS
C
      NAME = 'DGELP\0'
      IERROR = 0
      CALL FCPSB_ERRACTIONSAVE(ERR_ACT)

      isnotran = psb_lsame(trans,'N')
      istran   = psb_lsame(trans,'T') .or. psb_lsame(trans,'C')
      IF     (M.LT.0) THEN
         IERROR = 7
      ELSE IF(N.LT.0) THEN
         IERROR = 3
      ELSE IF(M.GT.LDB) THEN
         IERROR = -6
      ELSE IF (LWORK.LT.M) THEN
         IF (LWORK.EQ.0) THEN
C
C           Return minimum LWORK
C
            IERROR = 8
            WORK(1) = DBLE(M)
            GOTO 9998
         ELSE IF(LWORK.NE.0) THEN
            IERROR = -8
         ENDIF
      ELSE IF (.not.istran.and..not.isnotran) THEN
         IERROR = -1
      ENDIF
C
C     Error handling
C
      IF(IERROR.LT.0) THEN
         CALL FCPSB_ERRPUSH(IERROR,NAME,INT_VAL)
         GOTO 9999
      ENDIF
C
C     Check for M, N, P
C
      IF(M.LE.0 .OR. N.LE.0 .OR. P(1).EQ.0) THEN
         IERROR=5
         CALL FCPSB_ERRPUSH(IERROR,NAME,INT_VAL)
         GOTO 9999
      ENDIF
C
C     Switch on TRANS
C
      IF (isnotran) THEN
C
C        Permuting with P
C
         DO 60 J = 1, N
            DO 20 I = 1, M
               WORK(I) = B(P(I),J)
 20         CONTINUE
            DO 40 I = 1, M
               B(I,J) = WORK(I)
 40         CONTINUE
 60      CONTINUE
      ELSE IF (istran) THEN
C
C        Permuting with the transpose of P.
C
         DO 160 J = 1, N
            DO 120 I = 1, M
               WORK(P(I)) = B(I,J)
 120        CONTINUE
            DO 140 I = 1, M
               B(I,J) = WORK(I)
 140        CONTINUE
 160     CONTINUE
      END IF

 9998 CONTINUE
      CALL FCPSB_ERRACTIONRESTORE(ERR_ACT)
      RETURN

 9999 CONTINUE
      CALL FCPSB_ERRACTIONRESTORE(ERR_ACT)

      IF ( ERR_ACT .NE. 0 ) THEN 
         CALL FCPSB_SERROR()
         RETURN
      ENDIF

      RETURN
      END
