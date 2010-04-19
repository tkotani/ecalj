      SUBROUTINE HBLSTR(E,LE,STRX,NLMF,NLM,NLF,NF0,HL,
     .    NLMST,CG,INDXCG,JCG,VOL)
C  MAKES LATTICE STRUCTURE CONSTANTS OUT OF BLOCH-SUMS HL.
C     implicit none
      integer nf0,nlf,nlm,nlmf,nlmst,le,indxcg(1),jcg(1)
      double precision cg(1),strx(nlf,nlm),hl(1),e,vol
      integer lmxx,icg,icg1,icg2,ii,ilm,indx,ipow,klm,l,lk,ll,llm,lm,
     .  lmax,mlm
      parameter (lmxx=28)
      double precision efac(0:lmxx),sig(0:lmxx),fourpi,fpibv,sum

      FOURPI=16.D0*DATAN(1.D0)
      FPIBV=FOURPI/VOL
      LMAX=LL(NLMF)+LL(NLM)
      IF(LMAX.GT.12) CALL RX('CHANGE DIMENSIONS IN BLSTRX''')
      EFAC(0)=1.D0
      SIG(0)=1.D0
      DO 1 L=1,LMAX
      EFAC(L)=-E*EFAC(L-1)
  1   SIG(L)=-SIG(L-1)
C ------ ADD TOGETHER CG-SUMS -----------
      DO 11 MLM=1,NLMF
      LM=LL(MLM)
      DO 11 KLM=1,NLM
      LK=LL(KLM)
      SUM=0.D0
      II=MAX0(MLM,KLM)
      INDX=(II*(II-1))/2+MIN0(MLM,KLM)
      ICG1=INDXCG(INDX)
      ICG2=INDXCG(INDX+1)-1
      DO 21 ICG=ICG1,ICG2
      LLM=JCG(ICG)
      IPOW=(LM+LK-LL(LLM))/2
  21  SUM=SUM+CG(ICG)*EFAC(IPOW)*HL(LLM)
  11  STRX(MLM+NF0,KLM)=SUM*FOURPI*SIG(LK)
C -- THE FOLLOWING INCLUDES EXTRA P TERMS 'IMPLICITLY' -----
      IF(LE.EQ.0) THEN
        DO 35 ILM=2,MIN0(4,NLM,NLMF)
  35    STRX(ILM+NF0,ILM)=STRX(ILM+NF0,ILM)-FPIBV
        ENDIF
      RETURN
      END

