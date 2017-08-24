*DECLARE SFLINT7A
*I SFLINT7A.370
c IF SQRTCD_K too small, CDTEMP1 explodes and CDR10M goes NaN
            IF ( CDR10M(I).ne.CDR10M(I) ) CDR10M(I)=0.
            IF ( abs(CDR10M(I)).gt.1e6 ) CDR10M(I)=0.
*I SFLINT7A.410
c IF SQRTCD_K too small, CDTEMP1 explodes and CHR1P5M goes -Inf
            IF ( abs(CHR1P5M(I)).gt.1e6 ) CHR1P5M(I)=0.
