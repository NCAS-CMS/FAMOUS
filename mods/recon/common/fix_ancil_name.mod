*ID FIX_NAME
*/ Allow environment variables to be used in use prognostic file
*DECLARE CONTROL1
*I  CONTROL1.378
      character*80 aname,bname,cname,dname
      integer index,len

*D GPB1F305.17 
        aname=C_NAMELIST(J)
        if (aname(1:1).eq."$") then
           index=scan(aname,'/')
           bname=aname(2:index-1)
           call get_environment_variable(bname,cname)
           len=len_trim(cname)
           dname=cname(1:len)//aname(index:)
           CALL FILE_OPEN(NFT_UPROG,dname,80,0,1,JERR)
        else
           CALL FILE_OPEN(NFT_UPROG,C_NAMELIST(J),80,0,1,JERR)
        endif
