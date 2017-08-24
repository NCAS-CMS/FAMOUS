#include <stdio.h>                                                              
                                                                                
/* Program to reverse the order of bytes in each word of a file.       */       
/* The program can be used to swap UM data files between Big-endian    */       
/* and Little-endian IEEE data formats and vice versa.                 */       
/* Both 32-bit and 64-bit word lengths recognised                      */       
                                                                                
                 /* Written by A. Dickinson 9/9/93 */                          

static int precision;
                                                                                
main(int argc, char *argv[])                                                    
{                                                                               
    FILE *fp1,*fp2;                                                             
    void filecopy(FILE *, FILE *);                                              
    char *prog=argv[0]; 
    char *prec=argv[1];                                                        
                                                                          
    if(argc!=4){                                                                
        fprintf(stderr,"%s: usage: -32|-64 file1 file2\n",prog);                
        exit(1);                                                                
     }                                                               
    
    if(strcmp(prec,  "-32") == 0){
       precision=32;
    }
    else if(strcmp(prec, "-64") == 0){
       precision=64;
    }
    else{
       fprintf(stderr,"%s: usage: -32|-64 file1 file2\n",prog); 
       exit(2);
    }
                                                                            
    if((fp1 = fopen(argv[2],"r")) == NULL){                                     
        fprintf(stderr,"%s: can't open %s\n",prog,argv[2]);                     
        exit(3);                                                                
    }                                                                           
                                                                                
    if((fp2 = fopen(argv[3],"w")) == NULL){                                     
        fprintf(stderr,"%s: can't open %s\n",prog,argv[3]);                     
        exit(4);                                                                
    }                                                                           
                                                                                
    filecopy(fp1,fp2);                                                          
    fclose(fp1);                                                                
    fclose(fp2);                                                                
                                                                                
    if(ferror(fp2)){                                                            
        fprintf(stderr,"%s: error writing %s\n",prog,argv[2]);                  
        exit(4);                                                                
    }                                                                           
                                                                                
    exit(0);                                                                    
                                                                                
}                                                                               
                                                                                
void filecopy(FILE*ifp, FILE *ofp)                                              
                                                                                
{                                                                               
    int i,c,c1,c2,c3,c4;                                                                                                     
    int c5,c6,c7,c8;
                                                                                
    i=0;                                                                        
    while((c1 = getc(ifp)) != EOF){
	c2=getc(ifp); c3=getc(ifp); c4=getc(ifp);
        if(precision == 64){
           c5=getc(ifp); c6=getc(ifp); c7=getc(ifp); c8=getc(ifp);
	   putc(c8,ofp);putc(c7,ofp);putc(c6,ofp);putc(c5,ofp);
        }
	putc(c4,ofp);                                                  
	putc(c3,ofp);                                                  
	putc(c2,ofp);                                                  
	putc(c1,ofp);


                                                  
    }                                                                           
                                                                                
}                                                                               
