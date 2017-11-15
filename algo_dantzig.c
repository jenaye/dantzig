/*
   /\     /\
  {  `---'  }
  {  O   O  }
  ~~>  V  <~~
   \  \|/  /   Dantzig in C
    `-----'__
    /     \  `^\_
   {       }\ |\_\_   W
   |  \_/  |/ /  \_\_( )
    \__/  /(_E     \__/
      (  /
       MM
*/

#include <math.h>
#define NMAX 6
#define MMAX 6
#define VARMAX 12

// prototype des fonctions
int pl_aps_entrant(double a[MMAX][NMAX],int hb[NMAX],int m,int n,int phase);

int pl_aps_sortant(double a[MMAX][NMAX],int m,int k);

void pl_pivotage(double a[MMAX][NMAX],int db[MMAX],int hb[NMAX], int m,int n,int l,int k);

void pl_aps_affich(double a[MMAX][NMAX],int db[MMAX],int hb[NMAX], int m,int n,int phase);

int  (double a[MMAX][NMAX],double sol[VARMAX], int ineq1,int ineq2, int eq,int n);


main()
{
 int i,j,ineq1,ineq2,eq,n,err;
 double a[MMAX][NMAX],sol[VARMAX];
 /* nb de variables */
 n=2;
 /*
 exo :
      2X1 + X2 + X3 = 18
      2X1 + 3X2 + X4 = 42
      3X1 + 1X2 + X5 = 24
      et Z(X1 X2) = 3X1 + 2X2
 */

 /* fonction economique */
a[0][0]=0;a[0][1]=3;a[0][2]=2;
 /* inequations en <= (nb=ineq1) */
 ineq1=3;
 a[1][0]=18;a[1][1]=2;a[1][2]=1;
 a[2][0]=42;a[2][1]=2;a[2][2]=2;
 a[3][0]=24;a[3][1]=3;a[3][2]=2;
 /* inequations en >= (nb=ineq2) */
 ineq2=0;
 /* equations (nb=eq) */
 eq=0;
 printf("Algorithme dantzig\n");
 err=pl_simplexe_primal(a,sol,ineq1,ineq2,eq,n);
    if(err==1)printf("Solution infinie\n");
        else
    if(err==2)printf("Domaine vide\n");
    else {
        printf("Solution optimale:\n");
        for(i=1;i<=ineq1+ineq2+n;i++)
            printf("x%d =%23.15e\n",i,sol[i]);
            printf("Valeur optimale: z=%23.15f\n",-a[0][0]);
    }
 }


int pl_aps_entrant(double a[MMAX][NMAX],int hb[NMAX],int m,int n,int phase)
{
 int i,j,k,l;
 double d,s,max;
 k=0;
 max=0.0;
 if(phase==2)l=0;
    else l=m+1;
    for(j=1;j<=n;j++) {
        d=a[l][j];
        s=0.0;
    if((d>0)&&(hb[j]!=n+m)) {
        for(i=1;i<=m;i++)
        s+=fabs(a[i][j]);
        d/=s;
        if(d>max){
            max=d;
            k=j;
        }
    }
 }
 return(k);
}

int pl_aps_sortant(double a[MMAX][NMAX],int m,int k)
{
 int i,l;
 double rap,min;
 min=1e308;
 l=0;
 for(i=1;i<=m;i++)
    if(a[i][k]>0) {
        rap=a[i][0]/a[i][k];
        if(rap<min) {
            min=rap;
            l=i;
        }
    }
 return(l);
}

void pl_pivotage(double a[MMAX][NMAX],int db[MMAX],int hb[NMAX],int m,int n,int l,int k)
{
 int i,j;
 double pivot,coef;
 pivot=a[l][k];
 for(i=0;i<=m;i++)
    if(i!=l) {
        coef=a[i][k]/pivot;
        a[i][k]=-coef;
        for(j=0;j<=n;j++)
            if(j!=k)
                a[i][j]=a[i][j]-coef*a[l][j];
    }
    coef=1/pivot;
    a[l][k]=coef;
    for(j=0;j<=n;j++)
        if(j!=k)
            a[l][j]=coef*a[l][j];
            i=db[l];
            db[l]=hb[k];
            hb[k]=i;
 }

void pl_aps_affich(double a[MMAX][NMAX],int db[MMAX],int hb[NMAX], int m,int n,int phase)
{
 int i,j;
 printf("         ");
 for(j=1;j<=n;j++)
    if((phase==1)||(hb[j]!=n+m))
        printf("          x%d",hb[j]);
        printf("\n");
        if(phase==1) {
            printf("z'+");
            for(j=0;j<=n;j++)
                printf("%11.4e ",a[m+1][j]);
                printf("\n");
        }
 for(i=0;i<=m;i++)
 {
  if(i==0)printf("z+ ");
        else
        if(db[i]!=0)printf("x%d ",db[i]);else printf("   ");
            for(j=0;j<=n;j++)
                if((phase==1)||(hb[j]!=n+m))
                    printf("%11.4e ",a[i][j]);
                    printf("\n");
   }
}

int pl_simplexe_primal(double a[MMAX][NMAX],double sol[VARMAX], int ineq1,int ineq2, int eq,int n)
{
 int i,j,k,l,phase,m,m1;
 int db[MMAX],hb[NMAX];
 double min;
 m=ineq1+ineq2+eq;
    for(i=ineq1+1;i<=ineq1+ineq2;i++)
        for(j=0;j<=n;j++)
            a[i][j]=-a[i][j];
            for(i=1;i<=ineq1+ineq2;i++)
                db[i]=n+i;
                for(i=ineq1+ineq2+1;i<=m;i++)
                    db[i]=0;
                    for(j=1;j<=n;j++)
                        hb[j]=j;
                        if(eq!=0) {
  /* printf("Création d'une base de départ\n"); */
   for(i=ineq1+ineq2+1;i<=m;i++){
        l=i;
        k=0;
        for(j=1;j<=n;j++)
            if(a[i][j]!=0)k=j;
                if(k==0) {
                    if(a[i][0]!=0)return(2);
   }
                    else {
    /* pl_aps_affich(a,db,hb,m,n,2); */
    /* printf("var.entrante: x%d\n",hb[k]); */
                        pl_pivotage(a,db,hb,m,n,l,k);
                        hb[k]=hb[n];
                        for(j=0;j<=m;j++)
                        a[j][k]=a[j][n];
                        n-=1;
   }
  }
 }
 n+=1;
 m1=m;
 hb[n]=n+m;
 phase=2;
 l=0;
 min=0;
    for(i=1;i<=m;i++)
        if(a[i][0]<min) {
            min=a[i][0];
            l=i;
        }
        if(l!=0)phase=1;
             k=1;
        if(phase==1){
  /* printf("Phase 1\n"); */
            m1=m+1;
            for(j=0;j<n;j++)
                a[m1][j]=0;
                for(i=1;i<=m;i++)
                    if(a[i][0]<0)
                        a[i][n]=-1;
                    else a[i][n]=0;
                        a[0][n]=0;
                        a[m1][n]=-1;
  /* pl_aps_affich(a,db,hb,m,n,phase); */
  /* printf("var.entrante: x%d   var.sortante: x%d\n",hb[n],db[l]); */
  pl_pivotage(a,db,hb,m1,n,l,n);
 }
 while(phase<=2)
 {
  do
  {
   /* pl_aps_affich(a,db,hb,m,n,phase); */
   k=pl_aps_entrant(a,hb,m,n,phase);
   if(k!=0)
   {
    l=pl_aps_sortant(a,m,k);
    if(l==0)return(1);
    /* printf("var.entrante: x%d   var.sortante: x%d\n",hb[k],db[l]); */
    pl_pivotage(a,db,hb,m1,n,l,k);
   }
  }
  while(k!=0);
  if(phase==1)
  {
   l=0;
   for(i=1;i<=m;i++)
   if(db[i]==n+m)l=i;
   if(l!=0)
   {
    if(fabs(a[l][0])>1e-15)return(2);
    else
    {
     for(j=1;j<=n;j++)
     if(a[l][j]!=0)
     k=j;
     /* printf("var.entrante: x%d   var.sortante: x%d\n",hb[k],db[l]); */
     pl_pivotage(a,db,hb,m1,n,l,k);
    }
   }
  }
  phase+=1;
  m1-=1;
  /* if(phase==2)printf("Phase 2\n"); */
 }
 for(i=1;i<m+n;i++)
    sol[i]=0;
 for(i=1;i<=m;i++)
    sol[db[i]]=a[i][0];
 return(0);
}



