#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "nrsrc/nr.h"
#include "nrsrc/nrutil.h"

#include "prototypes.h"
#include "globvars.h"





void set_gashalo_velocities(void)
{
  int   i;
  double R,z,r,x,y;
  /*GB  added r */
  double vx,vy,vz,vr,vphi;
  double vch,vcg,vcd,vcb,vccg,mh,mg,md,mb,mcg;
  double upart,spinfactor,lz;
  long dum;
  
  double *rgh,*ugh,*ugh2r;
  const int Nugh = 500;

  double *mgh,*mgh2r;
  const int Nmgh = 500;

  if(N_GASHALO==0) return;

  dum = drand48() * 1e8;

  printf("set gashalo velocities...\t"); fflush(stdout);  


  rgh=dvector(1,Nugh);
  ugh=dvector(1,Nugh);
  ugh2r=dvector(1,Nugh);
  for(i=1;i<=Nugh;i++){
	rgh[i]= pow(10.0, log10(0.01) + i*5.0/(double)(Nugh));
    ugh[i]= G/(GAMMA-1.0)*M_HALO/beta_profile(rgh[i])*qromb(temp_beta,rgh[i],1e4);
    //printf("%f %f\n",rgh[i],ugh[i]);
  }
  spline(rgh,ugh,Nugh,1e40,1e40,ugh2r);

  mgh=dvector(1,Nmgh);
  mgh2r=dvector(1,Nmgh);
  for(i=1;i<=Nmgh;i++){
	rgh[i]= pow(10.0, log10(0.01) + i*5.0/(double)(Nmgh));
    mgh[i]= 4.0*M_PI*RHO_0*qromb(mass_beta,0.0,rgh[i]);
    //printf("%f %f\n",rgh[i],mgh[i]);
  }
  spline(rgh,mgh,Nmgh,1e40,1e40,mgh2r);

  lz=0.0;


  for(i=1;i<=N_GASHALO;i++)
    {
      x=xp_gashalo[i];
      y=yp_gashalo[i];
      z=zp_gashalo[i];
      R=sqrt(xp_gashalo[i]*xp_gashalo[i] + yp_gashalo[i]*yp_gashalo[i]);
      r = sqrt(R*R + z*z);


      mh = M_HALO*(R*R+z*z)/(sqrt(R*R+z*z)+RH)/(sqrt(R*R+z*z)+RH);
	  md = M_DISK*tanh(z/2/Z0)*(1.0-exp(-R/H)*(H+R)/H);
	  mcg= M_GAS*tanh(z/2/Z0)*(1.0-exp(-R/H)*(H+R)/H);
      mb = M_BULGE*(R*R+z*z)/(sqrt(R*R+z*z)+A)/(sqrt(R*R+z*z)+A);
      splint(rgh,mgh,mgh2r,Nmgh,sqrt(R*R+z*z),&mg);

      md = abs(md);
      mcg= abs(mcg);
      //printf("%f %f %f %f %f %f %f %f ",R,z,sqrt(R*R+z*z),mh,md,mcg,mb,mg);

      vr=0;
      vz=0;
	  
	  //printf("%f %f %f\n",sqrt(R*R+z*z),mh,mg);
	  
	  vch = sqrt(G*mh/sqrt(R*R+z*z)*R*R/(R*R+z*z));
	  vcg = sqrt(G*mg/sqrt(R*R+z*z)*R*R/(R*R+z*z));
	  vcd = sqrt(G*md/sqrt(R*R+z*z)*R*R/(R*R+z*z));
	  vcb = sqrt(G*mb/sqrt(R*R+z*z)*R*R/(R*R+z*z));
	  vccg= sqrt(G*mcg/sqrt(R*R+z*z)*R*R/(R*R+z*z));
	  
	  vphi= sqrt(vch*vch+vcg*vcg+vcd*vcd+vcb*vcb+vccg*vccg);

      //printf("%f %f %f %f %f %f\n",vch,vcg,vcd,vcb,vccg,vphi);
      vx=vr*xp_gashalo[i]/R - vphi*yp_gashalo[i]/R;
      vy=vr*yp_gashalo[i]/R + vphi*xp_gashalo[i]/R;

      vxp_gashalo[i]=vx;
      vyp_gashalo[i]=vy;
      vzp_gashalo[i]=vz;

      if (r<R200) lz+= mp_gashalo[i]*(xp_gashalo[i]*vyp_gashalo[i]-yp_gashalo[i]*vxp_gashalo[i]);

      splint(rgh,ugh,ugh2r,Nugh,r,&upart);

      u_gashalo[i] = upart; 

    }
	
  spinfactor = Fspin * JHALO / (M_HALO/(1.0 + RH/R200)/(1.0 + RH/R200)) * MR200_GH / lz;
  if (lz==0.0) spinfactor=0.0;
  
  printf("\n%f %f %f %f %f %f\n",Fspin,JHALO,(M_HALO/(1.0 + RH/R200)/(1.0 + RH/R200)),MR200_GH,lz,spinfactor);


  for(i=1;i<=N_GASHALO;i++)
    {
      vxp_gashalo[i]*=spinfactor;
      vyp_gashalo[i]*=spinfactor;
      vzp_gashalo[i]*=spinfactor;
    }





  //printf("\nJhothalo=%f\nFspin=%f\n",lz,Fspin);
  //printf("\nLx=%f\nLy=%f\nLz=%f\nJhothalo=%f\nmFspin=%f\nJdarkhalo=%f\nDMspin=%f\njhh=%f\njdm=%f\nSpinfactor=%f\nNumber of HH particles=%d\nrspin=%f\n",lx,ly,lz,sqrt(lx*lx+ly*ly+lz*lz),mFspin,jhalo,mDMspin,sqrt(lx*lx+ly*ly+lz*lz)/mFspin,jhalo/mDMspin,spinfactor,Ncount,rspin);
  //printf("\nLx=%f\nLy=%f\nLz=%f\nJdmhalo=%f\nmDMspin=%f\nJhalo=%f\nM200=%f\n",ldmx,ldmy,ldmz,sqrt(ldmx*ldmx+ldmy*ldmy+ldmz*ldmz),mDMspin,jhalo,M_HALO/(1.0 + RH/R200)/(1.0 + RH/R200));


  printf("done.\n"); fflush(stdout);
}




void set_gashalo_positions(void)
{
  int   i,countr,countz;
  double q,R,phi,theta;
  double gasmass;

  if(N_GASHALO==0) return;

  srand48(73);

  printf("set gashalo positions...\t"); fflush(stdout);

  for(i=1,countr=countz=0;i<=N_GASHALO;)
    {
      do
	{
	  q=drand48();
	  R=gashalo_q_to_r(q);
	}
      while(R > 10 * R200);

      phi=drand48()*PI*2;
      theta=acos(drand48()*2-1);
	  
      xp_gashalo[i]=R*sin(theta)*cos(phi);
      yp_gashalo[i]=R*sin(theta)*sin(phi);
      zp_gashalo[i]=R*cos(theta);
      /* hole nominal is 5*H 5*Z0  

      if((R < 5*H) && (fabs(zp_gashalo[i]) < 5*Z0))
	{
	  printf("\nrejected i=%d R=%g  z=%g .... ",i,R,zp_gashalo[i]); fflush(stdout);
	  continue;
	}
*/
       
      i++;

    }


  gasmass = M_GASHALO/N_GASHALO;
  if (N_BULGE!=0) gasmass = mp_bulge[1];
  if (N_DISK!=0) gasmass = mp_disk[1];
  if (N_GAS!=0) gasmass = mp_gas[1];


  for(i=1;i<=N_GASHALO;i++)
#ifdef PURE_HALO
    mp_gashalo[i]=M_GASHALO/N_GASHALO;
#else
    mp_gashalo[i]=gasmass;
#endif

  printf("done\n"); fflush(stdout);

}




