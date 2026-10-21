;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Example of how file format of Gadget can be read-in  ;;
;; in IDL                                               ;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


; Filename

fname="./ics_etg.dat"
;fname="./ics_disk_model_sf_20_hres_v160.dat"
;fname="./ics_disk.dat"
;fname="./ics_disk_exp_stream_35.dat"
BH_logfile  = 'blackholes.txt'

 
    hub=0.71
    rscale=3.11/hub

    redshift=0.0

    npart=lonarr(6)	
    massarr=dblarr(6)
    time=0.0D
    redshift=0.0D
    flag_sfr=0L
    flag_feedback=0L
    npartTotal=lonarr(6)	
    bytesleft=256-6*4 - 6*8 - 8 - 8 - 2*4-6*4
    la=intarr(bytesleft/2)


    print,fname
    openr,1,fname,/f77_unformatted;,/swap_endian
    readu,1,npart,massarr,time,redshift,flag_sfr,flag_feedback,npartTotal,la
  
    print,"Time= ",time

    N=total(npart)
    pos=fltarr(3,N)
    vel=fltarr(3,N)
    id=lonarr(N)
 
    mass=fltarr(N) 


    readu,1,pos
    readu,1,vel
    readu,1,id
 ;   readu,1,mass

    pos=pos/hub
    mass=mass/hub
    massarr=massarr/hub

    NGas=npart(0)
    NHalo=npart(1)
    NDisk=npart(2)
    NBulge=npart(3)
    NStars=npart(4)
    NBH=npart(5)

    if Ngas gt 0 then begin
        u=fltarr(Ngas)
        readu,1,u

 ;       rho=fltarr(Ngas)
 ;       readu,1,rho

    endif  

   close,1


   
    if Ngas gt 0 then begin
        xgas=fltarr(Ngas) &  ygas=fltarr(Ngas)  & zgas=fltarr(Ngas) & mgas=fltarr(Ngas)
        xgas(*)=pos(0,0:Ngas-1)
        ygas(*)=pos(1,0:Ngas-1)
        zgas(*)=pos(2,0:Ngas-1)
        vxgas=fltarr(Ngas) &  vygas=fltarr(Ngas)  & vzgas=fltarr(Ngas) 
        vxgas(*)=vel(0,0:Ngas-1)
        vygas(*)=vel(1,0:Ngas-1)
        vzgas(*)=vel(2,0:Ngas-1)
        if massarr(0) eq 0 then begin
            mgas(*)=mass(0:Ngas-1)	
        endif else begin
            mgas(*)= massarr(0)
            mass(0:Ngas-1)=massarr(0)
	endelse
    endif

    if Nhalo gt 0 then begin
        xhalo=fltarr(NHalo) &  yhalo=fltarr(Nhalo) & zhalo=fltarr(Nhalo) & mhalo=fltarr(Nhalo)
        xhalo(*)=pos(0,0+Ngas:Nhalo+Ngas-1)
        yhalo(*)=pos(1,0+Ngas:Nhalo+Ngas-1)
        zhalo(*)=pos(2,0+Ngas:Nhalo+Ngas-1)
        if massarr(1) eq 0 then begin
	    skip=0L
            for t=0,0 do begin	
                if (npart(t) gt 0) and (massarr(t) eq 0) then begin
			skip=skip + npart(t)
                endif
            endfor
            mhalo(*)=mass(0+skip:Nhalo-1+skip)	
        endif else begin
        skip=0L
            for t=0,0 do begin	
                skip=skip + npart(t)
            endfor
            mhalo(*)= massarr(1)
            mass(0+skip:Nhalo-1+skip)=massarr(1)
        endelse
    endif

    if Ndisk gt 0 then begin
        xdisk=fltarr(NDisk) &  ydisk=fltarr(NDisk) &  zdisk=fltarr(NDisk) & mdisk=fltarr(NDisk)
        vxdisk=fltarr(NDisk) &  vydisk=fltarr(NDisk) &  vzdisk=fltarr(NDisk) 
        xdisk(*)=pos(0,Nhalo+Ngas:Nhalo+Ndisk+Ngas-1)
        ydisk(*)=pos(1,Nhalo+Ngas:Nhalo+Ndisk+Ngas-1)
        zdisk(*)=pos(2,Nhalo+Ngas:Nhalo+Ndisk+Ngas-1)

        vxdisk(*)=vel(0,Nhalo+Ngas:Nhalo+Ndisk+Ngas-1)
        vydisk(*)=vel(1,Nhalo+Ngas:Nhalo+Ndisk+Ngas-1)
        vzdisk(*)=vel(2,Nhalo+Ngas:Nhalo+Ndisk+Ngas-1)



        if massarr(2) eq 0 then begin
	    skip=0L
            for t=0,1 do begin	
                if (npart(t) gt 0) and (massarr(t) eq 0) then begin
			skip=skip + npart(t)
                endif
            endfor
            mdisk(*)=mass(0+skip:Ndisk-1+skip)	
        endif else begin
            skip=0L
            for t=0,1 do begin	
                skip=skip + npart(t)
            endfor
            mdisk(*)= massarr(2)
            mass(0+skip:Ndisk-1+skip)=massarr(2)
        endelse
    endif

    if Nbulge gt 0 then begin
        xbulge=fltarr(NBulge) &  ybulge=fltarr(NBulge) &  zbulge=fltarr(NBulge) & mbulge=fltarr(NBulge)
        xbulge(*)=pos(0,Nhalo+Ndisk+Ngas:Nhalo+Ndisk+Ngas+Nbulge-1)
        ybulge(*)=pos(1,Nhalo+Ndisk+Ngas:Nhalo+Ndisk+Ngas+Nbulge-1)
        zbulge(*)=pos(2,Nhalo+Ndisk+Ngas:Nhalo+Ndisk+Ngas+Nbulge-1)
        if massarr(3) eq 0 then begin
	    skip=0L
            for t=0,2 do begin	
                if (npart(t) gt 0) and (massarr(t) eq 0) then begin
			skip=skip + npart(t)
                endif
            endfor
            mbulge(*)=mass(0+skip:Nbulge-1+skip)	
        endif else begin
            skip=0L
            for t=0,2 do begin	
                skip=skip + npart(t)
            endfor
            mbulge(*)= massarr(3)
            mass(0+skip:Nbulge-1+skip)=massarr(3)
        endelse
    endif


    if Nstars gt 0 then begin
        xstars=fltarr(Nstars) &  ystars=fltarr(Nstars)  & zstars=fltarr(Nstars)  & mstars=fltarr(Nstars)
        xstars(*)=pos(0,Nhalo+Ngas+Ndisk+Nbulge:Nhalo+Ndisk+Ngas+NBulge+Nstars-1)
        ystars(*)=pos(1,Nhalo+Ngas+Ndisk+Nbulge:Nhalo+Ndisk+Ngas+Nbulge+Nstars-1)
        zstars(*)=pos(2,Nhalo+Ngas+Ndisk+Nbulge:Nhalo+Ndisk+Ngas+NBulge+Nstars-1)
        vxstars=fltarr(Nstars) &  vystars=fltarr(Nstars)  & vzstars=fltarr(Nstars) 
        vxstars(*)=vel(0,Nhalo+Ngas+Ndisk+Nbulge:Nhalo+Ndisk+Ngas+NBulge+Nstars-1)
        vystars(*)=vel(1,Nhalo+Ngas+Ndisk+Nbulge:Nhalo+Ndisk+Ngas+Nbulge+Nstars-1)
        vzstars(*)=vel(2,Nhalo+Ngas+Ndisk+Nbulge:Nhalo+Ndisk+Ngas+NBulge+Nstars-1)
        if massarr(4) eq 0 then begin
	    skip=0L
            for t=0,2 do begin	
                if (npart(t) gt 0) and (massarr(t) eq 0) then begin
			skip=skip + npart(t)
                endif
            endfor
            mstars(*)=mass(0+skip:Nstars-1+skip)	
        endif else begin
            skip=0L
            for t=0,2 do begin	
                skip=skip + npart(t)
            endfor
            mstars(*)= massarr(4)
            mass(0+skip:Nstars-1+skip)=massarr(4)
	endelse
    endif

    if NBH gt 0 then begin


        xBH=fltarr(NBH) &  yBH=fltarr(NBH)  & zBH=fltarr(NBH)  & mBH=fltarr(NBH)
        xBH(*)=pos(0,Nhalo+Ngas+Ndisk+Nbulge+Nstars:Nhalo+Ndisk+Ngas+NBulge+Nstars+NBH-1)
        yBH(*)=pos(1,Nhalo+Ngas+Ndisk+Nbulge+Nstars:Nhalo+Ndisk+Ngas+Nbulge+Nstars+NBH-1)
        zBH(*)=pos(2,Nhalo+Ngas+Ndisk+Nbulge+Nstars:Nhalo+Ndisk+Ngas+NBulge+Nstars+NBH-1)
        
        if massarr(5) eq 0 then begin
      	  
        N_mass_tot=n_elements(mass) 
        mBH(*)=mass(N_mass_tot-1)
        
        endif
          
	
    endif
 
;----------------------------------------------------------------------------------------------
  
;  replace disk with bulge 

rdisk=sqrt(xdisk^2+ydisk^2+zdisk^2)
ic=where(rdisk lt 10000000.) 
cmx=total(xdisk(ic)*mdisk(ic))/total(mdisk(ic))
cmy=total(ydisk(ic)*mdisk(ic))/total(mdisk(ic))
cmz=total(zdisk(ic)*mdisk(ic))/total(mdisk(ic))

vcmx=total(vxdisk(ic)*mdisk(ic))/total(mdisk(ic))
vcmy=total(vydisk(ic)*mdisk(ic))/total(mdisk(ic))
vcmz=total(vzdisk(ic)*mdisk(ic))/total(mdisk(ic))

print,'center disk now ',cmx,cmy,cmz
print,'center vel disk now ',vcmx,vcmy,vcmz

rdisk=sqrt(pos(0,*)^2+pos(1,*)^2+pos(2,*)^2)
ic=where(rdisk lt 10000000.) 
cmx=total(pos(0,ic)*mass(ic))/total(mass(ic))
cmy=total(pos(1,ic)*mass(ic))/total(mass(ic))
cmz=total(pos(2,ic)*mass(ic))/total(mass(ic))

vcmx=total(vel(0,ic)*mass(ic))/total(mass(ic))
vcmy=total(vel(1,ic)*mass(ic))/total(mass(ic))
vcmz=total(vel(2,ic)*mass(ic))/total(mass(ic))


print,'center now ',cmx,cmy,cmz
print,'center now ',vcmx,vcmy,vcmz




    
set_plot,'x'
;stop
;device,/portrait,encapsulated=1,/color,filename='initial.ps',xsize=20,ysize=20, xoffset=1,yoffset=1


    loadct,2

    xlen=20.0

    !p.multi=[0,2,2]
    !x.style=1
    !y.style=1
    
;         window,0,xsize=500,ysize=500
; 	wset,0 & plot,(xhalo(*)),(yhalo(*)),psym=3,xrange=[-xlen,xlen],yrange=[-xlen,xlen],color=red
if  (ngas gt 0) then begin 
 	plot,(xgas(*)),(ygas(*)),psym=3,xrange=[-xlen,xlen],yrange=[-xlen,xlen],color=red,/iso
;        oplot,xdisk,ydisk,psym=3,color=100
        
;        plots,(xbh(*)),(ybh(*)),psym=4,color=100,symsize=2

	plot,(xgas(*)),(zgas(*)),psym=3,xrange=[-xlen,xlen],yrange=[-xlen,xlen],color=red,/iso
;        oplot,xdisk,zdisk,psym=3,color=100
;                  plots,(xbh(*)),(zbh(*)),psym=4,color=100,symsize=2
        endif else begin 

 	plot,(xdisk(*)),(ydisk(*)),psym=3,xrange=[-xlen,xlen],yrange=[-xlen,xlen],color=red,/iso
;        oplot,xdisk,ydisk,psym=3,color=100
        
;        plots,(xbh(*)),(ybh(*)),psym=4,color=100,symsize=2

	plot,(xdisk(*)),(zdisk(*)),psym=3,xrange=[-xlen,xlen],yrange=[-xlen,xlen],color=red,/iso
;        oplot,xdisk,zdisk,psym=3,color=100
;                  plots,(xbh(*)),(zbh(*)),psym=4,color=100,symsize=2
endelse


;print,mean(pos(0,*)),mean(pos(1,*)),mean(pos(2,*))
;print,sigma(pos(0,*)),sigma(pos(1,*)),sigma(pos(2,*))

;print,mean(vel(0,*)),mean(vel(1,*)),mean(vel(2,*)) 
;print,sigma(vel(0,*)),sigma(vel(1,*)),sigma(vel(2,*))

;print,'BH position',xBH,yBH,zBH
;print,'BH particle mass',mBH

;print,'BH mass',BH_mass
;print,'BH accretion rate',BH_mdot

vscale=sqrt(43007.)
G=43007.1
x=pos(0,*)
y=pos(1,*)
z=pos(2,*)
m=mass(*)
r=x^2.+y^2.+z^2.
r=sqrt(r)
id=sort(r)
r_sort=r(id)
mcum=total(mass(id),/cumulative)
plot,r(id),vscale*sqrt(mcum/r(id)),xr=[0,20],linestyle=3,yr=[0,300]
print,'total mass',max(mcum)
mcum_total=mcum
r_total=r_sort
vmax=max(vscale*sqrt(mcum/r(id)))
mdisktot=total(mdisk)

print,'Maximum rotation velocity',vmax
print,'Disk mass',mdisktot
print,'Scale length',rscale
alpha=1./rscale
stab=vmax/sqrt(alpha*mdisktot*G)
print,'Stability parameter',stab






if (ngas gt 0) then begin 
rg=xgas^2+ygas^2+zgas^2
rg=sqrt(rg)
id=sort(rg)
r_sort=rg(id)
mcum=total(mgas(id),/cumulative)
oplot,rg(id),vscale*sqrt(mcum/rg(id)),linestyle=2
print,'total gas mass',max(mcum)
mcum_gas=mcum
r_gas=r_sort
endif 

rd=xdisk^2+ydisk^2+zdisk^2
rd=sqrt(rd)
id=sort(rd)
r_sort=rd(id)
mcum=total(mdisk(id),/cumulative)
oplot,rd(id),vscale*sqrt(mcum/rd(id)),linestyle=0
mcum_disk=mcum
r_disk=r_sort
print,'total disk mass',max(mcum)

if (nbulge gt 0) then begin 
rb=xbulge^2+ybulge^2+zbulge^2
rb=sqrt(rb)
id=sort(rb)
r_sort=rb(id)
mcum=total(mbulge(id),/cumulative)
oplot,rb(id),vscale*sqrt(mcum/rb(id)),linestyle=1
mcum_bulge=mcum
r_bulge=r_sort
print,'total bulge mass',max(mcum)
endif

rh=xhalo^2+yhalo^2+zhalo^2
rh=sqrt(rh)
id=sort(rh)
r_sort=rh(id)
mcum=total(mhalo(id),/cumulative)
oplot,rh(id),vscale*sqrt(mcum/rh(id)),linestyle=3
mcum_halo=mcum
r_halo=r_sort
print,'total halo mass',max(mcum)


nrbins=41
rbin=findgen(nrbins)/2.
sd_disk=fltarr(nrbins)
sd_bulge=fltarr(nrbins)
sd_gas=fltarr(nrbins)
sd_halo=fltarr(nrbins)


for k=0,nrbins-2 do begin 
area=!pi*(rbin(k+1)^2-rbin(k)^2)

rd=sqrt(xdisk^2+ydisk^2)
id=where((rd gt rbin(k)) and (rd lt rbin(k+1)),count)
if(count gt 0) then begin  
m_in_area=total(mdisk(id))
sd_disk(k)=m_in_area/area
endif

if (nbulge gt 0) then begin 
rb=sqrt(xbulge^2+ybulge^2)
id=where((rb gt rbin(k)) and (rb lt rbin(k+1)),count)
if(count gt 0) then begin  
m_in_area=total(mbulge(id))
sd_bulge(k)=m_in_area/area
endif
endif

if(ngas gt 0) then begin 
rg=sqrt(xgas^2+ygas^2)
id=where((rg gt rbin(k)) and (rg lt rbin(k+1)),count)
;print,count
;print,area
if(count gt 0) then begin  
m_in_area=total(mgas(id))
sd_gas(k)=m_in_area/area
endif
endif 

rh=sqrt(xhalo^2+yhalo^2)
id=where((rh gt rbin(k)) and (rh lt rbin(k+1)),count)
if(count gt 0) then begin  
m_in_area=total(mhalo(id))
sd_halo(k)=m_in_area/area
endif

endfor 

;sscale=1.5e9/1.e6
sscale=1.e10/(1.d3)^2

plot,/ylog,rbin,sd_disk*sscale,linestyle=0,yr=[0.01,2000]
oplot,rbin,sd_bulge*sscale,linestyle=1
if(ngas gt 0) then begin 
oplot,rbin,sd_gas*sscale,linestyle=2
endif
oplot,rbin,sd_halo*sscale,linestyle=3


legend,['Disk','Bulge','Gas','Halo'],linestyle=[0,1,2,3],charsize=2,pos=[3,0.5],spacing=0.8

;device,/close_file
;set_plot,'x'

;stop
end








