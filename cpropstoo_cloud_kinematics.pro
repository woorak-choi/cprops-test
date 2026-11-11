Function parabolic_function, x, p

;  p - parameters of parabolic function

   yfit = p[0]*(x+p[1])^2.+p[2]

   return, yfit

END

PRO cpropstoo_cloud_kinematics, cloudnum = cloudnum  $
								, smfac = smfac $
								, fitmode = fitmode $
								, infile = infile $
								, inassign = inassign $
								, inprop = inprop $
							    , outfile = outfile $
							    , outchisq_line = outchisq_line $
							    , outchisq_para = outchisq_para 

	                 

;+
; NAME: 
;   cpropstoo_cloud_kinematics
;
; PURPOSE:
;   Perform the kinematics analysis for an individual cloud
;  
;   
; INPUTS:
;   cloudnum     --  the ID of the cloud
;   smfac        --  the smooth factor (in unit of beamsize) of mom1 image
;   fitmode      --  the fitting mode of rotation:
;					 'linear', 'parabolic'
;   infile       --  the fits file of the original (observed) data
;   inassign    --  the fits file of the GMCs assignment cube
;   inprop      --  the idl file of GMC's properties
;    
; OUTPUTS:
;   outfile        --  the '.eps' file for output figure
;   outchisq_line   --  the reduced chiq-square from the binned d_v to
;                      the straight line fit
;   [outchisq_para] --  the reduced chiq-square from the binned d_v to
;                      the parabolic cuver fit
; 
; HISTORY:
;   20170306  LJ   introduced

; HEAD MESSAGE
  print, '========================================='
  print, 'CPROPSTOO: GENERATE GMC KINEMATIC PLOT'
  print, '========================================='

if n_elements(smfac) eq 0 then begin
	smfac = 0.5
endif

if n_elements(fitmode) eq 0 then begin
	fitmode = 'linear'
endif

outchisq_line = !values.f_nan
outchisq_para = !values.f_nan


; READ ORIGIANL DATA CUBE
  data = readfits(infile, hdr)

  vectorify, data, x = x, y = y, v = v, $
	  t = t

; EXTRACT XIN/YIN/VIN FOR THE CLOUDNUM-TH GMC
  assign = readfits(inassign, hd)
  cloud_ind = where(assign eq cloudnum, ct)

  if ct lt 1 then begin
	  print, 'There is no GMC with this cloudnum!'
	  return
  endif

  xin = x[cloud_ind]
  yin = y[cloud_ind]
  vin = v[cloud_ind]
  tin = t[cloud_ind]

; THE FIRST-MOMENT MAP OF THE CLOUD
  pad = 2
  cubify $
	  , x=xin, y=yin, v=vin, t=tin $
	  , cube = minicube $
	  , pad = pad $
	  , location = location $
	  , mask = mask

  l3d = array_indices(minicube, location)

  sz = size(minicube)
  xaxis = (indgen(sz[1])+min(xin)-pad)
  yaxis = (indgen(sz[2])+min(yin)-pad)
  vaxis = (indgen(sz[3])+min(vin)-pad)
  xmap = (intarr(sz[2])+1) ## xaxis
  ymap = yaxis ## (intarr(sz[1])+1)
  vcube = rebin(reform(vaxis,1,1,sz[3]), sz[1:3])

  mom0 = total(minicube,3,/nan)
  mom1 = total(vcube*minicube,3,/nan) / mom0
  mask2d = total(mask,3, /nan)/(total(mask, 3, /nan) > 1)

  ; change velocity from pixel to km/s 
  rdhd, hdr, s = hdstruct
  vaxis_old = hdstruct.v
  if strmatch(sxpar(hdr,'CUNIT3'),'*KM/S*') or $
      strmatch(sxpar(hdr,'CUNIT3'),'*km/s*') then begin
      vaxis_old = vaxis_old*1e3    ; in [km/s]
  endif
  mom1 = mom1*(vaxis_old[1]-vaxis_old[0])+vaxis_old[0]   ; in [km/s]


; SHIFT V_CTR TO ZERO
  ; the cloud's central pixel coordinate
  v_ctr = total(tin*vaxis_old[vin], /NAN)/total(tin, /NAN)
  mom1 = mom1 - v_ctr

; SMOOTH THE FIRST-MOMENT MAP WITH A GAUSSIAN KERNEL WHOSE DISPERSION
; IS HALF THE TELESCOPE BEAMWIDTH
  restore, inprop
  ;ppbeam = !pi/(alog(2))*(props.beamfwhm_pix/2.)^2.
  ppbeam = props.pixperbeam
  ppbeam = ppbeam[uniq(ppbeam, sort(ppbeam))]
  beaminpix = sqrt(ppbeam/!pi*alog(2.))*2./2.354  ; sigma

  if smfac gt 0.0 then begin
	  smfwhm = beaminpix*2.354*smfac
	  psf=psf_gaussian(npixel=[sz[1],sz[2]], $
		  FWHM = smfwhm, /normal)
	  mom1_sm=convol(mom1, psf, /nan, /edge_mirror, /center, $
				missing=NaN, /normalize)
	  mom0_sm=convol(mom0, psf, /nan, /edge_mirror, /center, $
				missing=NaN, /normalize)
  endif else begin
	  mom1_sm=mom1
	  mom0_sm=mom0
  endelse

  mom1_sm = mask2d * mom1_sm
  non_ind=where(mask2d eq 0)
  mom1_sm[non_ind]=-!values.f_nan

  mom0_sm = mask2d * mom0_sm
  non_ind=where(mask2d eq 0)
  mom0_sm[non_ind]=-!values.f_nan

; FIT A PLANE TO THE MOMENT MAP (REQ FREUDENREICH ROUTINES)
  ; vfit = coef(0)+coef(1)*(xin)+coef(2)*(yin)
  mom1vec = reform(mom1_sm[l3d[0,*],l3d[1,*]])

  ; estimate uncertainties
  errcube = (vcube - rebin(mom1_sm, sz[1:3]))^2
  mom1err = (props[cloudnum-1]).noise/mom0_sm * sqrt(total(errcube,3,/nan))
  weights = 1.0/(mom1err^2)
  wt = reform(weights[l3d[0,*],l3d[1,*]])
  ;wt = 0.

  ; xin_pc, yin_pc
  daxis = hdstruct.dec
  distance = (props.dist_pc)[0]
  arctopc = !pi/(180.*3600.)*distance
  pcperpix = abs(daxis[1]-daxis[0])*!dtor*distance
  ;xin_pc = (xin+0.5)*pcperpix
  ;yin_pc = (yin+0.5)*pcperpix
  xin_pc = (xin)*pcperpix
  yin_pc = (yin)*pcperpix

  coef = PLANEFIT(xin_pc, yin_pc, mom1vec, wt, vfit)
  ;print, 'coef', coef
  
  ; refer to Equation (13) of Utomo+ 2015
  v0 = coef[0]
  a = coef[1]
  b = coef[2]

; THE CLOUD ROTATION AXIS
  ; the rotation axis: (y-Y0) = -a/b*(x-X0) ==> (a/b)x+y-(a/b * x0+y0)=0
  ; distance(a1x+b1y+c1=0, (x0,y0)) = abs(a1x0+b1y0+c1)/sqrt(a1^2+b1^2)
  ; the perpendicular distance of each pixel from the rotation axis
  X0 = total(tin*xin, /NAN)/total(tin, /NAN)+0.5
  Y0 = total(tin*yin, /NAN)/total(tin, /NAN)+0.5
  X0_pc = X0*pcperpix
  Y0_pc = Y0*pcperpix
  ;print, 'X0, Y0', X0, Y0
  
  a1 = a/b
  b1 = 1.
  c1 = -(a/b * X0_pc+Y0_pc)
  ; if dot-product of vector (xin-X0, yin-Y0) and (a, b) is positive
  ; then distance of (xin, yin) to rotation axis (-b, a) is positive
  sign = (xin_pc-X0_pc)*a+(yin_pc-Y0_pc)*b
  sign = sign/abs(sign)
  d_v = sign * abs(a1*xin_pc + b1*yin_pc + c1)/sqrt(a1^2. + b1^2.)  ;in [pc]

  
; PLOT THE FIGURE
  set_plot, 'ps'
  device, filename=outfile, bits_per_pixel=64, /encapsulated $
	  , xsize = 8, ysize=8, /inches, /color, xoffset=0., yoffset=0. 

  !p.multi=[0,2,1]


  ; ------------------------------------
  ; subplot - 1. plot smoothed mom1-map 
  ; ------------------------------------
  im = mom1_sm
  ; mom1 image
  ;cgloadct, 27, /silent
  sauron_colormap
  minvalue=min(im, /nan)
  maxvalue=max(im, /nan)
  if minvalue eq maxvalue then begin
	  minvalue=-1 & mavalue=1
  endif
  pos = [0.21, 0.42, 0.94, 0.80]
  CGIMAGE, im, position=pos, stretch=1, /noe $
	  ;, minvalue=minvalue-abs(minvalue)*0.1 $
	  , minvalue=minvalue  $
	  , maxvalue=maxvalue $
	  ;, xtitle='pixels', ytitle='pixels'  $
	  ;, title='Cloud No. '+str(cloudnum) $
	  ;, font=1.2, charsize=2.0 $
	  ;, /Keep_aspect_ratio 
	  , /axis, AXKEYWORDS={xthick:4,ythick:4,xstyle:4, ystyle:4}
  cgloadct, 0, /silent
  
  ; rotation axis
  ; + 0.5 is to locate each pixel by its center position 
  ; i.e., (0.5,0.5) is for (0,0) pixel
  xc = (X0-min(xin)) + pad
  yc = (Y0-min(yin)) + pad

  ;     plot rotation axis
  xarr=findgen(1e4)/1e4*sz(1)
  yarr=-a/b*(xarr-xc)+yc

  cgoplot,[xc],[yc], psym=16
  cgoplot, xarr, yarr, linestyle=0, thick=5., color=cgcolor('black')

  ; add axis
  arctopc = !pi/(180.*3600.)*distance
  CDELT1=SXPAR(hdr, 'CDELT1')
  CDELT2=SXPAR(hdr, 'CDELT2')
  ;xcrange = ([0-0.5,(size(im))[1]-1+0.5]*1. - (xc))*cdelt1*3600.
  ;ycrange = ([0-0.5,(size(im))[2]-1+0.5]*1. - (yc))*cdelt2*3600.
  xcrange = ([0,(size(im))[1]]*1. - (xc))*cdelt1*3600.
  ycrange = ([0,(size(im))[2]]*1. - (yc))*cdelt2*3600.
  
  cgAxis, xAxis=0, xrange=[xcrange[0],xcrange[1]], $
     xtitle=' ', font=1.95, charsize=1.9, $
	 xthick=5, charthick=5,xst=1,yst=1
  cgAxis, xAxis=1, xrange=[xcrange[0]*arctopc,xcrange[1]*arctopc], $
     xtitle=' ', font=1.95, charsize=1.9, $
	 xthick=5, charthick=5,xst=1,yst=1
  cgAxis, yAxis=0, yrange=[ycrange[0],ycrange[1]], $
     ytitle=' ', font=1.95, charsize=1.9, $
	 ythick=5, charthick=5,xst=1,yst=1
  cgAxis, yAxis=1, yrange=[ycrange[0]*arctopc,ycrange[1]*arctopc], $
     ytitle=' ', font=1.95, charsize=1.9, $
	 ythick=5, charthick=5,xst=1,yst=1

  cgtext, 0.29, 0.33, '$\Delta$$\alpha$ (arcsec)', alignment=0.5,$
      charsize = 1.7, charthick=6.0, orientation = 0,/normal
  cgtext, 0.29, 0.86, '$\Delta$$\alpha$ (pc)', alignment=0.5, $
      charsize = 1.7, charthick=6.0, orientation = 0,/normal
  cgtext, 0.03, 0.62, '$\Delta$$\delta$ (arcsec)', alignment=0.5,$
      charsize = 1.7, charthick=6.0, orientation = 90,/normal
  cgtext, 0.52, 0.63, '$\Delta$$\delta$ (pc)', alignment=0.5, $
      charsize = 1.7, charthick=6.0, orientation = -90,/normal

  title='Cloud No. '+str(cloudnum)
  cgtext, 0.29, 0.93, title, alignment=0.5, $
      charsize = 2.0, charthick=6.0, orientation = 0,/normal


  ; color bar
  ;loadct, 27, /silent
  sauron_colormap
  crange=[minvalue-abs(minvalue)*0.1,maxvalue]
  tickint=fix(crange[1]-crange[0])/5
  cgCOLORBAR, range=crange, POSITION=[0.11, 0.23, 0.47, 0.27] $
         , title='v (km s!U-1!N)', tlocation='bottom' $
		 , tickinterval=5, charsize=2.0, font=1.5, charthick=5., thick=5
		 
  loadct, 0, /silent

  ; ------------------------------------
  ; subplot - 2. plot d_v vs v
  ; ------------------------------------
  ; data points for each pixel
  d_v = d_v/arctopc   ; in [arcsec]
  xmin=min(d_v)-0.4*abs(min(d_v))
  xmax=max(d_v)+0.4*abs(max(d_v))
  ymin=min(mom1vec)-0.4*abs(min(mom1vec))
  ymax=max(mom1vec)+0.4*abs(max(mom1vec))

  cgplot, d_v, mom1vec $
	  , psym=16, symcolor='orange', symsize=0.5 $
	  , pos=[0.66, 0.15, 0.96, 0.90] $
	  , xrange=[xmin,xmax], xst=9 $
	  , yrange=[ymin,ymax], yst=1 $
	  , xtitle='Offset from axis (arcsec)' $
	  , xthick = 4, ythick=4 $
	  , xticklen = 0.02, yticklen = 0.05 $
	  ;, ytitle='v (km s!U-1!N)', $
	  , charsize=2.0, charthick=3., font=1.5

  cgtext, 0.60, 0.53, 'v (km s!U-1!N)', alignment=0.5, $
      charsize = 1.7, charthick=6.0, orientation = 90,/normal

  ; x-axis in [pc]
  xcrange = !x.crange
  cgAxis, /xAxis, xrange=[xcrange[0]*arctopc,xcrange[1]*arctopc], $
     xtitle='Offset from axis (pc)', font=1.5, charsize=2.0, xthick=4, $
     charthick=3, xst=1, yst=1




  ; binned d_v and v points
  if fitmode eq 'linear' then begin
	  nbins=10.
  endif else if fitmode eq 'parabolic' then begin
	  nbins=8.
  endif
  yhist=histogram(d_v, nbins=nbins, locations=xbins_hist, reverse_indice = R)
  binsize=(max(d_v)-min(d_v))/(nbins-1.)

  xbins=[]
  ybins=[]
  xerr=[]
  yerr=[]
  xprev=[]
  yprev=[]
  the_total = n_elements(uniq(mom1vec[sort(mom1vec)]))
  for i =0, nbins-1 do begin
    if R[i] ne (R[i+1]) then begin
	  the_array = [mom1vec[R[R[i]:R[i+1]-1]],yprev]
	  if (stddev([mom1vec[R[R[i]:R[i+1]-1]],yprev]) gt 0.0) and $
		  (n_elements(uniq(the_array[sort(the_array)])) ge $
		  (2. > the_total/nbins)) then begin
		;print,'i,[mom1vec[R[R[i]:R[i+1]-1]],prev]',i,[mom1vec[R[R[i]:R[i+1]-1]],prev]
		;xbins=[xbins, xbins_hist[i]+0.5*binsize]
		xbins=[xbins, mean([d_v[R[R[i]:R[i+1]-1]], xprev])]
		ybins=[ybins, mean([mom1vec[R[R[i]:R[i+1]-1]], yprev])]
		xerr=[xerr, 0.0]
		yerr=[yerr, stddev([mom1vec[R[R[i]:R[i+1]-1]], yprev])]

	    ;print,'i, xbins[i], ybins[i], yerr[i]',i,xbins_hist[i]+0.5*binsize,$
		;  mean(mom1vec[R[R[i]:R[i+1]-1]]), stddev([mom1vec[R[R[i]:R[i+1]-1]], prev])
		xprev=[]
		yprev=[]
	  endif else begin
		xprev = [xprev, d_v[R[R[i]:R[i+1]-1]]]
		yprev = [yprev, mom1vec[R[R[i]:R[i+1]-1]]]
	  endelse
	endif
  endfor
  

  oploterror, xbins, ybins, xerr, yerr $
	  , psym=15, color=cgcolor('blue'), errthick=2.5

  ; fitted line or curve
  ;  a) fit a straight line
  sixlin, d_v, mom1vec, a1, siga1, b1, sigb1
  xarr=findgen(1e5)/1e5*(max(d_v) - min(d_v))+min(d_v)
  yarr=a1[1]+xarr*b1[1]
 
  ;  b) refer to Andrae, Schulze-Hartung & Melchior 2010 for the 
  ;  calculation of reduced chi-square:
  if n_elements(ybins) gt 2. then begin
	  outchisq_line=total((ybins-(a1[1]+xbins*b1[1]))^2./yerr^2.)/ $
		  (n_elements(ybins)-2.)
  endif else begin
	  outchisq_line = !values.f_nan
  endelse
  ;print, 'ybins', ybins
  ;print, 'xbins', xbins
  ;print, 'yfit', a1[1]+xbins*b1[1]
  ;print, 'yerr',yerr
  print,'Omega', sqrt(a^2. + b^2.)
  print,'velchisq',outchisq_line

  ;  c) fit paraboic curve
  if fitmode eq 'parabolic' then begin
	  expr = 'parabolic_function(x, p[0:2])'
	  start = [-1., 1., 1.]
	  pi = replicate({fixed:0, limited:[0,0], limits:[0.D, 0.D]}, 3)


	  yerr = make_array(n_elements(mom1vec), value=1.)
	  tranc_result=mpfitexpr(expr, d_v, mom1vec, yerr, start, $
		  parinfo=pi, bestnorm=bestnorm, perror=perror, status=status, $
		  errmsg=errmsg, /quiet)
	  if status le 0 then begin
		  print, 'The fit to the parabolic curve failed!'
	  endif

	  xarr=findgen(1e5)/1e5*(xmax-xmin)+xmin
	  yarr=parabolic_function(xarr, tranc_result)

	  yfits=parabolic_function(xbins, tranc_result)
	  outchisq_para=total((ybins-yfits)^2./yerr^2.)/(nbins-3.)
	  print,'outchisq_para:',outchisq_para
  endif 

  oplot, xarr, yarr, color=cgcolor('blue') $
	  , linestyle=0, thick=3.0
     


  device, /close
  set_plot, 'X'

; CONVERT TO PDF
; https://www.onlineconverter.com/eps-to-pdf


end
