load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
begin
    
    OUTPUT_TYPE = "png"
    OUTDIR      = "/home/reproduce/sunc/enso/plot/"
    CESMDIR     = "/home/reproduce/sunc/enso/CESM1-FASTCHEM_picontrol_atm/"
    FGOALSDIR   = "/home/reproduce/sunc/enso/FGOALS-g2_picontrol_atm/"
    ERA40       = "/home/reproduce/sunc/enso/obs_data/"
    ;---------------------------------------------------------------------------------------

    wks_type = OUTPUT_TYPE
    wks_BJ_fd = gsn_open_wks(wks_type,OUTDIR+"atmopheric_Bjerkes_feedback")

    ;--------------------------------------------------------------------------------------
    data = addfile(ERA40+"ERA40_wind_ts.195801-200112.nc","r") 
    u_obs = short2flt(data->u10(:,:,:))
    v_obs = short2flt(data->v10(:,:,:))    
    t2m_obs = short2flt(data->t2m(:,:,:))
    lat_obs = data->latitude
    lon_obs = data->longitude
    time_obs = data->time

    pi=4.*atan(1.0)
    rad=(pi/180.)
    coswgt_obs=cos(rad*lat_obs)
    coswgt_obs!0 = "lat" 
    coswgt_obs&lat= lat_obs 

    rho_air=1.23;  kg/m^3 
    cd=1.2e-3; non-dimensional drag
    wmag = sqrt(u_obs^2+v_obs^2);% wind speed
    zws_obs=rho_air*cd*wmag*u_obs; %Zonal Wind Stress N/m^2
    ;tauy=rho_air*cd*wmag*v; %wind stress curl"/units="N/m^2
    copy_VarCoords(u_obs, zws_obs)

    zws_obs_anl   = rmMonAnnCycTLL(zws_obs)
    copy_VarCoords(zws_obs, zws_obs_anl)
    zws_obs_nino4 = wgt_areaave_Wrap(zws_obs_anl(:,{-5:5},{160:210}),coswgt_obs({-5:5}),1.0,1)

    d = addfile("$NCARG_ROOT/lib/ncarg/data/cdf/landsea.nc","r")   
    basemap = d->LSMASK
    lsm = landsea_mask(basemap,t2m_obs&latitude,t2m_obs&longitude)
    sst_obs = mask(t2m_obs,conform(t2m_obs,lsm,(/1,2/)).ge.1,False)
    delete(lsm)
    copy_VarCoords(t2m_obs, sst_obs)   
    sst_obs_anl = rmMonAnnCycTLL(sst_obs)
    copy_VarCoords(sst_obs, sst_obs_anl)
    sst_obs_nino3 = wgt_areaave_Wrap(sst_obs_anl(:,{-5:5},{210:270}),coswgt_obs({-5:5}),1.0,1)
    printVarSummary(zws_obs_anl)
    bj_obs = regCoef(sst_obs_nino3(:),zws_obs_nino4(:))
    printVarSummary(bj_obs)
    ;----------------------------------------------------------------------------------------------
    data = addfile(CESMDIR+"tauu_Amon_CESM1-FASTCHEM_piControl_r1i1p1_007001-029112.nc","r") 
    zws_cesm = data->tauu(0:12*200-1,:,:)    ;N/m^2
    lat_cesm = data->lat
    lon_cesm = data->lon
    time_cesm = data->time
    ;printVarSummary(zws_cesm)

    data = addfile(CESMDIR+"ts_data/"+"CESM1-FASTCHEM_piControl_r1i1p1.ts.007001-029112.nc","r")
    sst_cesm  = data->ts(0:200*12-1,:,:)
    sst_cesm =(/sst_cesm-273.15/)
    sst_cesm = where(sst_cesm.le.-1.8,-1.8,sst_cesm)
    d = addfile("$NCARG_ROOT/lib/ncarg/data/cdf/landsea.nc","r")   
    basemap = d->LSMASK
    lsm = landsea_mask(basemap,sst_cesm&lat,sst_cesm&lon)
    sst_cesm = mask(sst_cesm,conform(sst_cesm,lsm,(/1,2/)).ge.1,False)
    delete(lsm)

    pi=4.*atan(1.0)
    rad=(pi/180.)
    coswgt_cesm=cos(rad*lat_cesm)
    coswgt_cesm!0 = "lat" 
    coswgt_cesm&lat= lat_cesm 

    sst_cesm_anl = rmMonAnnCycTLL(sst_cesm)
    zws_cesm_anl = rmMonAnnCycTLL(zws_cesm)
    copy_VarCoords(zws_cesm, zws_cesm_anl)
    copy_VarCoords(sst_cesm, sst_cesm_anl)
    sst_cesm_nino3 = wgt_areaave_Wrap(sst_cesm_anl(:,{-5:5},{210:270}),coswgt_cesm({-5:5}),1.0,1)
    zws_cesm_nino4 = wgt_areaave_Wrap(zws_cesm_anl(:,{-5:5},{160:210}),coswgt_cesm({-5:5}),1.0,1)
    bj_cesm = regCoef(sst_cesm_nino3(:),zws_cesm_nino4(:))
    printVarSummary(zws_cesm_anl)

    ;-----------------------------------------------------------------------------------------------------------

    data = addfile(FGOALSDIR+"t0225_ocn_S.gamil.h0.TAUX.060001-079912.nc","r") 
    zws_fgoals = data->TAUX(0:12*200-1,:,:)    ;N/m^2
    lat_fgoals = data->lat
    lon_fgoals = data->lon
    time_fgoals = data->time
    ;printVarSummary(zws_fgoals)
    data = addfile(FGOALSDIR+"ts_data/"+"t0225_ocn_S.gamil.h0.ts.060001-079912.nc","r")
    sst_fgoals  = data->TS(:,:,:)
    sst_fgoals =(/sst_fgoals-273.15/)
    sst_fgoals = where(sst_fgoals.le.-1.8,-1.8,sst_fgoals)
    d = addfile("$NCARG_ROOT/lib/ncarg/data/cdf/landsea.nc","r")   
    basemap = d->LSMASK
    lsm = landsea_mask(basemap,sst_fgoals&lat,sst_fgoals&lon)
    sst_fgoals = mask(sst_fgoals,conform(sst_fgoals,lsm,(/1,2/)).ge.1,False)
    delete(lsm)
    pi=4.*atan(1.0)
    rad=(pi/180.)
    coswgt_fgoals=cos(rad*lat_fgoals)
    coswgt_fgoals!0 = "lat" 
    coswgt_fgoals&lat= lat_fgoals 
    
    zws_fgoals =(/-zws_fgoals/)
 
    sst_fgoals_anl = rmMonAnnCycTLL(sst_fgoals)
    zws_fgoals_anl = rmMonAnnCycTLL(zws_fgoals)
    copy_VarCoords(zws_fgoals, zws_fgoals_anl)
    copy_VarCoords(sst_fgoals, sst_fgoals_anl)
    sst_fgoals_nino3 = wgt_areaave_Wrap(sst_fgoals_anl(:,{-5:5},{210:270}),coswgt_fgoals({-5:5}),1.0,1)
    zws_fgoals_nino4 = wgt_areaave_Wrap(zws_fgoals_anl(:,{-5:5},{160:210}),coswgt_fgoals({-5:5}),1.0,1)
    bj_fgoals = regCoef(sst_fgoals_nino3(:),zws_fgoals_nino4(:))
    printVarSummary(zws_fgoals_anl)
;------------------------------------------------------------------------------------------------------------------
    x = (/sst_obs_nino3,sst_cesm_nino3(0:44*12-1),sst_fgoals_nino3(0:44*12-1)/)
    y = (/zws_obs_nino4,zws_cesm_nino4(0:44*12-1),zws_fgoals_nino4(0:44*12-1)/)
    xyres               = True 
    xyres@gsnDraw       = False
    xyres@gsnFrame      = False
    xyres@vpHeightF     = 0.3
    xyres@vpWidthF      = 0.7
    xyres@gsnLeftString =""
    xyres@gsnRightString =""
    xyres@gsnStringFontHeightF = 0.022
    xyres@gsnLeftStringFontHeightF = 0.022
    xyres@gsnRightStringFontHeightF = 0.022
    xyres@tmXBLabelFontHeightF  = 0.020
    ;xyres@tmLabelAutoStride = True
    xyres@tmYLLabelFontHeightF  = 0.020
    xyres@tmYLMode = "Manual"
    xyres@tmYLTickStartF = -0.04
    xyres@tmYLTickEndF = 0.08
    xyres@tmYLTickSpacingF = 0.02
    xyres@xyMarkLineModes     = (/"Markers","Markers","Markers"/)
    xyres@xyMarkers           = (/16,10,4/)                     ; choose type of marker 
    xyres@xyMarkerColors       = (/"black","red","blue"/)                  ; Marker color
    xyres@xyMarkerSizeF       = 0.005  
    ;xyres@tmXBMode     = "Explicit"
    ;xyres@tmXBValues   =(/120.,140.,160.,180.,200.,220.,240.,260.,280./)
    ;xyres@tmXBLabels   = (/"120~S~o~N~E","140~S~o~N~E","160~S~o~N~E","180","160~S~o~N~W","140~S~o~N~W","120~S~o~N~W","100~S~o~N~W","90~S~o~N~W"/)
    ;xyres@xyLineThicknesses   = (/5.,5.,5./)      ; line thicknesses
    ;xyres@xyLineColors      = (/"black","blue","red"/)
    
    xyres@xyDashPatterns  = (/0.,0.,0./)      ; make all lines solid    
  
  ;---Turn off some tickmarks
    xyres@tmXTOn                 = False    
    xyres@tmYROn                 = False 
    xyres@tmXBOn                 = True
    xyres@tmXBLabelsOn           = True    

      
    xyres@gsnLeftString    = ""
    xyres@gsnCenterString  = "Atmospheric Bjerkes feedback"
    xyres@tiYAxisString    = "nino4 zws anomaly"+"   (N/m~S~2~N~)"
    xyres@tiXAxisString    = "nino3 SST anomaly" +"   (~S~o~N~C)"
    plot3=gsn_csm_xy(wks_BJ_fd,x,y,xyres)

    res_text                    = True                  ; text mods desired
    res_text@txFontHeightF      = 0.015                 ; change text size
    res_text@txJust             = "BottomLeft"          ; text justification
  
    res_marker                   = True                  
    res_marker@gsMarkerIndex     = 10                     ; choose type of marker 
    res_marker@gsMarkerColor    = "red"                 ; Marker color
    res_marker@gsMarkerSizeF     = 0.005 
    xx = (/-2.8,-2.7,-2.6/)
    yy = (/0.05,0.05,0.05/)
    dum1 = gsn_add_polymarker(wks_BJ_fd,plot3,xx,yy,res_marker)              ; add polyline
    dum2 = gsn_add_text(wks_BJ_fd,plot3,"CESM"+"  ~F33~m="+sprintf("%5.3f",bj_cesm),-2.5,0.045,res_text); add text
    
    yy = (/0.06,0.06,0.06/)
    res_marker@gsMarkerIndex          = 4   
    res_marker@gsMarkerColor    = "blue"
    dum3 = gsn_add_polymarker(wks_BJ_fd,plot3,xx,yy,res_marker)                ; add polymarker
    dum4 = gsn_add_text(wks_BJ_fd,plot3,"FGOALS-g2"+"  ~F33~m="+sprintf("%5.3f",bj_fgoals),-2.5,0.055,res_text)       ; add text
    
    yy = (/0.07,0.07,0.07/)
    res_marker@gsMarkerIndex         = 16   
    res_marker@gsMarkerColor    = "black"
    dum5 = gsn_add_polymarker(wks_BJ_fd,plot3,xx,yy,res_marker)                ; add polymarker
    dum6 = gsn_add_text(wks_BJ_fd,plot3,"ERA40"+"  ~F33~m="+sprintf("%5.3f",bj_obs),-2.5,0.065,res_text)     ; add text
    res_lines = True
    res_lines@gsLineDashPattern = 1. 
    ;res_lines@gsLineThicknessF  = 5.     
    xx1 = (/-3.,4./)
    yy1 = (/0.0,0.0/)
    dum7 = gsn_add_polyline(wks_BJ_fd,plot3,xx1,yy1,res_lines)
    xx1 = (/0.,0./)
    yy1 = (/-0.04,0.08/)
    dum8 = gsn_add_polyline(wks_BJ_fd,plot3,xx1,yy1,res_lines)
    draw(wks_BJ_fd)
    frame(wks_BJ_fd)









end
   