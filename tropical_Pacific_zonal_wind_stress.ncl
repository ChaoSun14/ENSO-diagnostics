load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
begin
    
    OUTPUT_TYPE = "png"
    OUTDIR      = "/home/reproduce/sunc/enso/plot/"
    CESMDIR     = "/home/reproduce/sunc/enso/CESM1-FASTCHEM_picontrol_atm/"
    FGOALSDIR   = "/home/reproduce/sunc/enso/FGOALS-g2_picontrol_atm/"
    JRA55       = "/home/reproduce/sunc/AMWG/obs_data_20140804/"

    wks_type = OUTPUT_TYPE


    wks_zws = gsn_open_wks(wks_type,OUTDIR+"tropical_zonal_wind_stress")
    
    wks_eqt_zws = gsn_open_wks(wks_type,OUTDIR+"tropical_equatorial_zonal_wind_stress")

    data = addfile(JRA55+"JRA25_ANN_climo.nc","r") 
    zws_obs = data->TAUX(:,:,:)    ;N/m^2
    lat_obs = data->lat
    lon_obs = data->lon
    time_obs = data->time
    ;printVarSummary(zws_obs)

    pi=4.*atan(1.0)
    rad=(pi/180.)
    coswgt_obs=cos(rad*lat_obs)
    coswgt_obs!0 = "lat" 
    coswgt_obs&lat= lat_obs 

    zws_obs_mon = dim_avg_n_Wrap(zws_obs, 0)
    zws_obs_mon_eqt = dim_avg_n_Wrap(zws_obs_mon({-5:5},:), 0)
    zws_obs_mon_ave = wgt_areaave_Wrap(zws_obs_mon({-5:5},{120:290}),coswgt_obs({-5:5}),1.0,1)
    print(zws_obs_mon_ave)


    data = addfile(CESMDIR+"tauu_Amon_CESM1-FASTCHEM_piControl_r1i1p1_007001-029112.nc","r") 
    zws_cesm = data->tauu(0:12*200-1,:,:)    ;N/m^2
    lat_cesm = data->lat
    lon_cesm = data->lon
    time_cesm = data->time
    ;printVarSummary(zws_cesm)

    pi=4.*atan(1.0)
    rad=(pi/180.)
    coswgt_cesm=cos(rad*lat_cesm)
    coswgt_cesm!0 = "lat" 
    coswgt_cesm&lat= lat_cesm 

    zws_cesm_mon = dim_avg_n_Wrap(zws_cesm, 0)
    zws_cesm_mon_eqt = dim_avg_n_Wrap(zws_cesm_mon({-5:5},:), 0)
    zws_cesm_mon_ave = wgt_areaave_Wrap(zws_cesm_mon({-5:5},{120:290}),coswgt_cesm({-5:5}),1.0,1)
    print(zws_cesm_mon_ave)

    data = addfile(FGOALSDIR+"t0225_ocn_S.gamil.h0.TAUX.060001-079912.nc","r") 
    zws_fgoals = data->TAUX(0:12*200-1,:,:)    ;N/m^2
    lat_fgoals = data->lat
    lon_fgoals = data->lon
    time_fgoals = data->time
    ;printVarSummary(zws_fgoals)

    pi=4.*atan(1.0)
    rad=(pi/180.)
    coswgt_fgoals=cos(rad*lat_fgoals)
    coswgt_fgoals!0 = "lat" 
    coswgt_fgoals&lat= lat_fgoals 
    
    zws_fgoals =(/-zws_fgoals/)
    zws_fgoals_mon = dim_avg_n_Wrap(zws_fgoals, 0)
    zws_fgoals_mon_eqt = dim_avg_n_Wrap(zws_fgoals_mon({-5:5},:), 0)
    zws_fgoals_mon_ave = wgt_areaave_Wrap(zws_fgoals_mon({-5:5},{120:290}),coswgt_fgoals({-5:5}),1.0,1)
    print(zws_fgoals_mon_ave)

    plot1 = new(3,graphic)
    res1=True
    res1@gsnDraw        = False            ; don't draw yet
    res1@gsnFrame       = False
    ;res1@mpCenterLonF   = 180.
    res1@vpHeightF      = 0.3
    res1@vpWidthF       = 0.7
    res1@gsnStringFontHeightF      = 0.022
    res1@gsnLeftStringFontHeightF  = 0.022
    res1@gsnRightStringFontHeightF = 0.022
    res1@tmXBLabelFontHeightF  = 0.022
    res1@tmLabelAutoStride     = True
    res1@tmYLLabelFontHeightF  = 0.022
    res1@cnFillOn  = True
    res1@cnLinesOn = False
    res1@lbLabelBarOn = False
    res1@mpMinLatF  = -30.
    res1@mpMaxLatF  = 30.
    res1@mpMinLonF  = 120.
    res1@mpMaxLonF  = 290.
    res1@mpCenterLonF   = 205.
    ;res1@gsnSpreadColors = True
    ;res1@gsnSpreadColorStart = -60
    res1@lbOrientation       = "vertical"
    res1@cnLevelSelectionMode = "ManualLevels"
    res1@cnLevelSpacingF = 0.02
    res1@cnMinLevelValF = -0.16
    res1@cnMaxLevelValF = 0.16

    dum = new(24,graphic)
    dum2= new(24,graphic)
    ypts4 = (/5,5,-5,-5,5/)
    xpts4 = (/160,210,210,160,160/)
    ypts3 = (/5,5,-5,-5,5/)
    xpts3 = (/210,270,270,210,210/)

    resl4 = True
    resl4@gsLineColor = "red"
    resl4@gsLineThicknessF = 3.
;   resl4@gsLineLabelString = "Nino 4"
    resl3 = True
    resl3@gsLineColor = "black"
    resl3@gsLineThicknessF = 3.
;   resl3@gsLineLabelString = "Nino 3"
     
    
    res1@gsnLeftString  ="JRA55"
    res1@gsnRightString ="mean="+sprintf("%3.2f", zws_obs_mon_ave)+"N/m~S~2~N~"
    plot1(0) = gsn_csm_contour_map_ce(wks_zws,zws_obs_mon,res1)

    do i=0,3
        dum(i) = gsn_add_polyline(wks_zws,plot1(0),xpts4(i:i+1),ypts4(i:i+1),resl4)
        dum2(i) = gsn_add_polyline(wks_zws,plot1(0),xpts3(i:i+1),ypts3(i:i+1),resl3)
    end do

    res1@gsnLeftString  ="CESM"
    res1@gsnCenterString = "RMSE="+sprintf("%3.2f", abs(zws_cesm_mon_ave/zws_obs_mon_ave-1.))
    res1@gsnRightString ="mean="+sprintf("%3.2f", zws_cesm_mon_ave)+"N/m~S~2~N~"
    plot1(1) = gsn_csm_contour_map_ce(wks_zws,zws_cesm_mon,res1)
    do i=0,3
        dum(4+i) = gsn_add_polyline(wks_zws,plot1(1),xpts4(i:i+1),ypts4(i:i+1),resl4)
        dum2(4+i) = gsn_add_polyline(wks_zws,plot1(1),xpts3(i:i+1),ypts3(i:i+1),resl3)
    end do


    res1@gsnLeftString  ="FGOALS-g2"
    res1@gsnCenterString = "RMSE="+sprintf("%3.2f", abs(zws_fgoals_mon_ave/zws_obs_mon_ave-1.))
    res1@gsnRightString ="mean="+sprintf("%3.2f", zws_fgoals_mon_ave)+"N/m~S~2~N~"
    plot1(2) = gsn_csm_contour_map_ce(wks_zws,zws_fgoals_mon,res1)
    do i=0,3
        dum(8+i) = gsn_add_polyline(wks_zws,plot1(2),xpts4(i:i+1),ypts4(i:i+1),resl4)
        dum2(8+i) = gsn_add_polyline(wks_zws,plot1(2),xpts3(i:i+1),ypts3(i:i+1),resl3)
    end do
   
    
    pres=True
    pres@txString           ="Tropical Pacific zonal wind stress"
    pres@gsnPanelLabelBar   = True       ; common label bar
        ;pres@gsnFrame            = False      ; don't advance frame yet
        ;pres1@lbOrientation       = "vertical"
        ;pres@lbLabelFontHeightF  = 0.015
        ;pres@pmLabelBarWidthF = 0.85
    ;pres@cnLevelSelectionMode = "ExplicitLevels"
    ;pres@cnLevels = (/0.2,0.5,1,2,3,4,5,6,7,8,9,10,12,14,17/)
    gsn_panel(wks_zws,plot1,(/3,1/),pres)

;-------------------- equatorical zws ------------------------------------------------------------------------------------------------------

   zws_obs_mon_eqt_tofgoals = linint1_n_Wrap(zws_obs_mon_eqt&lon, zws_obs_mon_eqt, True, lon_fgoals, 0, 0)
   zws_cesm_mon_eqt_tofgoals = linint1_n_Wrap(lon_cesm, zws_cesm_mon_eqt, True, lon_fgoals, 0, 0)
   ;print(zws_obs_mon_eqt_tofgoals)
   x = lon_fgoals({120:290})
   ;print(x)
   y = (/zws_obs_mon_eqt_tofgoals({120:290}),zws_cesm_mon_eqt_tofgoals({120:290}),zws_fgoals_mon_eqt({120:290})/)
    
    xyres=True
    xyres               =True 
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
    ;xyres@tmXBMode = "Manual"
    ;xyres@tmXBTickStartF = 120.
    ;xyres@tmXBTickEndF = 280.
    ;xyres@tmXBTickSpacingF = 20.
    xyres@trXMinF = 120.
    xyres@trXMaxF = 280.
    xyres@tmXBMode = "Explicit"
    xyres@tmXBValues =(/120.,140.,160.,180.,200.,220.,240.,260.,280./)
    xyres@tmXBLabels = (/"120~S~o~N~E","140~S~o~N~E","160~S~o~N~E","180","160~S~o~N~W","140~S~o~N~W","120~S~o~N~W","100~S~o~N~W","90~S~o~N~W"/)
    xyres@xyLineThicknesses   = (/5.,5.,5./)      ; line thicknesses
    xyres@xyLineColors      = (/"black","red","blue"/)
    
    xyres@xyDashPatterns  = (/0.,0.,0./)      ; make all lines solid    
  
  ;---Turn off some tickmarks
    xyres@tmXTOn                 = False    
    xyres@tmYROn                 = False 
    xyres@tmXBOn                 = True
    xyres@tmXBLabelsOn           = True    

      
    xyres@gsnLeftString    = "b) Equatorial zonal wind stress"
    xyres@gsnCenterString  = ""
    xyres@tiYAxisString    = "Zonal wind stress"+"   (N/m~S~2~N~)"
    xyres@tiXAxisString    = "Longitude" 
    plot3=gsn_csm_xy(wks_eqt_zws,x,y,xyres)

    res_text                    = True                  ; text mods desired
    res_text@txFontHeightF      = 0.015                 ; change text size
    res_text@txJust             = "BottomLeft"          ; text justification
  
    res_lines                   = True                  ; polyline mods desired
    res_lines@gsLineDashPattern = 0.                    ; solid line
    res_lines@gsLineThicknessF  = 5.                    ; line thicker
    res_lines@gsLineColor       = "red"                 ; line color
    xx = (/125,140/)
    yy = (/-0.09,-0.09/)
    dum1 = gsn_add_polyline(wks_eqt_zws,plot3,xx,yy,res_lines)              ; add polyline
    dum2 = gsn_add_text(wks_eqt_zws,plot3,"CESM",145,-0.09,res_text); add text
    
    yy = (/-0.08,-0.08/)
    res_lines@gsLineColor       = "blue"                                 ; change to blue
    dum3 = gsn_add_polyline(wks_eqt_zws,plot3,xx,yy,res_lines)                ; add polyline
    dum4 = gsn_add_text(wks_eqt_zws,plot3,"FGOALS-g2",145,-0.08,res_text)       ; add text
    
    yy = (/-0.07,-0.07/)
    res_lines@gsLineColor       = "black"                                ; change to black
    dum5 = gsn_add_polyline(wks_eqt_zws,plot3,xx,yy,res_lines)                ; add polyline
    dum6 = gsn_add_text(wks_eqt_zws,plot3,"JRA55",145,-0.07,res_text)     ; add text
    res_lines@gsLineDashPattern = 1. 
    xx = (/120,280/)
    yy = (/0.0,0.0/)
    dum7 = gsn_add_polyline(wks_eqt_zws,plot3,xx,yy,res_lines)
    draw(wks_eqt_zws)
    frame(wks_eqt_zws)



end

