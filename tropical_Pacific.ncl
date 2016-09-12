load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
begin
    
    OUTPUT_TYPE = "png"
    OUTDIR      = "/home/reproduce/sunc/enso/plot/"
    CESMDIR     = "/home/reproduce/sunc/enso/CESM1-FASTCHEM_picontrol_atm/ts_data/"
    FGOALSDIR   = "/home/reproduce/sunc/enso/FGOALS-g2_picontrol_atm/ts_data/"
    HadISST     = "/home/reproduce/sunc/enso/obs_data/"

    wks_type = OUTPUT_TYPE

    wks_sst = gsn_open_wks(wks_type,OUTDIR+"tropical_sst")  
    wks_sst_std = gsn_open_wks(wks_type,OUTDIR+"tropical_sst_std")
    wks_pre = gsn_open_wks(wks_type,OUTDIR+"tropical_pre")
    wks_pre_rmse = gsn_open_wks(wks_type,OUTDIR+"tropical_pre_rmse")
    wks_zws = gsn_open_wks(wks_type,OUTDIR+"tropical_zonal_wind_stress")
    wks_zws_rmse = gsn_open_wks(wks_type,OUTDIR+"tropical_zonal_wind_stress_rmse")
    wks_eqt_sst = gsn_open_wks(wks_type,OUTDIR+"tropical_equatorial_sst")
    wks_eqt_ws = gsn_open_wks(wks_type,OUTDIR+"tropical_equatorial_wind_stress")

;-------------------------------------------- SST & SST standard deviation-------------------------------------------------------------------------------
    
    data = addfile(HadISST+"HadISST_sst.ts.187001-201503.nc","r")
    sst_obs  = data->sst(0:145*12-1,:,:)
    time_obs = data->time(0:145*12-1)
    lat_obs  = data->latitude
    sst_obs = where(sst_obs.le.-1.8,-1.8,sst_obs)

    d = addfile("$NCARG_ROOT/lib/ncarg/data/cdf/landsea.nc","r")   
    basemap = d->LSMASK
    lsm = landsea_mask(basemap,sst_obs&latitude,sst_obs&longitude)
    sst_obs = mask(sst_obs,conform(sst_obs,lsm,(/1,2/)).ge.1,False)
    delete(lsm)
    
    sst_obs = lonFlip(sst_obs)             ; longitude convert (-180:180) to (0:360) or reserve 
    ;sst_lon_conv = new((/dimsizes(time_obs),180,360/),"float")
    ;sst_lon_conv!2 = "longitude"
    ;sst_lon_conv&longitude=fspan(0.5,359.5,360)
    ;sst_lon_conv(:,:,0:179)=sst_obs(:,:,180:359)
    ;sst_lon_conv(:,:,180:359)=sst_obs(:,:,0:179)
    ;sst_sst_obs = sst_lon_conv
    ;sst_obs!0 = "time"
    ;sst_obs&time=time_obs
    ;sst_obs!1 = "latitude"
    ;sst_obs&latitude=lat_obs
    ;sst_obs!2 = "longitude"
    ;sst_obs&longitude=fspan(0.5,359.5,360)
    printVarSummary(sst_obs)


    pi=4.*atan(1.0)
    rad=(pi/180.)
    coswgt=cos(rad*sst_obs&latitude)
    coswgt!0 = "lat" 
    coswgt&lat= lat_obs 
    ;printVarSummary(sst_obs)
    ;printVarSummary(coswgt)
    sst_obs_anl = rmMonAnnCycTLL(sst_obs)
    sst_obs_trop = dim_avg_n_Wrap(sst_obs,0)
    sst_obs_trop_std = dim_stddev_n_Wrap(sst_obs_anl,0)
    ;printVarSummary(sst_obs_trop)
    sst_obs_trop_ave = wgt_areaave_Wrap(sst_obs_trop({-30:30},{120:290}),coswgt({-30:30}),1.0,1)
    sst_obs_trop_std_ave = wgt_areaave_Wrap(sst_obs_trop_std({-30:30},{120:290}),coswgt({-30:30}),1.0,1)
    sst_obs_trop_std_nino3 = wgt_areaave_Wrap(sst_obs_trop_std({-5:5},{210:270}),coswgt({-5:5}),1.0,1)
    sst_obs_trop_std_nino4 = wgt_areaave_Wrap(sst_obs_trop_std({-5:5},{160:210}),coswgt({-5:5}),1.0,1)
    delete(coswgt)



    data = addfile(CESMDIR+"CESM1-FASTCHEM_piControl_r1i1p1.ts.007001-029112.nc","r")
    sst_cesm  = data->ts(0:200*12-1,:,:)
    time_cesm = data->time(0:200*12-1)
    ;printVarSummary(sst_cesm)
    sst_cesm =(/sst_cesm-273.15/)
    sst_cesm = where(sst_cesm.le.-1.8,-1.8,sst_cesm)
    d = addfile("$NCARG_ROOT/lib/ncarg/data/cdf/landsea.nc","r")   
    basemap = d->LSMASK
    lsm = landsea_mask(basemap,sst_cesm&lat,sst_cesm&lon)
    sst_cesm = mask(sst_cesm,conform(sst_cesm,lsm,(/1,2/)).ge.1,False)
    delete(lsm)

    coswgt=cos(rad*sst_cesm&lat)
    coswgt!0 = "lat" 
    coswgt&lat= sst_cesm&lat 
    
    sst_cesm_anl = rmMonAnnCycTLL(sst_cesm)
    sst_cesm_trop = dim_avg_n_Wrap(sst_cesm,0)
    sst_cesm_trop_std = dim_stddev_n_Wrap(sst_cesm_anl,0)
    sst_cesm_trop_ave = wgt_areaave_Wrap(sst_cesm_trop({-30:30},{120:290}),coswgt({-30:30}),1.0,1)
    sst_cesm_trop_std_ave = wgt_areaave_Wrap(sst_cesm_trop_std({-30:30},{120:290}),coswgt({-30:30}),1.0,1)
    sst_cesm_trop_std_nino3 = wgt_areaave_Wrap(sst_cesm_trop_std({-5:5},{210:270}),coswgt({-5:5}),1.0,1)
    sst_cesm_trop_std_nino4 = wgt_areaave_Wrap(sst_cesm_trop_std({-5:5},{160:210}),coswgt({-5:5}),1.0,1)
    delete(coswgt)


    data = addfile(FGOALSDIR+"t0225_ocn_S.gamil.h0.ts.060001-079912.nc","r")
    sst_fgoals  = data->TS(:,:,:)
    time_fgoals = data->time(:)
    ;printVarSummary(sst_fgoals)
    sst_fgoals =(/sst_fgoals-273.15/)
    sst_fgoals = where(sst_fgoals.le.-1.8,-1.8,sst_fgoals)
    d = addfile("$NCARG_ROOT/lib/ncarg/data/cdf/landsea.nc","r")   
    basemap = d->LSMASK
    lsm = landsea_mask(basemap,sst_fgoals&lat,sst_fgoals&lon)
    sst_fgoals = mask(sst_fgoals,conform(sst_fgoals,lsm,(/1,2/)).ge.1,False)

    coswgt=cos(rad*sst_fgoals&lat)
    coswgt!0 = "lat" 
    coswgt&lat= sst_fgoals&lat 
 
    sst_fgoals_anl = rmMonAnnCycTLL(sst_fgoals)
    sst_fgoals_trop = dim_avg_n_Wrap(sst_fgoals,0)
    sst_fgoals_trop_std = dim_stddev_n_Wrap(sst_fgoals_anl,0)
    sst_fgoals_trop_ave = wgt_areaave_Wrap(sst_fgoals_trop({-30:30},{120:290}),coswgt({-30:30}),1.0,1)
    sst_fgoals_trop_std_nino3 = wgt_areaave_Wrap(sst_fgoals_trop_std({-5:5},{210:270}),coswgt({-5:5}),1.0,1)
    sst_fgoals_trop_std_nino4 = wgt_areaave_Wrap(sst_fgoals_trop_std({-5:5},{160:210}),coswgt({-5:5}),1.0,1)
    ;printVarSummary(sst_fgoals_trop)
    pi=4.*atan(1.0)
    rad=(pi/180.)

    plot1 = new(3,graphic)
    plot2 = new(3,graphic)
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
    ;res1@cnLevelSelectionMode = "ExplicitLevels"
    ;res1@cnLevels = (/0.2,0.5,1,2,3,4,5,6,7,8,9,10,12,14,17/)

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
     
    
    res1@gsnLeftString  ="HadISST"
    res1@gsnRightString ="mean="+sprintf("%5.3f", sst_obs_trop_ave)+"~S~o~N~C"
    plot1(0) = gsn_csm_contour_map_ce(wks_sst,sst_obs_trop,res1)

    do i=0,3
        dum(i) = gsn_add_polyline(wks_sst,plot1(0),xpts4(i:i+1),ypts4(i:i+1),resl4)
        dum2(i) = gsn_add_polyline(wks_sst,plot1(0),xpts3(i:i+1),ypts3(i:i+1),resl3)
    end do

    res1@gsnLeftString  ="CESM"
    res1@gsnRightString ="mean="+sprintf("%5.3f", sst_cesm_trop_ave)+"~S~o~N~C"
    plot1(1) = gsn_csm_contour_map_ce(wks_sst,sst_cesm_trop,res1)
    do i=0,3
        dum(4+i) = gsn_add_polyline(wks_sst,plot1(1),xpts4(i:i+1),ypts4(i:i+1),resl4)
        dum2(4+i) = gsn_add_polyline(wks_sst,plot1(1),xpts3(i:i+1),ypts3(i:i+1),resl3)
    end do


    res1@gsnLeftString  ="FGOALS-g2"
    res1@gsnRightString ="mean="+sprintf("%5.3f", sst_fgoals_trop_ave)+"~S~o~N~C"
    plot1(2) = gsn_csm_contour_map_ce(wks_sst,sst_fgoals_trop,res1)
    do i=0,3
        dum(8+i) = gsn_add_polyline(wks_sst,plot1(2),xpts4(i:i+1),ypts4(i:i+1),resl4)
        dum2(8+i) = gsn_add_polyline(wks_sst,plot1(2),xpts3(i:i+1),ypts3(i:i+1),resl3)
    end do
   
    
    pres=True
    pres@txString           ="Tropical Pacific SST"
    pres@gsnPanelLabelBar   = True       ; common label bar
    ;pres@gsnFrame            = False      ; don't advance frame yet
    ;pres1@lbOrientation       = "vertical"
    ;pres@lbLabelFontHeightF  = 0.015
    ;pres@pmLabelBarWidthF = 0.85
    ;pres@cnLevelSelectionMode = "ExplicitLevels"
    ;pres@cnLevels = (/0.2,0.5,1,2,3,4,5,6,7,8,9,10,12,14,17/)
    gsn_panel(wks_sst,plot1,(/3,1/),pres)

    res1@gsnLeftString   ="HadISST"
    res1@gsnCenterString ="nino4="+sprintf("%5.3f", sst_obs_trop_std_nino4)+"~S~o~N~C"
    res1@gsnRightString  ="nino3="+sprintf("%5.3f", sst_obs_trop_std_nino3)+"~S~o~N~C"
    plot2(0) = gsn_csm_contour_map_ce(wks_sst_std,sst_obs_trop_std,res1)
    do i=0,3
        dum(12+i) = gsn_add_polyline(wks_sst_std,plot2(0),xpts4(i:i+1),ypts4(i:i+1),resl4)
        dum2(12+i) = gsn_add_polyline(wks_sst_std,plot2(0),xpts3(i:i+1),ypts3(i:i+1),resl3)
    end do


    res1@gsnLeftString  ="CESM"
    res1@gsnCenterString ="nino4="+sprintf("%5.3f", sst_cesm_trop_std_nino4)+"~S~o~N~C"
    res1@gsnRightString  ="nino3="+sprintf("%5.3f", sst_cesm_trop_std_nino3)+"~S~o~N~C"
    plot2(1) = gsn_csm_contour_map_ce(wks_sst_std,sst_cesm_trop_std,res1)
    do i=0,3
        dum(16+i) = gsn_add_polyline(wks_sst_std,plot2(1),xpts4(i:i+1),ypts4(i:i+1),resl4)
        dum2(16+i) = gsn_add_polyline(wks_sst_std,plot2(1),xpts3(i:i+1),ypts3(i:i+1),resl3)
    end do

    res1@gsnLeftString  ="FGOALS-g2"
    res1@gsnCenterString ="nino4="+sprintf("%5.3f", sst_fgoals_trop_std_nino4)+"~S~o~N~C"
    res1@gsnRightString  ="nino3="+sprintf("%5.3f", sst_fgoals_trop_std_nino3)+"~S~o~N~C"
    plot2(2) = gsn_csm_contour_map_ce(wks_sst_std,sst_fgoals_trop_std,res1)
    do i=0,3
        dum(20+i) = gsn_add_polyline(wks_sst_std,plot2(2),xpts4(i:i+1),ypts4(i:i+1),resl4)
        dum2(20+i) = gsn_add_polyline(wks_sst_std,plot2(2),xpts3(i:i+1),ypts3(i:i+1),resl3)
    end do

    pres@txString           ="Tropical Pacific SST standard deviation"
    gsn_panel(wks_sst_std,plot2,(/3,1/),pres)
    

    delete(res1)
    delete(pres)
;------------------------------------------------------------------------------------------------------------------------------------
end  